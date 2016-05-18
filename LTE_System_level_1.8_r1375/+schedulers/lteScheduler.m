classdef lteScheduler < handle
    % Implements common methods needed by any implementation of an LTE
    % scheduler (eg. Round Robin, Best CQI...)
    % (c) Josep Colom Ikuno, INTHFT, 2009
    
    properties
        % Scheduler name
        name = 'Generic LTE scheduler superclass';
        
        % Where this scheduler is attached
        attached_eNodeB
        
        % Copy of the CQI tables. % In order to avoid getting this every TTI
        CQI_range
        CQI_tables
        CQIs_efficiency
        
        % The algorithm used to average the several SINRs into a single TB SINR. Typically MIESM
        SINR_averager
        
        % BLER data and CQI mapping
        BLER_curves
        CQI_mapper
        
        % The TX mode to be used by default
        default_tx_mode
        
        % To control power allocation
        max_power
        
        % Target system BLER
        target_BLER = 0.1;
        
        % Set this option to false if you are doing SL-LL validation: if not, this leads to problems with the simulated
        % BLER. For CQI 1 one could not simulate BLERs worse than 0.1 as those RBs are just skipped
        skip_null_CQIs = true;
        
        % Genie information (eg: all of the eNodeBs and UEs)
        genie
        
        % UE trace: from where the historical information (ie. UE throughput) is extracted
        UE_traces
        % In order to know how many TTIs to skip when looking for throughput
        % (in reality you would not know instantly whether the TB was received correctly or not)
        feedback_delay_TTIs
        
        % Network clock. Tells the scheduler in what TTI he is
        clock
        
        % RB grid attached to this eNodeB
        RB_grid
        
        % Other things
        av_const
        fairness
        k
        d
        MUMIMO
        overhead_ref
        overhead_sync
        
        % Whether the precoding is precalculated at a trace or done at run time
        runtime_precoding
        
        % Rel 8 codebook
        Rel8_codebook
        
        % To use with FFR (or similar) scheduling
        fractional_BW_allocation = false; % Whether the allocation is normal or is part of a fractional one (i.e. a parent scheduler with children schedulers)
        fractional_allocation             % The actual allocation
    end
    
    methods(Abstract)
        schedule_users(obj,attached_UEs,UE_feedback)
        add_UE(obj,UE_id)
        remove_UE(obj,UE_id)
    end
    
    methods
        
        % Class constructor
        function obj = lteScheduler(scheduler_params,attached_eNodeB)
            obj.attached_eNodeB = attached_eNodeB;
            the_CQI_range       = LTE_common_get_CQI_params(scheduler_params,'range');
            obj.CQI_tables      = LTE_common_get_CQI_params(scheduler_params,the_CQI_range(1)); % To initialize the struct
            
            for i_= the_CQI_range(1):the_CQI_range(2)
                obj.CQI_tables(i_) = LTE_common_get_CQI_params(scheduler_params,i_);
            end
            
            obj.CQI_range       = the_CQI_range(1):the_CQI_range(2);
            obj.CQIs_efficiency = [0 [obj.CQI_tables.efficiency]]; % Take note of CQI 0 also
            obj.default_tx_mode = scheduler_params.default_tx_mode;
            obj.max_power       = scheduler_params.max_power;
            obj.clock           = attached_eNodeB.parent_eNodeB.clock;
            obj.RB_grid         = attached_eNodeB.RB_grid;
            obj.av_const        = scheduler_params.av_window;
            obj.k               = scheduler_params.k;
            obj.d               = scheduler_params.d;
            obj.overhead_ref    = scheduler_params.overhead_ref;
            obj.overhead_sync   = scheduler_params.overhead_sync;
            
            % This is not used by every scheduler
            if isfield(scheduler_params,'fairness')
                obj.fairness = scheduler_params.fairness;
            else
                obj.fairness = [];
            end
        end
        
        % Print some info
        function print(obj)
            fprintf('%s\n',obj.name);
        end
        
        % Find the optimum CQI values for a set of N codewords. It is assumed that the CQI vector/matrix is NOT empty. Each column contains the SINRs of one codeword
        function [assigned_CQI, predicted_BLER, predicted_SINR_dB, predicted_Is, predicted_Is_min] =...
                get_optimum_CQIs(obj,CQIs_to_average_all)
            % Preallocation
            nCodewords        = size(CQIs_to_average_all,2);
            assigned_CQI      = zeros(nCodewords,1);
            predicted_BLER    = zeros(nCodewords,1);
            predicted_SINR_dB = zeros(nCodewords,1);
            
            UE_estimated_SINR_dB_all = obj.CQI_mapper.CQI_to_SINR(CQIs_to_average_all);
            UE_estimated_MI_all      = obj.SINR_averager.SINR_to_I(UE_estimated_SINR_dB_all,obj.CQI_range);
            sum_UE_estimated_MI_all  = sum(UE_estimated_MI_all,2);
            
            % Don't use the SINR value on the lower bound of a CQI interval
            % but something inbetween (otherwise this is too conservative)
            % - as in the Link Level Simulator (heuristically optimized)
            % SINR_temp1 = obj.CQI_mapper.CQI_to_SINR(CQIs_to_average);
            % SINR_temp2 = obj.CQI_mapper.CQI_to_SINR(CQIs_to_average+1);
             
            % UE_estimated_SINR_dB = SINR_temp1+(SINR_temp2-SINR_temp1)/2.5;
            % UE_estimated_SINR_linear = 10.^(0.1*UE_estimated_SINR_dB);
            
            SINR_MCS_dependent_dB = reshape(obj.SINR_averager.average_I(UE_estimated_MI_all,1,obj.CQI_range),[nCodewords length(obj.CQI_range)])';
            
            for cw_=1:nCodewords
                predicted_BLERs = obj.BLER_curves.get_BLER_CQI(obj.CQI_range,SINR_MCS_dependent_dB(:,cw_));
                
                % Objective is the closest smaller or equal to 10% BLER (BLER 0 is preferred to BLER 1)
                if predicted_BLERs(end) == 0
                    % Case of a very good channel
                    cw_assigned_CQI = 15;
                elseif predicted_BLERs(1) >= obj.target_BLER
                    % Case of a bad channel
                    cw_assigned_CQI = 1;
                else
                    % Case in the middle
                    abs_diffs = predicted_BLERs-obj.target_BLER;
                    abs_diffs = round(abs_diffs*1000)/1000; % To avoid small statistical mistakes in the BLER plots. No change assuming that the target BLER is in the order of 10%
                    cw_assigned_CQI = find(abs_diffs<=0,1,'last');
                end
                assigned_CQI(cw_)      = cw_assigned_CQI;
                predicted_BLER(cw_)    = predicted_BLERs(cw_assigned_CQI);
                predicted_SINR_dB(cw_) = SINR_MCS_dependent_dB(cw_assigned_CQI,cw_);     % effective logarithmic SINR
            end
            predicted_Is     = mean(sum_UE_estimated_MI_all(:,1,cw_assigned_CQI));
            predicted_Is_min = min(sum_UE_estimated_MI_all(:,1,cw_assigned_CQI));
        end
        
        function CQI_spectral_efficiency = get_spectral_efficiency(obj,CQI_matrix)
            % Returns the spectral efficiencies related to the CQIs in the
            % matrix. In case of unquantized feedback, it floors the CQI
            % values to the nearest usable CQI
            CQI_idx                 = floor(CQI_matrix)+1; % Get CQI indexes
            CQI_idx(CQI_idx<1)      = 1;       % Clip extremes
            CQI_idx(CQI_idx>16)     = 16;
            sizes_CQI_feedback      = size(CQI_idx);
            CQI_idx_vector          = uint16(CQI_idx(:)); % To avoid an error caused by Matlab sometimes saying those were not integers
            CQI_spectral_efficiency = reshape(obj.CQIs_efficiency(CQI_idx_vector),sizes_CQI_feedback);
        end
        
        function UE_resource_assignment = get_max_UEs(obj,metric_matrix,UE_assignment)
            % Scans a nUExM matrix and returns for each column the index of
            % the highest value. In case more than one UE has the same
            % metric value, one of these ones is randomly selected. Assume
            % all values to be non-negative and 0 as "out-of-range"
            max_metric = max(metric_matrix,[],1);
            UE_resource_assignment = zeros(length(max_metric),1);
            for rb_=1:length(max_metric)
                if max_metric(rb_)~=0
                    candidates = find(metric_matrix(:,rb_)==max_metric(rb_));
                    index = ceil(rand*length(candidates));
                    if index==0
                        index = 1;
                    end
                    UE_resource_assignment(rb_) = UE_assignment(candidates(index)); % Choose a random UE from the set
                end
            end
        end
        
        function [N_assigned_RBs, CQIs_to_average_all, UE_scheduled, new_UE_RB_map] = ...
                filter_out_zero_RBs_and_get_CQIs(obj,...
                RB_grid,...
                nCodewords,...
                UE_CQI_feedback,...
                current_UE)
            
            % Do not use RBs with a CQI of 0 (they are lost)
            UE_RBs = RB_grid.user_allocation==current_UE.id;
            if obj.skip_null_CQIs
                if nCodewords == 1
                    if (size(UE_CQI_feedback,1)==1) % SISO case, to ensure correct dimensions
                        UE_CQI_feedback = reshape(UE_CQI_feedback,[],1);
                    end
                    zero_CQIs     = (UE_CQI_feedback(:,1)<1);
                    non_valid_RBs = UE_RBs & zero_CQIs;  % RBs to filter out
                    valid_RBs     = UE_RBs & ~zero_CQIs; % CQIs that will be averaged
                    RB_grid.user_allocation(non_valid_RBs) = 0;
                    CQIs_to_average_all = UE_CQI_feedback(valid_RBs);
                else
                    % For the case where more than 1 codewords are sent, all CWs must have a CQI >0
                    zero_CQIs     = sum(UE_CQI_feedback<1,2)>=1;
                    non_valid_RBs = UE_RBs & zero_CQIs; % RBs to filter out
                    valid_RBs     = UE_RBs & ~zero_CQIs; % CQIs that will be averaged
                    RB_grid.user_allocation(non_valid_RBs) = 0;
                    CQIs_to_average_all = UE_CQI_feedback(valid_RBs,:);
                end
                new_UE_RB_map  = valid_RBs;
                N_assigned_RBs = sum(valid_RBs);
            else
                if nCodewords == 1
                    CQIs_to_average_all = UE_CQI_feedback(UE_RBs);
                else
                    CQIs_to_average_all = UE_CQI_feedback(UE_RBs,:);
                end
                new_UE_RB_map  = UE_RBs;
                N_assigned_RBs = sum(UE_RBs);
            end
            
            if isempty(CQIs_to_average_all)
                UE_scheduled   = false;
            else
                UE_scheduled = true;
            end
        end
        
        function schedule_users_common(obj,attached_UEs,UE_feedback,current_TTI,tx_mode)
            % NOTE: since for the zero-delay case, RI and nLayers of the
            % transmission may not match, nCodewords is not used here
            
            % Common tasks for all schedulers
            the_RB_grid        = obj.RB_grid;
            RB_grid_size_bits  = 0;
            max_codewords      = 2;
            predicted_UE_BLERs = NaN(max_codewords,length(attached_UEs));
            nRBs               = the_RB_grid.n_RB;
            
            % Homogeneous power allocation
            if ~isempty(attached_UEs)
                if isprop(the_RB_grid, 'max_ues')
                    power_per_rb = obj.max_power / the_RB_grid.n_RB;
                    scheduled_idxs = the_RB_grid.user_allocation>0;
                     the_RB_grid.power_allocation = scheduled_idxs.*power_per_rb./kron(ones(1, the_RB_grid.max_ues),sum(scheduled_idxs, 2));
                else
                    the_RB_grid.power_allocation(:) = obj.max_power / the_RB_grid.n_RB;
                end
            end
            
            % Clean up assigned precoders
            if obj.runtime_precoding
                empty_precoder.W        = complex(0);
                the_RB_grid.precoder(:) = empty_precoder;
            end
            
            % Add here a tx_mode-dependent feedback calculation
            % needs for mode 400 to translate the Hs into "standard"
            % feedback with CQI a precoder struct
            
            % Continue UE common scheduling procedures
            for u_=1:length(attached_UEs)
                current_UE = attached_UEs(u_);

                % If the feedback is present: process normally
                total_RB_power = the_RB_grid.power_allocation + the_RB_grid.power_allocation_signaling;
                if UE_feedback.feedback_received(u_)
                    tx_mode = UE_feedback.tx_mode(u_);
                    if tx_mode==1||tx_mode==2
                        % SIXO, TxD
                        nLayers    = 1;
                        nCodewords = 1;
                    elseif tx_mode==3||tx_mode==4||tx_mode==9 % OLSM, CLSM, and others
                        nLayers    = UE_feedback.RI(u_);
                        nCodewords = min(2,nLayers);
                    else
                        nLayers    = UE_feedback.RI(u_);
                        nCodewords = min(2,nLayers);
                    end
                    
                    UE_CQI_feedback = UE_feedback.CQI(:,1:nCodewords,u_);
                    
                    % Do not use RBs with a CQI of 0 (they are lost).
                    % This function also averages the CQIs and assigns an overall TB CQI with predicted BLER < 10%
                    [num_assigned_RB, CQIs_to_average_all, UE_scheduled, RB_map] = ...
                        obj.filter_out_zero_RBs_and_get_CQIs(...
                        the_RB_grid,...
                        nCodewords,...
                        UE_CQI_feedback,...
                        current_UE);
                    
                    if UE_scheduled
                        % Calculate the optimum CQI (MCS)
                        [assigned_CQI, predicted_UE_BLERs(1:nCodewords,u_), estimated_TB_SINR, Is, I_min] = ...
                            obj.get_optimum_CQIs(CQIs_to_average_all);

                        % UE RB allocation
                        if ~obj.fractional_BW_allocation
                            attached_UEs(u_).eNodeB_signaling.assigned_RB_map = RB_map;
                        else
                            UE_full_RB_map = obj.fractional_allocation(:);
                            UE_full_RB_map(obj.fractional_allocation) = RB_map;
                            attached_UEs(u_).eNodeB_signaling.assigned_RB_map = UE_full_RB_map;
                        end

                        attached_UEs(u_).eNodeB_signaling.assigned_power               = sum(total_RB_power(RB_map));
                        attached_UEs(u_).eNodeB_signaling.tx_mode                      = tx_mode;
                        attached_UEs(u_).eNodeB_signaling.TB_CQI                       = assigned_CQI;
                        attached_UEs(u_).eNodeB_signaling.nCodewords                   = nCodewords;
                        attached_UEs(u_).eNodeB_signaling.nLayers                      = nLayers;
                        attached_UEs(u_).eNodeB_signaling.genie_TB_SINR                = estimated_TB_SINR;
                        attached_UEs(u_).eNodeB_signaling.adaptive_RI.avg_MI           = Is;
                        attached_UEs(u_).eNodeB_signaling.adaptive_RI.min_MI           = I_min;
                        attached_UEs(u_).eNodeB_signaling.adaptive_RI.RBs_for_feedback = true(1,nRBs);
                        
                        % Assign precoder (if applicable)
                        if obj.runtime_precoding
                            % CLSM
                            if  tx_mode == 4 || tx_mode == 9
                                layer_codebook = obj.Rel8_codebook{nLayers,obj.attached_eNodeB.total_nTX,tx_mode};
                                RB_list        = find(RB_map);
                                for rb_idx=1:length(RB_list)
                                    rb_ = RB_list(rb_idx);
                                    the_RB_grid.precoder(rb_).W = layer_codebook.W(:,:,UE_feedback.PMI(rb_,u_));
                                    if isreal(the_RB_grid.precoder(rb_).W)
                                        the_RB_grid.precoder(rb_).W = complex(the_RB_grid.precoder(rb_).W);
                                    end
                                end
                            elseif  sum(tx_mode==[4 5 6])
                                layer_codebook = obj.Rel8_codebook{nLayers,obj.attached_eNodeB.total_nTX,tx_mode};
                                RB_list        = find(RB_map);
                                for rb_idx=1:length(RB_list)
                                    rb_ = RB_list(rb_idx);
                                    the_RB_grid.precoder(rb_).W = layer_codebook.W(:,:,UE_feedback.PMI(rb_,u_));
                                    if isreal(the_RB_grid.precoder(rb_).W)
                                        the_RB_grid.precoder(rb_).W = complex(the_RB_grid.precoder(rb_).W);
                                    end
                                end
                                
                            % In these modes, the precoder is calculated at
                            % the eNodeB in the feedback_calculation.process_UE_feedback_at_eNodeB
                            % function
                            elseif tx_mode>=100
                                RB_list        = find(RB_map);
                                for rb_idx=1:length(RB_list)
                                    rb_ = RB_list(rb_idx);
                                    the_RB_grid.precoder(rb_).W = UE_feedback.precoder(u_).W(:,:,rb_);
                                    if isreal(the_RB_grid.precoder(rb_).W)
                                        the_RB_grid.precoder(rb_).W = complex(the_RB_grid.precoder(rb_).W);
                                    end
                                end
                            end
                        end
                    else
                        attached_UEs(u_).eNodeB_signaling.assigned_RB_map              = [];
                        attached_UEs(u_).eNodeB_signaling.assigned_power               = 0;
                        attached_UEs(u_).eNodeB_signaling.tx_mode                      = 0;
                        attached_UEs(u_).eNodeB_signaling.TB_CQI                       = 0;
                        attached_UEs(u_).eNodeB_signaling.nCodewords                   = 0;
                        attached_UEs(u_).eNodeB_signaling.nLayers                      = 0;
                        attached_UEs(u_).eNodeB_signaling.genie_TB_SINR                = NaN;
                        attached_UEs(u_).eNodeB_signaling.adaptive_RI.avg_MI           = [];
                        attached_UEs(u_).eNodeB_signaling.adaptive_RI.std_MI           = [];
                        attached_UEs(u_).eNodeB_signaling.adaptive_RI.RBs_for_feedback = true(1,nRBs);
                    end
                    
                else
                    % If the feedback is not present: assign a default CQI value of 1 with rank one.
                    UE_scheduled    = true;
                    RB_map          = the_RB_grid.user_allocation==current_UE.id;
                    num_assigned_RB = sum(RB_map);
                    nCodewords      = 1;
                    nLayers         = 1;
                    
                    % Set default precoder
                    if obj.runtime_precoding
                        switch tx_mode
                            case {4,5,6,400} % CLSM
                                % Assign the first precoder of rank one 
                                default_W = obj.Rel8_codebook{nLayers,obj.attached_eNodeB.total_nTX,4}.W(:,:,1);
                            case 9  % Up-to-8-Layer transmission
                                default_W = obj.Rel8_codebook{nLayers,obj.attached_eNodeB.total_nTX,9}.W(:,:,1);
                            otherwise
                                % No precoder. Set as 0 so that "[]" would not give us unexpected outputs
                                default_W = 0;
                        end
                    end
                    
                    % Assign a default precoder
                    if obj.runtime_precoding
                        default_precoder.W = default_W;
                        if isreal(default_precoder.W)
                            default_precoder.W = complex(default_precoder.W);
                        end
                        the_RB_grid.precoder(RB_map') = default_precoder;
                    end
                    
                    % Signal down the user CQI assignment
                    if ~obj.fractional_BW_allocation
                        attached_UEs(u_).eNodeB_signaling.assigned_RB_map = RB_map;
                    else
                        UE_full_RB_map = obj.fractional_allocation(:);
                        UE_full_RB_map(obj.fractional_allocation) = RB_map;
                        attached_UEs(u_).eNodeB_signaling.assigned_RB_map = UE_full_RB_map;
                    end
                    
                    attached_UEs(u_).eNodeB_signaling.assigned_power               = sum(total_RB_power(RB_map));
                    attached_UEs(u_).eNodeB_signaling.tx_mode                      = tx_mode;
                    attached_UEs(u_).eNodeB_signaling.TB_CQI(1:nCodewords)         = 1;
                    attached_UEs(u_).eNodeB_signaling.nCodewords                   = nCodewords;
                    attached_UEs(u_).eNodeB_signaling.nLayers                      = nLayers;
                    attached_UEs(u_).eNodeB_signaling.genie_TB_SINR                = NaN;
                    attached_UEs(u_).eNodeB_signaling.adaptive_RI.avg_MI           = [];
                    attached_UEs(u_).eNodeB_signaling.adaptive_RI.std_MI           = [];
                    attached_UEs(u_).eNodeB_signaling.adaptive_RI.RBs_for_feedback = true(1,nRBs);
                    predicted_UE_BLERs(u_) = 0; % Dummy value to avoid a NaN
                end
                
                % Calculate TB size
                if UE_scheduled
                    TB_CQI_params    = obj.CQI_tables(attached_UEs(u_).eNodeB_signaling.TB_CQI);
                    modulation_order = [TB_CQI_params.modulation_order];
                    coding_rate      = [TB_CQI_params.coding_rate_x_1024]/1024;
                    
                    if mod(current_TTI-1,5) || strcmp(class(obj.attached_eNodeB.scheduler),'schedulers.ffrScheduler') % TB without sync symbols. Solution (for now) for FFR
                        % The factor of two is because there are two time slots per subframe, 24 CRC bits are attached
                        % Segmentation prior to channel coding NOT taken into account.
                        TB_size_bits = max(8*round(1/8*(the_RB_grid.sym_per_RB_nosync .* num_assigned_RB .* modulation_order .* coding_rate * 2))-24,0);
                    else % TB with sync symbols
                        sync_pos = false(size(the_RB_grid.user_allocation));
                        sync_pos(floor(length(sync_pos)/2)-2:floor(length(sync_pos)/2)+3) = true;
                        sync_RBs = sum((the_RB_grid.user_allocation==attached_UEs(u_).id) .* sync_pos);
                        non_sync_RBs = sum(the_RB_grid.user_allocation==attached_UEs(u_).id)-sync_RBs;
                        TB_size_bits = max(8*round(1/8*(the_RB_grid.sym_per_RB_sync .* sync_RBs + the_RB_grid.sym_per_RB_nosync * non_sync_RBs) .* modulation_order .* coding_rate * 2)-24,0);
                    end
                    
                    % Adapt TB size calculation for layer mapping (TS 36.211)
                    switch nLayers
                        case 3
                            TB_size_bits(2) = TB_size_bits(2)*2; % The second CW comprises the second and third layers. The first one just one layer
                        case 4
                            TB_size_bits = TB_size_bits*2;       % Each CW comprises two layers
                    end
                else
                    num_assigned_RB = 0;
                    TB_size_bits = 0;
                end
                
                % Update accumulated RB size (bits)
                RB_grid_size_bits = RB_grid_size_bits + TB_size_bits;
                
                % Fill the UE-specific signaling
                attached_UEs(u_).eNodeB_signaling.num_assigned_RBs = num_assigned_RB;
                attached_UEs(u_).eNodeB_signaling.TB_size          = TB_size_bits;
                attached_UEs(u_).eNodeB_signaling.rv_idx           = 0;
                % attached_UEs(u_).eNodeB_signaling.
                
                % Traffic model-related code (packet parts)
                for cw_ = 1:length(TB_size_bits)
                    if TB_size_bits(cw_) ~= 0
                        if (~attached_UEs(u_).traffic_model.is_fullbuffer) && (strcmp(attached_UEs(u_).traffic_model.type,'voip') || strcmp(attached_UEs(u_).traffic_model.type,'video') || strcmp(attached_UEs(u_).traffic_model.type,'gaming') || strcmp(attached_UEs(u_).traffic_model.type,'MLaner'))
                            packet_parts = attached_UEs(u_).traffic_model.decrease_packets(TB_size_bits(cw_));
                            if ~isempty(packet_parts)
                                attached_UEs(u_).eNodeB_signaling.N_used_bits(cw_) = sum(packet_parts.get_size);
                            else
                                attached_UEs(u_).eNodeB_signaling.N_used_bits(cw_) = 0;
                            end
                            attached_UEs(u_).eNodeB_signaling.packet_parts{cw_} = packet_parts;
                        else
                            attached_UEs(u_).traffic_model.decrease_packets(TB_size_bits(cw_));
                            attached_UEs(u_).eNodeB_signaling.N_used_bits(cw_) = min(attached_UEs(u_).traffic_model.get_buffer_length,TB_size_bits(cw_));
                        end
                    else
                        attached_UEs(u_).eNodeB_signaling.N_used_bits(cw_) =  0;
                        attached_UEs(u_).eNodeB_signaling.packet_parts{cw_} = [];
                    end
                end
            end
            
            % Total size in bits of the RB grid
            the_RB_grid.size_bits = RB_grid_size_bits;
            
            % TODO: HARQ handling, #streams decision and tx_mode decision.
        end
        
        function TP = compute_av_throughput(obj,u_,UE_feedback,TTI_to_read)
            UE_id = UE_feedback.UE_id(u_);
            if UE_id
                TP = sum(obj.UE_traces(UE_id).avg_throughput(:,TTI_to_read))*10^-3; % Mean throughput, averaged with an exponential window
            else
                TP = 0;
            end
        end
               
        function [c,user_ind] = get_efficiency(obj,N_UE,N_RB,UE_feedback)
            c                   = zeros(N_UE,N_RB);
            nCodewords_feedback = size(UE_feedback.CQI,2);
            
            CQIs_per_UE = reshape(UE_feedback.CQI,[],N_UE);
            max_UE_CQI  = max(CQIs_per_UE,[],1);
            CQI_bars    = max_UE_CQI+1;
            
            % Assumes that the CQIs in both codewords will be similar
            a = obj.k(CQI_bars)';
            b = obj.d(CQI_bars)';
            
            RI          = UE_feedback.RI;
            RI(~UE_feedback.feedback_received)         = 1;

            % Begin permuted part
            user_ind      = randperm(N_UE);
            a_vector      = a(user_ind);
            b_vector      = b(user_ind);
            a_matrix      = a_vector(:,ones(1,nCodewords_feedback));
            b_matrix      = b_vector(:,ones(1,nCodewords_feedback));
            permuted_CQIs = UE_feedback.CQI(:,:,user_ind);
            for u_=1:N_UE
                UE_idx = user_ind(u_);
                UE_RI  = RI(UE_idx);
                
                % Assuming always a maximum possible number of codewords
                % i.e., min(2,nLayers)
                switch UE_RI
                    case {0,1} % SISO sends RI of 0
                        a_matrix(u_,2) = 0;
                        b_matrix(u_,2) = 0;
                    case 2
                        % Do nothing
                    case 3
                        a_matrix(u_,2) = 2*a_matrix(u_,2);
                        b_matrix(u_,2) = 2*b_matrix(u_,2);
                    case 4
                        a_matrix(u_,:) = 2*a_matrix(u_,:);
                        b_matrix(u_,:) = 2*b_matrix(u_,:);
                end
            end
            for rb = 1:N_RB
                % c(:,rb) = a_vector.*reshape(permuted_CQIs(rb,1,:),[],1)+b_vector;
                c(:,rb) = sum(a_matrix .* reshape(reshape(permuted_CQIs(rb,:,:),[],1),2,[])' + b_matrix,2);
            end
        end
        
        % This wrapping functions are needed to be able to support FFR schedulers
        function set_SINR_averager(obj,SINR_averager)
            obj.SINR_averager = SINR_averager;
        end
        function set_CQI_mapper(obj,CQI_mapper)
            obj.CQI_mapper = CQI_mapper;
        end
        function set_BLER_curves(obj,BLER_curves)
            obj.BLER_curves = BLER_curves;
        end
        function set_genie_UEs(obj,UEs)
            obj.genie.UEs = UEs;
        end
        function set_genie_eNodeBs(obj,eNodeBs)
            obj.genie.eNodeBs = eNodeBs;
        end
        function set_feedback_delay_TTIs(obj,feedback_channel_delay)
            obj.feedback_delay_TTIs = feedback_channel_delay;
        end
        function set_UE_traces(obj,the_UE_traces)
            obj.UE_traces = the_UE_traces;
        end
    end
    
end

