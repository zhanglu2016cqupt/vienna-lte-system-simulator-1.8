classdef UE < handle
    % Class that represents an LTE UE (user)
    % (c) Josep Colom Ikuno, INTHFT, 2008
    
    properties
        
        id                    % Unique UE id
        pos                   % pos in meters (x,y)
        attached_site         % Site to where this user is attached
        attached_sector_idx   % sector index to which the UE is attached
        attached_eNodeB       % eNodeB to which the UE is attached to
        
        walking_model         % Walking model for this user
        downlink_channel      % Downlink channel model for this user
        RB_grid               % This links to obj.downlink_channel.RB_grid Stored again here for efficiency reasons
        
        uplink_channel        % Uplink channel model for this user
        receiver_noise_figure % Noise figure for this specific UE (dB)
        thermal_noise_W_RB    % Calculated based on the thermal noise density and the receiver noise figure in Watts/RB
        penetration_loss      % Penetration loss in dB for this user
        wall_loss = 0;        % wall loss that only affects interferers [dB]
        nRX                   % Number of receive antennas
        antenna_gain          % Antenna gain of the UE
        
        trace                 % Trace that stores info about what happened
        clock                 % Network clock. Tells the UE in what TTI he is
        CQI_mapper            % Performs the mapping between SINR and CQI
        
        % Output of the link quality (measurement) model.
        link_quality_model_output
        
        % Data to be fedbacked to the eNodeB. It is used to pass the feedback data to the send_feedback() function
        feedback
        
        % Whether the CQI feedback should be unquantized. Having this set
        % to true is equivalent to directly sending the post-equalization
        % SINR for each codeword (note that there is still a layer mapping)
        unquantized_CQI_feedback
        
        % Will decide whether a give TB made it or not
        BLER_curves
        
        % Gives the means to average the several Transport Block (TB) SINRs
        SINR_averager
        
        % Contains the LTE precoding codebook
        codebook
        
        % Signaling from the eNodeB to this UE. This is a direct channel
        % between the eNodeB and this UE, where it gets signaled
        % UE-specific signaling information. The signaled information and
        % where it is located is as follows:
        %   UEsignaling:
        %     - TB_CQI           % CQI used for the transmission of each codeword
        %     - TB_size          % size of the current TB, in bits
        %     - tx_mode          % transmission mode used (SISO, tx diversity, spatial multiplexing)
        %     - rv_idx           % redundancy version index for each codeword
        %   downlink_channel.RB_grid
        %     - user_allocation  % what UE every RB belongs to
        %     - power_allocation % how much power to allocate to each RB,
        %     - n_RB             % RB grid size (frequency)
        %     - sym_per_RB       % number of symbols per RB (12 subcarriers, 0.5ms)
        %     - size_bits        % total size of the RB grid in bits
        %     - numStreams       % maximum number of allowed streams. Resource allocation is described for all of them
        eNodeB_signaling
        
        % Extra tracing options (default options)
        trace_SINR
        
        % average preequalization SNR at current position (averaged over microscopic fading and noise)
        SNR_avg_preequal
        
        % This is an "overall SINR", calculated by summing up all of the
        % signal power and dividing it by the sum of all interfering and
        % noise power.
        wideband_SINR
        
        % Lets you use instead of the link quality and performance models
        % dummy funcitons that output dummy values. Useful if you would like
        % to "deactivate" some UEs to shorten simulation time. Note that the
        % feedback will also be deactivated.
        deactivate_UE
        
        % This variable is used for the default feedback calculation sent in
        % case the UE was not scheduled. For the old (v1) trace format, this
        % sets the actual transmit mode.
        default_tx_mode
        
        traffic_model
        lambda = 0;
        
        adaptive_RI
        
        % Helper variables for the very simple handover management. Just
        % meant as an example as to how that could be implemented.
        cell_change
        
        % Safe temporarily for tracing
        rx_power_tb_in_current_tti
        rx_power_interferers_in_current_tti
        
        % If the attached cell and pathloss (RSRP) come from a trace
        trace_UE = false;
        
        % If the simulation is set to employ runtime precoding (true) or
        % precalculated precoders (false)
        runtime_precoding = false;
        
        Zero_FB_delay_SINR_store  % auxiliary variable to store the calculated SINRs in case of 0 delay feedback
    
    end
    
    methods
        % Constructor with the default UE parameter values
        function obj = UE
            obj.unquantized_CQI_feedback  = false;
            obj.trace_SINR                = false;
            obj.deactivate_UE             = false;
            obj.cell_change.requested     = false;
            obj.cell_change.target_eNodeB = [];
        end
        
        function print(obj)
            if isempty(obj.attached_site)
                fprintf('User %d, (%d,%d), not attached to an eNodeB\n',obj.id,obj.pos(1),obj.pos(2));
            else
                fprintf('User %d, (%d,%d), Site %d, sector %d (eNodeB %d)\n',obj.id,obj.pos(1),obj.pos(2),obj.attached_site.id,obj.attached_sector_idx,obj.attached_eNodeB.eNodeB_id);
            end
            obj.walking_model.print;
        end
        
        % Clear variables
        function clear(obj)
            obj.attached_site             = [];
            obj.attached_eNodeB           = [];
            obj.walking_model             = [];
            obj.downlink_channel          = [];
            obj.RB_grid                   = [];
            obj.uplink_channel            = [];
            obj.trace                     = [];
            obj.clock                     = [];
            obj.CQI_mapper                = [];
            obj.link_quality_model_output = [];
            obj.feedback                  = [];
            obj.BLER_curves               = [];
            obj.SINR_averager             = [];
            obj.eNodeB_signaling          = [];
            obj.traffic_model             = [];
            obj.codebook                  = [];
        end
        
        % Move this user according to its settings
        function move(obj)
            new_pos = obj.walking_model.move(obj.pos);
            obj.pos = new_pos;
        end
        % Move this user to where it was the last TTI before according to
        % its settings
        function move_back(obj)
            old_pos = obj.walking_model.move(obj.pos);
            obj.pos = old_pos;
        end
        function UE_in_roi = is_in_roi(a_UE,roi_x_range,roi_y_range)
            % Tells you whether a user in in the Region of Interest (ROI) or not
            % (c) Josep Colom Ikuno, INTHFT, 2008
            % input:    a_UE         ... the UE in question
            %           roi_x_range  ... roi x range. minimum and maximum x coordinates
            %                            which are valid
            %           roi_y_range  ... roi y range. minimum and maximum y coordinates
            %                            which are valid
            % output:   UE_in_roi  ... true or false, whether the UE is inside or not
            
            UE_pos_x = a_UE.pos(1);
            UE_pos_y = a_UE.pos(2);
            
            if UE_pos_x<roi_x_range(1) || UE_pos_x>roi_x_range(2)
                UE_in_roi = false;
                return;
            end
            
            if UE_pos_y<roi_y_range(1) || UE_pos_y>roi_y_range(2)
                UE_in_roi = false;
                return;
            end
            UE_in_roi = true;
        end
        % Starts handover procedures from the currently attached eNodeB to
        % the specified target_eNodeB
        % for now... immediate handover. A proper implementation remains
        % pending.
        function start_handover(obj,new_eNodeB)
            % Remove the user from the eNodeB and its scheduler
            obj.attached_eNodeB.deattachUser(obj);
            
            % Add the user to the eNodeB and its scheduler
            new_eNodeB.attachUser(obj);
            
            % Set a new channel realization
            pregenerated_ff = obj.downlink_channel.fast_fading_model.ff_trace;
            N_eNodeBs       = length(obj.downlink_channel.fast_fading_model.interfering_starting_points);
            obj.downlink_channel.set_fast_fading_model_model(channel_gain_wrappers.fastFadingWrapper(pregenerated_ff,'random',N_eNodeBs));
        end
        
        % Measure whatever needs to be measured and send a feedback to the attached eNodeB
        function send_feedback(obj)
            obj.uplink_channel.send_feedback(obj.feedback);
        end
        
        function [...
                interfering_eNodeBs,...
                user_macroscopic_pathloss_dB,...
                user_shadow_fading_loss_dB,...
                there_are_interferers...
                ] = get_signal_macroscale_losses(obj)
            % Read macro scale pathloss values for the signal part
            if ~obj.trace_UE
                interfering_eNodeBs          = obj.attached_eNodeB.neighbors_eNodeB;
                user_macroscopic_pathloss_dB = obj.downlink_channel.macroscopic_pathloss + obj.penetration_loss - obj.antenna_gain;
                user_shadow_fading_loss_dB   = obj.downlink_channel.shadow_fading_pathloss;
                there_are_interferers        = ~isempty(interfering_eNodeBs);
            else
                UE_id       = obj.id;
                UE_trace    = obj.downlink_channel.macroscopic_pathloss_model.pathloss(UE_id);
                current_TTI = obj.clock.current_TTI;
                current_trace_time = floor((current_TTI-1)/UE_trace.TTIs_per_time_idx)+1;
                
                current_pathlosses = UE_trace.pathloss(current_trace_time,:);
                user_macroscopic_pathloss_dB = min(current_pathlosses);
                user_shadow_fading_loss_dB   = 0;
                there_are_interferers        = length(current_pathlosses)>1;
                
                attached_cell = UE_trace.attached_cell(current_trace_time);
                all_cells     = UE_trace.cellsIds(current_trace_time,:);
                if there_are_interferers
                    % is max_ues is a property, MU-MIMO is being done, in
                    % that case, the home eNodeB is also a potential
                    % interferer
                    if ~isprop(RB_grid, 'max_ues')
                        interfering_eNodeBs = obj.downlink_channel.eNodeBs(all_cells(all_cells~=attached_cell));
                    else
                        interfering_eNodeBs = obj.downlink_channel.eNodeBs(all_cells);
                    end
                else
                    interfering_eNodeBs = [];
                end
            end
        end
        
        function [interfering_macroscopic_pathloss_eNodeB_dB, interfering_shadow_fading_loss_dB]...
                = get_interfering_macroscale_losses(obj,interfering_eNodeB_ids,parent_sites_id)
            % Read macro scale pathloss values for the interference part
            if ~obj.trace_UE
                interfering_macroscopic_pathloss_eNodeB_dB = obj.downlink_channel.interfering_macroscopic_pathloss(interfering_eNodeB_ids) + obj.penetration_loss - obj.antenna_gain;
                interfering_shadow_fading_loss_dB          = obj.downlink_channel.interfering_shadow_fading_pathloss(parent_sites_id);
            else
                UE_id       = obj.id;
                UE_trace    = obj.downlink_channel.macroscopic_pathloss_model.pathloss(UE_id);
                current_TTI = obj.clock.current_TTI;
                current_trace_time = floor((current_TTI-1)/UE_trace.TTIs_per_time_idx)+1;
                
                % The ordering comes from the traces, so we do not need to
                % do a one-by-one mapping
                current_cells      = UE_trace.cellsIds(current_trace_time,:);
                current_pathlosses = UE_trace.pathloss(current_trace_time,:);
                interfering_macroscopic_pathloss_eNodeB_dB = current_pathlosses(current_cells~=obj.attached_eNodeB.eNodeB_id)'; % Must be a column vector
                interfering_shadow_fading_loss_dB          = 0;
            end
        end
        
        % Redirects the link quality function to different functions
        % depending on the configuration
        function link_quality_model(obj,config,before_FB)
            if obj.runtime_precoding
                obj.link_quality_model_v2(config,before_FB);
            else
                if ~obj.deactivate_UE 
                    obj.link_quality_model_v1(config);
                else
                    obj.dummy_link_quality_model(config);
                end
            end
        end
        
        % Calculates the receiver SINR (version with precalculated precoders)
        function link_quality_model_v1(obj,config)
            
            % Get current time
            t = obj.clock.time;

            % Get map-dependant parameters for the current user
            [...
                interfering_eNodeBs,...
                user_macroscopic_pathloss_lin,...
                user_shadow_fading_loss_lin,...
                there_are_interferers...
                ] = obj.get_signal_macroscale_losses;
            
            % Number of codewords, layers, power etc. assigned to this user
            DL_signaling = obj.eNodeB_signaling;
            tx_mode      = obj.default_tx_mode;         % Fixed tx mode according to LTE_config
            the_RB_grid  = obj.downlink_channel.RB_grid;
            nRB          = the_RB_grid.n_RB;
            nSC          = nRB*2;
            
            % Get the RX power (power allocation) from the target eNodeB
            TX_power_data      = the_RB_grid.power_allocation';
            TX_power_signaling = the_RB_grid.power_allocation_signaling';
            RX_total_RB        = (TX_power_data+TX_power_signaling)./user_macroscopic_pathloss_lin./user_shadow_fading_loss_lin;
            RX_total           = reshape([RX_total_RB; RX_total_RB],1,[])/(2);            
            
            % Get fast fading trace for this subframe
            [zeta,chi,psi] = obj.downlink_channel.fast_fading_model.generate_fast_fading_signal(t,tx_mode);
            
            % Emulate ICI power according to user velocity
            ICI_power = config.ICI_power.*repmat(RX_total,size(config.ICI_power,1),1);
            %% The SINR calculation is done under the following circumstances:
            % Power allocation is done on a per-subframe (1 ms) and RB basis
            % The fast fading trace is given for every 6 subcarriers (every
            % 90 KHz), so as to provide enough samples related to a
            % worst-case-scenario channel length
            
            % TX_power_signaling_half_RB =  TODO: add signaling interference in better-modeled way
            S_dims = size(zeta);
            S_dims(3) = 1; % All MATLAB variables have at least 2 dimensions, so not a problem.            
            % RX power
            switch tx_mode
                case 1     % SISO
                    RX_temp = repmat(RX_total,[S_dims(1),1]);
                    RX_power = RX_temp.*zeta;
                otherwise  % TxD, OLSM or CLSM
                    RX_total = shiftdim(RX_total,-1);
                    RX_power_half_RB_repmat = repmat(RX_total,S_dims);
                    RX_power = RX_power_half_RB_repmat.*zeta;
            end

            % Get interfering eNodeBs. Take only those with a power higher than 45 with respect to our signal
            if there_are_interferers % no interfering eNodeBs present (single eNodeB simulation)
                
                [...
                    interfering_eNodeB_ids,...
                    interfering_eNodeBs,...
                    interfering_macroscopic_pathloss_eNodeB_lin,...
                    interfering_power_allocations,...
                    interfering_shadow_fading_loss_lin...
                    ] = filter_out_interferers_below_SINR(obj,sum(RX_total_RB),interfering_eNodeBs,45);
                
                % Get interfering channel fading parameters
                theta = obj.downlink_channel.fast_fading_model.generate_fast_fading_interference(t,tx_mode,interfering_eNodeB_ids);
                SINR_interf_dims = size(theta);
                % Get assigned interfering power on a per-half-RB-basis
                if config.feedback_channel_delay~=0
                    interfering_power_allocations_temp = interfering_power_allocations/2;
                    interf_power_all_RB = reshape([interfering_power_allocations_temp(:) interfering_power_allocations_temp(:)]',2*size(interfering_power_allocations_temp,1),[]); % Take scheduled power
                else
                    if ndims(SINR_interf_dims)==3 % #ok<ISMAT>
                        TX_power_interferers = [interfering_eNodeBs.max_power]/SINR_interf_dims(3);
                        interf_power_all_RB  = TX_power_interferers(ones(1,SINR_interf_dims(3)),:);
                    else
                        interf_power_all_RB = repmat([interfering_eNodeBs.max_power]/SINR_interf_dims(3),[SINR_interf_dims(3) 1]); % Turn on all interferers
                    end
                end
                
                temp_macro_mat      = interfering_macroscopic_pathloss_eNodeB_lin';
                temp_macro_mat      = temp_macro_mat(ones(SINR_interf_dims(3),1),:);
                temp_shadow_mat     = interfering_shadow_fading_loss_lin';
                temp_shadow_mat     = temp_shadow_mat(ones(SINR_interf_dims(3),1),:);
                interf_power_all_RB = interf_power_all_RB./temp_macro_mat./temp_shadow_mat; % Add macro and shadow fading
                
                % Temporarily safe interference power on per RB block basis; For tracing and statistical evaluation
                % of interference.
                %                 if tx_mode==1 && true
                %                     if (obj.clock.current_TTI==1)&&(~obj.deactivate_UE)
                %                         % Assume unit power, i.e. transmit power per RB = 1
                %                         % Take only taps from first RB as representative - Reason: Shadow fading is constant for all RBs.
                %                         microscale_fading_taps_temp  = theta';
                %                         shadow_fading_taps_temp      = temp_shadow_mat;
                %                         composite_fading_taps_temp   = microscale_fading_taps_temp.*shadow_fading_taps_temp;
                %                         aggregated_interference_temp = sum(theta./temp_macro_mat'./temp_shadow_mat',1);
                %                         if (exist('Interference Statistics.mat','file'))
                %                             load('Interference Statistics.mat', 'microscale_fading_taps', 'shadow_fading_taps', 'composite_fading_taps', 'aggregated_interference_taps');
                %                             microscale_fading_taps       = [microscale_fading_taps;       microscale_fading_taps_temp(:)];
                %                             shadow_fading_taps           = [shadow_fading_taps;           shadow_fading_taps_temp(:)];
                %                             composite_fading_taps        = [composite_fading_taps;        composite_fading_taps_temp(:)];
                %                             aggregated_interference_taps = [aggregated_interference_taps; aggregated_interference_temp(:)];
                %                             save('Interference Statistics.mat', 'microscale_fading_taps', 'shadow_fading_taps','composite_fading_taps','aggregated_interference_taps');
                %                         % In case no file exists, create a new one
                %                         else
                %                             microscale_fading_taps       = microscale_fading_taps_temp(:);
                %                             shadow_fading_taps           = shadow_fading_taps_temp(:);
                %                             composite_fading_taps        = composite_fading_taps_temp(:);
                %                             aggregated_interference_taps = aggregated_interference_temp(:);
                %                             save('Interference Statistics.mat', 'microscale_fading_taps', 'shadow_fading_taps','composite_fading_taps','aggregated_interference_taps');
                %                         end
                %                     end
                %                 end
                obj.rx_power_tb_in_current_tti = mean(RX_total_RB,2); % linear scale !
                
                % To avoid errors. This trace is thought for SISO
                obj.rx_power_interferers_in_current_tti = zeros(2,size(interf_power_all_RB,2));
                if length(size(interf_power_all_RB))==2
                    obj.rx_power_interferers_in_current_tti(1,:) = 2*mean(interf_power_all_RB(:,:),1); % linear scale !
                    % Add eNodeB ID of interferer as second line to rx power (in order to identify interferer tiers)
                    obj.rx_power_interferers_in_current_tti(2,:) = interfering_eNodeB_ids;
                else
                    obj.rx_power_interferers_in_current_tti(1,:) = NaN;
                    obj.rx_power_interferers_in_current_tti(2,:) = NaN;
                end
                
                max_Layers = SINR_interf_dims(1);
                if length(SINR_interf_dims) > 4
                    N_RI = SINR_interf_dims(5);
                else
                    N_RI = 1;
                end
                
%                 switch tx_mode
%                     case 1
%                         interf_power_all_RB_repmat = interf_power_all_RB.';
%                     otherwise
                        interf_power_all_RB_repmat = zeros(SINR_interf_dims);
                        for nLayers = 1:max_Layers
                            for RI = 1:N_RI
                                interf_power_all_RB_repmat(nLayers,:,:,:,RI) = repmat(shiftdim(interf_power_all_RB,-2),1,SINR_interf_dims(2));
                            end
                        end
                        
                        % This line is totally incorrect!!!
                        % interf_power_all_RB_repmat = reshape(repmat(interf_power_all_RB,SINR_interf_dims_repmat),SINR_interf_dims); % Also valid for the case where more than one rank is used
%                 end
            else
                obj.rx_power_interferers_in_current_tti = 0;
            end
            
            % Calculate thermal noise
            thermal_noise_watts_per_half_RB = obj.thermal_noise_W_RB/2;
            % Calculate average preequalization SNR
            % This is a total SNR, the same as in the Link Level Simulator
            obj.SNR_avg_preequal = 10*log10(mean(RX_total)./thermal_noise_watts_per_half_RB); % mean over the subcarriers

            switch tx_mode
                case 1
                    % SINR calculation (SISO)
                    noise_plus_inter_layer_power = psi.*thermal_noise_watts_per_half_RB + psi.*ICI_power;
                    if there_are_interferers
                        % Also works for more than one rank (i.e. extra dimension)
                        size_vec = size(interf_power_all_RB_repmat);
                        interfering_rx_power = reshape(sum(interf_power_all_RB_repmat.*theta,4),[size_vec(2:3),size_vec(5:end)]);
%                         interfering_rx_power = squeeze(sum(interf_power_all_RB_repmat.*theta,4));
                        Interference_plus_noise_power = noise_plus_inter_layer_power + interfering_rx_power;
                    else
                        Interference_plus_noise_power = noise_plus_inter_layer_power;
                    end
                    SINR_linear = RX_power ./ Interference_plus_noise_power; % Divide thermal noise by 2: Half-RB frequency bins
                    
                    % Calculate SIR
                    % if there_are_interferers
                    %     SIR_linear = RX_power ./ interfering_rx_power.';
                    % else
                    %     SIR_linear = Inf(size(SINR_linear));
                    % end
                otherwise
                    % SINR calculation (TxD, OLSM, CLSM)
                    psi_size = size(psi);
                    psi_size = [psi_size,1]; % just to make sure that psi_size has four elements (not per se with one layer)
                    ICI_power = repmat(shiftdim(ICI_power,-1),[psi_size(1),1,1,psi_size(4)]);
%                     layer_norm = ones(size(ICI_power));
%                     for li = 1:size(layer_norm,4)
%                     size(ICI_power)
%                     error('j')
                    noise_plus_inter_layer_power = chi.*RX_power + psi.*thermal_noise_watts_per_half_RB + psi.*ICI_power; % Divide thermal noise by 2: Half-RB frequency bins

                    if there_are_interferers
                        % Also works for more than one rank (i.e. extra dimension)
%                         interfering_rx_power = squeeze(sum(interf_power_all_RB_repmat.*theta,4));
                        size_vec = size(interf_power_all_RB_repmat);
                        interfering_rx_power = reshape(sum(interf_power_all_RB_repmat.*theta,4),[size_vec(1:3),size_vec(5:end)]);
                        Interference_plus_noise_power = noise_plus_inter_layer_power + interfering_rx_power;
                    else
                        Interference_plus_noise_power = noise_plus_inter_layer_power;
                    end
                    SINR_linear = RX_power ./ Interference_plus_noise_power;
                    
                    % Calculate SIR
                    % if there_are_interferers
                    %     SIR_linear = RX_power ./ interfering_rx_power;
                    % else
                    %     SIR_linear = Inf(size(SINR_linear));
                    % end
            end

            % Calculation of the wideband SINR
            if there_are_interferers
                obj.wideband_SINR = 10*log10(sum(RX_total(:))/(sum(interf_power_all_RB(:))+thermal_noise_watts_per_half_RB*nSC));
            else
                obj.wideband_SINR = 10*log10(sum(RX_total(:))/(thermal_noise_watts_per_half_RB*nSC));
            end
            
            % Calculation of the post-equalization symbols SINR
            SINR_dB = 10*log10(SINR_linear);
            % SIR_dB  = 10*log10(SIR_linear);
            
            if tx_mode == 1
                SINR_dB = shiftdim(SINR_dB,-1);
            end
            
            % Store post-equalization SINR for the link performance model
            if config.feedback_channel_delay ~= 0 % in case of 0 delay feedback, the SINRs have to be stored after scheduling
                if DL_signaling.num_assigned_RBs>0
                    assigned_SCs = reshape([DL_signaling.assigned_RB_map DL_signaling.assigned_RB_map]',1,[]);
                    obj.store_SINR_for_link_performance_model(tx_mode,...
                        SINR_dB(1:DL_signaling.nLayers,:,assigned_SCs,DL_signaling.nLayers));
                else
                    obj.store_SINR_for_link_performance_model(tx_mode, -Inf);
                    obj.wideband_SINR       = -Inf;
                    obj.SNR_avg_preequal    = -Inf;
                end
            else
                obj.Zero_FB_delay_SINR_store      = SINR_dB;
            end          
                        
            % Calculate and save feedback, as well as the measured SINRs
            obj.feedback = feedback_calculation.CQI_and_rank(config,tx_mode,SINR_dB,nRB,[],DL_signaling,...
                obj.SINR_averager,obj.CQI_mapper,obj.unquantized_CQI_feedback,obj.adaptive_RI);

        end
        
        function set_SINRs_according_to_schedule(obj)
            DL_signaling = obj.eNodeB_signaling;
            tx_mode      = obj.default_tx_mode;         
            SINR_dB      = obj.Zero_FB_delay_SINR_store;
            if DL_signaling.num_assigned_RBs>0 && ~obj.deactivate_UE
                assigned_SCs = reshape([DL_signaling.assigned_RB_map DL_signaling.assigned_RB_map]',1,[]);
                obj.store_SINR_for_link_performance_model(tx_mode,...
                    SINR_dB(1:DL_signaling.nLayers,:,assigned_SCs,DL_signaling.nLayers));
            else
                obj.store_SINR_for_link_performance_model(tx_mode, -Inf);
                obj.wideband_SINR       = -Inf;
                obj.SNR_avg_preequal    = -Inf;
            end
        end
        
        
        % Next version of the Link Quality model, implemented with run-time
        % precoding. Only CLSM and similar precoded modes implemented
        function link_quality_model_v2(obj,config,before_FB)
            
            % To support zero delay feedback in case of runtime precoding
            FB_delay_zero_run = false;
            if before_FB && config.feedback_channel_delay == 0
                FB_delay_zero_run = true;
            end
            
            % Get current time
            t = obj.clock.time-1e-3;
            
            % Get map-dependant parameters for the current user
            [...
                interfering_eNodeBs,...
                pathloss_linear_0,...
                shadow_fading_linear_0,...
                there_are_interferers...
                ] = obj.get_signal_macroscale_losses;
            
            % Number of codewords, layers, power etc. assigned to this user
            DL_signaling = obj.eNodeB_signaling;
            tx_mode      = obj.default_tx_mode;         % Fixed tx mode according to LTE_config
            RB_grid_0    = obj.downlink_channel.RB_grid;
            nRB          = RB_grid_0.n_RB;
            UE_RBs_mat   = DL_signaling.assigned_RB_map; 
            if isempty(UE_RBs_mat)
                UE_RBs_mat = zeros(size(RB_grid_0.user_allocation));
            end
%             if config.feedback_channel_delay == 0 % assume full RB allocation in case of zero-delay feedback
%                 UE_RBs_mat = ones(size(RB_grid_0.user_allocation));
%             end
            % When not doing MU-MIMO, this is a logical index vector
            % When using MU-MIMO, this is a logical matrix given by (nRBxnTX). As it stands now, 
            % a UE is only allowed to be scheduled once per RB, i.e. there
            % is only one 'true' entry per row. Therefore, a index vector
            % can always be calculated like this:
            UE_RBs       = logical(sum(UE_RBs_mat, 2));
            
            % Get the RX power (power allocation) from the target eNodeB
            TX_W_data_0      = RB_grid_0.power_allocation';
            TX_W_signaling_0 = RB_grid_0.power_allocation_signaling';
            TX_W_RB_0        = (TX_W_data_0+TX_W_signaling_0);
            noise_W_half_RB  = obj.thermal_noise_W_RB/2;            
            H_0 = obj.downlink_channel.fast_fading_model.channel_matrix_signal(t)./ sqrt(pathloss_linear_0(1).*shadow_fading_linear_0(1));
            % Emulate ICI power according to user velocity
            TX_W_RB_0_half_RB = repmat(mean([TX_W_RB_0,TX_W_RB_0]/2,1),size(config.ICI_power,1),1);
            ICI_power = config.ICI_power./(pathloss_linear_0(1).*shadow_fading_linear_0(1))*size(H_0,2).*TX_W_RB_0_half_RB;

            % Add RRH H_0 channel matrix entries
            if ~isempty(obj.attached_eNodeB.RRHs)
                RRHs      = [obj.attached_eNodeB.RRHs];
                extra_nTX = [RRHs.nTX];
                %                 nTX_idxs  = size(H_0,2)+cumsum(extra_nTX);
                nTX_idxs = size(H_0,2) + (1:sum(extra_nTX));
%                 size(H_0,2)
%                 sum(extra_nTX)
%                 nTX_idxs
%                 extra_nTX
%                 size(H_0,2) + 1:sum(extra_nTX)
%                 error('j')
                % Preallocate complete channel matrix
                H_0(:,(end+1):(end+sum(extra_nTX)),:,:) = 0;
                % Fill in the extra positions in H_0
                for rrh_=1:length(RRHs)
                    H_0(:,nTX_idxs(rrh_):(nTX_idxs(rrh_)+extra_nTX(rrh_)-1),:,:) = ...                        
                        obj.downlink_channel.secondary_fast_fading_models(rrh_).channel_matrix_signal(t) ./ ...
                        sqrt(pathloss_linear_0(rrh_+1).*shadow_fading_linear_0(rrh_+1));
                    ICI_power = ICI_power + config.ICI_power./(pathloss_linear_0(rrh_+1).*shadow_fading_linear_0(rrh_+1))*extra_nTX(rrh_).*TX_W_RB_0_half_RB;
                end              
            end            
            % Get the correct precodig matrices
            nTX         = size(H_0,2);
            UE_half_RBs = logical(kron(UE_RBs, [1; 1]));
            
            % ICI power normalization (equal power allocation across
            % antennas assumed
            ICI_power = ICI_power/nTX;
            
            % This part of the calcualtion of the interfering power set
            % aside because it is reused for the feedback calculations so
            % that we don't have to repeat the same calls
            avg_pathloss_linear_0      = obj.attached_eNodeB.average_pathloss_value(pathloss_linear_0);
            avg_shadow_fading_linear_0 = obj.attached_eNodeB.average_pathloss_value(shadow_fading_linear_0);
            sum_macroscopic_RX_W       = sum(sum(TX_W_RB_0))/avg_pathloss_linear_0/avg_shadow_fading_linear_0;
            

            if there_are_interferers % no interfering eNodeBs present (single eNodeB simulation)
                % Reduce complexity: 
                [...
                    interfering_eNodeB_ids,...
                    interfering_eNodeBs,...
                    pathloss_linear_i,...
                    TX_W_RB_i,...
                    shadow_fading_linear_i,...
                    avg_pathloss_lin_i,...
                    avg_shadow_lin_i...
                    ] = obj.filter_out_interferers_below_SINR_v2(sum_macroscopic_RX_W,... 
                                                                UE_RBs_mat,...
                                                            interfering_eNodeBs,...
                                                            45);
                RB_grids_i = [interfering_eNodeBs.RB_grid];
                % Channel matrix coming from the main TX
                H_i_full = obj.downlink_channel.fast_fading_model.channel_matrix_interferers(t,interfering_eNodeB_ids);
                for i_=1:length(interfering_eNodeB_ids)
                    H_i_full(:,:,:,:,i_) = H_i_full(:,:,:,:,i_) ./ sqrt(pathloss_linear_i(i_,1).*shadow_fading_linear_i(i_,1));
                end
                
                % Channel matrix of the Remote Radio Heads (RRHs)
                % IMPORTANT: ASSUMES THE SAME NUMBER OF RRHs IN THE
                % INTERFERING ENODEBS AS IN THIS ENODEB!!!!
                
                % Add RRH H_i_full channel matrix entries
                if ~isempty(obj.attached_eNodeB.RRHs)
                    RRHs      = vertcat(interfering_eNodeBs.RRHs);
                    extra_nTX = [RRHs(1,:).nTX];
                    %                     nTX_idxs  = size(H_i_full,2)+cumsum(extra_nTX);
                    nTX_idxs  = size(H_i_full,2) + (1:sum(extra_nTX));
                    
                    % Preallocate complete channel matrix
                    H_i_full(:,(end+1):(end+sum(extra_nTX)),:,:,:) = 0;
                    
                    % Fill in the extra positions in H_i_full
                    for rrh_=1:size(RRHs,2)
                        H_i_RRH = obj.downlink_channel.secondary_fast_fading_models(rrh_).channel_matrix_interferers(t,interfering_eNodeB_ids);
                        for i_=1:length(interfering_eNodeB_ids)
                            H_i_RRH(:,:,:,:,i_) = H_i_RRH(:,:,:,:,i_) ./ sqrt(pathloss_linear_i(i_,rrh_+1).*shadow_fading_linear_i(i_,rrh_+1));
                        end
                        H_i_full(:,nTX_idxs(rrh_):(nTX_idxs(rrh_)+extra_nTX(rrh_)-1),:,:,:) = H_i_RRH;
                    end
                end
                
                % Interfering precoders
                W_i_temp = horzcat(RB_grids_i.precoder);
                W_i_full = repmat(W_i_temp, 1, 2);
                W_i_full(:, 1:2:end) = W_i_temp;
                W_i_full(:, 2:2:end) = W_i_temp;
                %W_i_full = reshape([W_i_full;W_i_full],[],size(H_i_full,4));
                
                % Reshape Interfering precoders into a per half_RB manner:
                % However, W_i is a matrix for each interferer. It now has
                % size: nTX x (nInterferers*nRBs)
                % Step one: Get a singleton dimension and transpose :
                % size: 1 x (nInt*nRBs) x nTX
                W_i_full = reshape(W_i_full.', 2*nRB, [], size(W_i_full, 1));
                
                %Step two: Double first dimension (for half_RBs) and
                % reshape first and second dimension in such a way that the
                % second dimension has size nInt --> first dimension has
                % size 2*nRB, third dimension stays the same
                % size(W_i_full) = 2*nRB x nInt x nTX
                %W_i_full = reshape([W_i_full;W_i_full],[],size(H_i_full,4), size(W_i_full, 3));
                
                % We are only interested on interference on the assigned RBs
                H_i          = H_i_full(:,:,UE_half_RBs,:,:);
                W_i          = W_i_full(UE_half_RBs,:, :);
                TX_W_RB_i_TB = TX_W_RB_i(UE_RBs,:, :);
            else
                TX_W_RB_i    = [];
                H_i_full     = [];
                W_i_full     = [];
                TX_W_RB_i_TB = [];
                H_i          = [];
                W_i          = [];
            end

            % Calculate receiver filter and effective channel matrix
             if sum(DL_signaling.num_assigned_RBs)>0 && ~FB_delay_zero_run
                % Assumes that all of the precoders for this UE are of the same size (i.e., same number of layers)
                % Get linear(!) indices of assigned RBs in RB matrix
                % (reshape is only to force column vector)
                
%                 if config.feedback_channel_delay == 0
%                     assigned_RBs_temp = nRB;
%                 else
                    assigned_RBs_temp = DL_signaling.num_assigned_RBs;
%                 end                
                
                RB_indexes = reshape(find(UE_RBs_mat'),[], 1);
                RB_indexes = kron(RB_indexes, [1; 1]);  % double entries for half RB
                
                % Reshape into size(precoding_vector) x assigned RBs, where
                % size(precoding_vector) itelf can be a matrix
                W_0        = reshape([RB_grid_0.precoder(RB_indexes).W],...
                    [size(RB_grid_0.precoder(RB_indexes(1)).W) sum(assigned_RBs_temp)*2]);
                % TB SINR
                [TB_SINR_lin, RX_total_half_RB_layers_UE, RX_W_half_RB_i_UE] = phy_modeling.post_equalization_SINR(...
                    TX_W_RB_0(UE_RBs_mat')',...
                    H_0(:,:,UE_half_RBs,:), W_0,...
                    noise_W_half_RB,...
                    TX_W_RB_i_TB,...
                    H_i, W_i,ICI_power(:,UE_half_RBs));
                
                if 0 && (obj.attached_eNodeB.scheduler.clock.current_TTI > 2) % feedback inspection
                    display('SINR inspection');
                end
                TB_SINR_dB  = 10*log10(TB_SINR_lin);
                % Store post-equalization SINR for the link performance model
                obj.store_SINR_for_link_performance_model(tx_mode,...
                    TB_SINR_dB);
                
                % Calculate average preequalization SNR: This is a total SNR, the same as in the Link Level Simulator
                obj.SNR_avg_preequal = 10*log10(sum_macroscopic_RX_W/(obj.thermal_noise_W_RB*RB_grid_0.n_RB));

                % Calculate the wideband SINR
                
                % Indexing matrix for the pathlosses for wideband SINR calculation
                pathloss_idxs = zeros(size(RX_W_half_RB_i_UE));
                for i_=1:size(RX_W_half_RB_i_UE,3)
                    pathloss_idxs(:,:,i_, :) = i_;
                end
                
                if there_are_interferers
                    if iscolumn(avg_pathloss_lin_i)
                        avg_pathloss_lin_i_calc = avg_pathloss_lin_i.' ;
                    else
                        avg_pathloss_lin_i_calc = avg_pathloss_lin_i ;
                    end

                    if iscolumn(avg_shadow_lin_i)                    
                        avg_pathloss_shadow_i_calc = avg_shadow_lin_i.';
                    else
                        avg_pathloss_shadow_i_calc = avg_shadow_lin_i ;
                    end
                else
                    avg_pathloss_lin_i_calc = ones(size(pathloss_idxs));
                    avg_pathloss_shadow_i_calc = ones(size(pathloss_idxs));
                end

                obj.wideband_SINR = 10*log10(sum(reshape(RX_total_half_RB_layers_UE/avg_pathloss_linear_0/avg_shadow_fading_linear_0,1,[])) /...
                    (sum(reshape(RX_W_half_RB_i_UE./avg_pathloss_lin_i_calc(pathloss_idxs)./avg_pathloss_shadow_i_calc(pathloss_idxs),1,[]))...
                    +noise_W_half_RB*numel(RX_total_half_RB_layers_UE)));
                if obj.attached_eNodeB.attached_CoMP_site ~= []
                for j_ = 1:length(interfering_eNodeBs)
                    ind = j_*ones(length(TX_W_RB_0), 1);
                    wideband_SIR(j_) = 10*log10(sum(reshape(TX_W_RB_0/avg_pathloss_linear_0/avg_shadow_fading_linear_0,1,[])) /...
                    (sum(reshape(TX_W_RB_i(:, j_)'./avg_pathloss_lin_i_calc(ind)./avg_pathloss_shadow_i_calc(ind),1,[]))));
                    
                end
                obj.attached_eNodeB.attached_CoMP_site.feedback_wideband_SIR(obj, obj.attached_eNodeB, wideband_SIR, interfering_eNodeBs);
                end
                %obj.wideband_SINR = 10*log10(sum(reshape(RX_total_half_RB_layers_UE/avg_pathloss_linear_0/avg_shadow_fading_linear_0,1,[])) /...
                %    (sum(reshape(RX_W_half_RB_i_UE./avg_pathloss_lin_i(pathloss_idxs)./avg_shadow_lin_i(pathloss_idxs),1,[]))...
                %    +noise_W_half_RB*numel(RX_total_half_RB_layers_UE)));
            else                
%                 if config.feedback_channel_delay ~= 0 % in case of 0 delay feedback, the SINRs have to be stored after scheduling
                    obj.store_SINR_for_link_performance_model(tx_mode, -Inf);
                    obj.wideband_SINR       = -Inf;
                    obj.SNR_avg_preequal    = -Inf;
%                 else
%                     obj.Zero_FB_delay_SINR_store      = -Inf;                    
%                 end                              
            end
            
            % Feedback calculation
            if FB_delay_zero_run || (~FB_delay_zero_run && config.feedback_channel_delay ~= 0)
                switch tx_mode
                    case {4, 5, 6, 9}
                        % CLSM, Rank1-CLSM, MU-MIMO
                        precoding_codebook = [obj.codebook{:,nTX,tx_mode}];
                        filtered_idx = (interfering_eNodeB_ids ~=obj.attached_eNodeB.eNodeB_id);
                        [PMI_feedback,TB_SINR_feedback_dB] = feedback_calculation.CLSM_PMI_feedback(...
                            precoding_codebook,...
                            H_0, TX_W_RB_0, noise_W_half_RB,...
                            H_i_full(:,:,:,:,filtered_idx), TX_W_RB_i(:, filtered_idx, :), W_i_full(:, filtered_idx, :), ICI_power );

            
                        obj.feedback = feedback_calculation.CQI_and_rank(config,tx_mode,TB_SINR_feedback_dB,nRB,PMI_feedback,DL_signaling,...
                            obj.SINR_averager,obj.CQI_mapper,obj.unquantized_CQI_feedback,obj.adaptive_RI);
                    case 400
                        % CLSM with precoder calculation in the eNodeB.
                        % Implemented as an example of eNodeB-computed precoding
                        obj.feedback.tx_mode = tx_mode;
                        obj.feedback.H_0     = H_0;
                        obj.feedback.H_i     = H_i_full;
                        obj.feedback.W_i     = W_i_full;
                        obj.feedback.TX_W_RB_0 = TX_W_RB_0;
                        obj.feedback.TX_W_RB_i = TX_W_RB_i;
                        obj.feedback.ICI_power = ICI_power;
                    otherwise
                        error('Not yet supported (tx mode %d)',tx_mode);
                end
            end
        end
        
        % From the whole set of interferers, filters out the ones with a power
        % level above a certain threshold (in dB) relative to the signal. If not,
        % the SINR calculation would need to take into account all of the
        % interferers, which could be innecessary computationally costly.
        % There is one version for v2, which is able to deal with MU-MIMO interference,
        % and the original version still used in v1
        function [...
                interfering_eNodeB_ids,...  
                interfering_eNodeBs,...      
                pathloss_lin_i,...
                TX_W_i_per_RB,...
                shadow_lin_i,...
                avg_pathloss_lin,...
                avg_shadow_lin...
                ] = filter_out_interferers_below_SINR_v2(obj,...
                                            RX_total_W_0,...
                                            UE_RBs_mat,...
                                            interfering_eNodeBs,...
                                            SINR_threshold_dB)
            % note that in this calculation, the home eNodeB is allowed to be interfering                            
            interfering_eNodeBs = [obj.attached_eNodeB interfering_eNodeBs];
            parent_sites           = [interfering_eNodeBs.parent_eNodeB];
            parent_sites_id        = [parent_sites.id];
            interfering_eNodeB_ids = [interfering_eNodeBs.eNodeB_id];
            interfering_RB_grids   = [interfering_eNodeBs.RB_grid];
            
            for i_ = 1:length(interfering_eNodeB_ids)
               
                % RB_filter defines which RB parts are interpreted as
                % interference. For the home RB, these are all parts not
                % assigned to the calling UE, for all other eNodeBs all RB
                % parts are considered. 
                if interfering_eNodeB_ids(i_) == obj.attached_eNodeB.eNodeB_id
                    RB_filter = (~UE_RBs_mat);
                else
                    RB_filter = ones(size(UE_RBs_mat));
                end
                % poewer allocation of interferers in one matrix
                TX_W_data_i(:, i_, :)            = RB_filter.*[interfering_RB_grids(i_).power_allocation];
                TX_W_signaling_i(:, i_, :)      = RB_filter.*[interfering_RB_grids(i_).power_allocation_signaling];
                interfering_RB_parts(:, i_, :) = RB_filter;
            end
            % Preallocate assuming all eNodeBs to have the same number of RRHs
            nRRHs            = length(interfering_eNodeBs(1).RRHs);
            pathloss_lin_i   = zeros(length(interfering_eNodeB_ids),1+nRRHs);
            shadow_lin_i     = zeros(length(interfering_eNodeB_ids),1+nRRHs);
            avg_pathloss_lin = zeros(length(interfering_eNodeB_ids),1);
            avg_shadow_lin   = zeros(length(interfering_eNodeB_ids),1);
            
            % Get macroscopic pathloss and shadow fading values
            for i_=1:length(interfering_eNodeB_ids)
                [pathloss_lin_i(i_,:),shadow_lin_i(i_,:)] = obj.get_interfering_macroscale_losses(...
                    interfering_eNodeB_ids(i_),parent_sites_id(i_));
                avg_pathloss_lin(i_) = interfering_eNodeBs(i_).average_pathloss_value(pathloss_lin_i(i_,:));
                avg_shadow_lin(i_)   = interfering_eNodeBs(i_).average_pathloss_value(shadow_lin_i(i_,:));
            end
            
            % Total power allocations
            TX_W_i_per_RB  = TX_W_data_i + TX_W_signaling_i;
            TX_W_i_sum_RBs = reshape(sum(sum(TX_W_i_per_RB,1), 3),[],1);
            
            RX_avg_dB_0 = 10*log10(RX_total_W_0);
            RX_avg_dB_i = 10*log10(TX_W_i_sum_RBs./avg_pathloss_lin./avg_shadow_lin);
            CI_dB       = RX_avg_dB_0-RX_avg_dB_i;
            
            % Overwrite variables to take into consideration just the interferers up to 45dB below our signal
            interfererIdxs = CI_dB < SINR_threshold_dB;
            if sum(interfererIdxs)==0 % Just to avoid a crash
                interfererIdxs(1) = true;
            end
            % Filter all vectors regarding interferers
            interfering_eNodeB_ids = interfering_eNodeB_ids(interfererIdxs);
            pathloss_lin_i         = pathloss_lin_i(interfererIdxs,:);
            TX_W_i_per_RB          = TX_W_i_per_RB(:,interfererIdxs, :);
            interfering_eNodeBs    = interfering_eNodeBs(interfererIdxs);
            avg_pathloss_lin = avg_pathloss_lin(interfererIdxs);
            avg_shadow_lin   = avg_shadow_lin(interfererIdxs);
            if ~isscalar(shadow_lin_i)
                % Only if the shadow fading is not a scalar (i.e., there is shadow fading)
                shadow_lin_i = shadow_lin_i(interfererIdxs,:);
            else
                shadow_lin_i = ones(size(pathloss_lin_i),1);
            end
        end
        
        function [...
                interfering_eNodeB_ids,...
                interfering_eNodeBs,...
                pathloss_lin_i,...
                TX_W_i_per_RB,...
                shadow_lin_i,...
                avg_pathloss_lin,...
                avg_shadow_lin...
                ] = filter_out_interferers_below_SINR(obj,RX_total_W_0,interfering_eNodeBs,SINR_threshold_dB)
            
            parent_sites           = [interfering_eNodeBs.parent_eNodeB];
            parent_sites_id        = [parent_sites.id];
            interfering_eNodeB_ids = [interfering_eNodeBs.eNodeB_id];
            interfering_RB_grids   = [interfering_eNodeBs.RB_grid];
            TX_W_data_i            = [interfering_RB_grids.power_allocation];
            TX_W_signaling_i       = [interfering_RB_grids.power_allocation_signaling];
            
            % Preallocate assuming all eNodeBs to have the same number of RRHs
            nRRHs            = length(interfering_eNodeBs(1).RRHs);
            pathloss_lin_i   = zeros(length(interfering_eNodeB_ids),1+nRRHs);
            shadow_lin_i     = zeros(length(interfering_eNodeB_ids),1+nRRHs);
            avg_pathloss_lin = zeros(length(interfering_eNodeB_ids),1);
            avg_shadow_lin   = zeros(length(interfering_eNodeB_ids),1);
            
            % Get macroscopic pathloss and shadow fading values
            for i_=1:length(interfering_eNodeB_ids)
                [pathloss_lin_i(i_,:),shadow_lin_i(i_,:)] = obj.get_interfering_macroscale_losses(...
                    interfering_eNodeB_ids(i_),parent_sites_id(i_));
                avg_pathloss_lin(i_) = interfering_eNodeBs(i_).average_pathloss_value(pathloss_lin_i(i_,:));
                avg_shadow_lin(i_)   = interfering_eNodeBs(i_).average_pathloss_value(shadow_lin_i(i_,:));
            end
            
            % Total power allocations
            TX_W_i_per_RB  = TX_W_data_i + TX_W_signaling_i;
            TX_W_i_sum_RBs = reshape(sum(TX_W_i_per_RB,1),[],1);
            
            RX_avg_dB_0 = 10*log10(RX_total_W_0);
            RX_avg_dB_i = 10*log10(TX_W_i_sum_RBs./avg_pathloss_lin./avg_shadow_lin);
            CI_dB       = RX_avg_dB_0-RX_avg_dB_i;
            
            % Overwrite variables to take into consideration just the interferers up to 45dB below our signal
            interfererIdxs = CI_dB < SINR_threshold_dB;
            if sum(interfererIdxs)==0 % Just to avoid a crash
                interfererIdxs(1) = true;
            end
            
            interfering_eNodeB_ids = interfering_eNodeB_ids(interfererIdxs);
            pathloss_lin_i         = pathloss_lin_i(interfererIdxs,:);
            TX_W_i_per_RB          = TX_W_i_per_RB(:,interfererIdxs);
            interfering_eNodeBs    = interfering_eNodeBs(interfererIdxs);
            avg_pathloss_lin       = avg_pathloss_lin(interfererIdxs);
            avg_shadow_lin         = avg_shadow_lin(interfererIdxs);
            if ~isscalar(shadow_lin_i)
                % Only if the shadow fading is not a scalar (i.e., there is shadow fading)
                shadow_lin_i = shadow_lin_i(interfererIdxs,:);
            else
                shadow_lin_i = ones(size(pathloss_lin_i),1);
            end
        end
        % Stores the calculated TB SINRs for use by the link performance model
        function store_SINR_for_link_performance_model(obj,tx_mode,SINR_dB)
            switch tx_mode
                case 2 % TxD
                    % Both layers have the same SINR
                    obj.link_quality_model_output.SINR_dB = SINR_dB(1,:);
                otherwise
                    % SISO, OLSM, CLSM {1,3,4}
                    % And any other mode. Add other cases if a specific
                    % transmit mode would need another mapping
                    obj.link_quality_model_output.SINR_dB = SINR_dB;
            end
        end
        
        % Evaluate whether this TB arrived correctly by using the data from
        % the link quality model and feeding it to the link performance
        % model (BLER curves). Just one version of it, so no need for a
        % gateway function as in the link quality model.
        function link_performance_model(obj)
            if ~obj.deactivate_UE || obj.runtime_precoding
                % Get SINRs from the link quality model. Only the dB (not linear) are needed.
                SINR_dB = obj.link_quality_model_output.SINR_dB;

                % Calculate TB SINR
                [...
                    TB_CQI,...
                    user_RBs,...
                    assigned_RBs,...
                    assigned_power,...
                    tx_mode,...
                    nLayers,...
                    nCodewords,...
                    rv_idxs,...
                    TB_size,...
                    N_used_bits,...
                    packet_parts...
                    ] = obj.eNodeB_signaling.get_TB_params;
                
                % Preallocate variables to store in trace
                TB_SINR_dB = zeros(1,nCodewords);
                BLER       = zeros(1,nCodewords);
                
                % Setting this "Kronecker product" to [1 0] would mean that you only take one SC per RB and repeat it. i.e. not use all of
                % your SCs in the trace.
                % Set feedback for all streams
                if assigned_RBs>0
                    % Layer mapping according to TS 36.212
                    for cw_=1:nCodewords
                        switch cw_
                            case 1
                                switch nLayers
                                    case {4 5}
                                        layers_cw = [1 2];
                                    case {6 7}
                                        layers_cw = [1 2 3]; 
                                    case 8
                                        layers_cw = [1 2 3 4];
                                    otherwise
                                        layers_cw = 1;
                                end
                            case 2
                                switch nLayers
                                    case 1
                                        error('2 codewords and 1 layers not allowed');
                                    case 2
                                        layers_cw = 2;
                                    case 3
                                        layers_cw = [2 3];
                                    case 4
                                        layers_cw = [3 4];
                                    case 5
                                        layers_cw = [3 4 5];
                                    case 6
                                        layers_cw = [4 5 6];
                                    case 7
                                        layers_cw = [4 5 6 7];
                                    case 8
                                        layers_cw = [5 6 7 8];
                                end
                        end
                        [TB_SINR_dB(cw_)] = obj.SINR_averager.average(SINR_dB(layers_cw,:),TB_CQI(cw_),true);
                        BLER(cw_)         = obj.BLER_curves.get_BLER(TB_CQI(cw_),TB_SINR_dB(cw_));
                    end
                    
                    % Receive
                    ACK = BLER<rand(1,nCodewords);
                else
                    % Dummy results
                    TB_SINR_dB = [];
                    ACK        = false(1,nCodewords);
                end
                
                % Needed for non-full-buffer simulations
                if ~obj.traffic_model.is_fullbuffer
                    obj.process_packet_parts(ACK, packet_parts,nCodewords);
                end
                
                % Add BLER/ACK feedback to the CQI and RI feedback
                if assigned_RBs~=0
                    obj.feedback.UE_scheduled = true;
                    obj.feedback.nCodewords   = nCodewords;
                    obj.feedback.TB_size      = TB_size;
                    obj.feedback.BLER         = BLER;
                    obj.feedback.ACK          = ACK;
                else
                    obj.feedback.UE_scheduled = false;
                    obj.feedback.nCodewords   = 0;
                    obj.feedback.TB_size      = 0;
                    obj.feedback.BLER         = NaN;
                    obj.feedback.ACK          = false;
                end
                
                % Optional traces
                if obj.trace_SINR
                    extra_traces{1} = obj.link_quality_model_output.SINR_dB;
                else
                    extra_traces{1} = [];
                end
                
                % Store trace of the relevant information
                tti_idx = obj.clock.current_TTI;
                % Store trace
                obj.trace.store(...
                    obj.feedback,...
                    nLayers,...
                    nCodewords,...
                    obj.attached_eNodeB,...
                    obj.pos,...
                    tti_idx,...
                    assigned_RBs,...
                    assigned_power,...
                    TB_CQI,...
                    TB_size,...
                    BLER,...
                    TB_SINR_dB,...
                    N_used_bits,...
                    obj.wideband_SINR,...
                    obj.deactivate_UE,...
                    obj.SNR_avg_preequal,...
                    obj.rx_power_tb_in_current_tti,...
                    obj.rx_power_interferers_in_current_tti,...
                    extra_traces);
            else
                obj.dummy_link_performance_model;
            end
        end
        
        % Dummy functions
        function dummy_link_quality_model(obj,config)
            obj.SNR_avg_preequal      = NaN;
            obj.wideband_SINR         = NaN;
            obj.feedback.CQI          = zeros(1,config.N_RB);
            obj.feedback.PMI          = ones(1,config.N_RB);
            obj.feedback.RI           = 1;
            obj.feedback.tx_mode      = 4;
            obj.feedback.UE_scheduled = false;
            obj.feedback.nCodewords   = 1;
            obj.feedback.TB_size      = 0;
            obj.feedback.BLER         = 0;
            obj.feedback.ACK          = false;
        end
        
        function dummy_link_performance_model(obj)
            % Store trace of the relevant information
            tti_idx = obj.clock.current_TTI;
            
            % Store trace
            extra_traces{1} = [];
            extra_traces{2} = [];
            obj.trace.store(...
                obj.feedback,...
                0,...
                0,...
                obj.attached_eNodeB,...
                obj.pos,...
                tti_idx,...
                0,...
                0,...
                0,...
                0,...
                0,...
                NaN,...
                0,...
                obj.wideband_SINR,...
                obj.deactivate_UE,...
                NaN,...
                0,...
                0,...
                extra_traces);
        end
        
        function distance_to_site = distance_to_attached_site(obj)
            % Return the distance to the attached site (function added for convenience)
            distance_to_site = sqrt(sum((obj.attached_site.pos -obj.pos).^2));
        end
        
        % Clear all non-basic info and leaves just basic information describing the UE
        function clear_non_basic_info(obj)
            obj.walking_model             = [];
            obj.downlink_channel          = [];
            obj.RB_grid                   = [];
            obj.uplink_channel            = [];
            obj.trace                     = [];
            obj.clock                     = [];
            obj.CQI_mapper                = [];
            obj.link_quality_model_output = [];
            obj.feedback                  = [];
            obj.unquantized_CQI_feedback  = [];
            obj.BLER_curves               = [];
            obj.SINR_averager             = [];
            obj.codebook                  = [];
            obj.eNodeB_signaling          = [];
            obj.trace_SINR                = [];
            obj.SNR_avg_preequal          = [];
            obj.wideband_SINR             = [];
            obj.deactivate_UE             = [];
            obj.default_tx_mode           = [];
            obj.traffic_model             = [];
            obj.lambda                    = [];
        end
        
        % Returns a struct containing the basic information (not deleted with the previous function) from the UE
        function struct_out = basic_information_in_struct(obj)
            struct_out.id                    = obj.id;
            struct_out.pos                   = obj.pos;
            if ~isempty(obj.attached_site)
                struct_out.attached_site = obj.attached_site.basic_information_in_struct();
            else
                struct_out.attached_site = [];
            end
            struct_out.attached_sector_idx   = obj.attached_sector_idx;
            if ~isempty(obj.attached_eNodeB)
                struct_out.attached_eNodeB = obj.attached_eNodeB.basic_information_in_struct();
            else
                struct_out.attached_eNodeB = [];
            end
            struct_out.receiver_noise_figure = obj.receiver_noise_figure;
            struct_out.thermal_noise_W_RB    = obj.thermal_noise_W_RB;
            struct_out.penetration_loss      = obj.penetration_loss;
            struct_out.nRX                   = obj.nRX;
            struct_out.antenna_gain          = obj.antenna_gain;
        end
        
        % Processes the packet parts
        function process_packet_parts(obj,ACK,packet_parts,nCodewords)
            for cw_ = 1:nCodewords
                if ACK(cw_)
                    if strcmp(obj.traffic_model.type,'voip') || strcmp(obj.traffic_model.type,'video') || strcmp(obj.traffic_model.type,'gaming')
                        for pp = 1:length(packet_parts{cw_}) % acknowledge all packet parts and remove them from the buffer
                            if packet_parts{cw_}(pp).data_packet_id
                                packet_ind = obj.traffic_model.get_packet_ids == packet_parts{cw_}(pp).data_packet_id;
                                if sum(packet_ind)
                                    [packet_done,packet_id] = obj.traffic_model.packet_buffer(packet_ind).acknowledge_packet_part(packet_parts{cw_}(pp).id,true);
                                    if packet_done && packet_id
                                        obj.traffic_model.remove_packet(packet_id,true);
                                    end
                                end
                            end
                        end
                    end
                else
                    % NOTE: activate me if HARQ exists
                    if strcmp(obj.traffic_model.type,'voip') || strcmp(obj.traffic_model.type,'video') || strcmp(obj.traffic_model.type,'gaming')
                        %                                             if current_rv_idx == max_rv_idx % if maximum of retransmissions is obtained --> delete the packet parts
                        for pp = 1:length(packet_parts{cw_})
                            if packet_parts{cw_}(pp).data_packet_id
                                packet_ind = obj.traffic_model.get_packet_ids == packet_parts{cw_}(pp).data_packet_id;
                                if sum(packet_ind)
                                    [packet_done,packet_id] = obj.traffic_model.packet_buffer(packet_ind).acknowledge_packet_part(packet_parts{cw_}(pp).id,false); % with the option "false" (last argument) the packet is deleted
                                    if packet_done && packet_id
                                        obj.traffic_model.remove_packet(packet_id,false); % with the option false, the packet is marked as non-successfully transmitted
                                    end
                                end
                            end
                        end
                        %                                             % else restore them for retransmission
                        %                                                  if strcmp(obj.traffic_model.type,'voip') || strcmp(obj.traffic_model.type,'video') || strcmp(obj.traffic_model.type,'gaming')
                        %                                                     for pp = 1:length(packet_parts{cw_})
                        %                                                         if packet_parts{cw_}(pp).data_packet_id
                        %                                                             packet_ind = obj.traffic_model.get_packet_ids == packet_parts{cw_}(pp).data_packet_id;
                        %                                                             if sum(packet_ind)
                        %                                                                 obj.traffic_model.packet_buffer(packet_ind).restore_packet_part(packet_parts{cw_}(pp).id);
                        %                                                             end
                        %                                                         end
                        %                                                     end
                        %                                                  end
                        %                                             end
                    end
                end
            end
        end
    end
end
