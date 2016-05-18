classdef roundRobinSchedulerMU < schedulers.MU_MIMO.lteSchedulerMU
    % An LTE round Robin scheduler.
    % (c) Josep Colom Ikuno, INTHFT, 2008
    
    properties
        % Where the scheduler will store which users to serve first (round robin fashion)
        UE_queue
        last_extracted
        length
        
        % See the lteScheduler class for a list of inherited attributes
    end
    
    methods
        
        % Class constructor. UE_queue size needs to be specified large
        % enough so it won't overflow
        function obj = roundRobinSchedulerMU(scheduler_params,attached_eNodeB_sector)
            % Fill in basic parameters (handled by the superclass constructor)
            obj      = obj@schedulers.MU_MIMO.lteSchedulerMU(scheduler_params,attached_eNodeB_sector);
            obj.name = 'Round Robin scheduler MU';
        end
        
        % Add a UE to the queue. It could be done so each TTI the scheduler
        % gets a UE list from the eNodeB, but such a query is not necessary.
        % Just updating when a UE attaches or drops is sufficient.
        function add_UE(obj,UE_id)
            % If not in queue: add
            if ~sum(obj.UE_queue==UE_id)
                if isempty(obj.UE_queue)
                    obj.last_extracted = 1;
                end
                obj.UE_queue = [obj.UE_queue UE_id];
                obj.length   = length(obj.UE_queue);
            end
        end
        
        % Delete an UE_id from the queue
        function remove_UE(obj,UE_id)
            % Remove the UE with this UE_id
            removed_set  = obj.UE_queue~=UE_id;
            obj.UE_queue = obj.UE_queue(removed_set);
            obj.length   = length(obj.UE_queue);
            
            % Adjust the last_extracted variable
            if ~isempty(find(removed_set,1))
                if obj.last_extracted>=find(removed_set,1)
                    obj.last_extracted = mod(obj.last_extracted-2,obj.length)+1; % One-indexed modulo-length adding of one
                end
            end
            
            if obj.length==0
                obj.last_extracted = [];
            end
        end
        
        % Next user to serve. If the queue is empty, returns 0
        function UE_id = get_next_users(obj,number)
            % Return the first item and shift the whole thing one position
            to_extract         = mod(obj.last_extracted:(obj.last_extracted+number-1),obj.length)+1;
            UE_id              = obj.UE_queue(to_extract);
            obj.last_extracted = to_extract(end);
        end
        
        % Schedule the users in the given RB grid
        function schedule_users(obj,attached_UEs,last_received_feedbacks)
            % Power allocation
            % Nothing here. Leave the default one (homogeneous)
            
            % For now use the static tx_mode assignment
            RB_grid = obj.RB_grid;
            RB_grid.size_bits = 0;
            tx_mode           = obj.default_tx_mode;
            current_TTI       = obj.clock.current_TTI;
            %             if ~sum(last_received_feedbacks.feedback_received == 0
            
            %                 %precoders_min = obj.get_precoder_map(1, obj.attached_eNodeB.total_nTX, tx_mode, ones(size(RB_grid.user_allocation, 1)), last_received_feedbacks.nonstandard.PMI_min);
            %             end
            RB_grid.user_allocation = zeros(size(RB_grid.user_allocation));
            if ~isempty(attached_UEs)
                % Assign RBs to UEs, RR fashion
                for rb_ = 1:size(RB_grid.user_allocation, 1)
                    current_ue =  obj.get_next_users(1);
                    RB_grid.user_allocation(rb_, 1) =   current_ue;
                    if current_TTI >3
                        precoders = obj.get_precoder_map(1, obj.attached_eNodeB.total_nTX, tx_mode, ones(size(RB_grid.user_allocation, 1)), last_received_feedbacks.PMI);
                        next_ind = 2;
                        schedule_type = 2;
                        % type 1: orthogonal scheduling
                        % type 2: Scheduling based on SINR degradation
                        switch schedule_type
                            case 1
                                rb_precoders = [precoders(rb_, :).W];
                                check_matrix = (0 == rb_precoders'*rb_precoders);
                                
                                possible_ues =  obj.get_orthogonal_UEs(check_matrix, attached_UEs, current_ue);
                                %while(~isempty(possible_ues))
                                if length(possible_ues)>0
                                    next_ue = randi(length(possible_ues));
                                    RB_grid.user_allocation(rb_, next_ind) = possible_ues(next_ue);
                                    next_ind = next_ind +1;
                                    current_ue = [current_ue  possible_ues(next_ue)];
                                    possible_ues =  obj.get_orthogonal_UEs(check_matrix, attached_UEs, current_ue);
                                end
                            case 2
                                % For this mode, insert the following text
                                % block into
                                % +feedback_calculation.CLSM_PMI_feedback.m:
                                %==============CUT HERE====================
%                                 if precoding_codebook.tx_mode == 5
%     % caluclate nonstandard genie information for mode 5
%     max_ues = size(precoding_codebook.W, 1);
%     for i_ =2:2
%         TX_W_RB_i_int = zeros ( size(TX_W_RB_i, 1), size(TX_W_RB_i, 2)+1, size(TX_W_RB_i, 3));
%         TX_W_RB_i_int(:, 1) = (i_-1)/i_ .* sum(TX_W_RB_i(:, 1, :),3);
%         TX_W_RB_i_int(:, 2:end, :) = TX_W_RB_i;
%         
%         for j_ = 1:length(precoding_codebook.codebook_index)
%             H_i_int = zeros(size(H_i,1),size(H_i,2),size(H_i,3),size(H_i,4)+1);
%             H_i_int(:,:,:,1) = H_0;
%             H_i_int(:,:,:,2:end) = H_i;
%             code.W = precoding_codebook(1).W(:,:,j_);
%             W_i_int = W_i;
%             W_i_int(:, end+1, :) = code;
%             W_i_int(:, 1, :) = W_i_int(:, end, :);
%             W_i_int(:, 2:end, :) = W_i;
%      
%             SINR_nonstandard(:, i_-1, j_) = phy_modeling.post_equalization_SINR(...
%         1/i_*TX_power_RB_0,...
%         H_0, W_0,...
%         noise_W_half_RB,...
%         TX_W_RB_i_int,...
%         H_i_int, W_i_int);
%         end
%     end
%     PMI_feedback.nonstandard.additional_SINR = 10*log10(SINR_nonstandard);
%     PMI_feedback.nonstandard.original_SINR = 10*log10(TB_SINR_feedback_lin);
% end
% ====================================CUT Here ===================
                                nonstandard = last_received_feedbacks.nonstandard;
                                for i_ = 1:length(nonstandard)
                                    nonstandard(i_).acceptable_additions = ...
                                        [(kron(ones(1, size(nonstandard(i_).additional_SINR, 3)), nonstandard(i_).original_SINR')...
                                        -squeeze(nonstandard(i_).additional_SINR)) < 5];
                                end
                                
                                ue_idx =[attached_UEs.id] ==current_ue;
                                prec = last_received_feedbacks.PMI(rb_,ue_idx);
                                allowed_precoders = nonstandard(ue_idx).acceptable_additions(2*rb_-1:2*rb_, 1, :);
                                allowed_precoders = find(allowed_precoders(1,:).* allowed_precoders(2,:));
                                if(~isempty(allowed_precoders))
                                    candidate_UEs = logical(sum(bsxfun(@eq, last_received_feedbacks.PMI (rb_, :), allowed_precoders'),1));
                                    candidate_UEs(ue_idx) = 0;
                                    allowed_UEs = [];
                                    candidate_idx = find(candidate_UEs);
                                    for i_ = (candidate_idx)
                                        allowed_UEs = [allowed_UEs nonstandard(i_).acceptable_additions(2*rb_-1:2*rb_, prec)];
                                        
                                    end
                                    if sum(sum(allowed_UEs)) > 0;
                                        allowed_UEs = allowed_UEs(1,:).*allowed_UEs(2,:);
                                        if sum(allowed_UEs)>0
                                            add_UEs = candidate_idx(logical(allowed_UEs));
                                            rand_idx = randi([1 length(add_UEs)]);
                                            RB_grid.user_allocation(rb_, 2) = attached_UEs(add_UEs(rand_idx)).id;
                                        end
                                    end
                                end
                            otherwise
                                error('not supported')
                        end
                    end
                end
                
                obj.schedule_users_common(attached_UEs,last_received_feedbacks,current_TTI,tx_mode);
            end
        end
    end
end
