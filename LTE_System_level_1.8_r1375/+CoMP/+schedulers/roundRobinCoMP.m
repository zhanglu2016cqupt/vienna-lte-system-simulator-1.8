classdef roundRobinCoMP < CoMP.schedulers.CoMP_scheduler
    % (c) Thomas Blazek, 2014,ITC,  thomas.blazek@nt.tuwien.ac.at
    % This is an example comp scheduler. It mirrors the behaviour of the
    % standard Round Robin scheduler (extraction from ringbuffer), but does
    % so simultaneously for all
    
    properties
        
        UE_matrix
        last_scheduled
    end
    
    methods
        function obj = roundRobinCoMP(site)
            obj = obj@CoMP.schedulers.CoMP_scheduler(site);
            obj.user_allocations = cell(size(site.attached_eNodeBs));
            
            obj.UE_matrix = [];
            obj.last_scheduled = ones(1, length(site.attached_eNodeBs));
        end
        
        
        
        function ues = extract_users(obj,nr_ues, eNodeB_idx)
            % Extract from ringbuffer
            mod_length = length(obj.UE_matrix{eNodeB_idx});
            idx_array = mod(((1:nr_ues)+obj.last_scheduled(eNodeB_idx)-1), mod_length)+1;
            ue_temp = obj.UE_matrix{eNodeB_idx}(idx_array);
            obj.last_scheduled(eNodeB_idx) = idx_array(end);
            ues = [ue_temp.id];
        end
        
        function schedule_users(obj, eNodeB,attached_UEs,last_received_feedbacks)
            eNodeB_array = obj.CoMP_site.attached_eNodeBs;
            for i_ = 1:length(eNodeB_array) % set the values
                sched(i_) = eNodeB_array(i_).scheduler;
                RB_grid(i_) = sched(i_).RB_grid;
                obj.user_allocations{i_}= zeros(size(RB_grid(i_).user_allocation));
                RB_grid(i_).user_allocation = zeros(size(RB_grid(i_).user_allocation));
                RB_grid(i_).size_bits = 0;
                tx_mode(i_)           = sched(i_).default_tx_mode;
                current_TTI(i_)       = sched(i_).clock.current_TTI;
            end
            
            site = obj.CoMP_site;
            for rb_ = 1:size(RB_grid(1).user_allocation,1)
                current_UEs = []; % array of UEs scheduled in same RB (in different eNodeBs)
                for i_ =randperm(length(eNodeB_array)) %randperm so no eNodeB is always first
                    forbidden_users = site.get_forbidden_UEs(eNodeB_array(i_));
                    % checks which users are not allowed alongside this
                    % eNodeB
                    
                    if isempty(current_UEs) ||isempty(forbidden_users)||...
                            sum(sum(bsxfun(@eq, forbidden_users, current_UEs)))==0
                        
                        next_ue =  obj.extract_users(1, i_);
                        dont_schedule = 0;
                        if(~isempty(site.cooperation_rules))
                            % check which eNodeBs are not allowed along the
                            % chosen UE
                            forbidden_eNodeBs_idx = [[site.cooperation_rules.UE_id] == next_ue];
                            forbidden_eNodeBs = [site.cooperation_rules(forbidden_eNodeBs_idx).dont_schedule];
                            for j_ = 1:length(forbidden_eNodeBs)
                                if obj.user_allocations{ [eNodeB_array.eNodeB_id] == forbidden_eNodeBs(j_).eNodeB_id}~=0
                                    % if a forbidden UE already has a RB
                                    % assigned, dont schedule
                                    dont_schedule =1;
                                    break;
                                end
                                
                            end
                        end
                        if ~dont_schedule
                            % Schedule user
                            obj.user_allocations{i_}(rb_) = next_ue;
                            current_UEs = [current_UEs next_ue];
                        else
                            % Push the ringbuffer back, so the user is not
                            % left out on the next run
                            obj.last_scheduled(i_) = mod(obj.last_scheduled(i_)-2, length(obj.UE_matrix{i_}))+1;
                        end
                    end
                end
            end
            
            
        end
        
        
    end
    
end

