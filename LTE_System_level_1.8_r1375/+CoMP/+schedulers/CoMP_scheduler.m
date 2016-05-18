classdef CoMP_scheduler < handle
    % (c) Thomas Blazek, 2014,ITC,  thomas.blazek@nt.tuwien.ac.at
    % The generic CoMP scheduler. It has to be able to tell which eNodeB
    % was already scheduled, and it has to be able to track set RB_grids
    
    properties
        already_scheduled
        user_allocations
        CoMP_site
    end
    
    methods
        function obj = CoMP_scheduler(site)
            obj.already_scheduled = zeros(1, length(site.attached_eNodeBs));
            obj.CoMP_site = site;
        end
        
        function clear_scheduler(obj)
            obj.already_scheduled = zeros(1, length( obj.CoMP_site.attached_eNodeBs));
        end
        function schedule_users(obj, eNodeB, attached_UEs, last_received_feedbacks)
        end
        
        function grid = fetch_grid(obj, eNodeB,attached_UEs,last_received_feedbacks)
            % This function is called by the LTE_scheduler 'CoMP'. Instead
            % of producing a RB_grid itself, it asks the CoMP-site for the
            % assigned grid. If a current grid is existent, it is assigned,
            % if not, it is generated. If all current grids are assigned
            % (every attached eNodeB has fetched their RB_grid), the grid
            % is deleted
            obj.update_UEs;
            eNodeB_array = obj.CoMP_site.attached_eNodeBs;
            eNodeB_idx = eNodeB_array == eNodeB;
            if sum(obj.already_scheduled) == 0
                obj.schedule_users(eNodeB,attached_UEs,last_received_feedbacks);
            end
            grid = obj.user_allocations{eNodeB_idx};
            obj.already_scheduled(eNodeB_idx) = 1;
            if prod(obj.already_scheduled) == 1
                obj.already_scheduled(:) = 0;
                obj.user_allocations = {};
            end
            
        end
        function update_UEs(obj)
            % This method looks for changes in the attached UEs at the
            % beginning of the TTI.
            eNodeB_array = obj.CoMP_site.attached_eNodeBs;
            for i_ = 1:length(eNodeB_array)
                UE_vector{i_} = eNodeB_array(i_).attached_UEs_vector;
            end
            if ~isequal(obj.UE_matrix, UE_vector)
                obj.add_and_remove_users(UE_vector);
            end
        end
        
        function add_and_remove_users(obj, UE_vector);
            obj.UE_matrix = UE_mat;
        end
        
        function feedback = adapt_feedback(obj, attached_UEs, received_feedback)
            % This function illustrates how one could adapt the CQI
            % feedback based on the CoMP scheduling. This function is
            % called from the LTE_scheduler, before passing to
            % schedule_users_common
            feedback = received_feedback;
            for i_ =  1:length(received_feedback.UE_id)
                if received_feedback.feedback_received(i_)
                    UE_id = received_feedback.UE_id(i_);
                    idx = [obj.CoMP_site.cooperation_rules.UE_id] == UE_id;
                    nr_canceled = length(obj.CoMP_site.cooperation_rules(idx).dont_schedule);
                    feedback.CQI(:,:,i_) = feedback.CQI(:,:,i_) + logical(nr_canceled)*ones(size(feedback.CQI(:,:,i_)));
                    
                    
                end
            end
        end
    end
    
end

