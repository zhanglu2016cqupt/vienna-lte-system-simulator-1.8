classdef sequential < CoMP.schedulers.CoMP_scheduler
    % This is the absolute easiest example of using the CoMP handles
    % everything
    
    properties
    end
    
    methods
        function obj =  sequential(site)
            obj = obj@CoMP.schedulers.CoMP_scheduler(site);
        end
        
        
       function schedule_users(obj, eNodeB, attached_UEs, last_received_feedbacks)
           eNodeB.scheduler.schedule_users(attached_UEs, last_received_feedbacks)
           obj.turn_off_interferers(eNodeB.scheduler);
           
       end
       function turn_off_interferers(obj, scheduler)
       end
    end
    
end

