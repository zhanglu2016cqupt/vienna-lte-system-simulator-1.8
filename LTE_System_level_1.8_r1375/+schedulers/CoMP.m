classdef CoMP < schedulers.lteScheduler
% CoMP wrapper
% Thomas Blazek
% (c) 2014 by ITC
% www.nt.tuwien.ac.at

   properties
       % See the lteScheduler class for a list of inherited attributes
       av_throughput % exponentially weighted throughputs
       CoMP_site
   end

   methods
       
       % Class constructor. Just specify where to attach the scheduler
       function obj = CoMP(scheduler_params,attached_eNodeB_sector)
           % Fill in basic parameters (handled by the superclass constructor)
           obj       = obj@schedulers.lteScheduler(scheduler_params,attached_eNodeB_sector);
           obj.CoMP_site = attached_eNodeB_sector.attached_CoMP_site;
%            obj.alpha = scheduler_params.alpha;
           obj.name  = 'CoMP wrapper scheduler';
       end
       
       % Dummy functions required by the lteScheduler Abstract class implementation
       % Add UE (no memory, so empty)
       function add_UE(obj,UE_id)
       end
       % Delete UE (no memory, so empty)
       function remove_UE(obj,UE_id)
       end
       
       % Schedule the users in the given RB grid
       function schedule_users(obj,attached_UEs,last_received_feedbacks)
           % This is a dummy scheduler. It fetches the Scheduled grid from
           % the CoMP site. The comp site in turn does the scheduling for
           % all attached eNodeBs simultaneously
           RB_grid = obj.RB_grid;
           RB_grid.size_bits = 0;
           tx_mode           = obj.default_tx_mode;
           current_TTI       = obj.clock.current_TTI; 
           
           N_RB = RB_grid.n_RB;
           UE_id_list = zeros(N_RB,1);
           
           if ~isempty(attached_UEs)
               UE_id_list = obj.CoMP_site.scheduler.fetch_grid(obj.attached_eNodeB,attached_UEs,last_received_feedbacks);
               %last_received_feedbacks = obj.CoMP_site.scheduler.adapt_feedback(attached_UEs, last_received_feedbacks);
               RB_grid.user_allocation(:) = UE_id_list;
               % CQI assignment. TODO: implement HARQ          
               obj.schedule_users_common(attached_UEs,last_received_feedbacks,current_TTI,tx_mode);
           end
       end
   end
end 
