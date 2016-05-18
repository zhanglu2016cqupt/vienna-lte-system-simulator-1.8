classdef bestCqiScheduler < schedulers.lteScheduler
% A best CQI LTE scheduler.
% (c) Josep Colom Ikuno, INTHFT, 2009

   properties
       % See the lteScheduler class for a list of inherited attributes
   end

   methods
       
       % Class constructor. Just specify where to attach the scheduler
       function obj = bestCqiScheduler(scheduler_params,attached_eNodeB_sector)
           % Fill in basic parameters (handled by the superclass constructor)
           obj      = obj@schedulers.lteScheduler(scheduler_params,attached_eNodeB_sector);
           obj.name = 'Best CQI scheduler';
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
           % Power allocation
           % Nothing here. Leave the default one (homogeneous)
           
           RB_grid = obj.RB_grid;
           RB_grid.size_bits = 0;
           
           % For now use the static tx_mode assignment
           RB_grid.size_bits = 0;
           tx_mode     = obj.default_tx_mode;
           current_TTI = obj.clock.current_TTI;
           
           if ~isempty(attached_UEs)
               % Sum-of-CQIs-wise maximization done via the calculation of the predicted sum spectral efficiency
               quantized_feedback_efficiency     = obj.get_spectral_efficiency(last_received_feedbacks.CQI);
               quantized_feedback_efficiency_sum = reshape(sum(quantized_feedback_efficiency,2),[RB_grid.n_RB length(attached_UEs)]);
               UE_id_list                        =  obj.get_max_UEs(quantized_feedback_efficiency_sum.',last_received_feedbacks.UE_id);
               
               % Fill in RB grid
               RB_grid.user_allocation(:) = UE_id_list;
               obj.schedule_users_common(attached_UEs,last_received_feedbacks,current_TTI,tx_mode);
           end
       end
   end
end 
