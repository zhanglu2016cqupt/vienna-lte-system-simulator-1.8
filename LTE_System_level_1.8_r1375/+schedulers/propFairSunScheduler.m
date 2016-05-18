classdef propFairSunScheduler < schedulers.lteScheduler
% A proportional fair LTE scheduler 
% (c) Stefan Schwarz, INTHFT, 2010

   properties
       % See the lteScheduler class for a list of inherited attributes
       av_throughput % exponentially weighted throughputs
   end

   methods
       
       % Class constructor. Just specify where to attach the scheduler
       function obj = propFairSunScheduler(scheduler_params,attached_eNodeB_sector)
           % Fill in basic parameters (handled by the superclass constructor)
           obj      = obj@schedulers.lteScheduler(scheduler_params,attached_eNodeB_sector);
           obj.name = 'Proportional fair Sun scheduler';
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
           N_UE        = length(attached_UEs);
           N_RB        = RB_grid.n_RB;
           UE_id_list  = zeros(N_RB,1);
           if ~isempty(attached_UEs)
                %% compute efficiency
                [c,user_ind] = obj.get_efficiency(N_UE,N_RB,last_received_feedbacks);          
                c = c';              
                             
               %% update average throughput
               TTI_to_read = max(current_TTI-obj.feedback_delay_TTIs-1,1); % Realistically read the ACKed throughput
               for uu = 1:N_UE
                    obj.av_throughput(uu) = obj.compute_av_throughput(uu,last_received_feedbacks,TTI_to_read);
               end
               
               %% PF scheduler
               RBs = obj.PF_scheduler(N_UE,N_RB,c,user_ind);
               
               for r_ = 1:N_RB
                   RB_tmp = RBs((r_-1)*N_UE+1:r_*N_UE);
                   ind = find(RB_tmp == 1);
                  if ~isempty(ind)
                    UE_id_list(r_) = attached_UEs(user_ind(ind)).id;
                  end
               end              
               RB_grid.user_allocation(:) = UE_id_list;
               % CQI assignment. TODO: implement HARQ          
               obj.schedule_users_common(attached_UEs,last_received_feedbacks,current_TTI,tx_mode);
           end
       end
       
       function RBs = PF_scheduler(obj,N_UE,N_RB,c,user_ind)
           % Core scheduling function (same in LL and SL)
           RB_set     = true(N_RB,1);
           RB_UEs     = false(N_RB,N_UE);
           alpha_temp = 1;
           % res        = find(RB_set);     % The RB set over which you will schedule (right now the same as the RB_set. Specified in linear indices. Not used in the newer implementation
           
           % Precalculated values taken out from the loop (speeds up simulations)
           cleaned_c_log_matrix = log10(max(c,eps)*12*7);
           avgd_UE_throughputs  = (obj.av_const-1)*obj.av_throughput(user_ind);
           
           % Calculate metric for each RB and attached user
           for rr = 1:N_RB
               res                    = find(RB_set);
               metric                 = -Inf(N_RB,N_UE);
               UE_avgd_pre_metric     = -alpha_temp*log10(max(avgd_UE_throughputs+sum(RB_UEs.*c,1)*12*7,eps));
               UE_avgd_pre_metric_mat = UE_avgd_pre_metric(ones(1,N_RB),:);
           
               metric(res(1:sum(RB_set)),:) = cleaned_c_log_matrix(res(1:sum(RB_set)),:)+UE_avgd_pre_metric_mat(res(1:sum(RB_set)),:);
               % for u_ = 1:N_UE
               % for r_ = 1:N_RB
               % metric(res(r_),u_) = c(res(r_),u_)*12*7/((obj.av_const-1)*obj.av_throughput(user_ind(u_))+RB_UEs(:,u_).'*c(:,u_)*12*7);      % 12*7 equals the number of elements in a RB
               % metric(res(r_),u_) = log10(max(c(res(r_),u_),eps)*12*7)-alpha_temp*log10(max((obj.av_const-1)*obj.av_throughput(user_ind(u_))+RB_UEs(:,u_).'*c(:,u_)*12*7,eps)); % Old implementation
               % metric(res(r_),u_) = log10(max(c(res(r_),u_)*12*7,eps))-alpha_temp*log10(max(obj.av_throughput(user_ind(u_)),eps));
               % end
               % end
               maxi            = max(metric(:));
               [RB_idx UE_idx] = find(metric == maxi);
               ind             = randi(length(RB_idx));
               
               RB_set(RB_idx(ind))             = false;
               RB_UEs(RB_idx(ind),UE_idx(ind)) = true;
           end
           RB_UEs = RB_UEs';
           RBs = RB_UEs(:);
       end
   end
end 
