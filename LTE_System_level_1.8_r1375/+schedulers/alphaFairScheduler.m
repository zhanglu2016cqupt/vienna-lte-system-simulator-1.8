classdef alphaFairScheduler < schedulers.lteScheduler
% A proportional fair LTE scheduler 
% (c) Stefan Schwarz, INTHFT, 2010

   properties
       % See the lteScheduler class for a list of inherited attributes
       av_throughput % exponentially weighted throughputs
       alpha
   end

   methods
       
       % Class constructor. Just specify where to attach the scheduler
       function obj = alphaFairScheduler(scheduler_params,attached_eNodeB_sector)
           % Fill in basic parameters (handled by the superclass constructor)
           obj       = obj@schedulers.lteScheduler(scheduler_params,attached_eNodeB_sector);
           obj.alpha = scheduler_params.alpha;
           obj.name  = 'Alpha fair scheduler';
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
               
                for u_ = 1:N_UE
                    attached_UEs(u_).traffic_model.check_TTI;
                end
                %% compute efficiency
                [c,user_ind] = obj.get_efficiency(N_UE,N_RB,last_received_feedbacks);          
                c = c';              
                             
               %% update average throughput
               TTI_to_read = max(current_TTI-obj.feedback_delay_TTIs-1,1); % Realistically read the ACKed throughput
               for uu = 1:N_UE
                    obj.av_throughput(uu) = obj.compute_av_throughput(uu,last_received_feedbacks,TTI_to_read);
               end
               
               %% PF scheduler
               RBs = obj.PF_scheduler(N_UE,N_RB,c,user_ind,attached_UEs);
               
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
       
       function RBs = PF_scheduler(obj,N_UE,N_RB,c,user_ind,attached_UEs)
           RB_set = ones(N_RB,1);
           RB_UEs = false(N_RB,N_UE);
           utility = zeros(size(c));
           for uu = 1:N_UE
               if ~strcmp(attached_UEs(user_ind(uu)).traffic_model.type,'fullbuffer')
                   buffer_length = attached_UEs(user_ind(uu)).traffic_model.get_buffer_length;
                   utility(:,uu) = min(c(:,uu)*12*7*2,buffer_length);
                   if strcmp(attached_UEs(user_ind(uu)).traffic_model.type,'MLaner')
                       if buffer_length <= 0 %% deactivate user if he has no traffic
                           attached_UEs(user_ind(uu)).deactivate_UE = true;
                       end
                       % reactivate user before traffic arrives, such that
                       % feedback is useful
                       next_traffic = attached_UEs(user_ind(uu)).traffic_model.generation_times - attached_UEs(user_ind(uu)).traffic_model.counter;
                       if attached_UEs(user_ind(uu)).deactivate_UE == true && min(abs(next_traffic)) < 3
                           attached_UEs(user_ind(uu)).deactivate_UE = false;
                       end
                   end
               else
                   utility(:,uu) = c(:,uu)*12*7*2;
               end
           end
           for rr = 1:N_RB
               if ~sum(sum(utility))
                   break
               end
               res = find(RB_set);
               metric = ones(N_RB,N_UE)*-Inf;
               for r_ = 1:sum(RB_set)
                   for u_ = 1:N_UE
                       metric(res(r_),u_) = real(log10(utility(res(r_),u_))-obj.alpha*log10(max((obj.av_const-1)*obj.av_throughput(user_ind(u_))+RB_UEs(:,u_).'*c(:,u_)*12*7*2,eps)));      % 12*7 equals the number of elements in a RB             
                   end
               end
               maxi = max(metric(:));
               indis = find(metric == maxi);
               ind = indis(randi(length(indis)));
               [temp_res,temp_ue] = ind2sub(size(metric),ind);
               if ~strcmp(attached_UEs(user_ind(temp_ue)).traffic_model.type,'fullbuffer')
                    attached_UEs(user_ind(temp_ue)).traffic_model.coarse_decrease(c(uu)*12*7*2);
                    bitnr = attached_UEs(user_ind(temp_ue)).traffic_model.bit_count;
                    utility(:,uu) = min(c(uu)*12*72*2,bitnr);
               end               
               RB_set(temp_res) = 0;
               RB_UEs(temp_res,temp_ue) = true;
           end
           RB_UEs = RB_UEs';
           RBs = RB_UEs(:);    
        end
       
   end
end 
