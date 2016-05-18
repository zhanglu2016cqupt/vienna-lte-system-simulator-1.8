classdef ConstrainedScheduler < schedulers.lteScheduler
% A proportional fair LTE scheduler 
% (c) Stefan Schwarz, INTHFT, 2010

   properties
       % See the lteScheduler class for a list of inherited attributes
       av_throughput % exponentially weighted throughputs
       lambda
       alpha
   end

   methods
       
       % Class constructor. Just specify where to attach the scheduler
       function obj = ConstrainedScheduler(scheduler_params,attached_eNodeB_sector)
           % Fill in basic parameters (handled by the superclass constructor)
           obj       = obj@schedulers.lteScheduler(scheduler_params,attached_eNodeB_sector);
           obj.alpha = scheduler_params.alpha;
           obj.name  = 'Constrained scheduler';
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
               
               %% update Lagrange multipliers
               for uu = 1:N_UE
                   TP = obj.av_throughput(uu);
%                    if ~isempty(UE_output(uu).rx_data_bits)
%                        for cc = 1:length(UE_output(uu).rx_data_bits)
%                              TP = TP + UE_output(uu).ACK(cc)*length(UE_output(uu).rx_data_bits{cc});
%                        end
%                        if ~isinf(obj.rate_constraints(uu))
%                             obj.lambda(uu) = max(0,obj.lambda(uu)-1/obj.av_const*(TP-obj.rate_constraints(uu)));
%                             obj.lambda_store{uu} = [obj.lambda_store{uu},obj.lambda(uu)];
%                        end
                       if strcmp(attached_UEs(uu).traffic_model.type,'voip')
                           attached_UEs(uu).lambda = max(0,attached_UEs(uu).lambda-1/obj.av_const*(attached_UEs(uu).traffic_model.delay_constraint-attached_UEs(uu).traffic_model.get_buffer_length/attached_UEs(uu).traffic_model.get_arrival_rate));
%                            obj.lambda_store{uu} = [obj.lambda_store{uu},obj.lambda(uu)];
                       end
%                    end
               end
               
               %% Constrained scheduler
               RBs = obj.C_scheduler(N_UE,N_RB,c,user_ind,attached_UEs);
               
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
       
       function RBs = C_scheduler(obj,N_UE,N_RB,c,user_ind,attached_UEs)
           % core scheduling function (same in LL and SL)
           RB_set = true(N_RB,1);
           RB_UEs = false(N_RB,N_UE);
           alpha_tmp = obj.alpha(end);
           if ~mod(obj.clock.current_TTI-1,5)
               overhead = obj.overhead_ref+obj.overhead_sync;
           else
               overhead = obj.overhead_ref;
           end
           for rr = 1:N_RB
               res = find(RB_set);
               metric = ones(N_RB,N_UE)*-Inf;
               for r_ = 1:sum(RB_set)
                   for u_ = 1:N_UE
                       if ~strcmp(attached_UEs(user_ind(u_)).traffic_model.type,'fullbuffer') 
                            bitnr = attached_UEs(user_ind(u_)).traffic_model.bit_count;
                            if  bitnr <= 0
                                metric(res(r_),u_) = -10^6;
                            else
                                metric(res(r_),u_) = min(c(res(r_),u_)*2*12*7,bitnr)*attached_UEs(user_ind(u_)).lambda;
                            end
                       else   
                            metric(res(r_),u_) = c(res(r_),u_)*2*12*7*(max(obj.av_throughput(user_ind(u_)),eps)^alpha_tmp+attached_UEs(user_ind(u_)).lambda);  
                       end
                   end
               end
               maxi = max(metric(:));
               indis = find(metric == maxi);
               ind = indis(randi(length(indis)));
               [temp_res,temp_ue] = ind2sub(size(metric),ind);
               RB_set(temp_res) = false;
               if maxi > -10^6
                   RB_UEs(temp_res,temp_ue) = true;
               else
                   RB_UEs(temp_res,temp_ue) = false;
               end
               if ~strcmp(attached_UEs(user_ind(temp_ue)).traffic_model.type,'fullbuffer')
                   if isinf(attached_UEs(user_ind(temp_ue)).traffic_model.bit_count)
                       attached_UEs(user_ind(temp_ue)).traffic_model.coarse_decrease(c(temp_res,temp_ue)*(2*12*7-overhead),1);
                   else
                       attached_UEs(user_ind(temp_ue)).traffic_model.coarse_decrease(c(temp_res,temp_ue)*(2*12*7-overhead));
                   end
               end
           end
           RB_UEs = RB_UEs';
           RBs = RB_UEs(:);
       end
   end
end 
