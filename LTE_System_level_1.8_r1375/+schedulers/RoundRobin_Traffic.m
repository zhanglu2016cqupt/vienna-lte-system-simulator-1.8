classdef RoundRobin_Traffic < schedulers.lteScheduler
% Round Robin scheduler that supports traffic models
% Martin Müller
% (c) 2011 by ITC
% www.nt.tuwien.ac.at

   properties
       % See the lteScheduler class for a list of inherited attributes
       av_throughput % exponentially weighted throughputs
       lambda
       alpha
   end

   methods
       
       % Class constructor. Just specify where to attach the scheduler
       function obj = RoundRobin_Traffic(scheduler_params,attached_eNodeB_sector)
           % Fill in basic parameters (handled by the superclass constructor)
           obj       = obj@schedulers.lteScheduler(scheduler_params,attached_eNodeB_sector);
%            obj.alpha = scheduler_params.alpha;
           obj.name  = 'round robin traffic scheduler';
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
           tx_mode           = obj.default_tx_mode;
           current_TTI       = obj.clock.current_TTI; 
           N_UE = length(attached_UEs);
           N_RB = RB_grid.n_RB;
           UE_id_list = zeros(N_RB,1);
           
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
               
               %% Round Robin traffic scheduler
               RBs = obj.RoundRobin_Traffic_scheduler(N_UE,N_RB,c,user_ind,attached_UEs);
               
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
       
       function RBs = RoundRobin_Traffic_scheduler(obj,N_UE,N_RB,c,user_ind,attached_UEs)
           % core scheduling function (same in LL and SL except for factor 2 --> 2 RBs are scheduled)
                      
           if ~mod(obj.clock.current_TTI-1,5)
               overhead = obj.overhead_ref+obj.overhead_sync;
           else
               overhead = obj.overhead_ref;
           end
           
           RBs = zeros(1,N_RB*N_UE);
           RB_UE = zeros(N_UE,N_RB);
           true_RB_UE = RB_UE;
           bits_left = zeros(N_UE,1);

           
           for ii = 1:length(user_ind)         % Check if data is available
                if strcmp(attached_UEs(user_ind(ii)).traffic_model.type,'fullbuffer')
                    bits_left(user_ind(ii)) = 1;
                else
                    bits_left(ii) = attached_UEs(user_ind(ii)).traffic_model.bit_count;
                end
                isbits = logical(bits_left);
            end
           
           for rr = 1:N_RB
               if sum(bits_left)
                    assigned_RBs = sum(RB_UE,2);
                    assigned_RBs(isbits == 0) = Inf;
                    [val,pos] = min(assigned_RBs(:,:,1));
                    RB_UE(pos,rr) = 1;
                    tmp_UE = user_ind(pos);
                    true_RB_UE(tmp_UE,rr) = 1;
                    if ~strcmp(attached_UEs(user_ind(ii)).traffic_model.type,'fullbuffer')
                        if  sum(RB_UE(pos,:)) <= 1      %coarse decrease with crc-bits subtracted (only for non-fullbuffer)
                            bits_left(pos) = attached_UEs(tmp_UE).traffic_model.coarse_decrease(c(rr,tmp_UE)*(2*12*7-overhead-24));
                        else    %crc is subtracted only once
                            bits_left(pos) = attached_UEs(tmp_UE).traffic_model.coarse_decrease(c(rr,tmp_UE)*(2*12*7-overhead));
                        end
                        isbits = logical(bits_left);
                    end
                end
            end
                           %putting the assigned_RBs in the right order, according to random order of user_ind
            RBs = reshape(true_RB_UE,[],1)';
       end          
   end
end 
