classdef PropFair_Traffic < schedulers.lteScheduler
% Proportional Fair scheduler that supports traffic models
% Martin Müller
% (c) 2011 by ITC
% www.nt.tuwien.ac.at

   properties
       % See the lteScheduler class for a list of inherited attributes
       av_throughput % exponentially weighted throughputs
       lambda
       alpha
       test_id
   end

   methods
       
       % Class constructor. Just specify where to attach the scheduler
       function obj = PropFair_Traffic(scheduler_params,attached_eNodeB_sector)
           % Fill in basic parameters (handled by the superclass constructor)
           obj       = obj@schedulers.lteScheduler(scheduler_params,attached_eNodeB_sector);
%            obj.alpha = scheduler_params.alpha;
           obj.name  = 'PropFair_Traffic scheduler';
           obj.test_id = attached_eNodeB_sector.id;
       end
       
       % Dummy functions required by the lteScheduler Abstract class implementation
       % Add UE (no memory, so empty)
       function add_UE(obj,UE_id)
       end
       % Delete UE (no memory, so empty)
       function remove_UE(obj,UE_id)
       end
       
       % Schedule the users in the given RB grid
%        function
%        schedule_users(obj,RB_grid,attached_UEs,last_received_feedbacks)       %RB_grid removed!!!
       function schedule_users(obj,attached_UEs,last_received_feedbacks)
           % Power allocation
           % Nothing here. Leave the default one (homogeneous)
           RB_grid = obj.RB_grid;
           RB_grid.size_bits = 0;
           
           % For now use the static tx_mode assignment
           
           RB_grid.size_bits = 0;
           tx_mode = obj.default_tx_mode;
           current_TTI = obj.clock.current_TTI;
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

               %% Proportional fair traffic scheduler
               RBs = obj.Propfair_Traffic_scheduler(N_UE,N_RB,c,user_ind,attached_UEs);
               
               for r_ = 1:N_RB
                   RB_tmp = RBs((r_-1)*N_UE+1:r_*N_UE);
                   ind = find(RB_tmp == 1);
                  if ~isempty(ind)
                    UE_id_list(r_) = attached_UEs(user_ind(ind)).id;
                  end
               end              
               RB_grid.user_allocation(:) = UE_id_list;
               % CQI assignment. TODO: implement HARQ          
%                obj.schedule_users_common(RB_grid,attached_UEs,last_received_feedbacks,current_TTI,tx_mode);   %RB_grid removed!!!
               obj.schedule_users_common(attached_UEs,last_received_feedbacks,current_TTI,tx_mode);
           end
       end
       
       function RBs = Propfair_Traffic_scheduler(obj,N_UE,N_RB,c,user_ind,attached_UEs)
           % core scheduling function (same in LL and SL except for factor 2 --> 2RBs are scheduled)
           
           if ~mod(obj.clock.current_TTI-1,5)
               overhead = obj.overhead_ref+obj.overhead_sync;
           else
               overhead = obj.overhead_ref;
           end
           alpha_temp   = 1;
           RBs          = zeros(N_RB*N_UE,1);
           bits_left    = zeros(N_UE,1);
           isbits       = false(N_UE,1);
           metric       = ones(N_RB,N_UE)*-Inf;
           RB_set     = true(N_RB,1);
           RB_UEs     = false(N_RB,N_UE);
           
           for ii = 1:N_UE         % Check if data is available
                if strcmp(attached_UEs(user_ind(ii)).traffic_model.type,'fullbuffer')
                    bits_left(user_ind(ii)) = 1;
                else
                    bits_left(user_ind(ii)) = attached_UEs(user_ind(ii)).traffic_model.bit_count;
                end
                isbits = logical(bits_left);
           end          
           
           % Precalculated values taken out from the loop (speeds up simulations)
           cleaned_c_log_matrix = log10(max(c,eps)*2*12*7);
           avgd_UE_throughputs  = (obj.av_const-1)*obj.av_throughput(user_ind);

           % Calculate metric for each RB and attached user
           for rr = 1:N_RB
               if sum(bits_left)
                   res                    = find(RB_set);
                   metric                 = -Inf(N_RB,N_UE);
                   UE_avgd_pre_metric     = -alpha_temp*log10(max(avgd_UE_throughputs+sum(RB_UEs.*c,1)*2*12*7,eps));
                   UE_avgd_pre_metric_mat = UE_avgd_pre_metric(ones(1,N_RB),:);
                   
                   metric(res(1:sum(RB_set)),:) = cleaned_c_log_matrix(res(1:sum(RB_set)),:)+UE_avgd_pre_metric_mat(res(1:sum(RB_set)),:);                  
                   metric(:,~isbits(user_ind)) = -Inf;          % if there are no bits left, set metric to -Inf
                   
                   maxi            = max(metric(:));
                   [RB_idx UE_idx] = find(metric == maxi);
                   ind             = randi(length(RB_idx));
                   
                   tmp_UE          = UE_idx(ind);
                   tmp_RB          = RB_idx(ind);
                   
                   RB_set(tmp_RB)               = false;
                   RB_UEs(tmp_RB,tmp_UE)        = true;
                   
                   % coarse decrease for UE who got the current RB and check if there are still bits left
                   if ~strcmp(attached_UEs(tmp_UE).traffic_model.type,'fullbuffer')                       
                       if sum(RB_UEs(:,tmp_UE)) <= 1       %coarse decrease with crc-bits subtracted (only for non-fullbuffer)
                           attached_UEs(tmp_UE).traffic_model.coarse_decrease(c(tmp_RB,tmp_UE)*(2*12*7-overhead-24));
                       else                                     %crc is subtracted only once
                           attached_UEs(tmp_UE).traffic_model.coarse_decrease(c(tmp_RB,tmp_UE)*(2*12*7-overhead));
                       end
                       bits_left(tmp_UE) = attached_UEs(tmp_UE).traffic_model.bit_count;
                       isbits(tmp_UE)    = logical(bits_left(tmp_UE));
                   end
               end
           end
           RB_UEs = RB_UEs';
           RBs = RB_UEs(:);           
       end
   end
end 
