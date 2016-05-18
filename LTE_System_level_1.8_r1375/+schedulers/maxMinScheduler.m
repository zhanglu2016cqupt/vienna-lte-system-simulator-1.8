classdef maxMinScheduler < schedulers.lteScheduler
%  A max. min. LTE scheduler (maximizes the minimum UE throughput)
% (c) Stefan Schwarz, INTHFT, 2010

   properties
       % See the lteScheduler class for a list of inherited attributes
       linprog_options % Options for the LP solver
   end

   methods
       
       % Class constructor. Just specify where to attach the scheduler
       function obj = maxMinScheduler(scheduler_params,attached_eNodeB_sector)
           % Fill in basic parameters (handled by the superclass constructor)
           obj      = obj@schedulers.lteScheduler(scheduler_params,attached_eNodeB_sector);
           obj.name = 'Max. Min. scheduler';
           obj.linprog_options = optimset('LargeScale','off','Simplex','on','Display','off');
           cd cvx
           cvx_setup
           cvx_quiet(true);
           cd ..
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
               C = zeros(N_UE,N_RB*N_UE);
               for rb = 1:N_RB
                    for uu = 1:N_UE
                        C(uu,N_UE*(rb-1)+uu) = c(uu,rb);
                    end
               end 
               
              %% Max. min. scheduler
              RBs = obj.Max_min_scheduler(N_UE,N_RB,C);
              RBs = RBs(1:end-1);
              for r_ = 1:N_RB
                  RB_tmp = RBs((r_-1)*N_UE+1:r_*N_UE);
                  ind = find(RB_tmp == 1);
                  if ~isempty(ind)
                    UE_id_list(r_) = attached_UEs(user_ind(ind)).id;
                  end
              end
              % Fill in RB grid
              RB_grid.user_allocation(:) = UE_id_list;
              % CQI assignment. TODO: implement HARQ
              obj.schedule_users_common(attached_UEs,last_received_feedbacks,current_TTI,tx_mode);
           end
       end
       
       function RBs = Max_min_scheduler(obj,N_UE,N_RB,C)
           % core scheduling function (same in LL and SL)
           A = kron(eye(N_RB),ones(1,N_UE));
           A = [A,zeros(size(A,1),1)];
           constraints = [-C,-ones(size(C,1),1)];
           constraints = [constraints;A];
%            RBs = linprog([zeros(N_UE*N_RB,1);1],constraints,[zeros(N_UE,1);ones(N_RB,1)],[],[],[zeros(N_RB*N_UE,1);-Inf],[ones(N_RB*N_UE,1);Inf],[],obj.linprog_options);
            c = [zeros(N_UE*N_RB,1);1];
            b = [zeros(N_UE,1);ones(N_RB,1)];
            max_val = [ones(N_RB*N_UE,1)];
            min_val = [zeros(N_RB*N_UE,1)];
            cvx_begin
                variable RBs(N_UE*N_RB+1);
                minimize(c'*RBs);
                subject to
                RBs(1:end-1) <= max_val;
                RBs(1:end-1) >= min_val;
                constraints*RBs <= b;
            cvx_end

            for rb = 1:N_RB % randomize UE choice if RBs is noninteger
                toss = rand;
                RB_temp = RBs((rb-1)*N_UE+1:rb*N_UE);
                temp_ind = cumsum(RB_temp) >= toss;
                RB_temp = zeros(size(RB_temp));
                RB_temp(find(temp_ind == 1,1)) = 1;
                RBs((rb-1)*N_UE+1:rb*N_UE) = RB_temp;
            end
%            RBs = IP([zeros(N_UE*N_RB,1);1],constraints,[zeros(N_UE,1);ones(N_RB,1)],[],[],[zeros(N_RB*N_UE,1);-Inf],[ones(N_RB*N_UE,1);Inf],1:24,10^-1);
       end
    end
end 
