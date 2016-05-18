classdef resourceFairScheduler < schedulers.lteScheduler
% A resource fair max. throughput LTE scheduler 
% (c) Stefan Schwarz, INTHFT, 2010

   properties
       % See the lteScheduler class for a list of inherited attributes
       linprog_options % Options for the LP solver
   end

   methods
       
       % Class constructor. Just specify where to attach the scheduler
       function obj = resourceFairScheduler(scheduler_params,attached_eNodeB_sector)
           % Fill in basic parameters (handled by the superclass constructor)
           obj                 = obj@schedulers.lteScheduler(scheduler_params,attached_eNodeB_sector);
           obj.linprog_options = optimset('LargeScale','off','Simplex','on','Display','off');
           obj.name            = 'Resource fair scheduler';
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
               c = c(:); 

               %% RF scheduler
               RBs = obj.RF_scheduler(N_UE,N_RB,c);
               
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
       
       function RBs = RF_scheduler(obj,N_UE,N_RB,c)
           % core scheduling function (same in LL and SL)
           A = kron(eye(N_RB),ones(1,N_UE));
           b = ones(N_RB,1);
           numb = zeros(N_UE,1);
           if ceil(N_RB/N_UE) == floor(N_RB/N_UE) % check wheter N_UE divides N_RB
               numb = N_RB/N_UE*ones(N_UE,1);
           else
               low = floor(N_RB/N_UE);
               high = ceil(N_RB/N_UE);
               x = (N_RB-low*N_UE)/(high-low);
               numb(1:x) = high;
               numb(x+1:end) = low;
               temp = randperm(size(numb,1)); % random permutation of the elements to be fair on average
               numb = numb(temp);
           end
           for i = 1:N_UE
               temp = zeros(1,size(A,2));
               temp(i:N_UE:end) = 1;
               A = [A;temp];
               b = [b;numb(i)];
           end
           RBs = linprog(-c,[],[],A,b,zeros(N_RB*N_UE,1),ones(N_RB*N_UE,1),[],obj.linprog_options);
       end
   end
end 
