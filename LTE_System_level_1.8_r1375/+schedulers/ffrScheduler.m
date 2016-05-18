classdef ffrScheduler < schedulers.lteScheduler
% An LTE scheduler implementing Fractional Frequency Reuse (FFR).
% (c) Josep Colom Ikuno, INTHFT, 2011

   properties
       FR_scheduler  % Schedules the FR part of the bandwidth
       PR_scheduler  % Schedules the PR part of the bandwidth
       
       UE_list       % list of UEs attached to this scheduler
       
       beta_FR       % ratio of the BW dedicated to the FR part
       FR_assignment % RBs assigned to the FR band
       PR_assignment % RBs assigned to the PR band
       
       UE_assignment % Object shared by all of the FFR schedulers that has the UE-to-PR/FR zone assignment
       
       reuse_factor = 3;
       PR_band      = 1; % Default value
       
       % Used for the extreme cases
       use_FR
       use_PR
       
       % See the lteScheduler class for a list of inherited attributes
       % (which in this case will not be used)
   end

   methods
       
       % Class constructor. UE_queue size needs to be specified large
       % enough so it won't overflow
       function obj = ffrScheduler(scheduler_params,attached_eNodeB_sector)
           % Fill in basic parameters (handled by the superclass constructor)
           obj      = obj@schedulers.lteScheduler(scheduler_params,attached_eNodeB_sector);
           obj.name = 'FFR scheduler';
           
           % Set PR band
           obj.PR_band = scheduler_params.frequency_assignment(attached_eNodeB_sector.eNodeB_id);
           
           % Create the two children schedulers: This structure allows for
           % any arbitrary scheduler to be added, as well as with any
           % arbitrary number/list of parameters
           scheduler_params_fieldlist = fieldnames(scheduler_params);
           child_scheduler_list = {'FR_scheduler' 'PR_scheduler'};
           for child_scheduler_idx = 1:length(child_scheduler_list)
               child_scheduler = child_scheduler_list{child_scheduler_idx};
               % Copy common parameters
               for field_idx = 1:length(scheduler_params_fieldlist)
                   current_fieldname = scheduler_params_fieldlist{field_idx};
                   switch current_fieldname
                       case 'FR_scheduler'
                           % Do not copy
                       case 'PR_scheduler'
                           % Do not copy
                       otherwise
                           full_scheduler_params.(child_scheduler).(current_fieldname) = scheduler_params.(current_fieldname);
                   end
               end
               
               % Copy scheduler-specific parameters (will overwrite the previous ones in case of conflict)
               current_child_scheduler_params = scheduler_params.(child_scheduler);
               current_scheduler_params_fieldlist = fieldnames(current_child_scheduler_params);
               for field_idx = 1:length(current_scheduler_params_fieldlist)
                   current_fieldname = current_scheduler_params_fieldlist{field_idx};
                   full_scheduler_params.(child_scheduler).(current_fieldname) = scheduler_params.(child_scheduler).(current_fieldname);
               end
           end
           
           % FFR-related config
           RBs_FR            = round(obj.RB_grid.n_RB*scheduler_params.beta_FR);
           if round(obj.RB_grid.n_RB*scheduler_params.beta_FR) ~= (obj.RB_grid.n_RB*scheduler_params.beta_FR)
               warning('Optimum B_FR set to %g, but rounded to %g (%g RBs). Non-integer assignment of RBs is not allowed.',scheduler_params.beta_FR,RBs_FR/obj.RB_grid.n_RB,RBs_FR);
           end
           RBs_PR_full       = obj.RB_grid.n_RB-RBs_FR;
           mod_reuse         = mod(RBs_PR_full,obj.reuse_factor);
           obj.UE_assignment = scheduler_params.FFR_UE_mapping;
           obj.beta_FR       = scheduler_params.beta_FR;
           
           % Handle extreme cases
           if obj.beta_FR == 0
               obj.use_FR = false;
           else
               obj.use_FR = true;
           end
           
           if obj.beta_FR==1 && mod(obj.RB_grid.n_RB-RBs_FR,obj.reuse_factor)==0
               obj.use_PR = false;
           else
               obj.use_PR = true;
           end
           
           % Create FR and PR schedulers
           obj.FR_scheduler = schedulers.schedulerFactory.create_scheduler(full_scheduler_params.FR_scheduler.scheduler,full_scheduler_params.FR_scheduler,attached_eNodeB_sector);
           obj.PR_scheduler = schedulers.schedulerFactory.create_scheduler(full_scheduler_params.FR_scheduler.scheduler,full_scheduler_params.FR_scheduler,attached_eNodeB_sector);
           
           if obj.beta_FR>1 || obj.beta_FR<0
               error('beta_FR has to be between 0 and 1: %d found',obj.beta_FR);
           end
           
           % Always take out from the FR zone: very slight (almost unnoticeable) degradation of peak and mean but increase in edge (according to theory)
           if mod_reuse~=0
               if RBs_FR>0
                   RBs_FR = RBs_FR - (obj.reuse_factor-mod_reuse);
               else
                   % Safeguard for the case where no FR are present
                   RBs_FR = RBs_FR + mod_reuse;
               end
           end
           
           RBs_PR_full = obj.RB_grid.n_RB-RBs_FR;
           RBs_PR      = RBs_PR_full/obj.reuse_factor;
           
           obj.FR_assignment           = false(1,obj.RB_grid.n_RB);
           obj.FR_assignment(1:RBs_FR) = true;
           
           PR_offset                                         = RBs_FR+(obj.PR_band-1)*RBs_PR+1;
           obj.PR_assignment                                 = false(1,obj.RB_grid.n_RB);
           obj.PR_assignment(PR_offset:(PR_offset+RBs_PR-1)) = true;
           
           obj.FR_scheduler.fractional_BW_allocation = true;
           obj.FR_scheduler.fractional_allocation    = obj.FR_assignment;
           obj.PR_scheduler.fractional_BW_allocation = true;
           obj.PR_scheduler.fractional_allocation    = obj.PR_assignment;
           
           % Substitute their RB_grid with a "fake" version
           obj.FR_scheduler.RB_grid = utils.ffrUtils.filter_RB_grid(obj.RB_grid,obj.FR_assignment);
           obj.PR_scheduler.RB_grid = utils.ffrUtils.filter_RB_grid(obj.RB_grid,obj.PR_assignment);
       end
       
       % Add a UE to the queue. It could be done so each TTI the scheduler
       % gets a UE list from the eNodeB, but such a query is not necessary.
       % Just updating when a UE attaches or drops is sufficient.
       function add_UE(obj,UE_id)
           if isempty(obj.UE_list)
               obj.UE_list = UE_id;
           else
               % If UE is not on the list, add it
               if isempty(find(obj.UE_list==UE_id,1))
                   obj.UE_list = [obj.UE_list UE_id];
               end
               % else: it is already there, so do nothing.
           end
           [UE_is_FR UE_is_PR] = obj.UE_assignment.get_assignment(UE_id);
           if UE_is_FR
               obj.FR_scheduler.add_UE(UE_id);
           end
           if UE_is_PR
               obj.PR_scheduler.add_UE(UE_id);
           end
       end
       
       % Delete an UE_id from the queue
       function remove_UE(obj,UE_id)
           obj.UE_list = obj.UE_list(obj.UE_list~=UE_id);
           % No need to check now.
           obj.FR_scheduler.remove_UE(UE_id);
           obj.PR_scheduler.remove_UE(UE_id);
       end
       
       % Schedule the users in the given RB grid
       function schedule_users(obj,attached_UEs,last_received_feedbacks)
           % Redirect the feedbacks to their correct scheduler and do the same with the attached UEs.
           attached_UEs_id = [attached_UEs.id];
           [FR_UEs_idx PR_UEs_idx] = obj.UE_assignment.get_assignment(attached_UEs_id);
           FR_UEs = attached_UEs(FR_UEs_idx);
           PR_UEs = attached_UEs(PR_UEs_idx);
           
           [FR_feedback PR_feedback] = utils.ffrUtils.separate_feedback_FFR(last_received_feedbacks,obj.FR_assignment,obj.PR_assignment,FR_UEs_idx,PR_UEs_idx);
           
           % Schedule FR and PR zones: The lteScheduler class has been modified to support a modified writing to the eNodeBSignaling object.
           if obj.use_FR
               obj.FR_scheduler.schedule_users(FR_UEs,FR_feedback);
           end
           if obj.use_PR
               obj.PR_scheduler.schedule_users(PR_UEs,PR_feedback);
           end
           
           % Normalize the powers of the FR and scheduled parts
           obj.normalize_partial_RB_power_allocation();
           
           % Fill in UE signaling telling the UEs to calculate any wideband feedback (RI, basically), based on only the following RBs
           if obj.use_FR
               the_FR_assignment = obj.FR_assignment;
               the_FR_RB_grid    = obj.FR_scheduler.RB_grid;
           else
               the_FR_assignment = [];
               the_FR_RB_grid    = [];
           end
           if obj.use_PR
               the_PR_assignment = obj.PR_assignment;
               the_PR_RB_grid    = obj.PR_scheduler.RB_grid;
           else
               the_PR_assignment = [];
               the_PR_RB_grid    = [];
           end
           
           for u_ = 1:length(FR_UEs)
               FR_UEs(u_).eNodeB_signaling.adaptive_RI.RBs_for_feedback = the_FR_assignment;
           end
           for u_ = 1:length(PR_UEs)
               PR_UEs(u_).eNodeB_signaling.adaptive_RI.RBs_for_feedback = the_PR_assignment;
           end
           
           % Merge the RB_grid object from the FR and PR parts
           utils.ffrUtils.merge_RB_grids(obj.RB_grid,the_FR_RB_grid,the_PR_RB_grid,obj.FR_assignment,obj.PR_assignment);
       end
       
       % Overloading of the functions defines in the superclass
       function set_SINR_averager(obj,SINR_averager)
           obj.SINR_averager = SINR_averager;
           obj.FR_scheduler.SINR_averager = SINR_averager;
           obj.PR_scheduler.SINR_averager = SINR_averager;
       end
       function set_CQI_mapper(obj,CQI_mapper)
           obj.CQI_mapper = CQI_mapper;
           obj.FR_scheduler.CQI_mapper = CQI_mapper;
           obj.PR_scheduler.CQI_mapper = CQI_mapper;
       end
       function set_BLER_curves(obj,BLER_curves)
           obj.BLER_curves = BLER_curves;
           obj.FR_scheduler.BLER_curves = BLER_curves;
           obj.PR_scheduler.BLER_curves = BLER_curves;
       end
       function set_genie_UEs(obj,UEs)
           obj.genie.UEs = UEs;
           obj.FR_scheduler.genie.UEs = UEs;
           obj.PR_scheduler.genie.UEs = UEs;
       end
       function set_genie_eNodeBs(obj,eNodeBs)
           obj.genie.eNodeBs = eNodeBs;
           obj.FR_scheduler.genie.eNodeBs = eNodeBs;
           obj.PR_scheduler.genie.eNodeBs = eNodeBs;
       end
       function set_feedback_delay_TTIs(obj,feedback_channel_delay)
           obj.feedback_delay_TTIs = feedback_channel_delay;
           obj.FR_scheduler.feedback_delay_TTIs = feedback_channel_delay;
           obj.PR_scheduler.feedback_delay_TTIs = feedback_channel_delay;
       end
       function set_UE_traces(obj,the_UE_traces)
            obj.UE_traces = the_UE_traces;
            obj.FR_scheduler.set_UE_traces(the_UE_traces);
            obj.PR_scheduler.set_UE_traces(the_UE_traces);
        end
       
       function normalize_partial_RB_power_allocation(obj)
           % Power correction factors
           FR_power_ratio = obj.FR_scheduler.RB_grid.n_RB / obj.RB_grid.n_RB;
           PR_power_ratio = obj.PR_scheduler.RB_grid.n_RB / obj.RB_grid.n_RB;
           
           % Adjust FR power
           obj.FR_scheduler.RB_grid.power_allocation           = obj.FR_scheduler.RB_grid.power_allocation           * FR_power_ratio;
           obj.FR_scheduler.RB_grid.power_allocation_signaling = obj.FR_scheduler.RB_grid.power_allocation_signaling * FR_power_ratio;
           
           % Adjust PR power
           obj.PR_scheduler.RB_grid.power_allocation           = obj.PR_scheduler.RB_grid.power_allocation           * PR_power_ratio;
           obj.PR_scheduler.RB_grid.power_allocation_signaling = obj.PR_scheduler.RB_grid.power_allocation_signaling * PR_power_ratio;
       end
   end
end
