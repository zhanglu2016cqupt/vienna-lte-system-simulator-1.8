classdef schedulerFactory
    % This class centralizes the creation of scheduler objects. If makes it easier to check whether a scheduler is valid, etc.
    % (c) Josep Colom Ikuno, INTHFT, 2011
    
    properties
    end
    
    methods(Static)
        function check_whether_scheduler_is_defined(scheduler_type_string)
            % Checks whether the name of this scheduler exists
            switch scheduler_type_string
                case 'round robin'
                    % Correct
                case 'best cqi'
                    % Correct
                case 'proportional fair'
                    error('"%s" scheduler not supported anymore. Please use instead "prop fair Sun"',scheduler_type_string);
                case 'max min'
                    % Correct
                case 'max TP'
                    % Correct
                case 'resource fair'
                    % Correct
                case 'prop fair Sun'
                    % Correct
                case 'constrained'
                    % Correct
                case 'alpha fair'
                    % Correct
                case 'FFR'
                    % Correct
                case 'prop fair traffic'
                    %Correct
                case 'round robin traffic'
                    %Correct
                case 'round robin MU'
                    %Correct
                    case 'round robin MUs'
                    %Correct
                case 'CoMP'
                    %Correct
                case 'CarScheduler'
                    %Correct
                otherwise
                    error('%s scheduler not supported',scheduler_type_string);
            end
        end
        
        function new_scheduler = create_scheduler(scheduler_type_string,scheduler_params,eNodeB)
            % Return a new scheduler
            switch scheduler_type_string
                case 'round robin'
                    new_scheduler = schedulers.roundRobinScheduler(scheduler_params,eNodeB);
                case 'CoMP'
                    new_scheduler = schedulers.CoMP(scheduler_params,eNodeB);
                case 'CarScheduler'
                    new_scheduler = schedulers.CarScheduler(scheduler_params,eNodeB);
                case 'round robin MU'
                    new_scheduler = schedulers.MU_MIMO.roundRobinSchedulerMU(scheduler_params,eNodeB);
                case 'round robin MUs'
                    new_scheduler = schedulers.MU_MIMO.roundRobinSchedulerMUsingle(scheduler_params,eNodeB);
                case 'round robin traffic'
                    new_scheduler = schedulers.RoundRobin_Traffic(scheduler_params,eNodeB);
                case 'best cqi'
                    new_scheduler = schedulers.bestCqiScheduler(scheduler_params,eNodeB);
                case 'max min'
                    new_scheduler = schedulers.maxMinScheduler(scheduler_params,eNodeB);
                case 'max TP'
                    new_scheduler = schedulers.maxThroughputScheduler(scheduler_params,eNodeB);
                case 'resource fair'
                    new_scheduler = schedulers.resourceFairScheduler(scheduler_params,eNodeB);
                case 'prop fair Sun'
                    new_scheduler = schedulers.propFairSunScheduler(scheduler_params,eNodeB);
                case 'constrained'
                    new_scheduler = schedulers.ConstrainedScheduler(scheduler_params,eNodeB);
                case 'alpha fair'
                    new_scheduler = schedulers.alphaFairScheduler(scheduler_params,eNodeB);
                case 'prop fair traffic'                    
                    new_scheduler = schedulers.PropFair_Traffic(scheduler_params,eNodeB);                    
                case 'proportional fair'
                    if isfield(scheduler_params,'alpha') && isfield(scheduler_params,'beta')
                        new_scheduler = schedulers.proportionalFairScheduler(scheduler_params,eNodeB,scheduler_params.alpha,scheduler_params.beta);
                    else
                        new_scheduler = schedulers.proportionalFairScheduler(scheduler_params,eNodeB);
                    end
                case 'FFR'
                    new_scheduler = schedulers.ffrScheduler(scheduler_params,eNodeB);
                otherwise
                    error('Scheduler %s not defined',scheduler_type_string);
            end
        end
    end
    
end

