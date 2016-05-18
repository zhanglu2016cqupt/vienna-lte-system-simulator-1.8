function [ CoMP_sites ] = initialize_CoMP_sites(LTE_config, eNodeBs )
% (c) Thomas Blazek, 2014, thomas.blazek@nt.tuwien.ac.at
% This script initializes the CoMP sites. Each site has eNodeBs assigned. 
% Right now, this is just a showcase. Possible configurations are: 
% 'trivial': Each eNodeB is a seperate site, no coordination
% 'global': All eNodeBs are connected, cooperation across all eNodeBs is
% possible

% One has to differentiate between the CoMP scheduler and the LTE
% scheduler. 
lenfor= length(eNodeBs);
for j_ = 1:lenfor
    % 2 Modes are possible: If this parameter is set to false, the eNodeB
    % calls the scheduler as always. The 'CoMP' LTE_scheduler illustrates how
    % this mode can easily allocate Resource Blocks in a CoMP way. This
    % requires less work, however it has little influence on the
    % postprocessing done in schedule_users_common. 
    % If set to true, the eNodeBs call the CoMP_scheduler directly. In this
    % mode, it supersedes the LTE scheduler completely, but one has to
    % replicate all functionalities of the LTE scheduler too!
    eNodeBs(j_).CoMP_handles_pre_and_postprocessing = 0;
end

if isfield(LTE_config, 'CoMP_configuration')
    switch LTE_config.CoMP_configuration
        case 'trivial'
            for i_= 1:lenfor
                CoMP_sites(i_) = CoMP.CoMP_site(LTE_config, eNodeBs(i_), i_);
                LTE_config.CoMP_scheduler = 'none';
            end       
        case 'global'
            CoMP_sites = CoMP.CoMP_site(LTE_config, eNodeBs, 1);
            
        otherwise
            error('Not implemented');
    end
    if ~isfield(LTE_config, 'CoMP_scheduler')
        LTE_config.CoMP_scheduler = 'none';
    end
    for i_ = 1:length(CoMP_sites)
        switch LTE_config.CoMP_scheduler
            case 'none'
                
                
            case 'sequential'
                for j_ = 1:lenfor
                    % Example for the second mode: the CoMP-scheduler will
                    % call the set LTE_schedulers, while being able to
                    % modify pre and post-processing. 
                    eNodeBs(j_).CoMP_handles_pre_and_postprocessing = 1;
                end
                CoMP_sites(i_).scheduler = CoMP.schedulers.sequential(CoMP_sites);
                
            case 'round robin CoMP'
                if ~strcmp(LTE_config.scheduler, 'CoMP')
                    error('LTE_config.scheduler has to be set to "CoMP" so the CoMP scheduler is correctly fed')
                end
                 for j_ = 1:lenfor
                    eNodeBs(j_).CoMP_handles_pre_and_postprocessing = 0;
                end
                CoMP_sites(i_).scheduler = CoMP.schedulers.roundRobinCoMP(CoMP_sites);
        end
        
    end
else
    CoMP_sites = [];
end
end
