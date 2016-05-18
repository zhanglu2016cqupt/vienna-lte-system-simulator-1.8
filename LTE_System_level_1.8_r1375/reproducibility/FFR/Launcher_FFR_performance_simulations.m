close all force;
clear
clc;
clear all
clear global;
clear classes;

cd ..
cd ..

beta_FR_vect          = 1:-0.03:0.01;  % 1=Reuse-1, 0=Reuse-3
SINR_threshold_vect   = -2:0.25:22.50; % If UE_SINR<SIRN_threshold, this UE is a PR UE (and viceversa). beta=1 and a very low is the "normal" reuse-1 LTE
antenna_configs       = {[4 4]};
scheduler_configs     = {'prop fair Sun'}; {'round robin' 'prop fair Sun'};
parallel_sim          = true;
simulations_per_loop  = 1; % Change to 100 for given simulations for the confidence interval simulations

LTE_configs           = cell(1,length(beta_FR_vect)*length(SINR_threshold_vect)*length(antenna_configs)*length(scheduler_configs));
i_                    = 1;

fprintf('Generating config files:\n');
for antenna_config_idx = 1:length(antenna_configs)
    for scheduler_config_idx = 1:length(scheduler_configs)
        % Get antenna and scheduler config
        TX_RX_antennas   = antenna_configs{antenna_config_idx};
        scheduler_to_use = scheduler_configs{scheduler_config_idx};
        
        % Load base parameters
        LTE_config                                    = LTE_load_params('tri_sector_tilted');
        LTE_config.nTX                                = TX_RX_antennas(1);
        LTE_config.nRX                                = TX_RX_antennas(2);
        LTE_config.compute_only_UEs_from_this_eNodeBs = [13 14 15 16 17 18 19 20 21 28 29 30 31 32 33 34 35 36 46 47 48];
        LTE_config.default_shown_GUI_cells            = [13 14 15 16 17 18 19 20 21 28 29 30 31 32 33 34 35 36 46 47 48];
        LTE_config.seedRandStream             = true; % Reproduces this exact simulation
        LTE_config.show_network               = 0;
        LTE_config.shadow_fading_type         = 'none';
        LTE_config.compact_results_file       = 3;
        LTE_config.delete_ff_trace_at_end     = true;
        LTE_config.UE_cache                   = false;
        LTE_config.adaptive_RI                = 2;
        LTE_config.simulation_time_tti        = 25;
        LTE_config.UE_per_eNodeB              = 30;
        LTE_config.keep_UEs_still             = true;
        LTE_config.channel_model.trace_length = 10;
        LTE_config.trace_version              = 'v1';    % 'v1' for pregenerated precoding. 'v2' for run-time-applied precoding
        
        % Generate baseline configs
        LTE_config.scheduler                  = scheduler_to_use;
        LTE_config.scheduler_params.av_window = 25;
        LTE_config.power_allocation           = 'homogeneous;'; % 'right now no power loading is implemented, so just leave it as 'homogeneous'
        
        for beta_FR_idx = 1:length(beta_FR_vect)
            beta_FR = beta_FR_vect(beta_FR_idx);
            for SINR_threshold_idx = 1:length(SINR_threshold_vect)
                for rep_idx = 1
                    
                    %% Scheduler options (FFR)
                    % LTE_config.scheduler = 'round robin';
                    LTE_config.scheduler                = 'FFR'; % 'round robin', 'best cqi', 'max min', 'max TP', 'resource fair', 'prop fair Sun'
                    % Scheduling options for the full reuse part scheduler
                    LTE_config.scheduler_params.FR_scheduler.scheduler = scheduler_to_use;
                    LTE_config.scheduler_params.FR_scheduler.av_window = 25;
                    % Scheduling options for the partial reuse part scheduler
                    LTE_config.scheduler_params.PR_scheduler.scheduler = scheduler_to_use;
                    LTE_config.scheduler_params.PR_scheduler.av_window = 25;
                    
                    % FFR-dependant parameters
                    LTE_config.FFR_params.beta_FR              = beta_FR;
                    LTE_config.FFR_params.SINR_threshold_value = SINR_threshold_vect(SINR_threshold_idx);
                    
                    %% Results folder
                    LTE_config.results_folder       = sprintf('./results/FFR/FFR_%s_%dx%d_%s',date,TX_RX_antennas(1),TX_RX_antennas(2),strrep(scheduler_to_use,' ','_'));
                    LTE_config.results_file         = sprintf('%gx%g_FFR_beta_%.2f_SINR_%.2fdB_%g.mat',TX_RX_antennas(1),TX_RX_antennas(2),LTE_config.FFR_params.beta_FR,LTE_config.FFR_params.SINR_threshold_value,rep_idx);
                    LTE_config.pregenerated_ff_file = fullfile('data_files/channel_traces',sprintf('4x4_%g.mat',rep_idx-1)); % As many channel traces as repetitions
                    
                    %% Write config file
                    LTE_configs{i_} = LTE_config;
                    fprintf('  - %g/%g\n',i_,length(beta_FR_vect)*length(SINR_threshold_vect)*length(antenna_configs)*length(scheduler_configs));
                    
                    % Increase index count
                    i_ = i_+1;
                end
            end
        end
    end
end
fprintf('\n');

if ~matlabpool('size') && parallel_sim
    matlabpool open
end

fprintf('Simulation begin:\n');
config_idxs = kron(1:length(LTE_configs),ones(1,simulations_per_loop));

% Call one simulation outside of the parfor loop just in case no cache files are present
parfor_idx = 1;
% Call simulator
i_ = config_idxs(parfor_idx);
currentLTE_config = LTE_configs{i_};
currentLTE_config.output_filename_suffix = sprintf('%s_%g',currentLTE_config.output_filename_suffix,parfor_idx);
allow_continuing_sims = true;
check_if_sim_is_done  = ~isempty(dir(fullfile(currentLTE_config.results_folder,currentLTE_config.results_file)));
if check_if_sim_is_done && allow_continuing_sims
    fprintf('Skipping simulation: "%s", FFR results found\n',currentLTE_config.results_file);
else
    generated_file = LTE_sim_main(currentLTE_config);
end

% Call the rest of the simulations
parfor parfor_idx=2:length(config_idxs)
    % Call simulator
    i_ = config_idxs(parfor_idx);
    currentLTE_config = LTE_configs{i_};
    currentLTE_config.output_filename_suffix = sprintf('%s_%g',currentLTE_config.output_filename_suffix,parfor_idx);
    check_if_sim_is_done  = ~isempty(dir(fullfile(currentLTE_config.results_folder,currentLTE_config.results_file)));
    if check_if_sim_is_done && allow_continuing_sims
        fprintf('Skipping simulation: "%s", FFR results found\n',currentLTE_config.results_file);
    else
        generated_file = LTE_sim_main(currentLTE_config);
    end
end

if matlabpool('size')
    matlabpool close
end

