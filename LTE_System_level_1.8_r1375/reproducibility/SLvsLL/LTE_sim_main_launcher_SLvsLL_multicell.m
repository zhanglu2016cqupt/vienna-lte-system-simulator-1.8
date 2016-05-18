close all force;
clc;
clear all
clear global;
clear classes;

simulation_type = 'LLvsSL';

cd ..
cd ..

simSet        = {[1 1 1] [2 2 2] [3 2 2] [4 2 2] [4 4 4]};
N_channels    = 25;
N_repetitions = 5;

%% Base configuration
LTE_config = LTE_load_params(simulation_type);
LTE_config.eNodeB_tx_power              = 46; % 46 dBm
LTE_config.bandwidth                    = 1.4e6;
LTE_config.simulation_time_tti          = 100;
LTE_config.network_source               = 'fixed pathloss';
LTE_config.pathlosses                   = [0 20 20 20]; % change accordingly depending on what scenario you want to reproduce
LTE_config.scheduler                    = 'prop fair Sun'; % prop fair Sun % round robin
LTE_config.channel_model.type           = 'TU';
LTE_config.UE_speed                     = 5/3.6;
LTE_config.UE_distribution              = 'constant UEs per ROI';
LTE_config.nUEs                         = 1;

% Misc options
LTE_config.non_parallel_channel_trace   = true;
LTE_config.show_network                 = 0;
LTE_config.channel_model.trace_length   = 1;
LTE_config.keep_UEs_still               = true;
LTE_config.compact_results_file         = 2;
LTE_config.compact_results_file         = 3;
LTE_config.delete_ff_trace_at_end       = true;
LTE_config.UE_cache                     = false;
LTE_config.pregenerated_ff_file         = 'auto';

try %#ok<TRYNC>
    matlabpool open
end

for s_=1:length(simSet)
    LTE_config.nTX     = simSet{s_}(2);
    LTE_config.nRX     = simSet{s_}(3);
    LTE_config.tx_mode = simSet{s_}(1);
    ticIdx = tic;
    parfor i_=1:N_channels
        for j_=1:N_repetitions
            LTE_config_current = LTE_config;
            LTE_config_current.RandStreamSeed = rand*(2^32-1);
            LTE_config_current.pregenerated_ff_file = sprintf('data_files/channel_traces/Fixed_pathloss_sims_%d_%dx%d_%d.mat',LTE_config.tx_mode,LTE_config.nTX,LTE_config.nRX,i_);
            LTE_config_current.results_file = sprintf('Fixed_pathloss_sims_mode_%d_%dx%d_%d_%d',LTE_config.tx_mode,LTE_config.nTX,LTE_config.nRX,i_,j_);
            LTE_sim_main(LTE_config_current);
        end
    end
    time{s_} = toc(ticIdx);
end

for s_=1:length(simSet)
    thr_all = [];
    for i_=1:N_channels
        for j_=1:N_repetitions
            results_file = sprintf('Fixed_pathloss_sims_mode_%d_%dx%d_%d_%d',simSet{s_}(1),simSet{s_}(2),simSet{s_}(3),i_,j_);
            simulation_data = load(fullfile(LTE_config.results_folder,[results_file '.mat']));
            thr = [simulation_data.the_UE_traces.average_throughput_Mbps];
            thr = thr(isfinite(thr));
            thr_all = [thr_all thr];
        end
    end
    sci = bootci(2000,{@mean,thr_all});
    fprintf('Mode %d, %dx%d: %.2f Mbit/s [%.2f %.2f], %.0fs\n',simSet{s_}(1),simSet{s_}(2),simSet{s_}(3),mean(thr_all),sci(1),sci(2),time{s_});
end