close all force;
clc;
% clear all
% clear global;
% clear classes;

% Reproduces a setting (to the best of the simulator's capabilities) as that in
% @inproceedings{kusume2010system,
%   title={System level performance of downlink {MU}-{MIMO} transmission for {3GPP} {LTE}-Advanced},
%   author={Kusume, K. and Dietl, G. and Abe, T. and Taoka, H. and Nagata, S.},
%   booktitle={IEEE 71st Vehicular Technology Conference (VTC2010-Spring)},
%   month=may,
%   address = {Taipei, Taiwan},
% }

example_class = 'tri_sector';
simulation_type = 'tri_sector_tilted';

% Possible simulation types now:
%   - 'tri_sector'
%   - 'tri_sector_tilted', 'tri_sector_tilted_4x2', 'tri_sector_tilted_4x4'
%   - 'tri_sector_plus_femtocells'
%   - 'six_sector_tilted'
%   - 'capesso_pathlossmaps'
%   - 'omnidirectional_eNodeBs'
%   - 'tri_sector_tilted_traffic'

LTE_config = LTE_load_params(simulation_type);

%% If you want to modify something taking as a base the configuration file, do it here: here an example is show that changes the inter-eNodeB distances based on the LTE_load_params_hex_grid_tilted config file.

% Some changes to the base configuration, in case you would need/want them

% Simulation length
LTE_config.simulation_time_tti          = 100;

% According to the simulation scenario 1
LTE_config.eNodeB_tx_power              = 46; % 46 dBm
LTE_config.frequency                    = 2e9;
LTE_config.bandwidth                    = 10e6;
LTE_config.nTX                          = 4;
LTE_config.nRX                          = 2;
LTE_config.tx_mode                      = 4;
LTE_config.UE_per_eNodeB                = 10;
LTE_config.shadow_fading_sd             = 8;
LTE_config.scheduler                    = 'prop fair Sun'; % prop fair Sun % round robin
LTE_config.antenna.antenna_gain_pattern = 'TS 36.942';
LTE_config.TS_36942_3dB_lobe            = 70;
LTE_config.antenna.electrical_downtilt  = 15;
LTE_config.antenna.max_antenna_gain     = 14;
LTE_config.tx_height                    = 32;
LTE_config.rx_height                    = 1.5;
LTE_config.calculate_3D_pathloss        = false;
LTE_config.network_geometry             = 'regular_hexagonal_grid';
LTE_config.inter_eNodeB_distance        = 500;
LTE_config.nr_eNodeB_rings              = 2;
LTE_config.channel_model.type           = 'TU';
LTE_config.feedback_channel_delay       = 6;
LTE_config.UE_speed                     = 3/3.6;

% Misc options
LTE_config.non_parallel_channel_trace   = false;
LTE_config.show_network                 = 0;
LTE_config.channel_model.trace_length   = 10;
LTE_config.keep_UEs_still               = true;
LTE_config.compact_results_file         = 2;
LTE_config.cache_network                = true;
LTE_config.compact_results_file         = 3;
LTE_config.delete_ff_trace_at_end       = true;
LTE_config.UE_cache                     = false;
LTE_config.pregenerated_ff_file         = 'auto';
LTE_config.network_cache                = sprintf('data_files/network_caches/Kusume_%ddeg.mat',LTE_config.TS_36942_3dB_lobe);

try
    matlabpool open
end
for i_=1:2
    LTE_config_current = LTE_config;
    LTE_config_current.pregenerated_ff_file = sprintf('data_files/channel_traces/Kusume_%d.mat',i_);
    LTE_config_current.results_file = sprintf('4x2_Kusume_%d',i_);
    LTE_sim_main(LTE_config_current);
end
try
    matlabpool close
end

thr_all = [];
for i_=1:2
    results_file = sprintf('4x2_Kusume_%d',i_);
    simulation_data = load(fullfile(LTE_config.results_folder,[results_file '.mat']));
    thr = [simulation_data.the_UE_traces.average_throughput_Mbps];
    thr = thr(isfinite(thr));
    thr_all = [thr_all thr];
end

[f,x]=ecdf(thr_all);
figure;
plot(x,f,'k');
xlim([0 6])
grid on