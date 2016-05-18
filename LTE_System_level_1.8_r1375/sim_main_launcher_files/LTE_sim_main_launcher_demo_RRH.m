close all force;
clc;
cd ..
simulation_type = 'tri_sector_tilted';
LTE_config = LTE_load_params(simulation_type);

% Example results (peak/mean/edge) [Mb/s] with a Rel'8 CLSM codebook and
% 3 UEs/eNodeB (round robin scheduling, for now no others possible). Rank
% distribution shown after the throughput results.
% 
% 2x2: 37.79/15.66/3.74 (82%/18%)
% 4x4: 60.97/25.20/6.62 (44%/43%/13%/<0.2%)
%
% 2x4 with 2xRRH (1TX each) (effective 4x4), independent RRH channels (with respect to
% the eNodeB's)
%  - high separation-250m radius, 50° arch: 45.21/21.86/5.62 (49%/50%/1%)
%  - medium separation-150m radius, 60° arch: 50.99/22.91/4.65 (0.5%/47%/46.5%/6%)
%  - low separation-100m radius, 60° arch: 48.42/20.56/6.13 (50%/43%/7%)
%
% 1x4 with 3xRRH (1TX each) (effective 4x4)
%  - 250m radius, 80° arch: 20.73/12.18/4.64 (1%/96/3%)

% Main config parameters
LTE_config.nTX     = 1;
LTE_config.nRX     = 4;
LTE_config.tx_mode = 4;

% RRH parameters
LTE_config.RRHs_enabled              = true;
LTE_config.RRH.nTX                   = 1;
LTE_config.RRH.antenna_type          = 'omnidirectional';
LTE_config.RRH.distribution_type     = 'on arch';
LTE_config.RRH.distribution_params.nRRH_per_eNodeB = 3;
LTE_config.RRH.distribution_params.arch_width_deg  = 80;  % Width of the arch around the antenna azimuth
                                                          % on which the RRHs are placed
LTE_config.RRH.distribution_params.arch_radius_m   = 150; % Radius of said arch
%% Base configuration
% LTE_config.eNodeB_tx_power = 40; % 40W=46dBm
LTE_config.bandwidth             = 20e6;
LTE_config.cyclic_prefix         = 'normal';
LTE_config.network_source        = 'generated';
LTE_config.network_geometry      = 'regular_hexagonal_grid';
LTE_config.inter_eNodeB_distance = 500;
LTE_config.nr_eNodeB_rings       = 2;
LTE_config.map_resolution        = 10; % Recommended 5m for 1 ring and 10m for two rings (unless you have lots of RAM)
LTE_config.channel_model.type    = 'TU';%'winner+';
LTE_config.feedback_channel_delay= 1; % In TTIs

LTE_config.scheduler          = 'round robin'; % prop fair Sun % round robin
LTE_config.shadow_fading_type = 'none'; % 'none' or 'claussen' (by default 8dB shadow fading)

LTE_config.UE_per_eNodeB       = 5;
LTE_config.simulation_time_tti = 10;
LTE_config.UE_speed            = 5/3.6; % in m/s

% Misc options
LTE_config.cache_network              = false;
LTE_config.UE_cache                   = true;
LTE_config.show_network               = 1;
LTE_config.channel_model.trace_length = 10;
LTE_config.keep_UEs_still             = true;
LTE_config.compact_results_file       = true;
LTE_config.delete_ff_trace_at_end     = true;
LTE_config.trace_version              = 'v2';

LTE_config.compute_only_UEs_from_this_eNodeBs = [13 14 15 16 17 18 19 20 21 28 29 30 31 32 33 34 35 36 46 47 48];
LTE_config.default_shown_GUI_cells            = [13 14 15 16 17 18 19 20 21 28 29 30 31 32 33 34 35 36 46 47 48];

ticIdx              = tic;
output_results_file = LTE_sim_main(LTE_config);
time                = toc(ticIdx);

simulation_data                   = load(output_results_file);
GUI_handles.aggregate_results_GUI = LTE_GUI_show_aggregate_results(simulation_data);
GUI_handles.positions_GUI         = LTE_GUI_show_UEs_and_cells(simulation_data,GUI_handles.aggregate_results_GUI);
