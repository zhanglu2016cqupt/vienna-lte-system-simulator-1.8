close all force;
clc;
cd ..

simulation_type = 'tri_sector_tilted';

% Possible simulation types now:
%   - 'tri_sector'
%   - 'tri_sector_tilted', 'tri_sector_tilted_4x2', 'tri_sector_tilted_4x4'
%   - 'tri_sector_plus_femtocells'
%   - 'six_sector_tilted'
%   - 'capesso_pathlossmaps'
%   - 'omnidirectional_eNodeBs'
%   - 'tri_sector_tilted_traffic'

simSet = [1 1 1];

%% Base configuration
LTE_config = LTE_load_params(simulation_type);
% LTE_config.eNodeB_tx_power              = 46; % 46 dBm
LTE_config.bandwidth                      = 20e6;
% LTE_config.simulation_time_tti          = 100;
% LTE_config.network_source               = 'capesso';
% LTE_config.pathlosses                   = [0 10 10];
% LTE_config.scheduler                    = 'prop fair Sun'; % prop fair Sun % round robin
% LTE_config.channel_model.type           = 'TU';
% LTE_config.UE_speed                     = 5/3.6;
% LTE_config.UE_distribution              = 'constant UEs per ROI';
LTE_config.network_geometry               = 'regular_hexagonal_grid';
LTE_config.shadow_fading_type             = 'none';
LTE_config.UE_per_eNodeB = 20;
LTE_config.simulation_time_tti = 10;
% 
% % Misc options
% LTE_config.non_parallel_channel_trace   = true;
LTE_config.show_network                 = 0;
% LTE_config.channel_model.trace_length   = 1;
LTE_config.keep_UEs_still               = true;
LTE_config.compact_results_file         = true;
% LTE_config.compact_results_file         = 3;
LTE_config.delete_ff_trace_at_end       = true;
% LTE_config.UE_cache                     = false;
% LTE_config.pregenerated_ff_file         = 'auto';
% LTE_config.cache_network = false;


LTE_config.nTX     = simSet(2);
LTE_config.nRX     = simSet(3);
LTE_config.tx_mode = simSet(1);
ticIdx = tic;
output_results_file = LTE_sim_main(LTE_config);
time = toc(ticIdx);

simulation_data                   = load(output_results_file);
GUI_handles.aggregate_results_GUI = LTE_GUI_show_aggregate_results(simulation_data);
GUI_handles.positions_GUI         = LTE_GUI_show_UEs_and_cells(simulation_data,GUI_handles.aggregate_results_GUI);
