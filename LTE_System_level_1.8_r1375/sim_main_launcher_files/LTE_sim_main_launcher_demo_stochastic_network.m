close all force;
clc;
cd ..
%clear all
%clear global;
%clear classes;

simulation_type = 'tri_sector_plus_femtocells';

% Possible simulation types now:
%   - 'tri_sector'
%   - 'tri_sector_tilted', 'tri_sector_tilted_4x2', 'tri_sector_tilted_4x4'
%   - 'tri_sector_plus_femtocells'
%   - 'six_sector_tilted'
%   - 'capesso_pathlossmaps'
%   - 'omnidirectional_eNodeBs'

LTE_config = LTE_load_params(simulation_type);

% Some changes to the base configuration, in case you would need/want them
LTE_config.show_network               = 0;
LTE_config.nTX                        = 1;
LTE_config.nRX                        = 1;
LTE_config.tx_mode                    = 1;
% stochastic macro network parameters
LTE_config.antenna.antenna_gain_pattern = 'omnidirectional';
LTE_config.sector_azimuths              = 0;
LTE_config.network_geometry             = 'hybrid';
LTE_config.average_eNodeB_distance      = 500; % [m]
LTE_config.network_size                 = 3; % determines height and width of region, in which eNodeBs are generated (can be rational, but > 0!)
LTE_config.inter_eNodeB_distance        = 600;
LTE_config.network_size                 = 7.5; % determines size of ROI
LTE_config.generate_strongest_interferer= false; % Generate strongest interferer at edge of exclusion region
LTE_config.eNodeB_tx_power              = 40; % Transmit power [W]
% User
LTE_config.UE_distribution              = 'predefined'; % In this setting, the position of the UEs has to be calculated before
UE_r                                    = 0.5 * LTE_config.inter_eNodeB_distance/2; % distance of UE to center
number_of_UEs                           = 2;
LTE_config.UE_positions                 = UE_r*[cos(2*pi*(1:number_of_UEs)/number_of_UEs)' sin(2*pi*(1:number_of_UEs)/number_of_UEs)'];
% Femto specific
LTE_config.add_femtocells                                                         = false;
LTE_config.femtocells_config.femtos_per_cell                                      = 2;
LTE_config.femtocells_config.tx_power_W                                           = 10^(20/10)*1/1000;  % Define as difference in [dB] to macrocell power
LTE_config.femtocells_config.mode                                                 = 'OSG'; 
LTE_config.femtocells_config.macroscopic_pathloss_model_settings.wall_loss        = 20; 
LTE_config.femtocells_config.macroscopic_pathloss_model_settings.penetration_loss = LTE_config.femtocells_config.macroscopic_pathloss_model_settings.wall_loss; % Desired signal experiences penetration loss
% General parameters
LTE_config.map_resolution             = 5;
LTE_config.shadow_fading_type         = 'claussen';
LTE_config.shadow_fading_sd           = 3; % according to 6 dB variance
LTE_config.simulation_time_tti        = 10;
LTE_config.compute_only_center_users  = true; % Inclusion radius set in LTE_init_determine_eNodeBs_to_compute.m
LTE_config.compact_results_file       = false;
LTE_config.keep_UEs_still             = true;


%%
output_results_file = LTE_sim_main(LTE_config);

simulation_data                   = load(output_results_file);
GUI_handles.aggregate_results_GUI = LTE_GUI_show_aggregate_results(simulation_data);
GUI_handles.positions_GUI         = LTE_GUI_show_UEs_and_cells(simulation_data,GUI_handles.aggregate_results_GUI);
