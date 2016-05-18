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
LTE_config.show_network               = 1;
LTE_config.nTX                        = 1;
LTE_config.nRX                        = 1;
LTE_config.tx_mode                    = 1;
LTE_config.nr_eNodeB_rings            = 2; % Number of eNodeB rings
LTE_config.UE_per_eNodeB              = [10 2]; % First number refers to 'macro', second to 'femto' ,i.e. [Nr_of_UEs_per_Macro Nr_of_UEs_per_Femto]   
LTE_config.simulation_time_tti        = 10;
LTE_config.map_resolution             = 5;
LTE_config.compute_only_center_users  = true; % Inclusion radius set in LTE_init_determine_eNodeBs_to_compute.m

LTE_config.compact_results_file       = false;
LTE_config.keep_UEs_still             = true;
% Femto specific
LTE_config.femtocells_config.femtos_per_cell    = 2;
LTE_config.femtocells_config.tx_power_W         = 10^(20/10)*1/1000;  
LTE_config.femtocells_config.mode               = 'CSG'; 
LTE_config.femtocells_config.macroscopic_pathloss_model_settings.wall_loss        = 20; 
LTE_config.femtocells_config.macroscopic_pathloss_model_settings.penetration_loss = LTE_config.femtocells_config.macroscopic_pathloss_model_settings.wall_loss; % Desired signal experiences penetration loss

%%
output_results_file = LTE_sim_main(LTE_config);

simulation_data                   = load(output_results_file);
GUI_handles.aggregate_results_GUI = LTE_GUI_show_aggregate_results(simulation_data);
GUI_handles.positions_GUI         = LTE_GUI_show_UEs_and_cells(simulation_data,GUI_handles.aggregate_results_GUI);
