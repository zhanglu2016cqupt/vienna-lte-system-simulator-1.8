close all force;
clc;
clear all
clear global;
clear classes;

results_folder                       = 'reproducibility/included_simulation_results';
just_plot_figures_from_existing_data = true;

% Move to the simulator folder
cd ..
cd ..

base_config = LTE_load_params;

% Some changes to the base configuration
if ~just_plot_figures_from_existing_data
    base_config.results_folder                     = results_folder;
    base_config.compute_only_UEs_from_this_eNodeBs = [13 14 15 16 17 18 19 20 21 28 29 30 31 32 33 34 35 36 46 47 48];
    base_config.default_shown_GUI_cells            = [13 14 15 16 17 18 19 20 21 28 29 30 31 32 33 34 35 36 46 47 48];
    base_config.show_network                       = 0;
    base_config.nTX                                = 2;
    base_config.nRX                                = 2;
    base_config.tx_mode                            = 4;
    base_config.compact_results_file               = true;
    base_config.delete_ff_trace_at_end             = true;
    base_config.UE_cache                           = true;
    base_config.UE_per_eNodeB                      = 20;
    base_config.simulation_time_tti                = 50;
    base_config.UE_cache_file                      = 'auto';
    base_config.keep_UEs_still                     = true;
    base_config.scheduler_params.av_window         = 20;
    base_config.map_resolution                     = 10;
    base_config.pregenerated_ff_file               = 'auto';
    
    % CLSM 2x2 RR throughput
    CLSM_2x2_config_RR               = base_config;
    CLSM_2x2_config_RR.scheduler     = 'round robin'; % prop fair Sun % round robin
    CLSM_2x2_config_RR.results_file  = 'CLSM_2x2_RR.mat';
    LTE_sim_main(CLSM_2x2_config_RR);
    
    % CLSM 2x2 Best CQI throughput
    CLSM_2x2_config_BESTCQI               = base_config;
    CLSM_2x2_config_BESTCQI.scheduler     = 'best cqi'; % prop fair Sun % round robin
    CLSM_2x2_config_BESTCQI.results_file  = 'CLSM_2x2_BESTCQI.mat';
    LTE_sim_main(CLSM_2x2_config_BESTCQI);
    
    % CLSM 2x2 Proportional fair throughput
    CLSM_2x2_config_PROPFAIR               = base_config;
    CLSM_2x2_config_PROPFAIR.scheduler     = 'prop fair Sun'; % prop fair Sun % round robin
    CLSM_2x2_config_PROPFAIR.results_file  = 'CLSM_2x2_PROPFAIR.mat';
    LTE_sim_main(CLSM_2x2_config_PROPFAIR);
end

zipfile = 'included_simulation_results.zip';
if exist(fullfile(results_folder,zipfile),'file') && ~exist(fullfile(results_folder,'CLSM_4x4_RR.mat'))
    fprintf('Unzipping results files\n');
    unzip(fullfile(results_folder,zipfile),results_folder);
end

RR_results                  = utils.resultsFileReader(fullfile(results_folder,'CLSM_2x2_RR.mat'));
UE_throughput_RR            = RR_results.get_UE_average_throughput_just_finite_values;
UE_throughput_RR_ecdf       = utils.miscUtils.ecdf(UE_throughput_RR);

BESTCQI_results             = utils.resultsFileReader(fullfile(results_folder,'CLSM_2x2_BESTCQI.mat'));
UE_throughput_BESTCQI       = BESTCQI_results.get_UE_average_throughput_just_finite_values;
UE_throughput_BESTCQI_ecdf  = utils.miscUtils.ecdf(UE_throughput_BESTCQI);

PROPFAIR_results            = utils.resultsFileReader(fullfile(results_folder,'CLSM_2x2_PROPFAIR.mat'));
UE_throughput_PROPFAIR      = PROPFAIR_results.get_UE_average_throughput_just_finite_values;
UE_throughput_PROPFAIR_ecdf = utils.miscUtils.ecdf(UE_throughput_PROPFAIR);

% Overlapped UE throughput CDFs
figure;
hold all;
plot(UE_throughput_RR_ecdf.x,      UE_throughput_RR_ecdf.f,      'DisplayName','2x2 CLSM: round robin scheduler');
plot(UE_throughput_BESTCQI_ecdf.x, UE_throughput_BESTCQI_ecdf.f, 'DisplayName','2x2 CLSM: best CQI scheduler');
plot(UE_throughput_PROPFAIR_ecdf.x,UE_throughput_PROPFAIR_ecdf.f,'DisplayName','2x2 CLSM: proportional fair scheduler');
legend('show','Location','SouthEast');
scatter(UE_throughput_RR_ecdf.mean_x,      UE_throughput_RR_ecdf.mean_f,      '.k','SizeData',200);
scatter(UE_throughput_BESTCQI_ecdf.mean_x, UE_throughput_BESTCQI_ecdf.mean_f, '.k','SizeData',200);
scatter(UE_throughput_PROPFAIR_ecdf.mean_x,UE_throughput_PROPFAIR_ecdf.mean_f,'.k','SizeData',200);
xlim([0 10]);
grid on;
title('2x2 CLSM: 20 MHz bandwidth, 5 km/h Winner II channel, 20 UEs/cell');
xlabel('UE throughput (Mbit/s)');
ylabel('UE throughput ECDF');

% Scheduler comparison figure
figure;
bar_data = [
    UE_throughput_RR_ecdf.mean_x   UE_throughput_BESTCQI_ecdf.mean_x   UE_throughput_PROPFAIR_ecdf.mean_x
    UE_throughput_RR_ecdf.p05      UE_throughput_BESTCQI_ecdf.p05      UE_throughput_PROPFAIR_ecdf.p05
    UE_throughput_RR_ecdf.p95      UE_throughput_BESTCQI_ecdf.p95      UE_throughput_PROPFAIR_ecdf.p95
    ];
the_axes = axes('Parent',gcf,...
    'XTickLabel',{'mean throughput','edge throughput','peak throughput'},...
    'XTick',[1 2 3]);
hold all;
bar(the_axes,bar_data);
grid on;
legend('show','Location','NorthWest',{sprintf('round robin scheduler. Fairness=%3.2f',UE_throughput_RR_ecdf.fairness) sprintf('best CQI scheduler. Fairness=%3.2f',UE_throughput_BESTCQI_ecdf.fairness) sprintf('proportional fair scheduler. Fairness=%3.2f',UE_throughput_PROPFAIR_ecdf.fairness)});
title('2x2 CLSM: 20 MHz bandwidth, 5 km/h Winner II channel, 20 UEs/cell');
ylabel('throughput (Mbit/s)');
