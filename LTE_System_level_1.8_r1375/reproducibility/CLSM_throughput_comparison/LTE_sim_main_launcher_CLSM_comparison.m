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
base_config.results_folder                     = results_folder;
base_config.compute_only_UEs_from_this_eNodeBs = [13 14 15 16 17 18 19 20 21 28 29 30 31 32 33 34 35 36 46 47 48];
base_config.default_shown_GUI_cells            = [13 14 15 16 17 18 19 20 21 28 29 30 31 32 33 34 35 36 46 47 48];
base_config.show_network                       = 0;
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

for i_=1:4
    configs{i_} = base_config;
    
    switch i_
        case 1
            % RR 2x2
            configs{i_}.nTX          = 2;
            configs{i_}.nRX          = 2;
            configs{i_}.scheduler    = 'round robin';
            configs{i_}.results_file = 'CLSM_2x2_RR.mat';
            configs{i_}.plot_title   = '2x2 CLSM:';
        case 2
            % RR 4x2
            configs{i_}.nTX          = 4;
            configs{i_}.nRX          = 2;
            configs{i_}.scheduler    = 'round robin';
            configs{i_}.results_file = 'CLSM_4x2_RR.mat';
            configs{i_}.plot_title   = '4x2 CLSM:';
        case 3
            % RR 4x4
            configs{i_}.nTX          = 4;
            configs{i_}.nRX          = 4;
            configs{i_}.scheduler    = 'round robin';
            configs{i_}.results_file = 'CLSM_4x4_RR.mat';
            configs{i_}.plot_title   = '4x4 CLSM:';
        case 4
            % RR SISO
            configs{i_}.nTX          = 1;
            configs{i_}.nRX          = 1;
            configs{i_}.scheduler    = 'round robin';
            configs{i_}.results_file = 'SISO_1x1_RR.mat';
            configs{i_}.plot_title   = 'SISO:';
            configs{i_}.tx_mode      = 1;
    end
    
    if ~just_plot_figures_from_existing_data
        LTE_sim_main(configs{i_});
    end
end

zipfile = 'included_simulation_results.zip';
if exist(fullfile(results_folder,zipfile),'file') && ~exist(fullfile(results_folder,'CLSM_4x4_RR.mat'),'file')
    fprintf('Unzipping results files\n');
    unzip(fullfile(results_folder,zipfile),results_folder);
end

results            = cell(1,length(configs));
UE_throughput      = cell(1,length(configs));
UE_SINRs           = cell(1,length(configs));
UE_throughput_ecdf = cell(1,length(configs));
for i_=1:length(configs)
    results{i_}            = utils.resultsFileReader(fullfile(results_folder,configs{i_}.results_file));
    UE_throughput{i_}      = results{i_}.get_UE_average_throughput_just_finite_values;
    UE_SINRs{i_}           = results{i_}.get_UE_average_wideband_SINR_just_finite_values;
    UE_throughput_ecdf{i_} = utils.miscUtils.ecdf(UE_throughput{i_});
end

% Overlapped UE throughput CDFs
figure;
hold all;
for i_=1:length(configs)
    plot(UE_throughput_ecdf{i_}.x,UE_throughput_ecdf{i_}.f,'DisplayName',sprintf('%s J=%3.2f',configs{i_}.plot_title,UE_throughput_ecdf{i_}.fairness));
end
legend('show','Location','SouthEast');
for i_=1:length(configs)
    scatter(UE_throughput_ecdf{i_}.mean_x,UE_throughput_ecdf{i_}.mean_f,'.k','SizeData',200);
end
xlim([0 10]);
grid on;
title('CLSM: 20 MHz bandwidth, 5 km/h Winner II channel, round robin, 20 UEs/cell');
xlabel('UE throughput (Mbit/s)');
ylabel('UE throughput ECDF');

% SIRN-to-throughput mapping
figure;
hold all;
for i_=1:length(configs)
    scatter(UE_SINRs{i_},UE_throughput{i_},'.','DisplayName',sprintf('%s J=%3.2f',configs{i_}.plot_title,UE_throughput_ecdf{i_}.fairness));
end
grid on;
legend('show','Location','NorthWest');
title('Wideband SINR-to-throughput mapping');
xlabel('UE Wideband SINR (dB)');
ylabel('UE throughput (Mbit/s)');

