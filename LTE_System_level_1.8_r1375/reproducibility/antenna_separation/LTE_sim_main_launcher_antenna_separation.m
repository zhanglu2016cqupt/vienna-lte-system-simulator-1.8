close all force;
clc;
clear all
clear global;
clear classes;

results_folder                       = 'reproducibility/included_simulation_results/antenna_separation';
just_plot_figures_from_existing_data = false;

frequencies_vector = [800e6 2.6e9]; % Set for LTE800 and LTE2600

% Move to the simulator folder
cd ..
cd ..

for frequencies_vector_idx = 1:length(frequencies_vector)
    base_config = LTE_load_params_hex_grid_tilted;
    
    % Some changes to the base configuration
    base_config.results_folder                     = results_folder;
    base_config.compute_only_UEs_from_this_eNodeBs = [13 14 15 16 17 18 19 20 21 28 29 30 31 32 33 34 35 36 46 47 48];
    base_config.default_shown_GUI_cells            = [13 14 15 16 17 18 19 20 21 28 29 30 31 32 33 34 35 36 46 47 48];
    base_config.frequency                          = frequencies_vector(frequencies_vector_idx);
    base_config.show_network                       = 0;
    base_config.tx_mode                            = 4;
    base_config.compact_results_file               = true;
    base_config.delete_ff_trace_at_end             = true;
    base_config.UE_cache                           = false;
    base_config.UE_per_eNodeB                      = 1;
    base_config.simulation_time_tti                = 50;
    base_config.UE_cache_file                      = 'auto';
    base_config.keep_UEs_still                     = true;
    base_config.scheduler_params.av_window         = 20;
    base_config.map_resolution                     = 10;
    base_config.pregenerated_ff_file               = 'auto';
    base_config.channel_model.trace_length         = 1;
    base_config.scheduler                          = 'round robin';
    
    simulation_repeat = 20;
    types_of_sims     = 11;  % Execute the 2x2 and 4x4 simulations. Change it to 11 if you also want 4x2 simulations
    
    repeat_index_offset = 0;
    
    switch base_config.frequency
        case 2.6e9
            freq_band = 'LTE2600';
        case 800e6
            freq_band = 'LTE800';
        otherwise
            freq_band = 'LTE';
    end
    
    %% Create config structs
    sim_idx=1;
    configs = cell(types_of_sims*simulation_repeat,1);
    for simulation_repeat_idx = 1:simulation_repeat
        for types_of_sims_idx = 1:types_of_sims
            configs{sim_idx}                            = base_config;
            configs{sim_idx}.channel_trace_id           = simulation_repeat_idx;
            configs{sim_idx}.non_parallel_channel_trace = true;
            
            switch types_of_sims_idx
                
                %% 4x4 setups
                case 1
                    % RR 4x4, TX XX-pol and RX XX-pol
                    configs{sim_idx}.nTX                        = 4;
                    configs{sim_idx}.nRX                        = 4;
                    configs{sim_idx}.winner_antenna_params.TX_antenna_polarization        = [45 -45 45 -45];
                    configs{sim_idx}.winner_antenna_params.TX_antenna_position_in_lambdas = [-1 -1 1 1];
                    configs{sim_idx}.winner_antenna_params.RX_antenna_polarization        = [45 -45 45 -45];
                    configs{sim_idx}.winner_antenna_params.RX_antenna_position_in_lambdas = [-1 -1 1 1];
                    configs{sim_idx}.results_file = sprintf('%s_CLSM_%gx%g_RR_TX_XX_RX_XX_%g.mat',freq_band,configs{sim_idx}.nTX,configs{sim_idx}.nRX,simulation_repeat_idx+repeat_index_offset);
                    configs{sim_idx}.plot_title   = sprintf('%gx%g CLSM: TX XX',configs{sim_idx}.nTX,configs{sim_idx}.nRX);
                case 2
                    % RR 4x4, TX X-X-pol and RX XX-pol
                    configs{sim_idx}.nTX                        = 4;
                    configs{sim_idx}.nRX                        = 4;
                    configs{sim_idx}.winner_antenna_params.TX_antenna_polarization        = [45 -45 45 -45];
                    configs{sim_idx}.winner_antenna_params.TX_antenna_position_in_lambdas = [-5 -5 5 5];
                    configs{sim_idx}.winner_antenna_params.RX_antenna_polarization        = [45 -45 45 -45];
                    configs{sim_idx}.winner_antenna_params.RX_antenna_position_in_lambdas = [-1 -1 1 1];
                    configs{sim_idx}.results_file = sprintf('%s_CLSM_%gx%g_RR_TX_X-X_RX_XX_%g.mat',freq_band,configs{sim_idx}.nTX,configs{sim_idx}.nRX,simulation_repeat_idx+repeat_index_offset);
                    configs{sim_idx}.plot_title   = sprintf('%gx%g CLSM: TX X-X',configs{sim_idx}.nTX,configs{sim_idx}.nRX);
                case 3
                    % RR 4x4, TX ||||-pol and RX XX-pol
                    configs{sim_idx}.nTX                        = 4;
                    configs{sim_idx}.nRX                        = 4;
                    configs{sim_idx}.winner_antenna_params.TX_antenna_polarization        = [90 90 90 90];
                    configs{sim_idx}.winner_antenna_params.TX_antenna_position_in_lambdas = [-3 -1 1 3];
                    configs{sim_idx}.winner_antenna_params.RX_antenna_polarization        = [45 -45 45 -45];
                    configs{sim_idx}.winner_antenna_params.RX_antenna_position_in_lambdas = [-1 -1 1 1];
                    configs{sim_idx}.results_file = sprintf('%s_CLSM_%gx%g_RR_TX_90909090_RX_XX_%g.mat',freq_band,configs{sim_idx}.nTX,configs{sim_idx}.nRX,simulation_repeat_idx+repeat_index_offset);
                    configs{sim_idx}.plot_title   = sprintf('%gx%g CLSM: TX ||||',configs{sim_idx}.nTX,configs{sim_idx}.nRX);
                case 4
                    % RR 4x4, TX |-|-|-|-pol and RX XX-pol
                    configs{sim_idx}.nTX                        = 4;
                    configs{sim_idx}.nRX                        = 4;
                    configs{sim_idx}.winner_antenna_params.TX_antenna_polarization        = [90 90 90 90];
                    configs{sim_idx}.winner_antenna_params.TX_antenna_position_in_lambdas = [-15 -5 5 15];
                    configs{sim_idx}.winner_antenna_params.RX_antenna_polarization        = [45 -45 45 -45];
                    configs{sim_idx}.winner_antenna_params.RX_antenna_position_in_lambdas = [-1 -1 1 1];
                    configs{sim_idx}.results_file = sprintf('%s_CLSM_%gx%g_RR_TX_90-90-90-90_RX_XX_%g.mat',freq_band,configs{sim_idx}.nTX,configs{sim_idx}.nRX,simulation_repeat_idx+repeat_index_offset);
                    configs{sim_idx}.plot_title   = sprintf('%gx%g CLSM: TX |-|-|-|',configs{sim_idx}.nTX,configs{sim_idx}.nRX);
                    
                    %% 2x2 setups
                case 5
                    % RR 2x2, TX X-pol and RX X-pol
                    configs{sim_idx}.nTX                        = 2;
                    configs{sim_idx}.nRX                        = 2;
                    configs{sim_idx}.winner_antenna_params.TX_antenna_polarization        = [45 -45];
                    configs{sim_idx}.winner_antenna_params.TX_antenna_position_in_lambdas = [0 0];
                    configs{sim_idx}.winner_antenna_params.RX_antenna_polarization        = [45 -45];
                    configs{sim_idx}.winner_antenna_params.RX_antenna_position_in_lambdas = [0 0];
                    configs{sim_idx}.results_file = sprintf('%s_CLSM_%gx%g_RR_TX_X_RX_X_%g.mat',freq_band,configs{sim_idx}.nTX,configs{sim_idx}.nRX,simulation_repeat_idx+repeat_index_offset);
                    configs{sim_idx}.plot_title   = sprintf('%gx%g CLSM: TX X',configs{sim_idx}.nTX,configs{sim_idx}.nRX);
                case 6
                    % RR 4x4, TX ||-pol and RX X-pol
                    configs{sim_idx}.nTX                        = 2;
                    configs{sim_idx}.nRX                        = 2;
                    configs{sim_idx}.winner_antenna_params.TX_antenna_polarization        = [90 90];
                    configs{sim_idx}.winner_antenna_params.TX_antenna_position_in_lambdas = [-1 1];
                    configs{sim_idx}.winner_antenna_params.RX_antenna_polarization        = [45 -45 45 -45];
                    configs{sim_idx}.winner_antenna_params.RX_antenna_position_in_lambdas = [0 0];
                    configs{sim_idx}.results_file = sprintf('%s_CLSM_%gx%g_RR_TX_9090_RX_XX_%g.mat',freq_band,configs{sim_idx}.nTX,configs{sim_idx}.nRX,simulation_repeat_idx+repeat_index_offset);
                    configs{sim_idx}.plot_title   = sprintf('%gx%g CLSM: TX ||',configs{sim_idx}.nTX,configs{sim_idx}.nRX);
                case 7
                    % RR 2x2, TX |-|-pol and RX X-pol
                    configs{sim_idx}.nTX                        = 2;
                    configs{sim_idx}.nRX                        = 2;
                    configs{sim_idx}.winner_antenna_params.TX_antenna_polarization        = [90 90];
                    configs{sim_idx}.winner_antenna_params.TX_antenna_position_in_lambdas = [-5 5];
                    configs{sim_idx}.winner_antenna_params.RX_antenna_polarization        = [45 -45];
                    configs{sim_idx}.winner_antenna_params.RX_antenna_position_in_lambdas = [0 0];
                    configs{sim_idx}.results_file = sprintf('%s_CLSM_%gx%g_RR_TX_90-90_RX_XX_%g.mat',freq_band,configs{sim_idx}.nTX,configs{sim_idx}.nRX,simulation_repeat_idx+repeat_index_offset);
                    configs{sim_idx}.plot_title   = sprintf('%gx%g CLSM: TX |-|',configs{sim_idx}.nTX,configs{sim_idx}.nRX);
                    
                    %% 4x2 setups
                case 8
                    % RR 4x2, TX XX-pol and RX X-pol
                    configs{sim_idx}.nTX                        = 4;
                    configs{sim_idx}.nRX                        = 2;
                    configs{sim_idx}.winner_antenna_params.TX_antenna_polarization        = [45 -45 45 -45];
                    configs{sim_idx}.winner_antenna_params.TX_antenna_position_in_lambdas = [-1 -1 1 1];
                    configs{sim_idx}.winner_antenna_params.RX_antenna_polarization        = [45 -45];
                    configs{sim_idx}.winner_antenna_params.RX_antenna_position_in_lambdas = [0 0];
                    configs{sim_idx}.results_file = sprintf('%s_CLSM_%gx%g_RR_TX_XX_RX_X_%g.mat',freq_band,configs{sim_idx}.nTX,configs{sim_idx}.nRX,simulation_repeat_idx+repeat_index_offset);
                    configs{sim_idx}.plot_title   = sprintf('%gx%g CLSM: TX XX',configs{sim_idx}.nTX,configs{sim_idx}.nRX);
                case 9
                    % RR 4x2, TX X-X-pol and RX X-pol
                    configs{sim_idx}.nTX                        = 4;
                    configs{sim_idx}.nRX                        = 2;
                    configs{sim_idx}.winner_antenna_params.TX_antenna_polarization        = [45 -45 45 -45];
                    configs{sim_idx}.winner_antenna_params.TX_antenna_position_in_lambdas = [-5 -5 5 5];
                    configs{sim_idx}.winner_antenna_params.RX_antenna_polarization        = [45 -45];
                    configs{sim_idx}.winner_antenna_params.RX_antenna_position_in_lambdas = [0 0];
                    configs{sim_idx}.results_file = sprintf('%s_CLSM_%gx%g_RR_TX_X-X_RX_X_%g.mat',freq_band,configs{sim_idx}.nTX,configs{sim_idx}.nRX,simulation_repeat_idx+repeat_index_offset);
                    configs{sim_idx}.plot_title   = sprintf('%gx%g CLSM: TX X-X',configs{sim_idx}.nTX,configs{sim_idx}.nRX);
                case 10
                    % RR 4x2, TX ||||-pol and RX X-pol
                    configs{sim_idx}.nTX                        = 4;
                    configs{sim_idx}.nRX                        = 2;
                    configs{sim_idx}.winner_antenna_params.TX_antenna_polarization        = [90 90 90 90];
                    configs{sim_idx}.winner_antenna_params.TX_antenna_position_in_lambdas = [-3 -1 1 3];
                    configs{sim_idx}.winner_antenna_params.RX_antenna_polarization        = [45 -45];
                    configs{sim_idx}.winner_antenna_params.RX_antenna_position_in_lambdas = [0 0];
                    configs{sim_idx}.results_file = sprintf('%s_CLSM_%gx%g_RR_TX_90909090_RX_X_%g.mat',freq_band,configs{sim_idx}.nTX,configs{sim_idx}.nRX,simulation_repeat_idx+repeat_index_offset);
                    configs{sim_idx}.plot_title   = sprintf('%gx%g CLSM: TX ||||',configs{sim_idx}.nTX,configs{sim_idx}.nRX);
                case 11
                    % RR 4x2, TX |-|-|-|-pol and RX XX-pol
                    configs{sim_idx}.nTX                        = 4;
                    configs{sim_idx}.nRX                        = 2;
                    configs{sim_idx}.winner_antenna_params.TX_antenna_polarization        = [90 90 90 90];
                    configs{sim_idx}.winner_antenna_params.TX_antenna_position_in_lambdas = [-15 -5 5 15];
                    configs{sim_idx}.winner_antenna_params.RX_antenna_polarization        = [45 -45 45 -45];
                    configs{sim_idx}.winner_antenna_params.RX_antenna_position_in_lambdas = [-1 -1 1 1];
                    configs{sim_idx}.results_file = sprintf('%s_CLSM_%gx%g_RR_TX_90-90-90-90_RX_X_%g.mat',freq_band,configs{sim_idx}.nTX,configs{sim_idx}.nRX,simulation_repeat_idx+repeat_index_offset);
                    configs{sim_idx}.plot_title   = sprintf('%gx%g CLSM: TX |-|-|-|',configs{sim_idx}.nTX,configs{sim_idx}.nRX);
            end
            sim_idx= sim_idx+1;
        end
    end
    
    %% Run simulations
    
    % Run the rest in parallel. Each simulation uses a different channel trace
    if ~matlabpool('size')
        matlabpool open
    end
    parfor sim_idx=1:(types_of_sims*simulation_repeat)
        if ~just_plot_figures_from_existing_data
            LTE_sim_main(configs{sim_idx});
        end
    end
    if ~matlabpool('size')
        matlabpool close
    end
    
    %% Plot sim results
    group_sims  = {[1 2 3 4] [5 6 7] [8 9 10 11]};
    group_names = {'4x4_CLSM' '2x2_CLSM' '4x2_CLSM' };
    for group_idx = 1:length(group_sims)
        types_of_sims_to_plot = group_sims{group_idx};
        results               = cell(1,types_of_sims);
        UE_throughput         = cell(1,types_of_sims);
        UE_SINRs              = cell(1,types_of_sims);
        UE_throughput_ecdf    = cell(1,types_of_sims);
        for sim_idx=types_of_sims_to_plot
            current_results_UE_throughput = cell(1,simulation_repeat);
            current_results_UE_SINRs      = cell(1,simulation_repeat);
            for repeat_idx=1:simulation_repeat
                file_idx                          = ((sim_idx-1) + (repeat_idx-1)*types_of_sims)+1;
                fprintf('Reading %s files\n',configs{file_idx}.plot_title);
                current_results                   = utils.resultsFileReader(fullfile(results_folder,configs{file_idx}.results_file));
                current_results_UE_throughput{repeat_idx} = current_results.get_UE_average_throughput_just_finite_values;
                current_results_UE_SINRs{repeat_idx}      = current_results.get_UE_average_wideband_SINR_just_finite_values;
            end
            UE_throughput{sim_idx}      = cell2mat(current_results_UE_throughput);
            UE_SINRs{sim_idx}           = cell2mat(current_results_UE_SINRs);
            UE_throughput_ecdf{sim_idx} = utils.miscUtils.ecdf(UE_throughput{sim_idx});
        end
        
        % Overlapped UE throughput CDFs
        figure;
        hold all;
        for sim_idx=types_of_sims_to_plot
            plot(UE_throughput_ecdf{sim_idx}.x,UE_throughput_ecdf{sim_idx}.f,'DisplayName',sprintf('%s',configs{sim_idx}.plot_title));
        end
        legend('show','Location','SouthEast');
        for sim_idx=types_of_sims_to_plot
            scatter(UE_throughput_ecdf{sim_idx}.mean_x,UE_throughput_ecdf{sim_idx}.mean_f,'.k','SizeData',200);
        end
        % xlim([0 10]);
        grid on;
        title('CLSM: 20 MHz bandwidth, 5 km/h Winner II channel, round robin');
        xlabel('UE throughput (Mbit/s)');
        ylabel('UE throughput ECDF');
        
        the_title = 'throughput_ECDF';
        print( fullfile(results_folder,sprintf('%s_%s_%s',freq_band,group_names{group_idx},the_title)),'-dpng');
        print( fullfile(results_folder,sprintf('%s_%s_%s',freq_band,group_names{group_idx},the_title)),'-depsc');
        hgsave(fullfile(results_folder,sprintf('%s_%s_%s',freq_band,group_names{group_idx},the_title)));
        
        % SIRN-to-throughput mapping
        figure;
        hold all;
        for sim_idx=types_of_sims_to_plot
            scatter(UE_SINRs{sim_idx},UE_throughput{sim_idx},'.','DisplayName',sprintf('%s',configs{sim_idx}.plot_title),'SizeData',150);
        end
        grid on;
        legend('show','Location','SouthEast');
        title('Wideband SINR-to-throughput mapping');
        xlabel('UE Wideband SINR (dB)');
        ylabel('UE throughput (Mbit/s)');
        
        the_title = 'SINR-to-throughput_mapping';
        print( fullfile(results_folder,sprintf('%s_%s_%s',freq_band,group_names{group_idx},the_title)),'-dpng');
        print( fullfile(results_folder,sprintf('%s_%s_%s',freq_band,group_names{group_idx},the_title)),'-depsc');
        hgsave(fullfile(results_folder,sprintf('%s_%s_%s',freq_band,group_names{group_idx},the_title)));
    end
end

