cd ..
cd ..

%% By default just re-generates the plots based on already post-processed results
just_redo_plots = true;

%% Write here where the simulation results are
results_folder = {
    './results/FFR/FFR_29-Mar-2012_4x4_round_robin_opt_fairness/'
    './results/FFR/FFR_29-Mar-2012_4x4_round_robin_no_mean_loss/'
    './results/FFR/FFR_29-Mar-2012_4x4_prop_fair_Sun_opt_fairness/'
    './results/FFR/FFR_29-Mar-2012_4x4_prop_fair_Sun_no_mean_loss/'
    './results/FFR/FFR_29-Mar-2012_4x4_round_robin_reuse1/'
    './results/FFR/FFR_29-Mar-2012_4x4_prop_fair_Sun_reuse1/'
    };

results_filter = {
    '*_beta_0.19_SINR_12.00dB_*.mat' % RR optimum fairness
    '*_beta_0.46_SINR_7.75dB_*.mat'  % RR no mean loss
    '*_beta_0.19_SINR_12.50dB_*.mat' % Prop fair optimum fairness
    '*_beta_0.37_SINR_8.25dB_*.mat'  % Prop fair no mean loss
    '*_beta_1.00_SINR_-2.00dB_*.mat' % Round Robin reuse-1 reference
    '*_beta_1.00_SINR_-2.00dB_*.mat' % Proportional fair reuse-1 reference
    };

results_output_folder = './results/FFR/conf_int';

for results_filter_idx = 1:length(results_filter)
    if ~just_redo_plots
        filelist                = dir(fullfile(results_folder{results_filter_idx},results_filter{results_filter_idx}));
        first_file              = true;
        SINR_vect               = [];
        beta_idx_to_val_mapping = [];
        
        fprintf('Processing %g files:\n',length(filelist));
        elapsedTime = zeros(1,length(filelist));
        
        for f_=1:length(filelist)
            fprintf('  File %g/%g',f_,length(filelist));
            tic;
            the_file      = fullfile(results_folder{results_filter_idx},filelist(f_).name);
            the_file_data = load(the_file);
            
            eNodeBs_to_use                   = the_file_data.LTE_config.compute_only_UEs_from_this_eNodeBs;
            [UEs_to_use cell_sum_throughput] = utils.resultsFileReader.get_UEs_in_given_cells(eNodeBs_to_use,the_file_data.the_UE_traces);
            
            beta_FR                = the_file_data.LTE_config.FFR_params.beta_FR;
            SINR_threshold_dB      = the_file_data.LTE_config.FFR_params.SINR_threshold_value;
            eNodeB_traces_filtered = the_file_data.the_eNodeB_traces(eNodeBs_to_use);
            
            FR_UEs_filtered = the_file_data.FFR_UE_mapping.FR_assignment(:) & UEs_to_use(:);
            PR_UEs_filtered = the_file_data.FFR_UE_mapping.PR_assignment(:) & UEs_to_use(:);
            
            UE_traces_FR = the_file_data.the_UE_traces(FR_UEs_filtered);
            UE_traces_PR = the_file_data.the_UE_traces(PR_UEs_filtered);
            
            UE_avg_throughput_Mbps_FR    = [UE_traces_FR.average_throughput_Mbps];
            UE_avg_throughput_Mbps_PR    = [UE_traces_PR.average_throughput_Mbps];
            
            wideband_SINR_FR = zeros(1,length(UE_traces_FR));
            wideband_SINR_PR = zeros(1,length(UE_traces_PR));
            avg_spectral_eff_FR = zeros(1,length(UE_traces_FR));
            avg_spectral_eff_PR = zeros(1,length(UE_traces_PR));
            for u_=1:length(UE_traces_FR)
                wideband_SINR_FR(u_)    = UE_traces_FR(u_).wideband_SINR(1);
                avg_spectral_eff_FR(u_) = UE_traces_FR(u_).average_spectral_efficiency_bit_per_cu;
            end
            for u_=1:length(UE_traces_PR)
                wideband_SINR_PR(u_)    = UE_traces_PR(u_).wideband_SINR(1);
                avg_spectral_eff_PR(u_) = UE_traces_PR(u_).average_spectral_efficiency_bit_per_cu;
            end
            
            results_struct(f_).SINR_threshold_dB    = SINR_threshold_dB;
            results_struct(f_).beta_FR              = beta_FR;
            results_struct(f_).FR_UE_avg_throughput = UE_avg_throughput_Mbps_FR;
            results_struct(f_).PR_UE_avg_throughput = UE_avg_throughput_Mbps_PR;
            results_struct(f_).wideband_SINR_FR     = wideband_SINR_FR;
            results_struct(f_).wideband_SINR_PR     = wideband_SINR_PR;
            results_struct(f_).avg_spectral_eff_FR  = avg_spectral_eff_FR;
            results_struct(f_).avg_spectral_eff_PR  = avg_spectral_eff_PR;
            
            elapsedTime(f_)       = toc;
            mean_time_per_file    = mean(elapsedTime(1,f_));
            estimated_remaining_s = (length(filelist)-f_)*mean_time_per_file;
            if estimated_remaining_s>60
                estimated_remaining_m = floor(estimated_remaining_s/60);
                estimated_remaining_s = mod(estimated_remaining_s,60);
            else
                estimated_remaining_m = 0;
            end
            fprintf('. %.0fm %02.0fs remaining\n',estimated_remaining_m,estimated_remaining_s);
            
            current_throughputs          = [results_struct(f_).FR_UE_avg_throughput results_struct(f_).PR_UE_avg_throughput];
            current_throughputs_ecdf(f_) = utils.miscUtils.ecdf(current_throughputs);
        end
        
        p05_values       = [current_throughputs_ecdf.p05];
        mean_values      = [current_throughputs_ecdf.mean_x];
        p95_values       = [current_throughputs_ecdf.p95];
        fairness_values  = [current_throughputs_ecdf.fairness];
        
        CI.fairness_mean = mean(fairness_values);
        CI.mean_mean     = mean(mean_values);
        CI.p05_mean      = mean(p05_values);
        CI.p95_mean      = mean(p95_values);
        
        [CI.fairness_95pp_conf CI.fairness_95pp_conf_relative] = utils.miscUtils.confidence_interval(fairness_values,95);
        [CI.mean_95pp_conf     CI.mean_95pp_conf_relative]     = utils.miscUtils.confidence_interval(mean_values,95);
        [CI.p05_95pp_conf      CI.p05_95pp_conf_relative]      = utils.miscUtils.confidence_interval(p05_values,95);
        [CI.p95_95pp_conf      CI.p95_95pp_conf_relative]      = utils.miscUtils.confidence_interval(p95_values,95);
        
        output_folder = fullfile(results_output_folder,'results');
        mkdir(results_output_folder);
        output_file   = fullfile(results_output_folder,sprintf('%d_conf_intervals_postprocessed.mat',results_filter_idx));
        save(output_file,'CI');
    else
        load(sprintf('./reproducibility/FFR/results/confidence_intervals/%d_conf_intervals_postprocessed.mat',results_filter_idx));
    end

    fprintf('Confidence interval results for folder %s\n.',results_folder{results_filter_idx});
	CI % Plot the confidence interval data
    fprintf('\n');
end