cd ..
cd ..

results_folder_all = {
     './results/FFR/FFR_16-Mar-2012_4x4_round_robin',...
     './results/FFR/FFR_23-May-2012_4x4_prop_fair_Sun'
    };

for results_folder_idx = 1:length(results_folder_all)
    results_folder                   = results_folder_all{results_folder_idx};
    [pathstr, results_filename, ext] = fileparts(results_folder);
    
    filelist                = dir(fullfile(results_folder,'*.mat'));
    first_file              = true;
    SINR_vect               = [];
    beta_idx_to_val_mapping = [];
    
    fprintf('Processing %g files:\n',length(filelist));
    elapsedTime = zeros(1,length(filelist));
    
    for f_=1:length(filelist)
        fprintf('  File %g/%g',f_,length(filelist));
        tic;
        the_file      = fullfile(results_folder,filelist(f_).name);
        the_file_data = load(the_file);
        
        beta_FR                 = the_file_data.LTE_config.FFR_params.beta_FR;
        SINR_threshold_dB       = the_file_data.LTE_config.FFR_params.SINR_threshold_value;
        UE_avg_throughput_Mbps  = [the_file_data.the_UE_traces.average_throughput_Mbps];
        UE_avg_spectral_eff_bcu = [the_file_data.the_UE_traces.average_spectral_efficiency_bit_per_cu];
        finite_UEs              = isfinite(UE_avg_throughput_Mbps);
        UE_avg_throughput_Mbps  = UE_avg_throughput_Mbps(finite_UEs);
        UE_avg_spectral_eff_bcu = UE_avg_spectral_eff_bcu(finite_UEs);
        
        results_struct(f_).filename                = the_file;
        results_struct(f_).SINR_threshold_dB       = SINR_threshold_dB;
        results_struct(f_).beta_FR                 = beta_FR;
        results_struct(f_).UE_avg_throughput_Mbps  = UE_avg_throughput_Mbps;
        results_struct(f_).UE_avg_spectral_eff_bcu = UE_avg_spectral_eff_bcu;
        
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
    end
    
    output_folder = fullfile(results_folder,'results');
    mkdir(output_folder);
    output_file   = fullfile(output_folder,sprintf('%s.mat',results_filename));
    save(output_file,'results_struct');
end