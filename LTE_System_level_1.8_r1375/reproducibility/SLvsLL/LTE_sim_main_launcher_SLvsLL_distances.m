close all force;
clc;
clear all
clear global;
clear classes;

cd ..
cd ..

simulation_type = 'LLvsSL';
LTE_config = LTE_load_params(simulation_type);
LTE_config.non_parallel_channel_trace = true;
output_results_file = LTE_sim_main(LTE_config);
simulation_data = load(output_results_file);

% Plot SNR as a function of distance
UE_dist = sqrt(sum([simulation_data.the_UE_traces.position].^2));
SNR_dB  = [simulation_data.the_UE_traces.SNR_dB];
figure;
plot(SNR_dB,UE_dist);
grid on;
title('SNR vs. distance to the eNodeB');
ylabel('distance [m]');
xlabel('SNR [dB]');
xlim([-10 40]);

[SNR_dB_sorted, IDXs] = sort(SNR_dB);
UE_dist_sorted = UE_dist(IDXs);
[SNR_dB_sorted, ia, ic] = unique(SNR_dB_sorted);
UE_dist_sorted = UE_dist_sorted(ia);
interp1(SNR_dB_sorted,UE_dist_sorted,-10:2:40)