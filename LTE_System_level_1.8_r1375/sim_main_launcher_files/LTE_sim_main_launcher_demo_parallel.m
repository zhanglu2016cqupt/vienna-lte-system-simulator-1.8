close all force;
clear classes
clc;
cd ..

simulation_type = 'tri_sector_tilted'; % see "LTE_load_params." for possible choices 

%% Load simulation parameters
% LTE_config = LTE_load_params(simulation_type);
LTE_config = simulation_config.MBSFN_test.apply_parameters;       
% LTE_config.keep_UEs_still     = true;
%% run the main simulation loop
ticIdx = tic;
parfor ii = 1:2
    output_results_file{ii} = LTE_sim_main(LTE_config);
end
time = toc(ticIdx);

%% start the GUIs for evaluation
simulation_data                   = load(output_results_file{1});
GUI_handles.aggregate_results_GUI = LTE_GUI_show_aggregate_results(simulation_data);
GUI_handles.positions_GUI         = LTE_GUI_show_UEs_and_cells(simulation_data,GUI_handles.aggregate_results_GUI);

% BLERs = zeros(length(simulation_data.simulation_traces.UE_traces),1);
% for k = 1:length(simulation_data.simulation_traces.UE_traces)
%     BLERs(k) = mean(simulation_data.simulation_traces.UE_traces(k).BLER(logical(simulation_data.simulation_traces.UE_traces(k).TB_size)));
% end
% mean(BLERs)

% av_count = 10;
% TP = zeros(length(simulation_data.simulation_traces.UE_traces.TB_size)/av_count,1);
% SNR_temp = zeros(size(TP));
% for nn = 1:length(TP)
%     TB_temp = simulation_data.simulation_traces.UE_traces.TB_size(:,(nn-1)*av_count+1:nn*av_count);
%     ACK_temp = simulation_data.simulation_traces.UE_traces.ACK(:,(nn-1)*av_count+1:nn*av_count);
%     TP(nn) = sum(TB_temp(ACK_temp))/av_count*1e-3;
%     SNR_temp(nn) = simulation_data.simulation_traces.UE_traces(1).SNR((nn-1)*av_count+1);
% end
% plot(SNR_temp,TP)
