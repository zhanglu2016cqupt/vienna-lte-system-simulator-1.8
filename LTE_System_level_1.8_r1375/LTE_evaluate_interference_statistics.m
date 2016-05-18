%% Evaluate interference statistics of users
%  Martin Taranetz, Dec 2012

clear all; close all; clc;
%clear classes;

%% Preallocate cumulation variables
eval_UE_traces = [];
macro_UE_idxs = [];
femto_UE_idxs = [];

%% Load all result files in list
filelist     = dir('./results/2.14GHz*');
for i_ = 1:length(filelist)
    % Load Files:
    file_to_load = fullfile('./results', filelist(i_).name);
    load(file_to_load);
    fprintf('%d: ' ,i_); fprintf(filelist(i_).name); fprintf('\n');
    
    % temporal variables
    n_RB                    = LTE_config.N_RB; % Number of RBs
    simulation_time_tti     = LTE_config.simulation_time_tti;
    k_                      = simulation_time_tti; % Choose one of the TTIs for evaluation (shouldn't be first due to delay of feedback!)
    
    all_UE_traces             =  the_UE_traces;
    enabled_UE_idxs           = ~[all_UE_traces(k_:simulation_time_tti:end).UE_was_disabled]; % Check whether UE was disabled at TTI of interest
    eval_UE_traces_temp        = all_UE_traces(enabled_UE_idxs); % The UE traces the evaluate
    
    % Indices of makro- and femtocell attached users
    attached_site_idx = [eval_UE_traces_temp.attached_site];
    attached_site_idx = attached_site_idx(k_:simulation_time_tti:end); %take attached site at TTI k_
    % Find out type of attached eNodeB: (makro / femto):
    macro_UE_idxs = logical([macro_UE_idxs strcmp({eNodeBs(attached_site_idx).site_type},'macro')]);
    femto_UE_idxs = logical([femto_UE_idxs strcmp({eNodeBs(attached_site_idx).site_type},'femto')]);
    
    % Cumulation variables
    eval_UE_traces = [eval_UE_traces eval_UE_traces_temp];
end

%% Evaluation of cumulation variables
for i_ = 1:length(eval_UE_traces)
    % Calculate cumulative interference for each user
    cumulative_interference(i_) = sum(eval_UE_traces(i_).rx_power_interferers{k_}(1,:));
    % Number of interferers (Consider delay when evaluating over various ttis!)
    number_of_interferers(i_) = size(eval_UE_traces(i_).rx_power_interferers{k_},2); % Number of interferers taken into account at the UE
    % Which part of interference is made up by the strongest interferer 
    [strongest_interferer_power strongest_interferer_idx] = max(eval_UE_traces(i_).rx_power_interferers{k_}(1,:),[],2);
    strongest_interferer_site_id = eNodeBs_sectors(eval_UE_traces(i_).rx_power_interferers{k_}(2,strongest_interferer_idx)).parent_eNodeB.id;
    strongest_vs_total_interference_power(i_) = strongest_interferer_power/cumulative_interference(i_);
    % Mean number of interferers making up for 90 % of the interference
    interferers_sorted = sort(eval_UE_traces(i_).rx_power_interferers{k_}(1,:),'descend');
    sorted_power_aggregation   = 0; % Sum up power, beginning with the strongest interferer
    nr_of_strongest_interferers = 0; % Number of strongest interferers to make up for 90% of the total interference
    while sorted_power_aggregation <= 0.9*cumulative_interference(i_)
       nr_of_strongest_interferers = nr_of_strongest_interferers + 1;
       sorted_power_aggregation = sorted_power_aggregation + interferers_sorted(nr_of_strongest_interferers);
    end
    number_of_strongest_interferers(i_) = nr_of_strongest_interferers; % nr_of_strongest_interferers is the temporal variable for this user
    % Distance to attached eNodeB
    attached_eNodeB_pos          = eNodeBs(eval_UE_traces(i_).attached_site(k_)).pos;
    strongest_interferer_pos     = eNodeBs(strongest_interferer_site_id).pos;
    UE_pos                       = eval_UE_traces(i_).position(:,k_)';
    distance_to_attached_eNodeB(i_)      = sqrt((UE_pos(1)-attached_eNodeB_pos(1))^2 + (UE_pos(2)-attached_eNodeB_pos(2))^2);
    distance_to_strongest_interferer(i_) = sqrt((UE_pos(1)-strongest_interferer_pos(1))^2 + (UE_pos(2)-strongest_interferer_pos(2))^2);
end
clc;

% Print output
fprintf('Mean / min / max number of interferers: %.1f / %d / %d.\n', mean(number_of_interferers), min(number_of_interferers), max(number_of_interferers));
fprintf('Standard deviation: %.1f\n', sqrt(1/length(number_of_interferers)*sum((number_of_interferers-mean(number_of_interferers)).^2)));
fprintf('Mean strongest vs. total interference power [perc]: %.1f\n', mean(strongest_vs_total_interference_power));
fprintf('Mean / min / max number of strongest interferers: %.1f / %d / %d\n', mean(number_of_strongest_interferers), min(number_of_strongest_interferers), max(number_of_strongest_interferers));
fprintf('Standard deviation: %.1f\n', sqrt(1/length(number_of_strongest_interferers)*sum((number_of_strongest_interferers-mean(number_of_strongest_interferers)).^2)));

% Figures
figure('name', 'Number of interferers'); hist(number_of_interferers,max(number_of_interferers)); grid on;% Number of bins equal to max number of interferers
figure('name', 'Number of strongest interferers'); hist(number_of_strongest_interferers,max(number_of_strongest_interferers)); grid on;
figure('name', 'Distance to attached eNodeB'); hist(distance_to_attached_eNodeB,100); grid on;
figure('name', 'Distance to strongest interferer'); hist(distance_to_strongest_interferer,100); grid on;

figure('name', 'Number of interferers'); ksdensity(number_of_interferers); grid on;% Number of bins equal to max number of interferers
figure('name', 'Number of strongest interferers'); ksdensity(number_of_strongest_interferers); grid on;
figure('name', 'Distance to attached eNodeB'); ksdensity(distance_to_attached_eNodeB); grid on;
figure('name', 'Distance to strongest interferer'); ksdensity(distance_to_strongest_interferer); grid on;