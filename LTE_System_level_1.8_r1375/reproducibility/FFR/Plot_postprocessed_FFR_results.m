function Plot_postprocessed_FFR_results
cd ..
cd ..

close all force

the_folder = './reproducibility/FFR/results/';
files = {
    'FFR_16-Mar-2012_4x4_round_robin.mat',...
    'FFR_23-May-2012_4x4_prop_fair_Sun.mat'
    };

output_plots_folder_base = './reproducibility/FFR/plots';
SINR_lims = [-2 22.5];
N_cells   = 21;

zipfile = 'results.zip';
if exist(fullfile(the_folder,zipfile),'file') && ~exist(fullfile(the_folder,files{1}))
    fprintf('Unzipping results files\n');
    unzip(fullfile(the_folder,zipfile),the_folder);
end

for file_idx = 1:length(files)
    clear results_struct_sorted2
    file = fullfile(the_folder,files{file_idx});
    [pathstr, name, ext] = fileparts(file);
    
    output_plots_folder = sprintf('%s%s%s_plots',output_plots_folder_base,filesep,name);
    if ~exist(output_plots_folder,'dir')
        mkdir(output_plots_folder);
    end
    
    FFR_data = load(file);
    
    results_struct  = FFR_data.results_struct;
    N_UEs(file_idx) = length(results_struct(1).UE_avg_throughput_Mbps);

    % Sort values
    beta_FR_values    = [results_struct.beta_FR];
    [B,IX]            = sort(beta_FR_values,'ascend');
    beta_FRs          = sort(unique(beta_FR_values),'ascend');
    N_beta_FRs        = length(beta_FRs);
    N_SINR_per_beta   = length(beta_FR_values) / N_beta_FRs;
    
    % Sort by B_FR first
    results_struct_sorted1 = results_struct(IX);
    
    % Now sort by SINR threshold
    SINR_threshold_dB = [results_struct_sorted1.SINR_threshold_dB];
    for b_idx=1:N_beta_FRs
        b_idx_offset = (b_idx-1)*N_SINR_per_beta;
        b_positions  = 1:N_SINR_per_beta;
        [B,IX]       = sort(SINR_threshold_dB(b_positions+b_idx_offset));
        results_struct_sorted2(b_positions+b_idx_offset) = results_struct_sorted1(IX+b_idx_offset);
    end
    
    UE_results.beta_FR        = zeros(N_beta_FRs,N_SINR_per_beta);
    UE_results.SINR_threshold = zeros(N_beta_FRs,N_SINR_per_beta);
    UE_results.mean           = zeros(N_beta_FRs,N_SINR_per_beta);
    UE_results.p95            = zeros(N_beta_FRs,N_SINR_per_beta);
    UE_results.p50            = zeros(N_beta_FRs,N_SINR_per_beta);
    UE_results.p05            = zeros(N_beta_FRs,N_SINR_per_beta);
    N_FR_over_N_tot    = zeros(N_beta_FRs,N_SINR_per_beta);
    scheduled_UEs      = zeros(N_beta_FRs,N_SINR_per_beta);
    beta_FR_mat        = zeros(N_beta_FRs,N_SINR_per_beta);
    SINR_threshold_mat = zeros(N_beta_FRs,N_SINR_per_beta);
    
    % figure;
    % hold all
    for beta_idx=1:N_beta_FRs
        for SINR_idx=1:N_SINR_per_beta
            
            % Continue processing the results
            results_idx         = (beta_idx-1)*N_SINR_per_beta + SINR_idx;
            current_results     = results_struct_sorted2(results_idx);
            current_throughputs = [current_results.UE_avg_throughput_Mbps];
            
            zero_throughput_UEs      = sum(current_throughputs==0);
            current_throughputs_ecdf = utils.miscUtils.ecdf(current_throughputs);
            
            % plot(current_throughputs_ecdf.x,current_throughputs_ecdf.f,'DisplayName',sprintf('%g dB',current_results.SINR_switching_point));
            
            beta_FR_mat(beta_idx,SINR_idx)        = current_results.beta_FR;
            SINR_threshold_mat(beta_idx,SINR_idx) = current_results.SINR_threshold_dB;
            
            UE_results.beta_FR(beta_idx,SINR_idx)        = current_results.beta_FR;
            UE_results.SINR_threshold(beta_idx,SINR_idx) = current_results.SINR_threshold_dB;
            UE_results.mean(beta_idx,SINR_idx)           = current_throughputs_ecdf.mean_x;
            UE_results.p95(beta_idx,SINR_idx)            = current_throughputs_ecdf.p95;
            UE_results.p50(beta_idx,SINR_idx)            = current_throughputs_ecdf.p50;
            UE_results.p05(beta_idx,SINR_idx)            = current_throughputs_ecdf.p05;
            
            UE_results.fairness(beta_idx,SINR_idx)           = sum(current_throughputs).^2 ./ (length(current_throughputs)*sum(current_throughputs.^2));
            scheduled_UEs(beta_idx,SINR_idx)      = length(current_throughputs)-zero_throughput_UEs;
        end
    end
    % legend('show','Location','Best');
    % grid on;
    % xlim([0 15]);
    
    R1_values.mean     = UE_results.mean(end,1);
    R1_values.p95      = UE_results.p95(end,1);
    R1_values.p50      = UE_results.p50(end,1);
    R1_values.p05      = UE_results.p05(end,1);
    R1_values.fairness = UE_results.fairness(end,1);
    
    UE_results.mean_diff_pp = UE_results.mean/R1_values.mean*100-100;
    UE_results.p95_diff_pp  = UE_results.p95/R1_values.p95*100-100;
    UE_results.p50_diff_pp  = UE_results.p50/R1_values.p50*100-100;
    UE_results.p05_diff_pp  = UE_results.p05/R1_values.p05*100-100;
    
    % Absolute throughput values
    plot_surf(beta_FR_mat,SINR_threshold_mat,UE_results.mean,'Mean throughput [Mbps]',SINR_lims,[],output_plots_folder);
    plot_surf(beta_FR_mat,SINR_threshold_mat,UE_results.p95,'Peak throughput [Mbps]',SINR_lims,[],output_plots_folder);
    plot_surf(beta_FR_mat,SINR_threshold_mat,UE_results.p05,'Edge throughput [Mbps]',SINR_lims,[],output_plots_folder);
    
    % Relative throughput values
    plot_surf(beta_FR_mat,SINR_threshold_mat,UE_results.mean_diff_pp, 'Mean throughput gain (%)',SINR_lims,[],output_plots_folder);
    plot_surf(beta_FR_mat,SINR_threshold_mat,UE_results.p95_diff_pp,  'Peak throughput gain (%)',SINR_lims,[],output_plots_folder);
    plot_surf(beta_FR_mat,SINR_threshold_mat,UE_results.p05_diff_pp,  'Edge throughput gain (%)',SINR_lims,[],output_plots_folder);
    
    % True/false relative gain zones
    UE_mean_throughput_gain_zone = (UE_results.mean-R1_values.mean) >= 0;
    UE_p95_throughput_gain_zone  = (UE_results.p95-R1_values.p95)   >= 0;
    UE_p50_throughput_gain_zone  = (UE_results.p50-R1_values.p50)   >= 0;
    UE_p05_throughput_gain_zone  = (UE_results.p05-R1_values.p05)   >= 0;
    fairness_gain_zone           = (UE_results.fairness>=R1_values.fairness);
    UE_mean_edge_throughput_gain_zone = UE_mean_throughput_gain_zone&UE_p05_throughput_gain_zone;
    
    %% Search for a specific value
    % search_values = [3.3 2 5.1]; % [mean edge peak]
    % search_zone = (UE_results.mean<(search_values(1)+0.1))&(UE_results.mean>(search_values(1)-0.1))&(UE_results.p05<(search_values(2)+0.1))&(UE_results.p05>(search_values(2)-0.1))&(UE_results.p95<(search_values(3)+0.1))&(UE_results.p95>(search_values(3)-0.1));
    % plot_surf(beta_FR_mat,SINR_threshold_mat,UE_results.p05_diff_pp, 'Positive Edge throughput search',SINR_lims,search_zone,output_plots_folder);
    
    %% Plot throughput gain plots (%) > 0%
    plot_surf(beta_FR_mat,SINR_threshold_mat,UE_results.p95_diff_pp, 'Positive Peak throughput gain',SINR_lims,UE_p95_throughput_gain_zone,output_plots_folder);
    plot_surf(beta_FR_mat,SINR_threshold_mat,UE_results.mean_diff_pp,'Positive Mean throughput gain',SINR_lims,UE_mean_throughput_gain_zone,output_plots_folder);
    plot_surf(beta_FR_mat,SINR_threshold_mat,UE_results.p05_diff_pp, 'Positive Edge throughput gain',SINR_lims,UE_p05_throughput_gain_zone,output_plots_folder);
    
    %% Plot throughput gain plots (%) > 0% everywhere
    plot_surf(beta_FR_mat,SINR_threshold_mat,UE_results.p95_diff_pp, 'Positive Peak throughput gain - gain zone',SINR_lims,UE_mean_edge_throughput_gain_zone,output_plots_folder);
    plot_surf(beta_FR_mat,SINR_threshold_mat,UE_results.mean_diff_pp,'Positive Mean throughput gain - gain zone', SINR_lims,UE_mean_edge_throughput_gain_zone,output_plots_folder);
    plot_surf(beta_FR_mat,SINR_threshold_mat,UE_results.p05_diff_pp, 'Positive Edge throughput gain - gain zone', SINR_lims,UE_mean_edge_throughput_gain_zone,output_plots_folder);
    
    %% Plot fairness-related plots
    plot_surf(beta_FR_mat,SINR_threshold_mat,UE_results.fairness,'Fairness',SINR_lims,[],output_plots_folder);
    plot_surf(beta_FR_mat,SINR_threshold_mat,UE_results.fairness,'Fairness gain zone',SINR_lims,fairness_gain_zone,output_plots_folder);
    plot_surf(beta_FR_mat,SINR_threshold_mat,UE_results.p95_diff_pp,'Peak throughput gain (%) in fairness gain zone', SINR_lims,fairness_gain_zone,output_plots_folder);
    plot_surf(beta_FR_mat,SINR_threshold_mat,UE_results.mean_diff_pp,'Mean throughput gain (%) in fairness gain zone',SINR_lims,fairness_gain_zone,output_plots_folder);
    plot_surf(beta_FR_mat,SINR_threshold_mat,UE_results.p05_diff_pp,'Edge throughput gain (%) in fairness gain zone', SINR_lims,fairness_gain_zone,output_plots_folder);
    
    %% Print lists
    % Fairness gain zone
    envelope(file_idx) = print_point_list(fairness_gain_zone,               UE_results,'mean',output_plots_folder,'fairness_gain',true);
    % Mean and edge throughput gain
    print_point_list(UE_mean_edge_throughput_gain_zone,UE_results,'mean',output_plots_folder,'throughput_mean_edge_gain',false);
    
    %% Plot fairness gain area
    plot_surf(beta_FR_mat,SINR_threshold_mat,fairness_gain_zone,'Fairness',SINR_lims,[],output_plots_folder);
    
    %% Additional plots
    % plot_surf(beta_FR_mat,SINR_threshold_mat,N_FR_over_N_tot,'#FR UEs vs. total',SINR_lims,[],output_plots_folder);
end

%% Aggregate plots
fairness_figure = figure;
fairness_axes   = axes('Parent',fairness_figure);
hold(fairness_axes,'all');
colors     = {'red' 'blue' 'green'};
markers    = {'+' '.'};
markerSize = [8 11];
min_max = [
    Inf -Inf
    Inf -Inf
    ];
for i_=1:length(envelope)
    switch i_
        case 1
            the_name = 'round robin scheduling';
        case 2
            the_name = 'proportional fair scheduling';
    end
    min_max = [
        min([min_max(1,1);envelope(i_).mean_throughput(:)]) max([min_max(1,2);envelope(i_).mean_throughput(:)])
        min([min_max(2,1);envelope(i_).fairness(:)])        max([min_max(2,2);envelope(i_).fairness(:)])
        ];
    plot(fairness_axes,envelope(i_).mean_throughput,envelope(i_).fairness,':','Marker',markers{i_},'MarkerSize',markerSize(i_),'DisplayName',the_name,'Color',colors{i_});
end
legend(fairness_axes,'show','Location','SouthWest');
grid(fairness_axes,'on');
xlim(fairness_axes,min_max(1,:));
ylim(fairness_axes,min_max(2,:));
xlabel(fairness_axes,'Mean throughput [Mbit/s]');
ylabel(fairness_axes,'Fairness');
title(fairness_axes,'Fairness-to-mean-throughput trade-off');
print( fullfile(output_plots_folder_base,'fairness_aggregate.png'),'-dpng');
print( fullfile(output_plots_folder_base,'fairness_aggregate.eps'),'-depsc');
hgsave(fullfile(output_plots_folder_base,'fairness_aggregate.fig'));

%% Aggregate plots (edge)
fairness_figure = figure;
fairness_axes   = axes('Parent',fairness_figure);
hold(fairness_axes,'all');
min_max = [
    Inf -Inf
    Inf -Inf
    ];
for i_=1:length(envelope)
    switch i_
        case 1
            the_name = 'round robin scheduling';
        case 2
            the_name = 'proportional fair scheduling';
    end
    min_max = [
        min([min_max(1,1);envelope(i_).mean_throughput(:)]) max([min_max(1,2);envelope(i_).mean_throughput(:)])
        min([min_max(2,1);envelope(i_).edge_throughput(:)])        max([min_max(2,2);envelope(i_).edge_throughput(:)])
        ];
    plot(fairness_axes,envelope(i_).mean_throughput,envelope(i_).edge_throughput,':','Marker',markers{i_},'MarkerSize',markerSize(i_),'DisplayName',the_name,'Color',colors{i_});
end
legend(fairness_axes,'show','Location','SouthWest');
grid(fairness_axes,'on');
xlim(fairness_axes,min_max(1,:));
ylim(fairness_axes,min_max(2,:));
xlabel(fairness_axes,'Mean throughput [Mbit/s]');
ylabel(fairness_axes,'Edge throughput [Mbit/s]');
title(fairness_axes,'Edge throughput vs. mean throughput');
print( fullfile(output_plots_folder_base,'mean_edge_aggregate.png'),'-dpng');
print( fullfile(output_plots_folder_base,'mean_edge_aggregate.eps'),'-depsc');
hgsave(fullfile(output_plots_folder_base,'mean_edge_aggregate.fig'));

%% Aggregate plots (peak)
fairness_figure = figure;
fairness_axes   = axes('Parent',fairness_figure);
hold(fairness_axes,'all');
min_max = [
    Inf -Inf
    Inf -Inf
    ];
for i_=1:length(envelope)
    switch i_
        case 1
            the_name = 'round robin scheduling';
        case 2
            the_name = 'proportional fair scheduling';
    end
    min_max = [
        min([min_max(1,1);envelope(i_).mean_throughput(:)]) max([min_max(1,2);envelope(i_).mean_throughput(:)])
        min([min_max(2,1);envelope(i_).peak_throughput(:)])        max([min_max(2,2);envelope(i_).peak_throughput(:)])
        ];
    plot(fairness_axes,envelope(i_).mean_throughput,envelope(i_).peak_throughput,':','Marker',markers{i_},'MarkerSize',markerSize(i_),'DisplayName',the_name,'Color',colors{i_});
end
legend(fairness_axes,'show','Location','SouthWest');
grid(fairness_axes,'on');
xlim(fairness_axes,min_max(1,:));
ylim(fairness_axes,min_max(2,:));
xlabel(fairness_axes,'Mean throughput [Mbit/s]');
ylabel(fairness_axes,'Peak throughput [Mbit/s]');
title(fairness_axes,'Peak throughput vs. mean throughput');
print( fullfile(output_plots_folder_base,'mean_peak_aggregate.png'),'-dpng');
print( fullfile(output_plots_folder_base,'mean_peak_aggregate.eps'),'-depsc');
hgsave(fullfile(output_plots_folder_base,'mean_peak_aggregate.fig'));

%% Change this if you want to plot other reference files
file_idx_to_plot = 1;

bad_fairness_example         = './reproducibility/FFR/results/references/bad_fairness_2.14GHz_freq_20fMHz_bw_winner+_5.0Kmph_50TTIs_20120219_170911_lab01_FFR_4x4CLSM_precomputed_precoding_r520_beta_0.31_SINR_18.00dB_2166.mat';
proportional_fair_reference  = './reproducibility/FFR/results/references/PF_reuse1_2.14GHz_freq_20fMHz_bw_winner+_5.0Kmph_50TTIs_20120315_181841_lab01_FFR_4x4CLSM_precomputed_precoding_r520_beta_1.00_SINR_-2.00dB_817.mat';
round_robin_reference        = './reproducibility/FFR/results/references/RR_reuse1_2.14GHz_freq_20fMHz_bw_winner+_5.0Kmph_50TTIs_20120315_144602_lab01_FFR_4x4CLSM_precomputed_precoding_r520_beta_1.00_SINR_-2.00dB_545.mat';
round_robin_optimum_fairness = './reproducibility/FFR/results/references/RR_opt_fairness_2.14GHz_freq_20fMHz_bw_winner+_5.0Kmph_50TTIs_20120217_131816_lab01_FFR_4x4CLSM_precomputed_precoding_r520_beta_0.19_SINR_12.00dB_2506.mat';
prop_fair_optimum_fairness   = './reproducibility/FFR/results/references/PF_opt_fairness_2.14GHz_freq_20fMHz_bw_winner+_5.0Kmph_50TTIs_20120218_222929_lab01_FFR_4x4CLSM_precomputed_precoding_r520_beta_0.19_SINR_12.50dB_5602.mat';

files = {
    bad_fairness_example
    proportional_fair_reference
    round_robin_reference
    round_robin_optimum_fairness
    prop_fair_optimum_fairness
    };

simulation_data                   = load(files{file_idx_to_plot});
GUI_handles.aggregate_results_GUI = LTE_GUI_show_aggregate_results(simulation_data);
GUI_handles.positions_GUI         = LTE_GUI_show_UEs_and_cells(simulation_data,GUI_handles.aggregate_results_GUI);

end

function plot_surf(Y,X,surf_data,the_title,X_lims,cut_area,output_plots_folder)

if ~isnumeric(surf_data)
    gray_colormap = true;
    surf_data     = double(surf_data);
else
    gray_colormap = false;
end

if isempty(cut_area)
    surf_data = surf_data;
else
    surf_data_in_zone    = surf_data(cut_area);
    surf_data_lims       = [min(surf_data_in_zone) max(surf_data_in_zone)];
    surf_data(~cut_area) = surf_data_lims(1);
end

a_figure = figure;
the_axes = axes('Parent',a_figure);
surf(X,Y,surf_data,'LineStyle','none');
% contour(X,Y,surf_data,20,'ShowText','on');
view(the_axes,[-90 90]);
if ~isempty(X_lims)
    xlim(the_axes,X_lims);
else
    xlim(the_axes,[min(X(:)) max(X(:))]);
end
ylim(the_axes,[min(Y(:)) max(Y(:))]);
xlabel(the_axes,'SINR threshold [dB]');
ylabel(the_axes,'\beta_{FR}');
title(the_axes,the_title);
colorbar;

if gray_colormap
    colormap gray
else
    colormap jet
end

print( fullfile(output_plots_folder,sprintf('%d_%s.png',gcf,the_title)),'-dpng');
print( fullfile(output_plots_folder,sprintf('%d_%s.eps',gcf,the_title)),'-depsc');
hgsave(fullfile(output_plots_folder,sprintf('%d_%s.fig',gcf,the_title)));
end

function envelope = print_point_list(points_idxs,UE_results,order,output_plots_folder,the_title,plot_results)

switch order
    case 'mean'
        order_idx = 1;
    case 'peak'
        order_idx = 3;
    case 'edge'
        order_idx = 2;
    otherwise
        error('Order not supported');
end
% columns:
% mean throughput diffference (%), edge (5%) throughput diffference (%), peak (95%) throughput diffference (%), beta_FR, SINR threshold (dB), fairness

if isempty(points_idxs)
    points_idxs = true(size(UE_results.mean_diff_pp));
end

matrix_to_sort  = [
    UE_results.mean_diff_pp(points_idxs),...
    UE_results.p05_diff_pp(points_idxs),...
    UE_results.p95_diff_pp(points_idxs),...
    UE_results.p50_diff_pp(points_idxs),...
    UE_results.beta_FR(points_idxs),...
    UE_results.SINR_threshold(points_idxs),...
    UE_results.fairness(points_idxs),...
    UE_results.mean(points_idxs),...
    UE_results.p05(points_idxs),...
    UE_results.p95(points_idxs)
    ];

all_gain_points = sortrows(matrix_to_sort,order_idx);

sorted_throughput_gain = all_gain_points(:,1);
sorted_throughput      = all_gain_points(:,8);
sorted_edge_gain       = all_gain_points(:,2);
sorted_edge            = all_gain_points(:,9);
sorted_peak_gain       = all_gain_points(:,3);
sorted_peak            = all_gain_points(:,10);
sorted_median_gain     = all_gain_points(:,4);
sorted_fairness        = all_gain_points(:,7);

fairness_ylims  = [min(sorted_fairness) max(sorted_fairness)];
edge_gain_ylims = [min(sorted_edge_gain) max(sorted_edge_gain)];
peak_gain_ylims = [min(sorted_peak_gain) max(sorted_peak_gain)];
edge_ylims      = [min(sorted_edge) max(sorted_edge)];
peak_ylims      = [min(sorted_peak) max(sorted_peak)];

gain_xlims        = [min(sorted_throughput_gain) max(sorted_throughput_gain)];
throughputs_xlims = [min(sorted_throughput)      max(sorted_throughput)];

if plot_results
    % Plots with the gain axis
    y_data_binned_max_idx = scatterplot_FFR(sorted_throughput_gain,sorted_fairness,   'Mean throughput gain (%)','fairness','fairness vs. throughput gain',output_plots_folder,gain_xlims,fairness_ylims);
    scatterplot_FFR(sorted_throughput_gain,sorted_edge_gain,  'Mean throughput gain (%)','Edge throughput gain (%)',sprintf('edge throughput vs. throughput gain, %s',the_title),output_plots_folder,gain_xlims,edge_gain_ylims,y_data_binned_max_idx);
    scatterplot_FFR(sorted_throughput_gain,sorted_peak_gain,  'Mean throughput gain (%)','Peak throughput gain (%)',sprintf('peak throughput vs. throughput gain, %s',the_title),output_plots_folder,gain_xlims,peak_gain_ylims,y_data_binned_max_idx);
    
    % Plots with the throughput axis
    scatterplot_FFR(sorted_throughput,sorted_fairness,   'Mean throughput (Mbps)','fairness','fairness vs. throughput',output_plots_folder,throughputs_xlims,fairness_ylims,y_data_binned_max_idx);
    scatterplot_FFR(sorted_throughput,sorted_edge,       'Mean throughput (Mbps)','Edge throughput (Mbps)',sprintf('edge throughput vs. throughput, %s',the_title),output_plots_folder,throughputs_xlims,edge_ylims,y_data_binned_max_idx);
    scatterplot_FFR(sorted_throughput,sorted_peak,       'Mean throughput (Mbps)','Peak throughput (Mbps)',sprintf('peak throughput vs. throughput, %s',the_title),output_plots_folder,throughputs_xlims,peak_ylims,y_data_binned_max_idx);
    
    envelope.mean_throughput_gain = sorted_throughput_gain(y_data_binned_max_idx);
    envelope.edge_throughput_gain = sorted_edge_gain(y_data_binned_max_idx);
    envelope.peak_throughput_gain = sorted_peak_gain(y_data_binned_max_idx);
    envelope.peak_throughput      = sorted_peak_gain(y_data_binned_max_idx);
    envelope.mean_throughput      = sorted_throughput(y_data_binned_max_idx);
    envelope.edge_throughput      = sorted_edge(y_data_binned_max_idx);
    envelope.peak_throughput      = sorted_peak(y_data_binned_max_idx);
    envelope.fairness             = sorted_fairness(y_data_binned_max_idx);
else
    envelope = [];
end

fid = fopen(fullfile(output_plots_folder,sprintf('%s.txt',the_title)),'w');
fprintf(fid,'Mean\tEdge\tPeak\tmedian\tbeta\tSINR\tFairness\n');
for i_=1:size(all_gain_points,1)
    fprintf(fid,'%03.2f\t%03.2f\t%03.2f\t%03.2f\t%03.2f\t%03.2f\t%03.2f\n',...
        sorted_throughput_gain(i_),...
        sorted_edge_gain(i_),...
        sorted_peak_gain(i_),...
        sorted_median_gain(i_),...
        all_gain_points(i_,5),... % beta_FR
        all_gain_points(i_,6),... % SINR threshold
        sorted_fairness(i_)...
        );
end
fclose(fid);
end

function bin_max_idx = scatterplot_FFR(x_data,y_data,x_label,y_label,the_title,output_folder,the_xlims,the_ylims,varargin)

if ~isempty(varargin)
    bin_max_idx = varargin{1};
else
    [bin_centers bin_means bin_count bin_max bin_min bin_max_idx bin_min_idx] = utils.miscUtils.fit_scatterplot_data(x_data,y_data,50);
end

figure;
scatter(x_data,y_data,'.b','SizeData',50);
hold all;

% Add the precalculated envelope pints
bin_max_idx_finite = bin_max_idx(isfinite(bin_max_idx));
scatter(x_data(bin_max_idx_finite),y_data(bin_max_idx_finite),'.r','SizeData',350);

grid on;
xlabel(x_label);
ylabel(y_label);

if ~isempty(the_xlims)
    xlim(the_xlims);
end

if ~isempty(the_ylims)
    ylim(the_ylims);
end

print( fullfile(output_folder,sprintf('%d_%s.png',gcf,the_title)),'-dpng');
print( fullfile(output_folder,sprintf('%d_%s.eps',gcf,the_title)),'-depsc');
hgsave(fullfile(output_folder,sprintf('%d_%s.fig',gcf,the_title)));
end