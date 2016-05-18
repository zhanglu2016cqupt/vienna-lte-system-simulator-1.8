function LTE_sim_main_launcher_SLvsLL
close all
simulation_type = 'LLvsSL';

cd ..
cd ..

LTE_config = LTE_load_params(simulation_type);
LTE_config.compact_results_file = 2;

% By default the distance-to-SNR mapping is plotted. Based on it, a loop
% over the SNR points is performed

results_file_template = 'SL_sv_LL_%dx%d_mode%d_%3.2fm_%s_%0.1fMHz';
results_folder = './reproducibility/included_simulation_results/SLvsLL';
BW = 1.4e6;
nRepeats = 10;
doSims   = false;
parallel = false;
useCache = true;

% Set this variable to false if it is the first time you run this script! (i.e., no cache file available)
usePathlossMapFromMemory = true;

if usePathlossMapFromMemory && doSims
    pathlossData = load('./data_files/network_caches/network_0_rings_1_sectors_0°_offset_1m_res_free_space_3_omnidirectional_antenna_2.14GHz_freq_no_shadow_fading.mat');
else
    pathlossData = [];
end

% -10:2:40 dB SNR, as in LTE_quick_test. Obtained from the SNR-to-distance mapping
UEdistances = [
    432.227789802208	370.720100249833	317.964929483700	272.717542163799	233.908576315514,...
    200.622859823039	172.073097721032	147.587043807230	126.585056951148	108.571789561582,...
    93.1211793596458	79.8698788107161	68.5052806274136	58.7567090422091	50.3964041957171,...
    43.2247572726074	37.0728803656206	31.7990446297751	27.2753351979060	23.3959215783363,...
    20.0636714238586	17.2120387015708	14.7648636198471	12.6673330518503	10.8627727750242,...
    9.32351655548161];

UEdistances_step = 1;

if doSims && parallel&& ~matlabpool('size');
    matlabpool open
end


% Simulations
antennaConfigs = {[1 1 1] [1 1 2] [2 2 2] [3 2 2] [3 4 2] [3 4 4] [4 2 2] [4 4 2] [4 4 4]};
i_=1;
LTE_config_parfor = cell(1,length(antennaConfigs)*nRepeats);
for antennaConfIdx=1:length(antennaConfigs)
    LTE_config.bandwidth = BW;
    LTE_config.tx_mode   = antennaConfigs{antennaConfIdx}(1);
    LTE_config.nTX       = antennaConfigs{antennaConfIdx}(2);
    LTE_config.nRX       = antennaConfigs{antennaConfIdx}(3);
    LTE_config.simulation_time_tti = 100;
    LTE_config.show_network = 0;
    LTE_config.recalculate_fast_fading    = false;
    LTE_config.non_parallel_channel_trace = true;
    LTE_config.usePathlossMapFromMemory   = usePathlossMapFromMemory;
    for nRepeat = 1:nRepeats
        LTE_config.pregenerated_ff_file = fullfile(sprintf('data_files/channel_traces/SLvsLL/%.1fMHz',BW/1e6),sprintf('SLvsLL_%dx%d_mode%d_%d_%.1fMHz',LTE_config.nTX,LTE_config.nRX,LTE_config.tx_mode,nRepeat,LTE_config.bandwidth/1e6));
        LTE_config.RandStreamSeed = rand*(2^32-1);
        for UE_distanceIdx = 1:UEdistances_step:length(UEdistances)
            UE_distance = UEdistances(UE_distanceIdx);
            LTE_config.UE_distribution_radii = UE_distance;
            LTE_config.results_folder = results_folder;
            LTE_config.results_file   = sprintf(results_file_template,LTE_config.nTX,LTE_config.nRX,LTE_config.tx_mode,UE_distance,sprintf('%d',nRepeat),LTE_config.bandwidth/1e6);
            LTE_config_parfor{i_} = LTE_config;
            i_=i_+1;
        end
    end
end

ticIdx = tic;
if doSims
    if parallel
        parfor i_=1:length(LTE_config_parfor)
            if ~exist(fullfile(LTE_config_parfor{i_}.results_folder,[LTE_config_parfor{i_}.results_file '.mat']),'file')
                LTE_sim_main(LTE_config_parfor{i_},[],pathlossData);
            else
                fprintf('%s exist. Skipping.\n',fullfile(LTE_config_parfor{i_}.results_folder,LTE_config_parfor{i_}.results_file));
            end
        end
    else
        for i_=1:length(LTE_config_parfor)
            if ~exist(fullfile(LTE_config_parfor{i_}.results_folder,[LTE_config_parfor{i_}.results_file '.mat']),'file')
                LTE_sim_main(LTE_config_parfor{i_},[],pathlossData);
            else
                fprintf('%s exist. Skipping.\n',fullfile(LTE_config_parfor{i_}.results_folder,LTE_config_parfor{i_}.results_file));
            end
        end
    end
end
elapsedTime_s = toc(ticIdx);

fprintf('Elapsed time: %3.0f s\n',elapsedTime_s);
save(sprintf('time_%.1fMHz',LTE_config.bandwidth/1e6),'elapsedTime_s');

% Load results
loaded_LL_results = load(sprintf('./reproducibility/included_simulation_results/SLvsLL/LL_PedA_%.1fMHz.mat',BW/1e6));

if ~exist(sprintf('./reproducibility/included_simulation_results/SLvsLL/SL_PedA_%.1fMHz.mat',BW/1e6),'file') || ~useCache
    for antennaConfIdx=1:length(antennaConfigs)
        UE_distance = [];
        UE_SNRs     = [];
        tx_mode     = antennaConfigs{antennaConfIdx}(1);
        nTX         = antennaConfigs{antennaConfIdx}(2);
        nRX         = antennaConfigs{antennaConfIdx}(3);
        LL_SNR      = loaded_LL_results.LL_results{100*tx_mode+10*nTX+1*nRX};
        UE_throughputs = [];
        UE_SNRs        = [];
        
        for UE_distanceIdx = 1:length(UEdistances)
            UE_distance = UEdistances(UE_distanceIdx);
            resultsFiles = dir(fullfile(LTE_config.results_folder,sprintf([results_file_template '.mat'],nTX,nRX,tx_mode,UE_distance,'*',BW/1e6)));
            for fileIdx=1:length(resultsFiles)
                resultsFile = fullfile(LTE_config.results_folder,resultsFiles(fileIdx).name);
                fprintf('Loading %s (%d/%d)\n',resultsFile,fileIdx,length(resultsFiles));
                simResults = load(resultsFile);
                UE_throughputs(UE_distanceIdx,fileIdx) = simResults.the_UE_traces.average_throughput_Mbps;
                UE_SNRs(UE_distanceIdx,fileIdx)        = simResults.the_UE_traces.SNR_dB(1);
            end
        end
        SL_results{100*tx_mode+10*nTX+1*nRX}.SNR = UE_SNRs(:,1);
        SL_results{100*tx_mode+10*nTX+1*nRX}.cell_throughput_over_SNR_Mbps = mean(UE_throughputs,2);
    end
    save(sprintf('./reproducibility/included_simulation_results/SLvsLL/SL_PedA_%.1fMHz.mat',BW/1e6),'SL_results');
else
    load(sprintf('./reproducibility/included_simulation_results/SLvsLL/SL_PedA_%.1fMHz.mat',BW/1e6));
end

export_plots_to_EPS = true;
plots_dir = './reproducibility/included_simulation_results/SLvsLL/plots';
% SIXO
plot_SLvsLL({[1 1 2] [1 1 1]},loaded_LL_results,SL_results,false);
if export_plots_to_EPS
    print(gcf,fullfile(plots_dir,sprintf('%d_%.1fMHz.eps',gcf,BW/1e6)),'-dpsc2');
end

% TxD
plot_SLvsLL({[2 2 2]},loaded_LL_results,SL_results,false);
if export_plots_to_EPS
    print(gcf,fullfile(plots_dir,sprintf('%d_%.1fMHz.eps',gcf,BW/1e6)),'-dpsc2');
end

% OLSM
plot_SLvsLL({[3 4 4] [3 4 2] [3 2 2]},loaded_LL_results,SL_results,false);
if export_plots_to_EPS
    print(gcf,fullfile(plots_dir,sprintf('%d_%.1fMHz.eps',gcf,BW/1e6)),'-dpsc2');
end

% CLSM
plot_SLvsLL({[4 4 4] [4 4 2] [4 2 2]},loaded_LL_results,SL_results,false);
if export_plots_to_EPS
    print(gcf,fullfile(plots_dir,sprintf('%d_%.1fMHz.eps',gcf,BW/1e6)),'-dpsc2');
end

% saved_data.SISO_SL_results = SL_results{111};
% saved_data.SISO_LL_results = loaded_LL_results.LL_results{111};
% save('scenario0_VUT','-struct','saved_data');

function plot_SLvsLL(antennaConfigs,LL_results,SL_results,plot_deviation)

f_throughput = figure;
a_throughput = axes('Parent',f_throughput);
hold(a_throughput,'on');

if plot_deviation
    f_difference = figure;
    a_difference = axes('Parent',f_difference);
    hold(a_difference,'on');
end

% Plot throughput plots
for antennaConfIdx=1:length(antennaConfigs)
    tx_mode    = antennaConfigs{antennaConfIdx}(1);
    nTX        = antennaConfigs{antennaConfIdx}(2);
    nRX        = antennaConfigs{antennaConfIdx}(3);
    
    LL_results_current = LL_results.LL_results{100*tx_mode+10*nTX+1*nRX};
    SL_results_current = SL_results{100*tx_mode+10*nTX+1*nRX};
    
    % Choose line color
    switch tx_mode
        case 1
            switch nRX
                case 1
                    plotName = 'SISO';
                    c        = 'b';
                    m        = 'x';
                otherwise
                    plotName = sprintf('1x%d MRC',nRX);
                    c        = 'b';
                    m        = 'o';
            end
        case 2
            plotName = 'TxD';
            c        = 'r';
        case 3
            plotName = 'OLSM';
            c        = 'k';
        case 4
            plotName = 'CLSM';
            c        = 'm';
    end
    
    switch tx_mode
        case {2,3,4}
            plotName = sprintf('%s %dx%d',plotName,nTX,nRX);
            if nTX==2&&nRX==2
                m = 'o';
            elseif nTX==4&&nRX==2
                m = '+';
            else % 4x4
                m = 's';
            end
    end
    
    LL_interp_results = interp1(LL_results_current.SNR,LL_results_current.cell_throughput_over_SNR_Mbps,SL_results_current.SNR);LL_interp_results(end) = LL_results_current.cell_throughput_over_SNR_Mbps(end);
    
    markerSize = 5;
    plot(a_throughput,SL_results_current.SNR,LL_interp_results,sprintf('%s%s%s',c,m,':'),'DisplayName',sprintf('%s (link level)',plotName),'MarkerSize',markerSize);
    plot(a_throughput,SL_results_current.SNR,SL_results_current.cell_throughput_over_SNR_Mbps,sprintf('%s%s%s',c,m,'-'),'DisplayName',sprintf('%s (system level)',plotName),'MarkerSize',markerSize);
    
    if plot_deviation
        deviation = LL_interp_results-SL_results_current.cell_throughput_over_SNR_Mbps;
        plot(a_difference,SL_results_current.SNR,deviation,sprintf('%s%s%s',c,m,'-'),'DisplayName',sprintf('%s',plotName));
    end
    
    if(antennaConfIdx==1)
        hold(a_throughput,'on');
        grid(a_throughput,'on');
        if plot_deviation
            hold(a_difference,'on');
            grid(a_difference,'on');
        end
    end
end
legend(a_throughput,'show','Location','NorthWest');
xlim(  a_throughput,[-10 40]);
xlabel(a_throughput,'SNR [dB]');
ylabel(a_throughput,'Throughput [Mbit/s]');

if plot_deviation
    legend(a_difference,'show','Location','NorthWest');
    xlim(  a_difference,[-10 40]);
    xlabel(a_difference,'SNR [dB]');
    ylabel(a_difference,'Throughput deviation [Mbit/s]');
end