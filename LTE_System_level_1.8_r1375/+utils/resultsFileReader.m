classdef resultsFileReader < handle
    % Encapsulates methods used for reading and plotting data from the results files
    % (c) Josep Colom Ikuno, INTHFT, 2011
    
    properties
        % The results data
        data
        % Whether the loaded data is in compact for or not
        compact_format
    end
    
    methods
        function obj = resultsFileReader(filename)
            fprintf('Loading results file: %s, ',filename);
            if exist(filename,'file')
                fprintf('exists, ');
                file_data = dir(filename);
                if length(file_data)==1
                    fprintf('%.1f MB, ',file_data.bytes/1024/1024);
                else
                    error('More than one file found. Non-unique filename');
                end
            else
                error('File does not exist');
            end
            obj.data = load(filename);
            if isfield(obj.data,'simulation_traces')
                obj.compact_format = false;
                fprintf('standard file format');
            else
                obj.compact_format = true;
                fprintf('compact file format');
            end
            fprintf('\n');
        end
        function clear(obj)
            if ~obj.compact_format
                utils.miscUtils.tidy_up_memory_before_closing(obj.data.UEs,obj.data.eNodeBs_sectors,obj.data.eNodeBs);
            end
            obj.data = [];
        end
        function [avg_UE_throughput_Mbps average_spectral_efficiency_bit_per_cu] = get_UE_average_throughput(obj)
            if obj.compact_format
                avg_UE_throughput_Mbps                 = [obj.data.the_UE_traces.average_throughput_Mbps];
                average_spectral_efficiency_bit_per_cu = [obj.data.the_UE_traces.average_spectral_efficiency_bit_per_cu];
            else
                avg_UE_throughput_Mbps                 = [obj.data.simulation_traces.UE_traces.average_throughput_Mbps];
                average_spectral_efficiency_bit_per_cu = [obj.data.simulation_traces.UE_traces.average_spectral_efficiency_bit_per_cu];
            end
        end
        function [avg_UE_throughput_Mbps average_spectral_efficiency_bit_per_cu] = get_UE_average_throughput_just_finite_values(obj)
            [avg_UE_throughput_Mbps, average_spectral_efficiency_bit_per_cu] = obj.get_UE_average_throughput;
            finite_idxs = isfinite(avg_UE_throughput_Mbps);
            avg_UE_throughput_Mbps                 = avg_UE_throughput_Mbps(finite_idxs);
            average_spectral_efficiency_bit_per_cu = average_spectral_efficiency_bit_per_cu(finite_idxs);
        end
        function avg_wideband_SINR = get_UE_average_wideband_SINR(obj)
            if obj.compact_format
                nUEs = length(obj.data.the_UE_traces);
            else
                nUEs = length(obj.data.simulation_traces.UE_traces);
            end
            avg_wideband_SINR = zeros(1,nUEs);
            for u_=1:nUEs
                if obj.compact_format
                    avg_wideband_SINR(u_) = mean(obj.data.the_UE_traces(u_).wideband_SINR);
                else
                    avg_wideband_SINR(u_) = median(obj.data.simulation_traces.UE_traces(u_).wideband_SINR);
                end
            end
        end
        function [cell_average_throughput_Mbps cell_fairness] = get_cell_average_throughput(obj)
            if obj.compact_format
                eNodeB_traces = [obj.data.the_eNodeB_traces];
                UE_traces     = [obj.data.the_UE_traces];
            else
                eNodeB_traces = [obj.data.simulation_traces.eNodeB_tx_traces];
                UE_traces     = [obj.data.simulation_traces.UE_traces];
            end
            cell_average_throughput_Mbps = zeros(1,length(eNodeB_traces));
            cell_throughputs             = cell(1,length(eNodeB_traces));
            cell_fairness                = NaN(1,length(eNodeB_traces));
            for u_=1:length(UE_traces)
                avg_throughput_Mbps = UE_traces(u_).average_throughput_Mbps;
                if isfinite(avg_throughput_Mbps)
                    attached_eNodeB = UE_traces(u_).attached_eNodeB(1);
                    cell_throughputs{attached_eNodeB}          = [cell_throughputs{attached_eNodeB} avg_throughput_Mbps];
                    cell_average_throughput_Mbps(attached_eNodeB) = cell_average_throughput_Mbps(attached_eNodeB) + avg_throughput_Mbps;
                end
            end
            for b_=1:length(cell_fairness)
                if ~isempty(cell_throughputs{b_})
                    throughput_vector = cell_throughputs{b_};
                    cell_fairness(b_) = (sum(throughput_vector)^2)/(length(throughput_vector)*sum(throughput_vector.^2));
                end
            end
        end
        function avg_wideband_SINR = get_UE_average_wideband_SINR_just_finite_values(obj)
            avg_wideband_SINR = obj.get_UE_average_wideband_SINR;
            finite_idxs = isfinite(avg_wideband_SINR);
            avg_wideband_SINR = avg_wideband_SINR(finite_idxs);
        end
        function LTE_config = get_simulation_config(obj)
            LTE_config = obj.data.LTE_config;
        end
        function simulation_basic_data = get_config_params_summary(obj)
            LTE_config = obj.data.LTE_config;
            simulation_basic_data.bandwidth = sprintf('%.1f MHz',LTE_config.bandwidth/1e6);
            switch LTE_config.tx_mode
                case 1
                simulation_basic_data.tx_mode = 'SISO';
                case 2
                    simulation_basic_data.tx_mode = 'TxD';
                case 3
                    simulation_basic_data.tx_mode = 'OLSM';
                case 4
                    simulation_basic_data.tx_mode = 'CLSM';
            end
            simulation_basic_data.antenna_config    = sprintf('%.0fx%.0f',LTE_config.nTX,LTE_config.nRX);
            simulation_basic_data.channel           = LTE_config.channel_model.type;
            simulation_basic_data.UE_speed          = sprintf('%.1f Km/h',LTE_config.UE_speed*3.6);
            simulation_basic_data.simulation_length = LTE_config.simulation_time_tti;
            simulation_basic_data.feedback_delay    = LTE_config.feedback_channel_delay;
            simulation_basic_data.scheduler         = LTE_config.scheduler;
            simulation_basic_data.site_count        = LTE_config.site_count;
            simulation_basic_data.eNodeB_count      = LTE_config.eNodeB_count;
            simulation_basic_data.UE_count          = LTE_config.UE_count;
        end
    end
    
    methods(Static)
        function [common_fields output_fieldnames] = find_common_params(sim_summarized_parameters)
            input_fieldnames = fieldnames(sim_summarized_parameters(1));
            for f_=1:length(input_fieldnames)
                current_field = sim_summarized_parameters(1).(input_fieldnames{f_});
                if ischar(current_field)
                    current_unique_fields = unique({sim_summarized_parameters.(input_fieldnames{f_})});
                else
                    current_unique_fields = unique([sim_summarized_parameters.(input_fieldnames{f_})]);
                end
                
                if length(current_unique_fields)==1
                    if ischar(current_field)
                        common_fields.(input_fieldnames{f_}) = current_unique_fields{1};
                    else
                        common_fields.(input_fieldnames{f_}) = current_unique_fields(1);
                    end
                end
            end
            output_fieldnames = fieldnames(common_fields);
        end
        function [UE_avg_throughput_ecdfs UE_avg_spectral_eff_ecdfs UE_wideband_SINR_ecdfs sim_summarized_parameters results_object] = plot_files(filenames,varargin)
            if ~iscell(filenames)
                filenames = {filenames};
            end
            if ~isempty(varargin)
                folder_set = true;
                folder = varargin{1};
                if length(varargin)>1
                    plot_figures = varargin{2};
                else
                    plot_figures = true;
                end
            else
                folder_set   = false;
                plot_figures = true;
            end
            [nSimulations nRepetitions] = size(filenames);
            UE_avg_throughput_Mbps = cell(1,nSimulations);
            
            %% Load results files
            for f_=1:nSimulations
                UE_avg_throughput_Mbps{f_} = [];
                UE_avg_spectral_efficiency_bit_per_cu{f_} = [];
                UE_wideband_SINR_dB{f_}  = [];
                clear sim_summary_repetition
                for r_=1:nRepetitions
                    if ~isempty(filenames{f_,r_})
                        if folder_set
                            current_filename = fullfile(folder,filenames{f_,r_});
                        else
                            current_filename = filenames{f_,r_};
                        end
                        if ~strcmp(current_filename(end-3:end),'.mat')
                            current_filename = [current_filename '.mat'];
                        end
                        current_results                           = utils.resultsFileReader(current_filename);
                        sim_summary_repetition(r_)                = current_results.get_config_params_summary;
                        [ current_avg_thr_Mbps,...
                          current_avg_spectral_eff_bit_per_cu ]   = current_results.get_UE_average_throughput;
                        UE_avg_throughput_Mbps{f_}                = [UE_avg_throughput_Mbps{f_} current_avg_thr_Mbps];
                        UE_avg_spectral_efficiency_bit_per_cu{f_} = [UE_avg_spectral_efficiency_bit_per_cu{f_} current_avg_spectral_eff_bit_per_cu];
                        UE_wideband_SINR_dB{f_}                   = [UE_wideband_SINR_dB{f_}  current_results.get_UE_average_wideband_SINR];
                    end
                end
                UE_avg_throughput_ecdfs(f_)        = utils.miscUtils.ecdf(UE_avg_throughput_Mbps{f_});
                UE_avg_throughput_ecdfs(f_).what   = 'UE average Throughput';
                UE_avg_throughput_ecdfs(f_).unit   = 'Mb/s';
                UE_avg_spectral_eff_ecdfs(f_)      = utils.miscUtils.ecdf(UE_avg_spectral_efficiency_bit_per_cu{f_});
                UE_avg_spectral_eff_ecdfs(f_).what = 'UE average Spectral efficiency';
                UE_avg_spectral_eff_ecdfs(f_).unit = 'bit/cu';
                UE_wideband_SINR_ecdfs(f_)          = utils.miscUtils.ecdf(UE_wideband_SINR_dB{f_});
                UE_wideband_SINR_ecdfs(f_).what     = 'UE median Wideband SINR';
                UE_wideband_SINR_ecdfs(f_).unit     = 'dB';
                [sim_summarized_parameters(f_) ~]  = utils.resultsFileReader.find_common_params(sim_summary_repetition);
                
                if nRepetitions==1 && nSimulations==1
                    results_object = current_results;
                else
                    results_object = [];
                end
            end
            
            if plot_figures
                %% If some data is equal for all simulations, print it on the title
                after_title   = [];
                title_suffix  = [];
                [common_fields common_fields_list] = utils.resultsFileReader.find_common_params(sim_summarized_parameters);
                if strmatch('bandwidth',common_fields_list)
                    after_title = ': ';
                    if isempty(title_suffix)
                        comma = [];
                    else
                        comma = ', ';
                    end
                    title_suffix = sprintf('%s%s%s',title_suffix,comma,common_fields.bandwidth);
                end
                if strmatch('channel',common_fields_list)
                    after_title = ': ';
                    if isempty(title_suffix)
                        comma = [];
                    else
                        comma = ', ';
                    end
                    title_suffix = sprintf('%s%s%s',title_suffix,comma,common_fields.channel);
                end
                if strmatch('site_count',common_fields_list)
                    after_title = ': ';
                    if isempty(title_suffix)
                        comma = [];
                    else
                        comma = ', ';
                    end
                    title_suffix = sprintf('%s%s%.0f sites',title_suffix,comma,common_fields.site_count);
                end
                if strmatch('UE_count',common_fields_list)
                    after_title = ': ';
                    if isempty(title_suffix)
                        comma = [];
                    else
                        comma = ', ';
                    end
                    title_suffix = sprintf('%s%s%.0f UEs',title_suffix,comma,common_fields.UE_count);
                end
                
                %% Plot UE throughput ECDFs
                figure;
                for f_=1:nSimulations
                    if isfield(sim_summarized_parameters(f_),'UE_speed')
                        speed_string = sprintf(', %s',sim_summarized_parameters(f_).UE_speed);
                    else
                        speed_string = [];
                    end
                    current_display_name = sprintf('%s %s%s',sim_summarized_parameters(f_).antenna_config,sim_summarized_parameters(f_).tx_mode,speed_string);
                    plot(UE_avg_throughput_ecdfs(f_).x,UE_avg_throughput_ecdfs(f_).f,'DisplayName',current_display_name);
                    if f_==1
                        hold all
                    end
                end
                grid on;
                legend('show','Location','SouthEast');
                xlabel(sprintf('Throughput [%s]',UE_avg_throughput_ecdfs(end).unit)); % Assume that the units are all the same
                ylabel('F(x)');
                title(sprintf('Throughput ECDF%s%s',after_title,title_suffix));
                
                %% Plot spectral efficiency ECDF
                figure;
                for f_=1:nSimulations
                    if isfield(sim_summarized_parameters(f_),'UE_speed')
                        speed_string = sprintf(', %s',sim_summarized_parameters(f_).UE_speed);
                    else
                        speed_string = [];
                    end
                    current_display_name = sprintf('%s %s%s',sim_summarized_parameters(f_).antenna_config,sim_summarized_parameters(f_).tx_mode,speed_string);
                    plot(UE_avg_spectral_eff_ecdfs(f_).x,UE_avg_spectral_eff_ecdfs(f_).f,'DisplayName',current_display_name);
                    if f_==1
                        hold all
                    end
                end
                grid on;
                legend('show','Location','SouthEast');
                xlabel(sprintf('Spectral efficiency [%s]',UE_avg_spectral_eff_ecdfs(end).unit)); % Assume that the units are all the same
                ylabel('F(x)');
                title(sprintf('Spectral efficiency ECDF%s%s',after_title,title_suffix));
            end
        end
        
        % Parsed a folder with results files and returns the antenna configurations found there
        function results = scan_results_filenames(folder)
            filenames = ls(fullfile(folder,'*.mat'));
            nFiles = size(filenames,1);
            
            antenna_config     = cell(1,nFiles);
            antenna_config_idx = zeros(1,nFiles);
            
            for f_=1:nFiles
                token_idx = 1;
                file = filenames(f_,:);
                remain = file;
                while ~isempty(remain)
                    [token{token_idx}, remain] = strtok(remain,'_');
                    
                    % Check antenna config
                    if regexp(token{token_idx},'.*[0-9]x[0-9].*')
                        antenna_config{f_} = token{token_idx};
                    end
                    token_idx = token_idx + 1;
                end
            end
            
            results.antenna_configs = unique(antenna_config);
        end
        
        % Plots after the end of one simulation
        function plot_simulation_results(filename)
            [UE_avg_throughput_ecdfs UE_avg_spectral_eff_ecdfs UE_wideband_SINR_ecdfs sim_summarized_parameters results_object] = utils.resultsFileReader.plot_files(filename,[],false);
            figure_name_aggregates = sprintf('%s %s, %s %s: aggregate plots',sim_summarized_parameters.antenna_config,sim_summarized_parameters.tx_mode,sim_summarized_parameters.bandwidth,sim_summarized_parameters.channel);
            
            figure('Name',figure_name_aggregates);
            the_axis = subplot(2,2,1);
            utils.resultsFileReader.plot_ecdf_on_axis(the_axis,UE_avg_throughput_ecdfs);
            the_axis = subplot(2,2,2);
            utils.resultsFileReader.plot_ecdf_on_axis(the_axis,UE_avg_spectral_eff_ecdfs);
            the_axis = subplot(2,2,3);
            utils.resultsFileReader.plot_ecdf_on_axis(the_axis,UE_wideband_SINR_ecdfs);
            the_axis = subplot(2,2,4);
            non_zero_throughput = UE_avg_throughput_ecdfs.input_data>0;
            scatter(UE_wideband_SINR_ecdfs.input_data(non_zero_throughput),UE_avg_throughput_ecdfs.input_data(non_zero_throughput),'.');
            grid on;
            xlabel(sprintf('%s [%s]',UE_wideband_SINR_ecdfs.what,UE_wideband_SINR_ecdfs.unit));
            ylabel(sprintf('%s [%s]',UE_avg_throughput_ecdfs.what,UE_avg_throughput_ecdfs.unit));
            title(sprintf('%s-to-%s',UE_wideband_SINR_ecdfs.what,UE_avg_throughput_ecdfs.what));
            
            % Plot eNodeb and UE locations
            % if ~results_object.compact_format
                % utils.plotUtils.plot_eNodeBs_and_UEs(results_object.data.LTE_config.plots.user_positions,results_object.data.eNodeBs,results_object.data.UEs,results_object.data.networkPathlossMap,false);
            % else
                % print_log(1,'More extensive results plots can only be shown from simulation results with LTE_config.compact_results_file=false\n');
            % end
        end
        
        function plot_ecdf_on_axis(the_axis,ecdf_data)
            plot(the_axis,ecdf_data.x,ecdf_data.f);
            grid on;
            xlabel(sprintf('%s [%s]',ecdf_data.what,ecdf_data.unit));
            ylabel('F(x)');
            title(sprintf('%s ECDF',ecdf_data.what));
        end
        
        function [UE_ids cell_sum_throughput] = get_UEs_in_given_cells(cell_ids,the_UE_traces)
            N_UEs  = length(the_UE_traces);
            cell_sum_throughput = zeros(1,length(cell_ids));
            UE_ids = false(1,N_UEs);
            for c_idx=1:length(cell_ids)
                c_ = cell_ids(c_idx);
                for u_=1:N_UEs
                    if ~isempty(find(the_UE_traces(u_).attached_eNodeB==c_,1,'first'))
                        % This UE is in the given cell
                        UE_ids(u_) = true;
                        cell_sum_throughput(c_idx) = cell_sum_throughput(c_idx) + the_UE_traces(u_).average_throughput_Mbps;
                    end
                end
            end
        end
    end
end

