classdef throughputMapper < handle
    % Performas a mapping between SNR and throughput according to data from LL simulations
    % (c) Josep Colom Ikuno, INTHFT, 2011
    
    properties
        mapping
    end
    
    methods
        function obj=throughputMapper(filename)
            the_data    = load(filename);
            obj.mapping = the_data.mapping; % Throughput is in Mb/s
        end
        function bits_per_channel_use = SNR_to_spectral_efficiency(obj,SNR_dB,tx_mode,nTX,nRX)
            % Error checking
            switch tx_mode
                case 1
                    if ~(nTX==1 && nRX==1)
                        error('One Tx and Rx antennas required for TX mode 1');
                    end
                    data_idx = 1;
                case 4
                    switch nTX
                        case 2
                            switch nRX
                                case 2
                                    data_idx = 2;
                                otherwise
                                    error('CLSM %.0fx%.0f not supported',nTX,nRX);
                            end
                        case 4
                            switch nRX
                                case 2
                                    data_idx = 3;
                                case 4
                                    data_idx = 4;
                                otherwise
                                    error('CLSM %.0fx%.0f not supported',nTX,nRX);
                            end
                        otherwise
                            error('CLSM %.0fx%.0f not supported',nTX,nRX);
                    end
                otherwise
                    error('Mode %.0f (%.0fx%.0f) not yet implemented',tx_mode,nTX,nRX);
            end
            
            sim_summary = obj.mapping(data_idx).sim_summary;
            
            switch sim_summary.bandwidth
                case '1.4 MHz'
                    N_RBs = 6;
                case '3 MHz'
                    N_RBs = 15;
                case '5 MHz'
                    N_RBs = 25;
                case '10 MHz'
                    N_RBs = 50;
                case '15 MHz'
                    N_RBs = 75;
                case '20 MHz'
                    N_RBs = 100;
                otherwise
                    error('Bandwidth not supported. Found: %.1 MHz',sim_summary.bandwidth/1e6);
            end
            
            % Get the proper data
            TTI_time_s = 1e-3;
            channel_uses_per_TTI = N_RBs*12*14;
            bits_per_channel_use = interp1(obj.mapping(data_idx).SNR_vec,obj.mapping(data_idx).throughput,SNR_dB)*1e6*TTI_time_s/channel_uses_per_TTI;
        end
        function plot(obj)
            % Plot throughput and gain
            throughput_figure = figure;
            throughput_axes = axes('Parent',throughput_figure);
            hold(throughput_axes,'all');
            
            for f_=1:length(obj.mapping)
                % Plot absolute throughput
                SNR_vec = obj.mapping(f_).SNR_vec;
                plot(throughput_axes,SNR_vec,obj.mapping(f_).throughput,'DisplayName',sprintf('%s (%s, %s)',obj.mapping(f_).sim_summary.tx_mode,obj.mapping(f_).sim_summary.antenna_config,obj.mapping(f_).sim_summary.bandwidth));
            end
            
            grid(throughput_axes,'on');
            legend(throughput_axes,'show','Location','NorthWest');
            xlabel(throughput_axes,'SNR [dB]');
            ylabel(throughput_axes,'throughput [Mb/s]');
            title(throughput_axes,'Throughput');
        end
        
        function plot_with_cell_layout(obj,max_SINR_dB_all,shadow_fading_used,networkPathlossMap,figure_handle,figure_handle2)
            bits_per_cu{1} = obj.SNR_to_spectral_efficiency(max_SINR_dB_all,1,1,1);
            bits_per_cu{2} = obj.SNR_to_spectral_efficiency(max_SINR_dB_all,4,2,2);
            bits_per_cu{3} = obj.SNR_to_spectral_efficiency(max_SINR_dB_all,4,4,2);
            bits_per_cu{4} = obj.SNR_to_spectral_efficiency(max_SINR_dB_all,4,4,4);

            for f_=1:length(bits_per_cu)
                ecdfs(f_) = utils.miscUtils.ecdf(bits_per_cu{f_}(:));
            end
            
            if isempty(figure_handle)
                figure_handle = figure;
            else
                figure(figure_handle);
            end
            
            if ~shadow_fading_used
                extra = '(no shadow fading)';
            else
                extra = '(with shadow fading)';
            end
            
            % Plot throughput and gain
            subplot(1,2,1);
            hold all
            displayname = cell(1,length(obj.mapping));
            all_summaries = [obj.mapping.sim_summary];
            for f_=1:length(obj.mapping)
                % Plot absolute throughput
                SNR_vec = obj.mapping(f_).SNR_vec;
                displayname{f_} = sprintf('%s %s',all_summaries(f_).tx_mode,all_summaries(f_).antenna_config);
                plot(SNR_vec,obj.mapping(f_).throughput,'DisplayName',displayname{f_});
            end
            grid on;
            legend('show','Location','NorthWest');
            xlabel('SNR [dB]');
            ylabel('throughput [Mb/s]');
            title(sprintf('Throughput, %s',all_summaries(f_).bandwidth));
            subplot(1,2,2);
            hold all;
            for f_=1:length(bits_per_cu)
                switch f_
                    case 1
                        legend_string = 'SISO';
                    case 2
                        legend_string = 'CLSM (2x2)';
                    case 3
                        legend_string = 'CLSM (4x2)';
                    case 4
                        legend_string = 'CLSM (4x4)';
                end
                plot(ecdfs(f_).x,ecdfs(f_).f,'DisplayName',legend_string);
            end
            grid on;
            legend('show','Location','SouthEast');
            xlabel('Spectral efficiency [bits/cu]');
            ylabel('F(x)');
            title(sprintf('Spectral efficiency (ROI), %s',extra));
            
            if isempty(figure_handle2)
                figure_handle2 = figure;
            else
                figure(figure_handle2);
            end
            
            for f_=1:length(obj.mapping)
                switch f_
                    case 1
                        c_lims = [min(bits_per_cu{1}(:)) max(bits_per_cu{1}(:))];
                    case 2
                        c_lims = [min([bits_per_cu{2}(:);bits_per_cu{3}(:)]) max([bits_per_cu{2}(:);bits_per_cu{3}(:)])];
                    case 3
                        c_lims = [min([bits_per_cu{2}(:);bits_per_cu{3}(:)]) max([bits_per_cu{2}(:);bits_per_cu{3}(:)])];
                    case 4
                        c_lims = [min(bits_per_cu{4}(:)) max(bits_per_cu{4}(:))];
                end
                % c_lims = [min([bits_per_cu{1}(:);bits_per_cu{2}(:);bits_per_cu{3}(:);bits_per_cu{4}(:)]) max([bits_per_cu{1}(:);bits_per_cu{2}(:);bits_per_cu{3}(:);bits_per_cu{4}(:)])];
                subplot(2,2,f_);
                imagesc(networkPathlossMap.roi_x,networkPathlossMap.roi_y,bits_per_cu{f_},c_lims);
                set(gca,'YDir','normal');
                xlabel('x pos [m]');
                ylabel('y pos [m]');
                title(sprintf('Spectral efficiency bits/cu: %s %s',displayname{f_},extra));
                colorbar;
            end
        end
    end
    
end

