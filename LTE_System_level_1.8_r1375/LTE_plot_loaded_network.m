function LTE_plot_loaded_network(LTE_config,...
    sites,eNodeBs,...
    networkPathlossMap,...
    CQI_mapper,use_subplots,...
    networkShadowFadingMap)
% Function that shows a few plots after a network has been loaded from a file
% (c) Josep Colom Ikuno, INTHFT, 2008
% www.nt.tuwien.ac.at

if ~isempty(networkShadowFadingMap) && ~isa(networkShadowFadingMap,'channel_gain_wrappers.shadowFadingDummyMap')
    plot_shadow_fading = true;
else
    plot_shadow_fading = false;
end

% Not very clean, but will do the trick
roi_to_map_x = networkPathlossMap.roi_x;
roi_to_map_y = networkPathlossMap.roi_y;

%% Plot the antenna gain pattern
antenna_types =  cell(1,length(eNodeBs));
if LTE_config.show_network>0
    for i_=1:length(eNodeBs) % Needed because Matlab does not have and built-in class-casting yet
        antenna_types{i_} = eNodeBs(i_).antenna.antenna_type;
        if isnumeric(antenna_types{i_})
            antenna_types{i_} = num2str(antenna_types{i_});
        end
    end
    [anntenna_types, m, n] = unique(antenna_types); %#ok
    number_of_subplots_row_col = ceil(sqrt(length(anntenna_types)));
    
    % Plot antenna gain pattern for each different antenna in the simulation
    figure(LTE_config.plots.antenna_gain_pattern);
    for ant_idx = 1:length(anntenna_types)
        example_ant_idx = find(n==ant_idx,1,'first');
        an_antenna = eNodeBs(example_ant_idx).antenna;
        
        % Choose correct antenna pattern
        the_axes = subplot(number_of_subplots_row_col,number_of_subplots_row_col,ant_idx,'replace');
        
        if ~an_antenna.pattern_is_3D
            % 2D antenna
            angle = -180:0.1:180;
            gain = zeros(1,length(angle));
            for i_=1:length(angle)
                gain(i_) = an_antenna.gain(angle(i_));
            end
            cla(the_axes);
            plot(the_axes,angle,gain);
            ylim(the_axes,ylim*1.1);
            title(the_axes,{['Antenna gain, ' an_antenna.antenna_type ' antenna']});
            xlabel(the_axes,{'\theta [°]'});
            ylabel(the_axes,{'gain [dB]'});
            box(the_axes,'on');
            grid(the_axes,'on');
        else
            % 3D antenna
            plot_tilt = 0; % Plot with 0° electrical (and mechanical) tilt
            data_limits = [-15 0 3];
            switch class(an_antenna)
                case 'antennas.TS36942_3DAntenna'
                    hor_degrees  = -179:180;
                    ver_degrees  = -179:180;
                    hor_gain     = an_antenna.gain(hor_degrees,0,plot_tilt);
                    ver_gain     = an_antenna.gain(0,ver_degrees,plot_tilt);
                    max_gain     = an_antenna.max_antenna_gain;
                    title_string = sprintf('%s antenna, blue: hor, red: vert [dBi], %.0f° electrical tilt',an_antenna.antenna_type,plot_tilt);
                otherwise
                    % Kathrein antenna
                    [hor_degrees, hor_gain, ver_degrees, ver_gain, max_gain] = ...
                        an_antenna.gain_patterns(plot_tilt);
                    title_string = sprintf('%d antenna, blue: hor, red: vert [dBi], %.0f° electrical tilt',an_antenna.antenna_type,plot_tilt);
            end
            utils.miscUtils.polar2(hor_degrees/180*pi, hor_gain-max_gain, data_limits,'blue');
            hold(the_axes,'all');
            utils.miscUtils.polar2(ver_degrees/180*pi, ver_gain-max_gain, data_limits,'red');
            title(the_axes,title_string);
            hold(the_axes,'off');
        end
    end
end

%% Plot of how the macroscopic pathloss looks like
if LTE_config.show_network>0 && LTE_config.macroscopic_pathloss_is_model
    sIdx = 1;
    all_pathloss_model_names = [];
    for s_=1:length(eNodeBs)
        if ~isempty(eNodeBs(s_).macroscopic_pathloss_model)
            all_pathloss_model_names{sIdx} = eNodeBs(s_).macroscopic_pathloss_model.name; %#ok<AGROW>
            sIdx = sIdx + 1;
        end
    end
    if ~isempty(all_pathloss_model_names)
        [unique_pathlosses m n] = unique(all_pathloss_model_names); %#ok
        % Will set the maximum distance as the diagonal that crosses the ROI
        range = sqrt((roi_to_map_x(2)-roi_to_map_x(1))^2+(roi_to_map_y(2)-roi_to_map_y(1))^2);
        distances = 0:LTE_config.map_resolution:range;
        
        if use_subplots
            figure(LTE_config.plots.macroscopic_pathloss);
        end
        for i_=1:length(m)
            if use_subplots
                subplot(1,length(m),i_);
            else
                figure;
            end
            macroscopic_pathloss_model = eNodeBs(m(i_)).macroscopic_pathloss_model;
            pathlosses = macroscopic_pathloss_model.pathloss(distances);
            plot(distances,pathlosses);
            title(['Macroscopic pathloss, using ' macroscopic_pathloss_model.name ' model']);
            xlabel('Distance [m]');
            ylabel('Pathloss [dB]');
            box on;
            grid on;
        end
    end
end

%% Plot of sector macroscopic pathlosses
if LTE_config.show_network>1
    number_rows = 3;
    number_cols = 4;
    num_figures = ceil(length(eNodeBs)/(number_rows*number_cols));
    
    b_ = 1;
    for f_=0:(num_figures-1)
        if use_subplots
            figure(LTE_config.plots.macroscopic_pathloss_sector+f_);
        end
        for subplot_idx = 1:(number_rows*number_cols)
            if b_ <= length(eNodeBs)
                parent_eNodeB = eNodeBs(b_).parent_eNodeB;
                
                if use_subplots
                    subplot(number_rows,number_cols,subplot_idx);
                else
                    figure;
                end
                imagesc(...
                    networkPathlossMap.roi_x,networkPathlossMap.roi_y,...
                    10*log10(networkPathlossMap.pathloss(:,:,b_)));
                set(gca,'YDir','normal');
                title(sprintf('eNodeB %d sector %d',parent_eNodeB.id,eNodeBs(b_).id));
                colorbar;
                hold on;
                scatter(parent_eNodeB.pos(1),parent_eNodeB.pos(2),'MarkerEdgeColor','white','MarkerFaceColor','black');
                text(parent_eNodeB.pos(1)+5*networkPathlossMap.data_res,parent_eNodeB.pos(2),num2str(b_),'Color','w');
                for bIdx=1:length(sites)
                    if sites(bIdx).id~=b_
                        scatter(sites(bIdx).pos(1),sites(bIdx).pos(2),'.','MarkerEdgeColor',0.75*[1 1 1]);
                    end
                end
                b_ = b_ + 1;
            end
        end
    end
end

%% Plot shadow fading
if LTE_config.show_network>1 && plot_shadow_fading
    num_eNodeBs = length(sites);
    N_cols = 3;
    N_rows = ceil(num_eNodeBs/N_cols);
    if use_subplots
        figure(LTE_config.plots.shadow_fading_loss);
    end
    for i_=1:num_eNodeBs
        if use_subplots
            subplot(N_rows,N_cols,i_);
        else
            figure;
        end
        imagesc(...
            networkShadowFadingMap.roi_x,networkShadowFadingMap.roi_y,...
            10*log10(networkShadowFadingMap.pathloss(:,:,i_)));
        set(gca,'YDir','normal');
        title(['Shadow fading, site ' num2str(i_)]);
        colorbar;
    end
    
    % Histogram, to see if they are really gaussian or not
    figure;
    hold all;
    histBins = 50;
    for i_=1:num_eNodeBs
        if i_==1
            Nrandn = 200000;
            [n_r, xout_r] = hist(networkShadowFadingMap.std*randn(1,Nrandn),histBins);
            plot(xout_r,n_r/Nrandn/(xout_r(2)-xout_r(1)),'r','DisplayName','Normal distribution');
        end
        reshaped_map = 10*log10(reshape(networkShadowFadingMap.pathloss(3:(end-3),3:(end-3),i_),1,[]));
        map_mean = mean(reshaped_map);
        map_sd = std(reshaped_map);
        [n, xout] = hist(reshaped_map,histBins);
        tot_area = numel(reshaped_map)*(xout(2)-xout(1));
        plot(xout,n/tot_area,'k','DisplayName','Correlated shadow fading maps');
        if i_==1
            legend('show','Location','SouthEast');
        end
    end
    grid on;
    plot(xout_r,n_r/Nrandn/(xout_r(2)-xout_r(1)),'r','DisplayName','Normal distribution');
    
    % Plot shadow fading cross-correlation
    nonOneCorrs = networkShadowFadingMap.sn_ccorr(networkShadowFadingMap.sn_ccorr~=1);
    figure;
    imagesc(networkShadowFadingMap.sn_ccorr);
    colorbar;
    title(sprintf('mean: %g, std: %g',mean(nonOneCorrs),std(nonOneCorrs)));
end

%% Plots previously called from the LTE_common_calculate_cell_capacity function
if LTE_config.show_network>0
    if plot_shadow_fading
        plot_SINR_plots(LTE_config,sites,eNodeBs,networkPathlossMap,CQI_mapper,1,use_subplots);
        plot_SINR_plots(LTE_config,sites,eNodeBs,networkPathlossMap,CQI_mapper,2,use_subplots);
    else
        plot_SINR_plots(LTE_config,sites,eNodeBs,networkPathlossMap,CQI_mapper,3,use_subplots);
    end
end

function plot_SINR_plots(LTE_config,eNodeBs,eNodeBs_sectors,networkPathlossMap,CQI_mapper,mode,use_subplots)

%% Plots previously in the LTE_common_calculate_cell_capacity function
switch mode
    case 1
        SINR_dB_to_plot   = networkPathlossMap.SINR;
        sector_assignment = networkPathlossMap.sector_assignment;
        cell_centers      = networkPathlossMap.sector_centers;
        diff_SINR_dB      = networkPathlossMap.diff_SINR_dB;
        plot_type         = 'macroscopic and shadow fading';
        figure(LTE_config.plots.sector_SINR);
    case 2
        SINR_dB_to_plot   = networkPathlossMap.SINR2;
        sector_assignment = networkPathlossMap.sector_assignment2;
        cell_centers      = networkPathlossMap.sector_centers2;
        diff_SINR_dB      = networkPathlossMap.diff_SINR_dB2;
        plot_type         = 'macroscopic fading';
        figure(LTE_config.plots.sector_SINR_no_shadowing);
    case 3
        SINR_dB_to_plot   = networkPathlossMap.SINR;
        sector_assignment = networkPathlossMap.sector_assignment;
        cell_centers      = networkPathlossMap.sector_centers;
        diff_SINR_dB      = networkPathlossMap.diff_SINR_dB;
        plot_type         = 'macroscopic fading';
        figure(LTE_config.plots.sector_SINR);
    otherwise
        error('Mode %d not valid',mode);
end

% Edge calculation
edge_function_exists = exist('edge','file');
if edge_function_exists
    all_edges = edge(sector_assignment,'sobel',0);
else
    all_edges = false(size(sector_assignment));
end
[all_edges_list(:,2), all_edges_list(:,1)] = find(all_edges);
all_edges_list_pos                         = LTE_common_pixel_to_pos(all_edges_list,networkPathlossMap.coordinate_origin,networkPathlossMap.data_res);

%% Plot SINR (all sectors)
if use_subplots
    subplot(2,2,1);
else
    figure;
end
SINR_limits = [-5 20];
imagesc(networkPathlossMap.roi_x,networkPathlossMap.roi_y,SINR_dB_to_plot); % Easily modifiable to plot the SNR
set(gca,'YDir','normal');
title(sprintf('ROI max SINR (SISO, %s)',plot_type));
colorbar;
caxis(SINR_limits);
xlabel('x pos [m]');
ylabel('y pos [m]');
hold on;

% Plot target sector boundary
if mode==2
    scatter(all_edges_list_pos(:,1),all_edges_list_pos(:,2),'.','MarkerEdgeColor','w','SizeData', 1);
end

% Plot where the BTs are and add a text legend to know where are the BTSs
for b_ = 1:length(eNodeBs)
    scatter(eNodeBs(b_).pos(1),eNodeBs(b_).pos(2),'MarkerEdgeColor','k','MarkerFaceColor','w');
    text(eNodeBs(b_).pos(1)+7*networkPathlossMap.data_res,eNodeBs(b_).pos(2),num2str(b_),'Color','k');
end

%% Plot CQIs
if use_subplots
    subplot(2,2,2);
else
    figure;
end
mapped_CQIs = CQI_mapper.SINR_to_CQI(SINR_dB_to_plot);
imagesc(networkPathlossMap.roi_x,networkPathlossMap.roi_y,mapped_CQIs);
set(gca,'YDir','normal');
title(sprintf('SISO CQI mapping (%s).',plot_type));
colorbar('YTick',0:15);
caxis([0 15]);
hold on;

% Plot target sector boundary
if mode~=1
    scatter(all_edges_list_pos(:,1),all_edges_list_pos(:,2),'.','MarkerEdgeColor','w','SizeData', 1);
end

% Plot where the BTs are and add a text legend to know where are the BTSs
for b_ = 1:length(eNodeBs)
    scatter(eNodeBs(b_).pos(1),eNodeBs(b_).pos(2),'MarkerEdgeColor','k','MarkerFaceColor','w');
    text(eNodeBs(b_).pos(1)+7*networkPathlossMap.data_res,eNodeBs(b_).pos(2),num2str(b_),'Color','k');
end

xlabel('x pos [m]');
ylabel('y pos [m]');

%% Plot difference between max and 2nd strongest SINR (visualizes cell edge)
if use_subplots
    subplot(2,2,3);
else
    figure;
end
caxis_max = 15;
imagesc(networkPathlossMap.roi_x,networkPathlossMap.roi_y,diff_SINR_dB);

set(gca,'YDir','normal');
title(sprintf('SINR difference (%s). caxis limited to %ddB',plot_type,caxis_max));
colorbar;
caxis([0 caxis_max]);
xlabel('x pos [m]');
ylabel('y pos [m]');
hold on;

% Plot where the BTs are and add a text legend to know where are the BTSs
for b_ = 1:length(eNodeBs)
    scatter(eNodeBs(b_).pos(1),eNodeBs(b_).pos(2),'MarkerEdgeColor','w','MarkerFaceColor','w');
    text(eNodeBs(b_).pos(1)+7*networkPathlossMap.data_res,eNodeBs(b_).pos(2),num2str(b_),'Color','w');
end

%% eNodeB assignment
if use_subplots
    subplot(2,2,4);
else
    figure;
end
imagesc(networkPathlossMap.roi_x,networkPathlossMap.roi_y,sector_assignment);
set(gca,'YDir','normal');
title(sprintf('eNodeB assignment (%s)',plot_type));
colormap('jet');
colorbar;
hold on;

% Plot where the BTs are and add a text legend to know where are the BTSs
for b_ = 1:length(eNodeBs)
    scatter(eNodeBs(b_).pos(1),eNodeBs(b_).pos(2),'MarkerEdgeColor','k','MarkerFaceColor','w');
    text(eNodeBs(b_).pos(1)+6*networkPathlossMap.data_res,eNodeBs(b_).pos(2),num2str(b_));
end

% Plot the center of the cells
for s_ = 1:length(eNodeBs_sectors)
    text(cell_centers(s_,1),cell_centers(s_,2),num2str(s_),'HorizontalAlignment','center','Verticalalignment','middle','Color',0.75*[1 1 1]);
end

xlabel('x pos [m]');
ylabel('y pos [m]');
