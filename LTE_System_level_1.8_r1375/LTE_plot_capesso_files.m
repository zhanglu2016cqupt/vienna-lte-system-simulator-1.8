function LTE_plot_capesso_files(LTE_config,eNodeBs,pathlossmaps, networkMacroscopicPathlossMap, elevation_map, capesso_params)
% Plots the loaded Capesso files based on the data stored on the struct
%
% (c) Josep Colom Ikuno, Martin Taranetz INTHFT, 2010

map_res      = capesso_params.maps_resolution;
first_figure = LTE_config.plots.capesso_maps_begin+0;

figure(first_figure);
hold on;
% Set a common pathloss scale and plot the pathloss files in the overall grid
for i_=1:length(pathlossmaps)
    if i_==1
        c_axis  = [min(pathlossmaps(i_).pathloss_map(:))  max(pathlossmaps(i_).pathloss_map(:))];
        c_axis2 = [min(elevation_map.data(:))             max(elevation_map.data(:))];
    else
        c_axis  = [min([c_axis(1);  pathlossmaps(i_).pathloss_map(:)])  max([c_axis(2);  pathlossmaps(i_).pathloss_map(:)])];
        c_axis2 = [min([c_axis2(1); elevation_map.data(:)])             max([c_axis2(2); elevation_map.data(:)])];
    end
    rectangle_pos = [pathlossmaps(i_).description.roi_x(1) pathlossmaps(i_).description.roi_y(1) diff(pathlossmaps(i_).description.roi_x) diff(pathlossmaps(i_).description.roi_y)];
    rectangle('Position',rectangle_pos);
    text(rectangle_pos(1),rectangle_pos(2),num2str(i_),'FontSize',10);
end

xlimits = xlim;
ylimits = ylim;
ROI_increase_factor = 0.05;
xlimits = ROI_increase_factor*abs(xlimits(1)-xlimits(2))*[-1 1] + xlimits;
ylimits = ROI_increase_factor*abs(ylimits(1)-ylimits(2))*[-1 1] + ylimits;
xlim(xlimits);
ylim(ylimits);
grid on

% Draw a suggested ROI
sug_roi_min    = [networkMacroscopicPathlossMap.roi_x(1) networkMacroscopicPathlossMap.roi_y(1)];
sug_roi_width  = diff(networkMacroscopicPathlossMap.roi_x);
sug_roi_height = diff(networkMacroscopicPathlossMap.roi_y);
rectangle_pos  = [sug_roi_min(1), sug_roi_min(2), sug_roi_width, sug_roi_height];
rectangle('Position',rectangle_pos,'EdgeColor','red');
title('Loaded pathloss maps: black, ROI: red');
hold off;

% Plot pathloss for each capesso file
pathloss_figure_idx = LTE_config.plots.capesso_maps_begin+1;
for i_=1:length(pathlossmaps)
    pathlossmap = pathlossmaps(i_);
    cleaned_filename = strrep(pathlossmap.description.filename,'#2F','/');
    
    % Plot pathloss
    figure_capesso_los_file = figure(pathloss_figure_idx+i_-1);
    colormap('default');
    
    % Pathloss map (1st subplot)
    axes1 = subplot(3,length(eNodeBs(i_).sectors),1,'Parent',figure_capesso_los_file);
    %imagesc(x_axis,y_axis,-pathlossmap.pathloss_map);
    imagesc(pathlossmap.description.roi_x,pathlossmap.description.roi_y,-pathlossmap.pathloss_map);
    set(gca,'YDir','normal');
    caxis(axes1,-[c_axis(2) c_axis(1)]);
    colorbar;
    grid(axes1,'on');
    title(axes1,sprintf('%d-Capesso: %s ',i_,strrep(cleaned_filename,'_','\_')));
    hold(axes1,'on');
    for b_=1:length(eNodeBs)
        [faces verts] = get_patch_data(eNodeBs(b_),pathlossmap.description.xdim);
        patch('Faces',faces,'Vertices',verts,'EdgeColor','k','LineWidth',2);
        %set(p,'EdgeColor','none');
    end
    xlim(axes1,networkMacroscopicPathlossMap.roi_x);
    ylim(axes1,networkMacroscopicPathlossMap.roi_y);
    hold(axes1,'off');
    
    % Elevation map (2nd subplot)
    axes2 = subplot(3,length(eNodeBs(i_).sectors),2,'Parent',figure_capesso_los_file);
    imagesc(elevation_map.description.roi_x,elevation_map.description.roi_y,elevation_map.data);
    set(gca,'YDir','normal');
    caxis(axes2,[c_axis2(1) c_axis2(2)]);
    grid(axes2,'on');
    hold(axes2,'on');
    for b_=1:length(eNodeBs)
        [faces verts] = get_patch_data(eNodeBs(b_),pathlossmap.description.xdim);
        patch('Faces',faces,'Vertices',verts,'EdgeColor','k','LineWidth',2);
        %set(p,'EdgeColor','none');
    end
    hold(axes2,'off');
    xlim(axes2,networkMacroscopicPathlossMap.roi_x);
    ylim(axes2,networkMacroscopicPathlossMap.roi_y);
    colorbar;
    title(axes2,sprintf('%d-Elevation map (m). %s area',i_,strrep(cleaned_filename,'_','\_')));
    
    % eNodeB position (3rd subplot)
    axes3 = subplot(3,length(eNodeBs(i_).sectors),3,'Parent',figure_capesso_los_file);
    hold(axes3,'on');
    for b_=1:length(eNodeBs)
        % Plot a line that tells where the antennas are pointing
        vector_length = 80;
        origin = eNodeBs(b_).pos;
        for s_=1:length(eNodeBs(b_).sectors)
            angle = wrapTo360(-eNodeBs(b_).sectors(s_).azimuth+90);
            vector = vector_length*[ cosd(angle) sind(angle) ];
            destiny = vector + origin;
            
            plot([origin(1) destiny(1)],[origin(2) destiny(2)]);
        end
        % Plot the eNodeBs
        scatter(axes3,eNodeBs(b_).pos(1),eNodeBs(b_).pos(2),'Marker','o','MarkerFaceColor','red','MarkerEdgeColor','black');
        text(eNodeBs(b_).pos(1)+map_res*5,eNodeBs(b_).pos(2),num2str(eNodeBs(b_).id));
    end
    hold(axes3,'off');
    grid(axes3,'on');
    xlim(axes3,networkMacroscopicPathlossMap.roi_x);
    ylim(axes3,networkMacroscopicPathlossMap.roi_y);
    colorbar; % This is just so the size of the plot is the same as the one above
    title(axes3,sprintf('eNodeB posititions'));
    
    % Plot sector antenna gain (assuming 3 sectors!!)
    % and pathlossmap with applied antenna gain in ROI for every sector

    for s_ = 1:length(eNodeBs(b_).sectors)
        s_idx = networkMacroscopicPathlossMap.site_sector_mapping(i_,s_);
        
        % Plot pathlossmaps with applied antenna gain in ROI
        axes1 = subplot(3,length(eNodeBs(i_).sectors),length(eNodeBs(i_).sectors)+s_,'Parent',figure_capesso_los_file);
        imagesc(networkMacroscopicPathlossMap.roi_x, networkMacroscopicPathlossMap.roi_y, -networkMacroscopicPathlossMap.pathloss(:,:,s_idx));
        set(gca,'YDir','normal');
        caxis(axes1,-[c_axis(2) c_axis(1)]);
        colorbar;
        title(strrep(sprintf('Pathloss map (Sec%d)',s_),'_','\_'));
        
        % Plot antenna gain patterns using polar plot for each sector without respecting the mechanical downtilt
        if capesso_params.plot_antenna_gain_patterns > 0
            sec_electrical_downtilt = eNodeBs(i_).sectors(s_).electrical_downtilt;
            sec_antenna             = eNodeBs(i_).sectors(s_).antenna;
            
            if find(sec_antenna.electrical_tilt == sec_electrical_downtilt)
                index_ = find(sec_antenna.electrical_tilt == sec_electrical_downtilt, 1, 'first');
            else
                error('Gain pattern for electrical tilt of %f° not available !\n', sec_electrical_downtilt);
                index_ = 0;
            end
            
            if index_  > 0
                axes4 = subplot(length(eNodeBs(i_).sectors),length(eNodeBs(i_).sectors),2*length(eNodeBs(i_).sectors)+s_,'Parent',figure_capesso_los_file);
                plot_tilt = sec_electrical_downtilt; % Plot with the electrical tilt from the sector. The effect of the mechanical tilt is not shown in the polar plot, but it is applied when calculating the pathloss map.
                data_limits = [-15 0 3];
                [hor_degrees hor_gain ver_degrees ver_gain max_gain] = sec_antenna.gain_patterns(plot_tilt);
                utils.miscUtils.polar2(hor_degrees/180*pi, hor_gain-max_gain, data_limits,'blue');
                hold(axes4,'all');
                utils.miscUtils.polar2(ver_degrees/180*pi, ver_gain-max_gain, data_limits,'red');
                title(axes4,sprintf('%d antenna. %3.0f° electrical tilt',sec_antenna.antenna_type,plot_tilt));
                hold(axes4,'off');
            end
            
        end
    end
end

function [faces verts] = get_patch_data(an_eNodeB,data_res)
% Processes the position data from an eNodeB to generate data which can be fed to the "patch" function.
num_sectors = length(an_eNodeB.sectors);
verts = zeros(num_sectors,2);
faces = zeros(num_sectors,2);
pos = an_eNodeB.pos;
vector_length = 1*data_res;

for s_=1:num_sectors
    angle = wrapTo360(-an_eNodeB.sectors(s_).azimuth+90);
    vector = vector_length*[ cosd(angle) sind(angle) ];
    verts(s_,:) = vector + pos;
    if s_>1
        faces(s_-1,:) = [s_-1 s_];
    end
end
faces(s_,:) = [s_ 1];
