classdef capessoUtils < handle
    % Wraps utilities needed to handle Capesso pathloss maps
    % (c) Josep Colom Ikuno , INTHFT, 2012
    
    properties
    end
    
    methods (Static)
        function [networkMacroscopicPathlossMap DTM_cut] = cut_pathloss_maps_to_sites_ROI(M_Capesso,sites,digital_terrain_model,LTE_config)
            capesso_params        = LTE_config.capesso_params;
            map_resolution        = capesso_params.maps_resolution/capesso_params.rescale_factor;
            enable_debug_plotting = capesso_params.enable_debug_plotting;
            eNodeBs               = [sites.sectors];
            sites_pos             = reshape([sites.pos],2,[]);
            
            if ~LTE_config.manually_set_ROI
                if LTE_config.debug_level>=1
                    fprintf('ROI bounded by the eNodeBs plus increase factor');
                end
                sites_roi_x    = [min(sites_pos(1,:)) max(sites_pos(1,:))];
                sites_roi_y    = [min(sites_pos(2,:)) max(sites_pos(2,:))];
                
                % Increase the ROI according to the ROI increase factor
                sites_roi_x = sites_roi_x + [-1 1].*diff(sites_roi_x)*capesso_params.eNodeB_ROI_increase_factor;
                sites_roi_y = sites_roi_y + [-1 1].*diff(sites_roi_y)*capesso_params.eNodeB_ROI_increase_factor;
            else
                if LTE_config.debug_level>=1
                    fprintf('ROI set manually');
                end
                sites_roi_x = LTE_config.roi_x;
                sites_roi_y = LTE_config.roi_y;
            end
            
            networkMacroscopicPathlossMap          = channel_gain_wrappers.macroscopicPathlossMap;
            networkMacroscopicPathlossMap.data_res = map_resolution ;
            networkMacroscopicPathlossMap.roi_x    = sites_roi_x;
            networkMacroscopicPathlossMap.roi_y    = sites_roi_y;
            networkMacroscopicPathlossMap.name     = 'capesso';

            roi_maximum_pixels = LTE_common_pos_to_pixel( [sites_roi_x(2) sites_roi_y(2)], [sites_roi_x(1) sites_roi_y(1)], map_resolution);
            roi_height_pixels  = roi_maximum_pixels(2);
            roi_width_pixels   = roi_maximum_pixels(1);
            
            networkMacroscopicPathlossMap.pathloss = zeros(roi_height_pixels,roi_width_pixels,length(eNodeBs));
            
            capesso_map_in_roi    = false(1,length(M_Capesso));
            for s_=1:length(M_Capesso)
                roi_margin_x           = (sites_roi_x-M_Capesso(s_).description.roi_x).*[1 -1];
                roi_margin_y           = (sites_roi_y-M_Capesso(s_).description.roi_y).*[1 -1];
                capesso_map_in_roi(s_) = sum(sign([roi_margin_x roi_margin_y])>=0)==4; % Check if the capesso map is in the ROI

                if ~capesso_map_in_roi(s_)
                    error('Capesso map %s not in the ROI defined by the site positions.',M_Capesso(s_).description.filename);
                end
                
                [pathloss_current_site padding_ratio] = capesso.capessoUtils.get_map_in_roi(sites_roi_x,sites_roi_y, M_Capesso(s_), map_resolution, enable_debug_plotting);
                for eNodeB_id = [sites(s_).sectors.eNodeB_id]
                    networkMacroscopicPathlossMap.pathloss(:,:,eNodeB_id) = pathloss_current_site;
                    eNodeBs(eNodeB_id).macroscopic_pathloss_model         = [];
                    % Mapping between s_idx and b_/s_ pair
                    networkMacroscopicPathlossMap.sector_idx_mapping(eNodeB_id,:) = [sites(s_).id eNodeBs(eNodeB_id).id];
                    networkMacroscopicPathlossMap.site_sector_mapping(sites(s_).id,eNodeBs(eNodeB_id).id)  = eNodeB_id;
                end
            end
            [DTM_cut padding_ratio] = capesso.capessoUtils.get_map_in_roi(sites_roi_x,sites_roi_y, digital_terrain_model, map_resolution, enable_debug_plotting);
        end
        
        function [cut_full padding_ratio_p] = get_map_in_roi(roi_x,roi_y, map,map_resolution, enable_debug_plotting)
            % Returns the map in the region of interest of the current pathlossmap.
            the_map_description = map.description;
            out_map_roi_min     = [roi_x(1) roi_y(1)];
            out_map_roi_max     = [roi_x(2) roi_y(2)];
            in_map_roi_min      = [map.description.roi_x(1) map.description.roi_y(1)];
            in_map_roi_max      = [map.description.roi_x(2) map.description.roi_y(2)];
            
            % make it work for bo th pathloss maps and elevation maps
            if isfield(map,'pathloss_map')
                map_data                 = map.pathloss_map;
                map_description_filename = the_map_description.filename;
            else
                map_data = map.data;
                map_description_filename = 'ROI';
            end
            
            % The objective here is to cut a piece of the map that matches the input ROI maps. IF the map would not be big enough, a
            % warning will be shown and just part of the map will be returned, with the rest being padded with the average map value for the cut area
            
            target_roi_min_pix = LTE_common_pos_to_pixel(out_map_roi_min,out_map_roi_min,map_resolution);
            target_roi_max_pix = LTE_common_pos_to_pixel(out_map_roi_max,out_map_roi_min,map_resolution);
            
            out_map_roi_min_pix = LTE_common_pos_to_pixel(out_map_roi_min,in_map_roi_min,map.description.xdim);
            out_map_roi_max_pix = LTE_common_pos_to_pixel(out_map_roi_max,in_map_roi_min,map.description.xdim);
            
            map_pixels_to_cut_x = [out_map_roi_min_pix(1) out_map_roi_max_pix(1)];
            map_pixels_to_cut_y = [out_map_roi_min_pix(2) out_map_roi_max_pix(2)];
            
            % Filter out the values out of the matrix
            map_pixels_to_cut_x_2 = [max(map_pixels_to_cut_x(1),1) min(map_pixels_to_cut_x(2),map.description.ncols)];
            map_pixels_to_cut_y_2 = [max(map_pixels_to_cut_y(1),1) min(map_pixels_to_cut_y(2),map.description.nrows)];
            
            % Place the cut DTM into the matrix
            map_diff_x = map_pixels_to_cut_x_2-map_pixels_to_cut_x;
            map_diff_y = map_pixels_to_cut_y_2-map_pixels_to_cut_y;
            
            % Postions in the cut_full matrix
            map_cut  = map_data(map_pixels_to_cut_y_2(1):map_pixels_to_cut_y_2(2),map_pixels_to_cut_x_2(1):map_pixels_to_cut_x_2(2));
            cut_full = zeros([(diff(map_pixels_to_cut_y)+1) (diff(map_pixels_to_cut_x)+1)]) + mean(map_cut(:));
            cut_full((1+map_diff_y(1)):(end+map_diff_y(2)),(1+map_diff_x(1)):(end+map_diff_x(2))) = map_cut;
            
            padding_ratio_p = (numel(cut_full)-numel(map_cut))/(numel(cut_full))*100;
            
            % Resize output
            cut_full = imresize(cut_full,[target_roi_max_pix(2) target_roi_max_pix(1)]);
            
            % clims = [min(elevation_map.data(:)) max(elevation_map.data(:))];
            % figure; imagesc(DTM_cut_full,clims); set(gca,'YDir','normal'); axis equal
            % figure; imagesc(DTM_cut,clims); set(gca,'YDir','normal'); axis equal
            % figure; imagesc(elevation_map.data,clims); set(gca,'YDir','normal'); axis equal
            
            if padding_ratio_p>0
                warning('Input map did not cover the whole output map. Padding was added. %3.2f%% of the area corresponding to the output map was padded at %3.1fm\n',padding_ratio_p,mean(map_cut(:)));
            end
            
            % figure;imagesc(pm_description.roi_x,pm_description.roi_y,cut_elevation_map);
            
            if enable_debug_plotting
                % Plot whole elevation map
                figure;
                imagesc(map.description.roi_x,map.description.roi_y,map_data);
                color_axis = caxis;
                grid on;
                hold on;
                set(gca,'YDir','normal');
                xlabel('x pos [m]');
                ylabel('y pos [m]');
                % Mark the part of the elevation map that corresponds to the pathloss file
                roi_width  = the_map_description.roi_x(2)-the_map_description.roi_x(1);
                roi_height = the_map_description.roi_y(2)-the_map_description.roi_y(1);
                rectangle('Position',[out_map_roi_min(1), out_map_roi_min(2), out_map_roi_max(1)-out_map_roi_min(1), out_map_roi_max(2)-out_map_roi_min(2)]);
                hold off;
                colorbar;
                title('Original map');
                % Plot only area corresponding to this pathloss map
                figure;
                imagesc([out_map_roi_min(1) out_map_roi_max(1)],[out_map_roi_min(2) out_map_roi_max(2)],cut_full);
                caxis(color_axis);
                grid on;
                set(gca,'YDir','normal');
                colorbar;
                title(sprintf('Map cut. %s area',strrep(map_description_filename,'_','\_')));
                xlabel('x pos [m]');
                ylabel('y pos [m]');
            end
        end
    end
    
end

