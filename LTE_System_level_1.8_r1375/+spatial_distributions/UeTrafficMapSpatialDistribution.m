classdef UeTrafficMapSpatialDistribution < spatial_distributions.networkElementSpatialDistribution
    % Generate spatial user distribution from stored traffic maps
    % Traffic maps contain traffic density - user positions are generated
    % by coin tossing.
    % (c) Martin Taranetz, Josep Colom Ikuno, ITC, 2012
    
    properties
        traffic_map_config
        rescale_factor
    end
    
    methods
        % Class constructor.
        function obj = UeTrafficMapSpatialDistribution(networkPathlossMap,rescale_factor,traffic_map_config)
            % Fill in basic parameters (handled by the superclass constructor)
            obj                    = obj@spatial_distributions.networkElementSpatialDistribution(networkPathlossMap);
            obj.rescale_factor     = rescale_factor;
            obj.traffic_map_config = traffic_map_config;
        end
        
        function user_positions = generate_positions(obj,varargin)
            % Use traffic maps for user generation
            networkPathlossMap = obj.networkPathlossMap;
            
            % OUR TESTCASE :
            % user density per pixel
            % user_density_ = 0.02;
            % traffic_map_in_roi = repmat(user_density_, size(networkMacroscopicPathlossMap.sector_assignment(:,:,1)));
            
            % Using User Density Traffic Maps from MKA
            user_density_traffic_maps = obj.init_user_density_traffic_map;
            % Set undefined values to zero
            user_density_traffic_maps.data(user_density_traffic_maps.data<0) = 0;
            % Now rescale the user density traffic map and calculate a user
            % density/pixel map
            
            % Resize by rescale factor also used for capesso pathlossmaps
            user_density_traffic_maps.data = imresize(user_density_traffic_maps.data , obj.rescale_factor);
            % Cut out the ROI
            pm.roi_x   = networkPathlossMap.roi_x;
            pm.roi_y   = networkPathlossMap.roi_y;
            udtm.roi_x = user_density_traffic_maps.description.roi_x;
            udtm.roi_y = user_density_traffic_maps.description.roi_y;
            % Calculate SW position of Pathlossmaps in pixel within Traffic maps
            pm_pos = LTE_common_pos_to_pixel([pm.roi_x(:,1) pm.roi_y(:,1)],[udtm.roi_x(:,1) udtm.roi_y(:,1)], networkPathlossMap.data_res);
            n_rows = size(networkPathlossMap.pathloss,1);
            n_cols = size(networkPathlossMap.pathloss,2);
            traffic_map_in_roi_pkm = user_density_traffic_maps.data(pm_pos(2):pm_pos(2)+n_rows-1, pm_pos(1):pm_pos(1)+n_cols-1);
            
            % Recalculate the traffic map in the ROI to a users/pixel map
            % Consider both dimensions to be equal:
            x_dim_n                 = user_density_traffic_maps.description.xdim / obj.rescale_factor;
            square_meters_per_pixel = x_dim_n^2;
            traffic_map_in_roi_pp   = traffic_map_in_roi_pkm/(10^6) * square_meters_per_pixel;
            
            % Place users in ROI via coin toss - use users/pixel density as
            % probability for positive outcome of random experiment
            
            % using traffic_map_in_roi as threshold and uniform distribution gives
            % logical matrix of pixel user positions in ROI
            user_matrix_in_roi    = (rand(size(traffic_map_in_roi_pp)) <= traffic_map_in_roi_pp);
            [row,col]             = find(user_matrix_in_roi);
            user_positions_pixels = [col, row];
            user_positions        = LTE_common_pixel_to_pos(user_positions_pixels, obj.networkPathlossMap.coordinate_origin, obj.networkPathlossMap.data_res);
        end
        
        function user_density_traffic_maps = init_user_density_traffic_map(obj)
            % Martin Taranetz, INTHFT 2010
            
            % Function returns user density traffic maps and descriptive data
            % 4 traffic maps (environments) with different penetration losses are available:
            %                           Penetration Loss[dB]
            % o)  DI    Deep Indoor     23
            % o)  ID    Indoor          17
            % o)  IC    Incar           7
            % o)  OD    Outdoor         0
            
            %% Read out UDTM from .bil file and additional information from .HDR file
            % udtm ... user traffic density map
            % MK specific :
            % filename_ = strcat(udtm_file_name,'_',udtm_environment,'_');
            % data_file = strcat(filename_,'focus.bil');
            % hdr_file  = strcat(filename_,'focus.HDR');
            
            traffic_map_upscaling = obj.traffic_map_config.traffic_map_upscaling;
            udtm_folder           = obj.traffic_map_config.udtm_folder;
            udtm_file_name        = obj.traffic_map_config.udtm_filename;
            enable_plotting       = false;
            
            filename_     = udtm_file_name;
            data_file     = strcat(filename_,'.bil');
            hdr_file      = strcat(filename_,'.HDR');
            udtm_file_hdr = fopen(fullfile(udtm_folder, hdr_file) ,'r');
            udtm_file     = fopen(fullfile(udtm_folder, data_file),'r');
            
            % Header file information
            udtm_hdr = textscan(udtm_file_hdr, '%s %f');
            
            udtm.description.NWxmap = udtm_hdr{2}(1);
            udtm.description.NWymap = udtm_hdr{2}(2);
            udtm.description.xdim   = udtm_hdr{2}(3);
            udtm.description.ydim   = udtm_hdr{2}(4);
            udtm.description.ncols  = udtm_hdr{2}(5);
            udtm.description.nrows  = udtm_hdr{2}(6);
            udtm.description.nbits  = udtm_hdr{2}(7);
            udtm.description.nbands = udtm_hdr{2}(8);
            
            udtm.description.SWxmap = udtm.description.NWxmap;                                                % Lower-leftmost corner (x) -> SW
            udtm.description.SWymap = udtm.description.NWymap - udtm.description.ydim*(udtm.description.nrows-1/100); % Lower-leftmost corner (y) -> SW
            
            udtm.description.NExmap = udtm.description.NWxmap + udtm.description.xdim*(udtm.description.ncols-1/100);
            udtm.description.NEymap = udtm.description.NWymap;
            
            udtm.description.SExmap = udtm.description.NExmap;
            udtm.description.SEymap = udtm.description.SWymap;
            
            udtm.description.roi_x = [udtm.description.SWxmap udtm.description.SExmap];
            udtm.description.roi_y = [udtm.description.SWymap udtm.description.NWymap];
            
            % Read the .bil file with the appropriate coding
            if udtm.description.nbits == 32
                udtm_map_t = fread(udtm_file, [udtm.description.ncols, udtm.description.nrows], 'single','l'); % Directly specify to use little endian, just in case...
                udtm_map_t(udtm_map_t==-realmax('single')) = 0; % In order to fill the 'holes' in the traffic map
            else
                error('Precision of User Traffic Denisty Map not supported');
            end
            
            udtm_map_t     = udtm_map_t*traffic_map_upscaling;
            udtm.data      = flipud(udtm_map_t');
            udtm.upscaling = traffic_map_upscaling;
            
            % Plotting
            if enable_plotting
                figure;
                imagesc(udtm.description.roi_x,udtm.description.roi_y, udtm.data);
                hold on;
                set(gca,'YDir','normal');
                xlabel('x pos [m]');
                ylabel('y pos [m]');
                colorbar;
                caxis([0 10]); % Using 'jet' colormap this fits perfectly to MKA legend
                title('User density traffic map (m) [Users/km^2] ');
            end
            
            fclose(udtm_file);
            fclose(udtm_file_hdr);
            
            user_density_traffic_maps = udtm;
        end
    end
    
end