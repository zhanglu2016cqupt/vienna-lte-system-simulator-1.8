classdef kmlUtils
    % Encapsulates a couple of methods needed for the manipulatio and
    % creation of Google-Earth-compatible KML files containing simulation
    % results
    %
    % (c) Josep Colom Ikuno, INTHFT, 2011
    properties
        map_markers_folder = './+utils/map_markers/';
        marker_radius = 12;
        radius = 12; % It should be an even number
        
        colors 
    end
    
    methods
        function obj = kmlUtils
            obj.colors.yellow = [1 1 0];
            obj.colors.magente = [1 0 1];
            obj.colors.cyan = [0 1 1];
            obj.colors.red = [1 0 0];
            obj.colors.green = [0 1 0];
            obj.colors.blue = [0 0 1];
            obj.colors.white = [1 1 1];
            obj.colors.black = [0 0 0];
        end
        
        function writePositionsToKmlFile(obj,latitude,longitude,file)
            % Check existence
            dot_filename = fullfile(obj.map_markers_folder,'dot.png');
            if ~exist(dot_filename,'file')
                obj.write_dot('dot',obj.colors.black);
            end
            % Write KML file
            geostruct_pos = struct('Lat',num2cell(latitude),'Lon',num2cell(longitude),'Geometry','Point');
            kmlwrite(file, geostruct_pos,'Icon',dot_filename,'IconScale',0.4);
        end
        
        function writeValuesToKmlFile(obj,latitude,longitude,colormap_idxs_to_plot,file,dataset_name,colormap_rgb)
            % Overwrite dots (for safety)
            markers_folder = obj.write_dots_colormap(dataset_name,colormap_rgb);
            
            icon_list = cell(length(colormap_idxs_to_plot),1);
            for i_=1:length(colormap_idxs_to_plot)
                icon_list{i_} = fullfile(markers_folder,sprintf('%d.png',colormap_idxs_to_plot(i_)));
            end

            % Write KML file
            geostruct_pos = struct('Lat',num2cell(latitude),'Lon',num2cell(longitude),'Geometry','Point');
            kmlwrite(file, geostruct_pos,'Icon',icon_list,'IconScale',0.2);
        end
        
        function write_dot(obj,filename,color_RGB)
            circle = obj.create_circle;
            alpha  = double(circle);
            alpha(alpha~=0) = 1;
            dot    = zeros([size(circle) 3]);
            for i_=1:3
                dot(:,:,i_)  = color_RGB(i_);
            end
            imwrite(dot,fullfile(obj.map_markers_folder,[filename,'.png']),'png','Alpha',alpha);
        end
        
        function markers_folder = write_dots_colormap(obj,colormap_name,the_colormap_rgb)
            markers_folder = fullfile(obj.map_markers_folder,colormap_name);
            if ~exist(markers_folder,'dir')
                mkdir(markers_folder);
            end
            
            % Overwrite any markers there to play it safe
            N_markers = size(the_colormap_rgb,1);
            for i_=1:N_markers
                obj.write_dot(fullfile(colormap_name,sprintf('%d',i_)),the_colormap_rgb(i_,:));
            end
        end
        
        function circle = create_circle(obj)
            % Returns a circle
            radius = obj.radius;
            if mod(radius,2)
                radius = radius-1; % Force even radius
            end
            center = [[radius+1 radius+1]];
            
            idxs   = 1:(radius*2+1);
            circle_pos_x = repmat(idxs,[length(idxs) 1]);
            circle_pos_y = repmat(idxs',[1 length(idxs)]);
            circle = (circle_pos_x-center(1)).^2 + (circle_pos_y-center(2)).^2 <= (radius+0.25)^2;
        end
    end
    
end