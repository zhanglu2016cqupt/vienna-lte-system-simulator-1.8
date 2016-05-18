classdef constantElementsPerCellSpatialDistribution < spatial_distributions.networkElementSpatialDistribution
    % A constant number of network elements (eNodeBs, HeNBs, UEs) is
    % distributed uniformly in each cell. The coverage area of each cell is
    % given by the calculated sector assignment.
    % (c) Martin Taranetz, Josep Colom Ikuno, ITC, 2012
   
   properties
       elementsPerCell
   end
   
   methods
    function obj = constantElementsPerCellSpatialDistribution(networkPathlossMap, elementsPerCell)
        obj                    = obj@spatial_distributions.networkElementSpatialDistribution(networkPathlossMap);
        obj.elementsPerCell    = elementsPerCell;
    end
   
    function elements_positions = generate_positions(obj,varargin)
            networkPathlossMap = obj.networkPathlossMap;
            sector_surfaces    = networkPathlossMap.sector_sizes;
            number_of_sectors  = size(networkPathlossMap.sector_sizes,2);
            
            % Check whether some sector surfaces are zero. This could be caused by
            % choosing a very small inter-eNodeB distance in relation to the map
            % resolution.
            if sum(sector_surfaces(:)==0) > 0
                warning('Some sector sizes are zero. No UEs will be generated there. Maybe a too big map resolution value?');
            end
            
            % Calculate how many elements per sector
            norm_sector_surface = ones(size(sector_surfaces));
            % Different UE densities per eNodeB sector
            if length(obj.elementsPerCell) > 1
                elements_sector = (strcmp(networkPathlossMap.site_type,'macro')*obj.elementsPerCell(1) + ...
                                   strcmp(networkPathlossMap.site_type,'femto')*obj.elementsPerCell(2))'.*...
                                  (sector_surfaces>0);
            else
                elements_sector = round(norm_sector_surface*obj.elementsPerCell).*(sector_surfaces>0);
            end
            
            % Alternatively use elementsPerCell as threshold probability that single
            % user is generated in sector: elementsPerCell in [0 1]
            % user_sector = rand(size(norm_sector_surface))<elementsPerCell;
            
            % All the possible positions for a given sector
            sector_positions = cell(size(networkPathlossMap.sector_sizes));
            % Where our elements will be positions (for each eNodeB and sector)
            positions_pixels_per_cell = cell(size(networkPathlossMap.sector_sizes));
            
            % Assign random positions to each UE
            for s_idx = 1:number_of_sectors
                if elements_sector(s_idx)~=0
                    [row,col] = find(networkPathlossMap.sector_assignment(:,:)==s_idx);
                    sector_positions{s_idx} = [col,row];
                    positions_pixels_per_cell{s_idx} = sector_positions{s_idx}(randi(size(sector_positions{s_idx},1),[1 elements_sector(s_idx)]),:);
                else
                    sector_positions{s_idx}      = [];
                    positions_pixels_per_cell{s_idx} = [];
                end
            end
            
            elements_positions_pixels_total_Elements = 0;
            for i_=1:length(positions_pixels_per_cell)
                if ~isempty(positions_pixels_per_cell{i_})
                    elements_positions_pixels_total_Elements = elements_positions_pixels_total_Elements + size(positions_pixels_per_cell{i_},1);
                end
            end
            position_pixels = zeros(elements_positions_pixels_total_Elements,2);
            l_=1;
            for i_=1:length(positions_pixels_per_cell)
                for j_=1:size(positions_pixels_per_cell{i_},1)
                    position_pixels(l_,:) = positions_pixels_per_cell{i_}(j_,:);
                    l_=l_+1;
                end
            end
            elements_positions = LTE_common_pixel_to_pos(position_pixels, obj.networkPathlossMap.coordinate_origin, obj.networkPathlossMap.data_res);
    end
   end
end