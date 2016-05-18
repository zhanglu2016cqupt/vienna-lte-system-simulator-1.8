classdef homogeneousSmallcellSpatialDistribution < spatial_distributions.networkElementSpatialDistribution
    % Homogeneous spatial distribution of network elements (HeNBs, Small Cells, Microcells, etc.). Assumes an already existing macrocell network and acts as overlay.
    % The ROI is the size of the networkPathlossMap in this case.
    % (c) Martin Taranetz, Josep Colom Ikuno, ITC, 2012
    
    properties
        elements_per_km2
    end
    
    methods
         % Class constructor.
       function obj = homogeneousSmallcellSpatialDistribution(networkPathlossMap, elements_per_km2)
           obj                    = obj@spatial_distributions.networkElementSpatialDistribution(networkPathlossMap);
           obj.elements_per_km2   = elements_per_km2;
       end
       
       % Generate position [m]
       function elements_positions = generate_positions(obj,varargin)
           networkPathlossMap = obj.networkPathlossMap;
           pixel_size_km2     = (networkPathlossMap.data_res^2)/1e6;
           avg_elements_pixel = obj.elements_per_km2*pixel_size_km2;
           element_placement  = rand([size(networkPathlossMap.pathloss,1) size(networkPathlossMap.pathloss,2)])<=avg_elements_pixel;
           [row,col]          = find(element_placement);
           positions_pixels   = [col,row];
           elements_positions = LTE_common_pixel_to_pos(positions_pixels, obj.networkPathlossMap.coordinate_origin, obj.networkPathlossMap.data_res);
       end
    end 
end