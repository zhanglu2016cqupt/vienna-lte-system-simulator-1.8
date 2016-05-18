classdef constantElementsPerRoiSpatialDistribution < spatial_distributions.networkElementSpatialDistribution
    % A constant number of network elements (eNodeBs, HeNBs, UEs) is
    % distributed uniformly in each cell. The coverage area of each cell is
    % given by the calculated sector assignment.
    % (c) Martin Taranetz, Josep Colom Ikuno, ITC, 2012
   
   properties
       nUEs
   end
   
   methods
       function obj = constantElementsPerRoiSpatialDistribution(networkPathlossMap, nUEs)
           obj      = obj@spatial_distributions.networkElementSpatialDistribution(networkPathlossMap);
           obj.nUEs = nUEs;
       end
       
       function elements_positions = generate_positions(obj,varargin)
           roi_x = obj.networkPathlossMap.roi_x;
           roi_y = obj.networkPathlossMap.roi_y;
           
           ROI_span = [(roi_x(2)-roi_x(1)) (roi_y(2)-roi_y(1))];
           elements_positions = rand(obj.nUEs,2).*ROI_span(ones(1,obj.nUEs),:)+ [roi_x(ones(obj.nUEs,1),1) roi_y(ones(obj.nUEs,1),1)];
       end
   end
end