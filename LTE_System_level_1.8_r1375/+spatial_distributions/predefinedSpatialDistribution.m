classdef predefinedSpatialDistribution < spatial_distributions.networkElementSpatialDistribution 
    % Places positions according to a predefined list of positions (e.g.,
    % from  trace).
    % (c) Josep Colom Ikuno, ITC, 2012
    
    properties
        UE_positions
    end
    
    methods
        % Class constructor.
        function obj = predefinedSpatialDistribution(networkPathlossMap, UE_positions)
            obj              = obj@spatial_distributions.networkElementSpatialDistribution(networkPathlossMap);
            obj.UE_positions = UE_positions;
        end
        
        function user_position = generate_positions(obj)
            user_position = obj.UE_positions;
        end
    end
    
end

