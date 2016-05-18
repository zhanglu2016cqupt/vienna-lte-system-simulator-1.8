classdef networkElementSpatialDistribution < handle
    % This class generates positions for network elements (eNodeBs, HeNBs, UEs)
    % according to a given spatial distribution
    % (c) Martin Taranetz, Josep Colom Ikuno, ITC, 2012
    
    properties
        networkPathlossMap
    end
    
    methods
        function obj = networkElementSpatialDistribution(networkPathlossMap)
            obj.networkPathlossMap = networkPathlossMap;
        end
    end
    
    methods(Abstract)
        elements_positions = generate_positions(obj,varargin)
    end
end