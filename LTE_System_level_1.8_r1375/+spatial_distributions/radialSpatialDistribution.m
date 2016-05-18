classdef radialSpatialDistribution < spatial_distributions.networkElementSpatialDistribution 
    % Places positions around each site
    % (c) Josep Colom Ikuno, ITC, 2012
    
    properties
        radii
        nUEs
        sites
        overlapUEs
    end
    
    methods
        % Class constructor.
        function obj = radialSpatialDistribution(networkPathlossMap,sites,radii,nUEs,overlapUEs)
            obj = obj@spatial_distributions.networkElementSpatialDistribution(networkPathlossMap);
            if length(nUEs)==1 && length(radii)>1
                nUEs = nUEs(ones(1,length(radii)));
            end
            obj.radii = radii;
            obj.nUEs  = nUEs;
            obj.sites = sites;
            obj.overlapUEs = overlapUEs;
        end
        
        function user_position = generate_positions(obj)
            site_pos           = reshape([obj.sites.pos],2,[]).';
            totalUEs           = sum(obj.nUEs);
            user_position_diff = zeros(totalUEs,2);
            cumsum_nUEs        = cumsum(obj.nUEs);
            for radiusIdx=1:length(obj.radii)
                if ~obj.overlapUEs
                    degrees     = 0:360/obj.nUEs(radiusIdx):359;
                else
                    degrees = zeros(1,obj.nUEs(radiusIdx));
                end
                UE_pos_diff = obj.radii(radiusIdx)*[cosd(degrees.') sind(degrees.')];
                user_position_diff((cumsum_nUEs(radiusIdx)-obj.nUEs(radiusIdx)+1):cumsum_nUEs(radiusIdx),:) = UE_pos_diff;
            end
            user_position = kron(site_pos,ones(totalUEs,1))+repmat(user_position_diff,[length(obj.sites) 1]);
        end
    end
    
end

