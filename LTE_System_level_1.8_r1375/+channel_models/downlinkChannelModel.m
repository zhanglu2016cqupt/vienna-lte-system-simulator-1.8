classdef downlinkChannelModel < handle
% Represents the downlink channel model that a specific user possesses. Each UE
% instance will have its own specific channel model.
% (c) Josep Colom Ikuno, INTHFT, 2008

   properties
       % These could actually be maps or an implementation that directly
       % calculates every time it is invoked.
       macroscopic_pathloss_model_is_set = false;
       macroscopic_pathloss_model
       shadow_fading_model_is_set = false;
       shadow_fading_model

       % Primary and secondary fast fading models (traces)
       fast_fading_model
       secondary_fast_fading_models
       
       % All eNodeBs
       eNodeBs
       
       % User to which this channel model is attached to
       attached_UE
       
       % This variables model the downlink signaling and the data that was
       % transmitted
       % RB_grid (retrieved via a function call)
       
   end

   methods
       % class constructor
       function obj = downlinkChannelModel(aUE)
           obj.attached_UE = aUE;
       end
       % Returns the macroscopic pathloss in dB between the given user's
       % position and his eNodeB sector. Returns 0 in case no model is specified.
       function [pathloss, is_set] = macroscopic_pathloss(obj)
           % Get eNodeB id
           attached_eNodeB = obj.attached_UE.attached_eNodeB;
           eNodeB_id       = attached_eNodeB.eNodeB_id;
           pos             = obj.attached_UE.pos;
           
           % Add RRHs (if applicable)
           if ~isempty(attached_eNodeB.RRHs)
               eNodeB_id = [eNodeB_id [attached_eNodeB.RRHs.id]];
           end
           
           % Now get the pathloss
           pathloss = reshape(obj.macroscopic_pathloss_model.get_pathloss_eNodeB(pos,eNodeB_id),1,[]);
           is_set   = true;
       end
       
       % Returns the macroscopic pathloss in dB between the given user's
       % position and a given eNodeB.
       function [pathloss, is_set] = interfering_macroscopic_pathloss(obj,eNodeB_id)
           % Get eNodeB id
           attached_eNodeB = obj.eNodeBs(eNodeB_id);
           eNodeB_id       = attached_eNodeB.eNodeB_id;
           pos             = obj.attached_UE.pos;
           
           % Add RRHs (if applicable)
           if ~isempty(attached_eNodeB.RRHs)
               eNodeB_id = [eNodeB_id [attached_eNodeB.RRHs.id]];
           end
           
           % Now get the pathloss
           pathloss = reshape(obj.macroscopic_pathloss_model.get_pathloss_eNodeB(pos,eNodeB_id),1,[]);
           is_set   = true;
       end
       
       % Returns the shadow fading pathloss in dB between the given user's
       % position and his eNodeB sector. Returns 0 in case no model is specified (for
       % example, when Odyssey data would be used).
       function [pathloss, is_set] = shadow_fading_pathloss(obj)
           if ~obj.shadow_fading_model_is_set
               pathloss = 1;
               is_set   = false;
               return
           else
               % Get eNodeB id
               attached_eNodeB = obj.attached_UE.attached_eNodeB;

               % Get eNodeB id and sector number
               if isa(obj.shadow_fading_model,'channel_gain_wrappers.shadowFadingMapClaussen')&&~obj.shadow_fading_model.oneMapPerSite
                   % Take the cell ID
                   shadow_map_id = attached_eNodeB.eNodeB_id;
               else
                   % Default
                   shadow_map_id = obj.attached_UE.attached_site.id;
               end
               
               % Add RRHs (if applicable)
               if ~isempty(attached_eNodeB.RRHs)
                   RRH_maps_ids  = [attached_eNodeB.RRHs.site_id];
                   if isempty(RRH_maps_ids) % in the case of dummy maps
                       RRH_maps_ids = zeros(1,length(attached_eNodeB.RRHs));
                   end
                   shadow_map_id = [shadow_map_id RRH_maps_ids];
               end
               
               pos      = obj.attached_UE.pos;
               pathloss = reshape(obj.shadow_fading_model.get_pathloss(pos,shadow_map_id),1,[]);
               is_set   = true;
           end
       end
       
       % Returns the shadow fading pathloss in dB between the given user's 
       % position and a given eNodeB. Returns 0 in case no model is specified (for
       % example, when Odyssey data would be used).
       function [pathloss, is_set] = interfering_shadow_fading_pathloss(obj,interferingSiteIds,interfering_eNodeB_ids)
           if ~obj.shadow_fading_model_is_set
               pathloss = ones([length(interferingSiteIds) 1]); % Column vector
               is_set   = false;
               return
           else
               % Get eNodeB id and sector number
               if isa(obj.shadow_fading_model,'channel_gain_wrappers.shadowFadingMapClaussen')&&~obj.shadow_fading_model.oneMapPerSite
                   % Take the cell ID
                   interferingIds = interfering_eNodeB_ids;
               else
                   % Default
                   interferingIds = interferingSiteIds;
               end
               pos      = obj.attached_UE.pos;
               pathloss = reshape(obj.shadow_fading_model.get_pathloss(pos,interferingIds),[],1);
               is_set   = true;
           end
       end
       
       % Returns the RB_grid so this UE can know what belongs to him
       function the_RB_grid = RB_grid(obj)
           if ~isempty(obj.attached_UE.attached_eNodeB)
               the_RB_grid = obj.attached_UE.attached_eNodeB.RB_grid;
           else
               the_RB_grid = [];
           end
       end
       % Set a macroscopic pathloss model
       function set_macroscopic_pathloss_model(obj,macroscopic_pathloss_model)
           obj.macroscopic_pathloss_model        = macroscopic_pathloss_model;
           obj.macroscopic_pathloss_model_is_set = true;
       end
       % Set a shadow fading model
       function set_shadow_fading_model(obj,shadow_fading_model)
           obj.shadow_fading_model        = shadow_fading_model;
           obj.shadow_fading_model_is_set = true;
       end
   end
end 
