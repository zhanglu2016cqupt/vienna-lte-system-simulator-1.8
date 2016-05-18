function on_arch_distribution(RRH_config,eNodeBs)
    % A distribution of RRHs on arch in the radiating direction of the eNodeB
    %
    % (c) Josep Colom Ikuno ITC, 2013
    % www.nt.tuwien.ac.at
    
    RRHs_per_eNodeB = RRH_config.distribution_params.nRRH_per_eNodeB;
    radius_m        = RRH_config.distribution_params.arch_radius_m;
    arch_deg        = RRH_config.distribution_params.arch_width_deg;
    pos_RRHs_deg    = linspace(0,arch_deg-1,RRHs_per_eNodeB);
    
    for b_=1:length(eNodeBs)
        % The azimuth angle in the eNodeBs is as in a clock, 0° being 90°
        % in cartesian and 30° 90-30° (and so on). We need to translate it
        % to the "normal" nomenclature
        azimuth_deg = 90-eNodeBs(b_).azimuth;
        eNodeB_pos  = eNodeBs(b_).parent_eNodeB.pos;
        
        RRH_azimuths      = azimuth_deg-arch_deg/2+pos_RRHs_deg;
        RRH_pos_offsets_m = [cosd(RRH_azimuths(:)) sind(RRH_azimuths(:))] * radius_m;
        RRH_pos_m         = eNodeB_pos(ones(1,RRHs_per_eNodeB),:) + RRH_pos_offsets_m;
        
        for rrh_=1:RRHs_per_eNodeB
            if rrh_==1
                eNodeBs(b_).RRHs       = network_elements.RRH;
            else
                eNodeBs(b_).RRHs(rrh_) = network_elements.RRH;
            end
            eNodeBs(b_).RRHs(rrh_).pos           = RRH_pos_m(rrh_,:);
            eNodeBs(b_).RRHs(rrh_).parent_eNodeB = eNodeBs(b_);
            eNodeBs(b_).RRHs(rrh_).azimuth       = RRH_azimuths(rrh_);
            eNodeBs(b_).RRHs(rrh_).antenna_type  = RRH_config.antenna_type;
            eNodeBs(b_).RRHs(rrh_).nTX           = RRH_config.nTX;
            
            antennas.antenna.attach_antenna_to_eNodeB(eNodeBs(b_).RRHs(rrh_),RRH_config);
            
            % Set the same macroscopic pathloss as in the parent eNodeB
            % (this could be changed in the future if needed, though)
            eNodeBs(b_).RRHs(rrh_).macroscopic_pathloss_model = ...
                eNodeBs(b_).RRHs(rrh_).parent_eNodeB.macroscopic_pathloss_model;
        end
    end
end
