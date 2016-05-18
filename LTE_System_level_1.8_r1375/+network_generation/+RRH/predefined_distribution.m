function predefined_distribution(RRH_config,eNodeBs)
% Distribute RRHs according to predefined positions and azimuths.
% (c) Martin Taranetz, 2013
% www.nt.tuwien.ac.at

RRH_positions = RRH_config.distribution_params.RRH_positions;
RRH_azimuths  = RRH_config.distribution_params.RRH_azimuths;

for b_=1:length(eNodeBs)
    for rrh_=1:size(RRH_positions(:,:,b_),1)
        if rrh_==1
            eNodeBs(b_).RRHs       = network_elements.RRH;
        else
            eNodeBs(b_).RRHs(rrh_) = network_elements.RRH;
        end
        eNodeBs(b_).RRHs(rrh_).pos           = RRH_positions(rrh_,:,b_);
        eNodeBs(b_).RRHs(rrh_).parent_eNodeB = eNodeBs(b_);
        eNodeBs(b_).RRHs(rrh_).azimuth       = RRH_azimuths(rrh_,:,b_);
        eNodeBs(b_).RRHs(rrh_).antenna_type  = RRH_config.antenna_type;
        eNodeBs(b_).RRHs(rrh_).nTX           = RRH_config.nTX;
        
        antennas.antenna.attach_antenna_to_eNodeB(eNodeBs(b_).RRHs(rrh_),RRH_config);
        
        % Set the same macroscopic pathloss as in the parent eNodeB
        % (this could be changed in the future if needed, though)
        eNodeBs(b_).RRHs(rrh_).macroscopic_pathloss_model = ...
            eNodeBs(b_).RRHs(rrh_).parent_eNodeB.macroscopic_pathloss_model;
    end
end