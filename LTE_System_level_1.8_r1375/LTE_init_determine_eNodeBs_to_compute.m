function eNodeBs_to_compute = LTE_init_determine_eNodeBs_to_compute(LTE_config, eNodeBs)
% Determine the eNodeBs in which users are computed in the simulations
% In a stochastic network, the eNodeB positions are unknown in the
% beginning.
% (c) Martin Taranetz, ITC, 2012

inclusion_radius                   = LTE_config.inter_eNodeB_distance*LTE_config.compute_only_center_users;
eNodeB_sites_positions             = reshape([eNodeBs.pos]',2,[])'; 
eNodeB_sites_distances_from_origin = sqrt(eNodeB_sites_positions(:,1).^2 + eNodeB_sites_positions(:,2).^2);
eNodeB_sites_to_compute_idxs       = eNodeB_sites_distances_from_origin < inclusion_radius;
% Check whether there are eNodeBs within the inclusion region.
if sum(eNodeB_sites_to_compute_idxs > 0)
    eNodeB_sectors_to_compute = [eNodeBs(eNodeB_sites_to_compute_idxs).sectors];
else
    % if no eNodeBs inside the conclusion region -> take all eNodeBs
    eNodeB_sectors_to_compute = [eNodeBs.sectors];
end
eNodeBs_to_compute = [eNodeB_sectors_to_compute.eNodeB_id];