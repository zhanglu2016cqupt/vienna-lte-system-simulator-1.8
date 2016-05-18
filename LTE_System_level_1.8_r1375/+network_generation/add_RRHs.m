function add_RRHs(LTE_config,eNodeBs,networkPathlossMap )
% Adds Remote Radio Heads to the actual simulation
% (c) Josep Colom Ikuno ITC, 2013
% www.nt.tuwien.ac.at

RRH_config = LTE_config.RRH;

% Create the RRHs attached to their respective eNodeBs
switch RRH_config.distribution_type
    case 'on arch'
        network_generation.RRH.on_arch_distribution(RRH_config,eNodeBs);
    case 'predefined'
        network_generation.RRH.predefined_distribution(RRH_config,eNodeBs);
    otherwise
        error('RRH distribution not defined');
end

RRHs           = [eNodeBs.RRHs];
last_eNodeB_id =  eNodeBs(end).eNodeB_id;
for rrh_=1:length(RRHs)
    RRHs(rrh_).id = rrh_+last_eNodeB_id;
end

%% Calculate final pathloss. With and without minimum coupling loss
if LTE_config.debug_level>=1
    fprintf('Creating RRH pathloss map\n');
end

% Fill in pathloss information
macroscopic_pathloss_models.generalPathlossModel.calculate_pathloss_maps(LTE_config,RRHs,networkPathlossMap);
