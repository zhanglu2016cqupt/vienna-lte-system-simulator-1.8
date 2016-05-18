function [eNodeBs_plus_femtocells, eNodeBs_sectors_plus_femtocells ] = ...
    add_femtocells(LTE_config,sites,eNodeBs,networkPathlossMap )
% Adds femtocells to the actual pathloss map
% (c) Josep Colom Ikuno, Martin Taranetz, ITC, 2012
% www.nt.tuwien.ac.at

% To place fixed number of femtos within eNodeB sectors, eNodeB sector assignment has to be calculated a priori
% Note: Assignment does not take into account shadow fading.
aPriori_sector_assignment = 0; 
switch LTE_config.femtocells_config.spatial_distribution
    case 'homogenous density'
        spatial_distribution_model = spatial_distributions.homogeneousSmallcellSpatialDistribution(networkPathlossMap, LTE_config.femtocells_config.femtocells_per_km2);
    case 'constant femtos per cell'
        % For femtocell placement calculate sector-assignment a-priori.
        [~,~,networkPathlossMap.sector_assignment,~,~,networkPathlossMap.sector_sizes,~] = LTE_common_calculate_cell_capacity(LTE_config, networkPathlossMap, sites, eNodeBs);
        spatial_distribution_model = spatial_distributions.constantElementsPerCellSpatialDistribution(networkPathlossMap, LTE_config.femtocells_config.femtos_per_cell);
        aPriori_sector_assignment = 1;
    otherwise
        error('Spatial distribution "%s" not supported',LTE_config.femtocells_config.spatial_distribution);
end

femtocell_pos = spatial_distribution_model.generate_positions;

% Clean up a-priori calculated eNodeB sector assignment
if aPriori_sector_assignment
    networkPathlossMap.sector_assignment = [];
    networkPathlossMap.sector_sizes      = [];
end
    
for f_=1:size(femtocell_pos,1)
    % Create the femtocell sites
    femto_sites(f_)           = network_elements.eNodeB;
    femto_sites(f_).id        = length(sites)+f_;
    femto_sites(f_).pos       = femtocell_pos(f_,:);
    femto_sites(f_).site_type = 'femto';
    
    % Create the femtocell eNodeBs
    femto_sites(f_).sectors               = network_elements.eNodeB_sector;
    femto_sites(f_).sectors.parent_eNodeB = femto_sites(f_);
    femto_sites(f_).sectors.id            = 1;
    femto_sites(f_).sectors.azimuth       = 0; % Omnidirectional antenna
    femto_sites(f_).sectors.max_power     = LTE_config.femtocells_config.tx_power_W;
    femto_sites(f_).sectors.antenna_type  = 'omnidirectional';
    femto_sites(f_).sectors.nTX           = LTE_config.nTX;
    femto_sites(f_).sectors.eNodeB_id     = length(eNodeBs) + f_;
    femto_sites(f_).sectors.antenna       = antennas.omnidirectionalAntenna;
    
    femto_sectors(f_) = femto_sites(f_).sectors;
    
    % Create the macroscopic pahloss model that will be used. Use the same one as in the macro sites (may change in the future)
    if LTE_config.debug_level>=1
        fprintf('Femto Site %d: ',f_);
    end
    
    femto_sectors(f_).macroscopic_pathloss_model = macroscopic_pathloss_models.generalPathlossModel.generateMacroscopicPathlossModel(...
            LTE_config,...
            LTE_config.femtocells_config.macroscopic_pathloss_model,...
            LTE_config.frequency,...
            LTE_config.femtocells_config.macroscopic_pathloss_model_settings);
end

if size(femtocell_pos,1)>0
    % There are femtocells
    eNodeBs_plus_femtocells         = [ sites femto_sites ];
    eNodeBs_sectors_plus_femtocells = [ eNodeBs femto_sectors ];
else
    % No femtocells. Return the same values
    eNodeBs_plus_femtocells            = sites;
    eNodeBs_sectors_plus_femtocells    = eNodeBs;
    return
end

%% Calculate final pathloss
if LTE_config.debug_level>=1
    fprintf('Creating femtocell pathloss map\n');
end
% Preallocate to make it faster for the following function (but not 100% necessary)
networkPathlossMap.pathloss(:,:,(length(eNodeBs)+1):(length(eNodeBs)+length(femto_sites))) = 1;
macroscopic_pathloss_models.generalPathlossModel.calculate_pathloss_maps(LTE_config,femto_sites,networkPathlossMap);

%% Fill in pathloss data in the pathlossMap
all_pathloss_models      = [femto_sectors.macroscopic_pathloss_model];
all_pathloss_model_names = {all_pathloss_models.name};

networkPathlossMap.name = {networkPathlossMap.name{1:end} all_pathloss_model_names{:}};
