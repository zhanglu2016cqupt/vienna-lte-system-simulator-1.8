function [eNodeBs_sites,eNodeBs,networkMacroscopicPathlossMap] = fixed_pathloss_network(LTE_config)
% A set of eNodeBs with a fixed, predetermined pathloss independent of
% distance. Useful for reproduction of scenarios where the conditions are
% fixed.
% (c) Josep Colom Ikuno, INTHFT, 2013
% www.nt.tuwien.ac.at

%% Calculate ROI
roi_x = [-10 10];
roi_y = [-10 10];

% Set dummy ROI of 10mx10m (arbitrarily chosen, just so that the plotting does not crash)
for b_ = 1:length(LTE_config.pathlosses)
    eNodeBs_sites(b_)           = network_elements.eNodeB;
    eNodeBs_sites(b_).id        = b_;
    eNodeBs_sites(b_).pos       = rand(1,2).*20 - 10;
    eNodeBs_sites(b_).site_type = 'macro';
end

%% Necessary data loaded in the LTE_init_config config
data_res               = 1;                          % Meters/pixel, resolution of the map
eNodeB_sector_tx_power = LTE_config.eNodeB_tx_power; % eNodeB tx power (Watts/sector)

% Add the Antennas to the eNodeBs
s_idx   = 1;
eNodeBs = network_elements.eNodeB_sector; % Initialization
for b_ = 1:length(eNodeBs_sites)
    % Create the eNodeB_sector objects
    eNodeBs_sites(b_).sectors               = network_elements.eNodeB_sector;
    eNodeBs_sites(b_).sectors.parent_eNodeB = eNodeBs_sites(b_);
    eNodeBs_sites(b_).sectors.id            = 1;
    eNodeBs_sites(b_).sectors.azimuth       = 0;
    eNodeBs_sites(b_).sectors.max_power     = eNodeB_sector_tx_power;
    eNodeBs_sites(b_).sectors.antenna_type  = 'omnidirectional';
    eNodeBs_sites(b_).sectors.nTX           = LTE_config.nTX;
    eNodeBs_sites(b_).sectors.tx_height     = 0;
    
    eNodeBs_sites(b_).sectors.eNodeB_id     = s_idx;
    eNodeBs(s_idx)                          = eNodeBs_sites(b_).sectors;
    
    % Attach the correct antenna to the eNodeB
    antennas.antenna.attach_antenna_to_eNodeB(eNodeBs_sites(b_).sectors,LTE_config);
    
    % Create the macroscopic pahloss model that will be used
    if LTE_config.debug_level>=1
        fprintf('Site %d: ',b_);
    end
    
    macroscopic_pathloss_model_settings.pathloss = LTE_config.pathlosses(b_);
    eNodeBs_sites(b_).sectors.macroscopic_pathloss_model = macroscopic_pathloss_models.generalPathlossModel.generateMacroscopicPathlossModel(...
        LTE_config,...
        LTE_config.network_source,...
        [],...
        macroscopic_pathloss_model_settings);
    
    s_idx = s_idx + 1;
end

%% Create pathlossMap
networkMacroscopicPathlossMap          = channel_gain_wrappers.macroscopicPathlossMap;
networkMacroscopicPathlossMap.data_res = data_res;
networkMacroscopicPathlossMap.roi_x    = roi_x;
networkMacroscopicPathlossMap.roi_y    = roi_y;

%% Put the pathloss for every pixel in the ROI in a 3D matrix.
% The equivalent "real-life" size of each pixel is determined by the
% data_res variable.

%% Calculate final pathloss. With and without minimum coupling loss
if LTE_config.debug_level>=1
    fprintf('Creating cell pathloss map\n');
end
macroscopic_pathloss_models.generalPathlossModel.calculate_pathloss_maps(LTE_config,eNodeBs_sites,networkMacroscopicPathlossMap);

all_sectors = [eNodeBs_sites.sectors];
all_pathloss_models = [all_sectors.macroscopic_pathloss_model];
all_pathloss_model_names = {all_pathloss_models.name};

networkMacroscopicPathlossMap.name = all_pathloss_model_names;
