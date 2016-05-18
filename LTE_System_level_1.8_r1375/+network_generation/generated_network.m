function [eNodeB_sites,eNodeBs,networkMacroscopicPathlossMap] = generated_network(LTE_config)
% Generate an hexagonal network with a free space pathloss model
% (c) Josep Colom Ikuno, INTHFT, 2008
% www.nt.tuwien.ac.at
% output:   eNodeBs            ... contains info reagarding the BTSs and its
%                                  sectors
%           pathloss_data      ... [heightxwidthx3xnBTS] double
%                                  Pathloss data for each sector (including
%                                  antenna gain).
%                                  [y,x,sector_num,brts_num]

%% Necessary data loaded in the LTE_init_config config
data_res               = LTE_config.map_resolution;           % Meters/pixel, resolution of the map
eNodeB_sector_tx_power = LTE_config.eNodeB_tx_power;          % eNodeB tx power (Watts/sector)

%% Some initialisation params
ROI_increase_factor = 0.1;

%% Create the eNodeBs
if LTE_config.debug_level>=1
    fprintf(['Creating eNodeBs, (geometry: "' LTE_config.network_geometry '")\n']);
end

if ~isfield(LTE_config,'network_geometry') % in case an older config file is used - hexagonal grid by default
   LTE_config.network_geometry = 'regular_hexagonal_grid'; 
end

switch LTE_config.network_geometry
    case 'regular_hexagonal_grid'
        eNodeB_positions = network_geometry.hexagonal_eNodeB_grid(LTE_config);
    case 'stochastic' % stochastic distribution of eNodeBs
        eNodeB_positions = network_geometry.stochastic_eNodeB_distribution(LTE_config);
    case 'hybrid' % stochastic and deterministic placement of eNodeBs
        eNodeB_positions = network_geometry.hybrid_eNodeB_distribution(LTE_config);
    case 'circular' % eNodeBs placed equidistantly on a circle.
        eNodeB_positions = network_geometry.circular_eNodeB_distribution(LTE_config);
    case 'predefined' % Predefined eNodeB distribution
        eNodeB_positions = LTE_config.eNodeB_positions;
    otherwise
        error(['Network geometry "' LTE_config.network_geometry '" is not supported']);
end

%% Create the eNodeB array
for b_ = 1:size(eNodeB_positions,1)
    eNodeB_sites(b_)           = network_elements.eNodeB;
    eNodeB_sites(b_).id        = b_;
    eNodeB_sites(b_).pos       = [eNodeB_positions(b_,1) eNodeB_positions(b_,2)];
    eNodeB_sites(b_).site_type = 'macro';
end

% Add the Antennas to the eNodeBs
s_idx   = 1;
eNodeBs = network_elements.eNodeB_sector; % Initialization
for b_ = 1:length(eNodeB_sites)
    % Create the eNodeB_sector objects
    % Writing eNodeBs(b_).sectors(1) gave me an error. Maybe a bug??
    eNodeB_sites(b_).sectors    = network_elements.eNodeB_sector;
    for s_ = 1:length(LTE_config.sector_azimuths)
        eNodeB_sites(b_).sectors(s_)               = network_elements.eNodeB_sector;
        eNodeB_sites(b_).sectors(s_).parent_eNodeB = eNodeB_sites(b_);
        eNodeB_sites(b_).sectors(s_).id            = s_;
        eNodeB_sites(b_).sectors(s_).azimuth       = utils.miscUtils.wrapTo359(LTE_config.antenna_azimuth_offsett + LTE_config.sector_azimuths(s_));
        eNodeB_sites(b_).sectors(s_).max_power     = eNodeB_sector_tx_power;
        eNodeB_sites(b_).sectors(s_).antenna_type  = LTE_config.antenna.antenna_gain_pattern;
        eNodeB_sites(b_).sectors(s_).nTX           = LTE_config.nTX;
        eNodeB_sites(b_).sectors(s_).tx_height     = LTE_config.tx_height;
        
        eNodeB_sites(b_).sectors(s_).eNodeB_id     = s_idx;
        eNodeBs(s_idx)                              = eNodeB_sites(b_).sectors(s_);
        
        % Attach the correct antenna to the eNodeB
        antennas.antenna.attach_antenna_to_eNodeB(eNodeB_sites(b_).sectors(s_),LTE_config);
        
        % Create the macroscopic pahloss model that will be used
        if LTE_config.debug_level>=1
            fprintf('Site %d, eNodeB %d: ',b_,s_);
        end
        
        eNodeB_sites(b_).sectors(s_).macroscopic_pathloss_model = macroscopic_pathloss_models.generalPathlossModel.generateMacroscopicPathlossModel(...
            LTE_config,...
            LTE_config.macroscopic_pathloss_model,...
            LTE_config.frequency,...
            LTE_config.macroscopic_pathloss_model_settings);

        s_idx = s_idx + 1;
    end
end

%% Calculate ROI

number_of_eNodeB_sites = length(eNodeB_sites);

if (number_of_eNodeB_sites > 1)
    tx_pos = reshape([eNodeB_sites.pos],2,[])';
    % Calculate ROI border points in ABSOLUTE coordinates
    roi_x = [min(tx_pos(:,1)),max(tx_pos(:,1))];
    roi_y = [min(tx_pos(:,2)),max(tx_pos(:,2))];
else % set the ROI equal to eNodeB distance if just one base station is used
    roi_x = [-LTE_config.inter_eNodeB_distance,LTE_config.inter_eNodeB_distance];
    roi_y = [-LTE_config.inter_eNodeB_distance,LTE_config.inter_eNodeB_distance];
end


%% Define an area of the ROI to map
% roi_reduction_factor times smaller and draw it. ABSOLUTE COORDINATES
if strcmp(LTE_config.network_geometry,'predefined') % Manually set roi
    roi_x = LTE_config.predefined_roi_x;
    roi_y = LTE_config.predefined_roi_y;
else
    roi_x = roi_x + ROI_increase_factor*abs(roi_x(2)-roi_x(1))*[-1,1];
    roi_y = roi_y + ROI_increase_factor*abs(roi_y(2)-roi_y(1))*[-1,1];
end

%% Create pathlossMap
networkMacroscopicPathlossMap                        = channel_gain_wrappers.macroscopicPathlossMap;
networkMacroscopicPathlossMap.data_res               = data_res;
networkMacroscopicPathlossMap.roi_x                  = roi_x;
networkMacroscopicPathlossMap.roi_y                  = roi_y;

%% Put the pathloss for every pixel in the ROI in a 3D matrix.
% The equivalent "real-life" size of each pixel is determined by the
% data_res variable.

%% Calculate final pathloss. With and without minimum coupling loss
if LTE_config.debug_level>=1
    fprintf('Creating cell pathloss map\n');
end
macroscopic_pathloss_models.generalPathlossModel.calculate_pathloss_maps(LTE_config,eNodeB_sites,networkMacroscopicPathlossMap);

all_sectors = [eNodeB_sites.sectors];
all_pathloss_models = [all_sectors.macroscopic_pathloss_model];
all_pathloss_model_names = {all_pathloss_models.name};

networkMacroscopicPathlossMap.name = all_pathloss_model_names;