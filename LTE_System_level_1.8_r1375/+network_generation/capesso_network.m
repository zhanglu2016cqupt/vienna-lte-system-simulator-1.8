function [sites, eNodeBs, networkMacroscopicPathlossMap] = capesso_network(LTE_config)
% (c) Martin Taranetz , INTHFT, 2010

% output:   eNodeBs            ... contains info reagarding the BTSs and its
%                                  sectors
%           pathloss_data      ... [heightxwidthx3xnBTS] double
%                                  Pathloss data for each sector (including
%                                  antenna gain).
%                                  [y,x,sector_num,brts_num]

%% Configuration parameters
capesso_params                            = LTE_config.capesso_params;
capesso_params.manually_set_ROI           = LTE_config.manually_set_ROI;
capesso_params.frequency                  = LTE_config.frequency;
capesso_params.maps_resolution            = LTE_config.map_resolution;
capesso_params.rescale_factor             = LTE_config.rescale_factor;
capesso_params.enable_plotting            = LTE_config.show_network;
capesso_params.debug_plotting             = false;
LTE_config.capesso_params                 = capesso_params;

capesso_files_reader = capesso.capessoFilesReader(capesso_params);

%% Create the eNodeBs
[sites,eNodeBs] = capesso_files_reader.read_cell_data;

%% Read the elevation map
digital_terrain_model = capesso_files_reader.read_DTM;

%% Read the pathloss maps
if LTE_config.debug_level>=1
    fprintf('Creating cell pathloss map from Capesso data\n');
end
M_Capesso = capesso_files_reader.read_pathloss_data(sites); % Read the pathloss data of the relevan sites

%% Cut the capesso pathloss maps to the ROI relevant to the site's positions
[networkMacroscopicPathlossMap, DTM_cut] = capesso.capessoUtils.cut_pathloss_maps_to_sites_ROI(M_Capesso,sites,digital_terrain_model,LTE_config);

%% Apply the antenna gain patterns
if LTE_config.debug_level>=1
    fprintf('Calculating sector antenna gains\n');
end
macroscopic_pathloss_models.generalPathlossModel.calculate_pathloss_maps(LTE_config,sites,networkMacroscopicPathlossMap,DTM_cut);

%% Plotting
% Creates Figures for each eNodeB containing
%   o  Capesso pathloss map
%   o  Elevation map
%   o  eNodeB positions
%   o  Antenna gain map for each sector
if capesso_params.enable_plotting
    LTE_plot_capesso_files(LTE_config,sites, M_Capesso, networkMacroscopicPathlossMap, digital_terrain_model, capesso_params);
end

% Add number of antennas information, as well as sector info
for b_=1:length(sites)
    for s_=1:length(sites(b_).sectors)
        sites(b_).sectors(s_).nTX = LTE_config.nTX;
    end
end

% Convert pathloss to linear
networkMacroscopicPathlossMap.pathloss = 10.^(networkMacroscopicPathlossMap.pathloss./10);

end
