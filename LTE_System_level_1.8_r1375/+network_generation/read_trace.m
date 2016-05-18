function [eNodeBs_sites,eNodeBs,networkMacroscopicPathlossMap] = read_trace(LTE_config,varargin)
% Generates the object structures necessary to have the simulator running
% while taking the values from a trace file. Most of the values are filled
% with dummy values so as not to crash the simulator.
% (c) Josep Colom Ikuno, INTHFT, 2013
% www.nt.tuwien.ac.at

if isempty(varargin) || isempty(varargin{1}) || isempty(varargin{1}{1}{3})
    UE_traces = utils.read_simulation_trace(LTE_config);
else
    UE_traces = varargin{1}{1}{3};
end

%% Create the eNodeBs
if LTE_config.debug_level>=1
    fprintf('Creating eNodeBs, (trace)\n');
end
N_cells = length(UE_traces(1).all_cell_ids);

%% Create the eNodeB array
for b_ = 1:N_cells
    eNodeBs_sites(b_)           = network_elements.eNodeB;
    eNodeBs_sites(b_).id        = b_;
    eNodeBs_sites(b_).pos       = [0 0];
    eNodeBs_sites(b_).site_type = 'trace';
end

% Add the Antennas to the eNodeBs
s_idx   = 1;
eNodeBs = network_elements.eNodeB_sector; % Initialization
for b_ = 1:length(eNodeBs_sites)
    % Create the eNodeB_sector objects
    % Writing eNodeBs(b_).sectors(1) gave me an error. Maybe a bug??
    eNodeBs_sites(b_).sectors    = network_elements.eNodeB_sector;
    eNodeBs_sites(b_).sectors               = network_elements.eNodeB_sector;
    eNodeBs_sites(b_).sectors.parent_eNodeB = eNodeBs_sites(b_);
    eNodeBs_sites(b_).sectors.id            = 1;
    eNodeBs_sites(b_).sectors.azimuth       = 0;
    eNodeBs_sites(b_).sectors.max_power     = LTE_config.eNodeB_tx_power;
    eNodeBs_sites(b_).sectors.nTX           = LTE_config.nTX;
    eNodeBs_sites(b_).sectors.tx_height     = LTE_config.tx_height;
    
    eNodeBs_sites(b_).sectors.eNodeB_id     = s_idx;
    eNodeBs(s_idx)                          = eNodeBs_sites(b_).sectors;
    
    s_idx = s_idx + 1;
end

% Add the TTI time step
for u_=1:length(UE_traces)
    UE_traces(u_).TTIs_per_time_idx = LTE_config.TTI_per_trace_step;
end

%% Create pathlossMap
networkMacroscopicPathlossMap          = channel_gain_wrappers.macroscopicPathlossMap;
networkMacroscopicPathlossMap.name     = 'trace';
networkMacroscopicPathlossMap.pathloss = UE_traces;
