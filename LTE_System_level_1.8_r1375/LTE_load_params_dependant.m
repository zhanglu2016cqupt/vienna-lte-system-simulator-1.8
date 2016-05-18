function LTE_config = LTE_load_params_dependant(LTE_config)
%% Automagically filled params and config params that will probably not be changed

% Hardcode the release number
LTE_config.release             = 'r1119';
matlab_release                 = version('-release');
matlab_year                    = str2double(matlab_release(1:(end-1)));
LTE_config.matlab_release_year = matlab_year;

% Check if the pathloss maps will need to be re-run
LTE_config.minimum_pathloss_map_version                      = 5;
LTE_config.maximum_supported_secondary_channel_traces_per_UE = 5;

% Default trace folders
LTE_config.default_channel_trace_folder = './data_fileschannel_traces';
LTE_config.default_UE_cache_folder      = './data_filesUE_caches';
LTE_config.default_network_cache_folder = './data_filesnetwork_caches';

%% Fill in some scheduler params

% Load the linear approximations to the BICM curves (necessary for scheduling)
load('./+utils/BICM_k_d_MSE.mat','k','d');
LTE_config.MI_data.k = k;
LTE_config.MI_data.d = d;
LTE_config.scheduler_params.k = LTE_config.MI_data.k; % to map from CQI to spectral efficiency
LTE_config.scheduler_params.d = LTE_config.MI_data.d; % to map from CQI to spectral efficiency

% Exponential averaging window:
% Number of TTIs used to calculate the average throughput (exponential filtering)
% See: P. Viswanath, D. Tse and R. Laroia, "Opportunistic Beamforming using Dumb Antennas", IEEE Transactions on Information Theory, vol. 48(6), June, 2002.
if ~isfield(LTE_config.scheduler_params,'av_window')
    LTE_config.scheduler_params.av_window = 25; % Default value
end

%% Throughput mapping for spectral efficiency calculations
LTE_config.SNR_to_throughput_mapping = './data_filesSNR_to_throughput_mapping.mat';
                                                                       
%% Check whether the parallel computing toolbox is installed
v = ver;
if isempty(strfind([v.Name],'Parallel Computing Toolbox'))
    LTE_config.parallel_toolbox_installed = false;
else
    LTE_config.parallel_toolbox_installed = true;
end

%% Check whether the results folder exists, and if not, create it
if ~exist(LTE_config.results_folder,'dir')
    mkdir(LTE_config.results_folder);
end

%% support of parfor to simulate multiple realizations of the network in parallel
if ~isfield(LTE_config,'parallel_network')
    LTE_config.parallel_network = false;    
elseif LTE_config.parallel_network
    if ~strcmp(LTE_config.network_cache,'auto')
        LTE_config.network_cache = 'auto';
        warning('LTE_config.network_cache is set to "auto" when multiple networks are simulated in parallel (LTE_config.parallel_network = true)')
    end
    if ~strcmp(LTE_config.results_file,'auto')
        LTE_config.results_file = 'auto';
        warning('LTE_config.results_file is set to "auto" when multiple networks are simulated in parallel (LTE_config.parallel_network = true)')
    end
    if ~strcmp(LTE_config.UE_cache_file,'auto')
        LTE_config.UE_cache_file = 'auto';
        warning('LTE_config.results_file is set to "auto" when multiple networks are simulated in parallel (LTE_config.parallel_network = true)')
    end
    if LTE_config.seedRandStream
        warning('Random number seeding (LTE_config.seedRandStream) not supported in parallel mode (LTE_config.parallel_network = true)')
    end
    LTE_config.seedRandStream = true; % make sure that we do not simulate the same stuff multiple times
    t = [];
    if LTE_config.parallel_toolbox_installed
        t = getCurrentTask();
    end
    if ~isempty(t)
        LTE_config.RandStreamSeed = num2str(t.ID) + rand(1)*10^4 + cputime;
    else
        LTE_config.RandStreamSeed = rand(1)*10^4 + cputime;
    end    
    if ~LTE_config.parallel_toolbox_installed
        error('Parallel computing toolbox not available! Set LTE_config.parallel_network = false.')
    end
end
%% Random number generation
% Number of sectors
if ~isfield(LTE_config,'seedRandStream')
    LTE_config.seedRandStream = false; % Mimic old behavior
end

if LTE_config.seedRandStream
    if LTE_config.matlab_release_year >= 2011
        RandStream.setGlobalStream(RandStream('mt19937ar','Seed',LTE_config.RandStreamSeed));
    else
        RandStream.setDefaultStream(RandStream('mt19937ar','Seed',LTE_config.RandStreamSeed)); %#ok<SETRS> Support for old versions
    end
else
    if LTE_config.matlab_release_year >= 2011
        RandStream.setGlobalStream(RandStream('mt19937ar','Seed',rand*intmax('uint32')));
    else
        RandStream.setDefaultStream(RandStream('mt19937ar','Seed',rand*intmax('uint32'))); %#ok<SETRS> Support for old versions
    end
end

%% Backwards-compatibility safeguards
% Trace version: v1 is the old one (precoding included), the new one implements run-time precoding
if ~isfield(LTE_config,'trace_version')
    LTE_config.trace_version = 'v1'; % Default to the old version (for now)
end

% Modes employing runtime precoding
if LTE_config.tx_mode >= 100
    LTE_config.trace_version = 'v2';
end

% Configure the runtime precoding
switch LTE_config.trace_version
    case 'v1'
        LTE_config.runtime_precoding = false;
    otherwise
        LTE_config.runtime_precoding = true;
end

% Whether Distributed Antenna Systems (DAS) is enabled
if ~isfield(LTE_config,'RRHs_enabled')
    LTE_config.RRHs_enabled = false;
end

if LTE_config.RRHs_enabled && ~LTE_config.runtime_precoding
    error('To use Remote Radio Heads (RRHs), run-time precoding mus be employed');
end

% Check for trace mode
if ~isfield(LTE_config,'trace_simulation_mode')
    LTE_config.trace_simulation_mode = false; % Mimic old behavior
end

if LTE_config.trace_simulation_mode
    LTE_config.cache_network      = false;
    LTE_config.UE_cache           = false;
    LTE_config.network_source     = 'trace';
    LTE_config.shadow_fading_type = 'none';
    LTE_config.UE_distribution    = 'trace';
    LTE_config.UE_cache_file      = [];
    LTE_config.keep_UEs_still     = true; % To avoid calls to the walking model
end

% Whether less feedback is to be stored in the result traces (for the case
% where a lot of UEs are present)
if ~isfield(LTE_config,'reduced_feedback_logs')
    LTE_config.reduced_feedback_logs = false; % Mimic old behavior
end

% Number of sectors
if ~isfield(LTE_config,'sector_azimuths')
    if strcmp(LTE_config.antenna.antenna_gain_pattern,'omnidirectional')
        LTE_config.sector_azimuths = 0; % only one eNodeB per site should be generated
    else
        LTE_config.sector_azimuths = 0:360/3:359; % Mimic old behavior     
    end
end

% Macroscopic pathloss settings
if ~isfield(LTE_config,'macroscopic_pathloss_model_settings')
    LTE_config.macroscopic_pathloss_model_settings = []; % Mimic old behavior
end

% Add femtocells
if ~isfield(LTE_config,'add_femtocells')
    LTE_config.add_femtocells = false; % Mimic old behavior
end

if LTE_config.add_femtocells
    if ~isfield(LTE_config.femtocells_config,'macroscopic_pathloss_model')
        LTE_config.femtocells_config.macroscopic_pathloss_model          = LTE_config.macroscopic_pathloss_model; % Mimic old behavior
        LTE_config.femtocells_config.macroscopic_pathloss_model_settings = LTE_config.macroscopic_pathloss_model_settings;
    end
    if ~isfield(LTE_config.femtocells_config,'macroscopic_pathloss_model_settings')
        LTE_config.femtocells_config.macroscopic_pathloss_model_settings = []; % Mimic old behavior
    end
end

% Minimum coupling loss
if ~isfield(LTE_config,'minimum_coupling_loss')
    LTE_config.minimum_coupling_loss = 0; % Mimic old behavior
end

% Interference
if ~isfield(LTE_config,'always_on')
    LTE_config.always_on = true; % Mimic old behavior
end

% Manual ROI setting
if ~isfield(LTE_config,'manually_set_ROI')
    LTE_config.manually_set_ROI = false; % Mimic old behavior
else
    if LTE_config.manually_set_ROI && (~isfield(LTE_config,'roi_x') || ~isfield(LTE_config,'roi_y'))
        error('When manually setting the ROI you need to specify the roi_x and roi_y variables');
    end
end

% Default antenna azimuth offsett of 30 (hexagonal sectors)
if ~isfield(LTE_config,'antenna_azimuth_offsett')
    LTE_config.antenna_azimuth_offsett = 30; % Mimic old behavior
end

% UE antenna gain (assumed omnidirectional)
if ~isfield(LTE_config.UE,'antenna_gain')
    LTE_config.UE.antenna_gain = 0; % Mimic old behavior
end

% eNodeB (macro) Tx height
if ~isfield(LTE_config,'tx_height')
    LTE_config.tx_height = 32; % [m] According to 3GPP 36.814 Table A.2.1.1-2. The transmitter height
end

% UE receiver height
if ~isfield(LTE_config,'rx_height')
    LTE_config.rx_height = 1.5; % [m] Mimic old behavior
end

% Whether TX and RX height are taken into account for the pathloss calculation
if ~isfield(LTE_config,'calculate_3D_pathloss')
    LTE_config.calculate_3D_pathloss = false; % [m] Mimic old behavior
end

% Whether to save results in a very compact way or the "traditional way" (results object)
if ~isfield(LTE_config,'compact_results_file') % Save just the UE throughput
    LTE_config.compact_results_file = false; % Mimic old behavior
end

% Traffic map upscaling factor
if ~isfield(LTE_config,'traffic_map_upscaling')
    LTE_config.traffic_map_upscaling = 1; % Mimic old behavior
end

%Wheter to utilize car user
if ~isfield(LTE_config,'support_car_user')
    LTE_config.support_car_user = false;
end;

if LTE_config.support_car_user == true
    
   if(~strcmp(LTE_config.scheduler,'CarScheduler')) 
        LTE_config.scheduler = 'CarScheduler'; 
   end
   
   if ~isfield(LTE_config,'traffic_models')
        LTE_config.traffic_models.type = 'car';
   end
   
   if(~strcmp(LTE_config.traffic_models.type,'car'))
       LTE_config.traffic_models.type = 'car';
   end
   
   if(~isfield(LTE_config.traffic_models,'usetraffic_model'))
      LTE_config.traffic_models.usetraffic_model = true; 
   end
   
   if(LTE_config.traffic_models.usetraffic_model == false)
        LTE_config.traffic_models.usetraffic_model = true; 
   end
   
   if ~isfield(LTE_config,'car_per_cent')
       LTE_config.car_per_cent = 0.1; 
   end    
    
end

% Whether to delete the microscale fading trace data at the end (saves some space)
if ~isfield(LTE_config,'delete_ff_trace_at_end')
    LTE_config.delete_ff_trace_at_end = false; % Mimic old behavior
end

% Whether to delete the pathloss data at the end (saves some space)
if ~isfield(LTE_config,'delete_pathloss_at_end')
    LTE_config.delete_pathloss_at_end = false; % Mimic old behavior
end

% Unquantized feedback
if ~isfield(LTE_config,'unquantized_CQI_feedback')
    LTE_config.unquantized_CQI_feedback = false; % Mimic old behavior
end

% Channel trace id
if ~isfield(LTE_config,'channel_trace_id')
    LTE_config.channel_trace_id = []; % Mimic old behavior
end

% Force non-parallel channel trace
if ~isfield(LTE_config,'non_parallel_channel_trace')
    LTE_config.non_parallel_channel_trace = false; % Mimic old behavior
end

% Extra tracing
if ~isfield(LTE_config,'trace_SINR')
    LTE_config.trace_SINR = false; % Mimic old behavior
end

% Adaptive RI
if ~isfield(LTE_config,'adaptive_RI')
    LTE_config.adaptive_RI = 0; % Mimic old behavior
end

% Whether the UEs should move or stay still
if ~isfield(LTE_config,'keep_UEs_still')
    LTE_config.keep_UEs_still = false; % Mimic old behavior
end

% Some GUI-related parameters
if ~isfield(LTE_config,'default_shown_GUI_cells')
    LTE_config.default_shown_GUI_cells = []; % Mimic old behavior
end

% Whether some UEs should be deactivated (faster simulation) if they are attached to any of the following eNodeBs
if ~isfield(LTE_config,'compute_only_UEs_from_this_eNodeBs') || strcmp(LTE_config.network_source,'fixed pathloss')
    LTE_config.compute_only_UEs_from_this_eNodeBs = []; % Mimic old behavior and disable for fixed pathloss mode
end

if strcmp(LTE_config.network_source,'fixed pathloss')
    LTE_config.default_shown_GUI_cells = []; % Disable for fixed pathloss mode
    LTE_config.network_cache           = false;
    LTE_config.shadow_fading_type      = 'none';
end

% Whether the MCL is re-applied at the end regardless of what is loaded from
% cache. If set to 'false', this extra step will not be performed when
% loading from cache or planning tool data (slightly shorter startup time)
if ~isfield(LTE_config,'reapply_MCL_to_cache_or_planning_tool')
    LTE_config.reapply_MCL_to_cache_or_planning_tool = false; % Mimic old behavior
end

% Whether to use wideband precoding
if LTE_config.tx_mode==4 || LTE_config.tx_mode == 9
    if ~isfield(LTE_config,'wideband_precoding')
        LTE_config.wideband_precoding = false; % Mimic old behavior
    end
end
if LTE_config.tx_mode==6 || LTE_config.tx_mode == 5
    if ~isfield(LTE_config,'wideband_precoding')
        LTE_config.wideband_precoding = false; % Mimic old behavior
    end
end
% By default the values specified in TS 36.942
if ~isfield(LTE_config,'TS_36942_3dB_lobe')
    LTE_config.TS_36942_3dB_lobe = 65;
end

% Winner extra antenna gain parameters
if strcmp(LTE_config.channel_model.type,'winner+')
    % Check if there are extra winner config params
    if~ isfield(LTE_config,'winner_antenna_params')
        LTE_config.winner_antenna_params = []; % Mimic old behavior
    end
end

% How the UEs are generated
if ~isfield(LTE_config,'UE_distribution')
    LTE_config.UE_distribution = 'constant UEs per cell'; % Mimic old behavior. And just in case...
end
switch LTE_config.UE_distribution
    case 'traffic map'
        LTE_config.UE.use_traffic_map = true;
    otherwise
        LTE_config.UE.use_traffic_map = false; % Mimic old behavior.
end

% FFR mode activation
switch LTE_config.scheduler
    case 'FFR'
        LTE_config.FFR_active = true;  % Activate FFR
    otherwise
        LTE_config.FFR_active = false; % Mimic old behavior
end

% FFR parameters that override the automatic setting of the thresholds
if ~isfield(LTE_config,'FFR_override')
    LTE_config.FFR_override = true; % Do not use the automatic mapping of beta and the SINR threshold.
end

% Configurable additional penetration loss (e.g. 7dB extra for "in car" losses)
if ~isfield(LTE_config,'additional_penetration_loss') || isempty(LTE_config.additional_penetration_loss)
    LTE_config.additional_penetration_loss = 0; % Mimic old behavior
end
if ~isnumeric(LTE_config.additional_penetration_loss)
    switch LTE_config.additional_penetration_loss
        case 'train'
            LTE_config.additional_penetration_loss = 30;
        case 'deep indoor'
            LTE_config.additional_penetration_loss = 23;
        case 'indoor'
            LTE_config.additional_penetration_loss = 17;
        case 'incar'
            LTE_config.additional_penetration_loss = 7;
        case 'outdoor'
            LTE_config.additional_penetration_loss = 0;
        otherwise
            error('"deep indoor", "indor", "incar", and "outdoor" values or a penetration loss value in dB are valid');
    end
end

% Preserve the previously-loaded pregenerated_ff object to avoid continuous reading of large files during long batches of short simulations
if ~isfield(LTE_config,'reuse_pregenerated_ff_trace_from_last_run')
    LTE_config.reuse_pregenerated_ff_trace_from_last_run = false; % Mimic old behavior
end

%Whether the pathloss map is to be read from memory (variable) or from a
%file. It can save considerable time when repeatedly reading the same
%network map.
if ~isfield(LTE_config,'usePathlossMapFromMemory')
    LTE_config.usePathlossMapFromMemory = false;
end

% Suffix to the output results file
if ~isfield(LTE_config,'output_filename_suffix')
    LTE_config.output_filename_suffix = []; % Mimic old behavior
end

% Ratio of the whole power dedicated to signaling
if ~isfield(LTE_config,'signaling_ratio')
    LTE_config.signaling_ratio = 0; % Mimic old behavior
else
    if LTE_config.signaling_ratio>1 || LTE_config.signaling_ratio<0
        error('Signaling ratio must be between 0 and 1.');
    end
end

% Pathloss options
switch LTE_config.network_source
    case 'generated'
        switch LTE_config.network_geometry
            case 'regular_hexagonal_grid'
                if ~isfield(LTE_config,'inter_eNodeB_distance')
                    error('LTE_config.inter_eNodeB_distance must be defined for the generated regular_hexagonal_grid network case');
                end
                if ~isfield(LTE_config,'nr_eNodeB_rings')
                    error('LTE_config.nr_eNodeB_rings must be defined for the generated regular_hexagonal_grid network case');
                end
            case 'stochastic'
                if ~isfield(LTE_config,'average_eNodeB_distance')
                    error('LTE_config.averag_eNodeB_distance must be defined for the generated stochastic network');
                end
                if ~isfield(LTE_config,'network_size')
                    error('LTE_config.network_size must be defined for the generated regular_hexagonal_grid network case');
                end
            case 'hybrid'
                if ~isfield(LTE_config,'inter_eNodeB_distance')
                    error('LTE_config.inter_eNodeB_distance must be defined for the generated hybrid network');
                end
            case 'circular'
                if ~isfield(LTE_config,'eNodeB_circle_radius')
                    error('LTE_config.eNodeB_circle_radius must be defined for the generated circular network');
                end
                if ~isfield(LTE_config,'eNodeB_circle_nBSs')
                    error('LTE_config.eNodeB_circle_nBSs must be defined for the generated circular network');
                end
            case 'predefined'
                if ~isfield(LTE_config,'eNodeB_positions')
                    error('LTE_config.eNodeB_positions must be specified in predefined setup');
                end
            otherwise
                error(['Network geometry "' LTE_config.network_geometry '" is not supported'])
        end
    case {'capesso' 'trace'}
    otherwise
        error([LTE_config.network_source ' network source not supported']);
end

% Whether there is one shadow fading map per site (default) or one map per
% sector (set it to true, then)
if ~isfield(LTE_config,'decouple_site_shadow_fading_maps')
   LTE_config.decouple_site_shadow_fading_maps = false; % Mimic old behavior
end

switch LTE_config.shadow_fading_type
    case 'claussen'
        if ~isfield(LTE_config,'deactivate_claussen_spatial_correlation')
            LTE_config.deactivate_claussen_spatial_correlation = false;
        end
        
        if ~isfield(LTE_config,'shadow_fading_map_resolution') ||...
                ~isfield(LTE_config,'shadow_fading_n_neighbors') ||...
                ~isfield(LTE_config,'shadow_fading_mean') ||...
                ~isfield(LTE_config,'shadow_fading_sd') ||...
                ~isfield(LTE_config,'r_eNodeBs')
                error('When specifying a "claussen" shadow fading map, the following parameters need to be specified:\n  %s\n  %s\n  %s\n  %s\n  %s\n',...
                'LTE_config.shadow_fading_map_resolution',...
                'LTE_config.shadow_fading_n_neighbors',...
                'LTE_config.shadow_fading_mean',...
                'LTE_config.shadow_fading_sd',...
                'LTE_config.r_eNodeBs');
        end
end

%% Some other checks

%% Check whether the macroscopic pathloss comes from a model or from data
switch LTE_config.network_source
    case 'capesso'
        LTE_config.shadow_fading_type = 'none';
        LTE_config.macroscopic_pathloss_is_model = false;
        if ~isfield(LTE_config.capesso_params,'enable_dtm')
            LTE_config.capesso_params.enable_dtm = false; % Mimic old behavior
        end
        LTE_config.antenna.antenna_gain_pattern = 'kathreinTSAntenna';
        LTE_config.antenna.kathrein_antenna_folder = LTE_config.capesso_params.kathrein_antenna_folder;
    otherwise
        LTE_config.macroscopic_pathloss_is_model = true;
end

% 3D antenna radiation pattern
switch LTE_config.network_source
    case {'capesso','trace'}
        % Do nothing
    otherwise
        if strcmp(LTE_config.antenna.antenna_gain_pattern, 'kathreinTSAntenna')
            LTE_config.AntennaPattern3d = true;
        end
end

%% Moved from the "do not touch" section
LTE_config.RB_bandwidth    = 180e3;         % Frequency in Hz
LTE_config.TTI_length      = 1e-3;          % Length of a TTI (subframe), in seconds.
LTE_config.cyclic_prefix   = 'normal';      % 'normal' or 'extended' cyclic prefix. Not working for values other that 'normal'
LTE_config.maxStreams      = 2;             % Maximum number of codewords per TTI. For LTE, that's 2. If you change it, the simulator will 100% sure crash.
                                            % I took the name from HSDPA. In the LTE standard is actually referred as 'codewords'

if isfield(LTE_config,'BF_short') % check if field exists
    if LTE_config.BF_short  % if short block-fading is activated, selected the time instants for which channel samples will be generated 
        if strcmp(LTE_config.cyclic_prefix,'normal')
%             LTE_config.BF_samples = [2,5,8,11,14];
            LTE_config.BF_samples = [1,5,9,13];
%             LTE_config.BF_samples = 1:14;
        else
            LTE_config.BF_samples = [2,5,8,11];
        end        
    else
        LTE_config.BF_samples =  [1];
    end    
else
    LTE_config.BF_short = false;
    LTE_config.BF_samples =  [1];
end

                                            
%% Warning about 0-delay feedback
if LTE_config.feedback_channel_delay==0
    warning('The assumed allocated power when calculating the SINR and performing the simulation will always be assumed to be %dW. The assigned power will affect the measured SINR, but in order to measure it, the scheduling has to be done first, which cannot be done without proper feedback. Therefore, power assignment is not used for 0-delay.',LTE_config.eNodeB_tx_power);
    LTE_config.always_on = true;
end

%% Check SISO and SIMO mode settings
if LTE_config.tx_mode==1 && LTE_config.nTX>1
    error('For SIxO mode, nTX must be set to 1 (or possibly 0 in RRH setup).');
end


%% Naming utility object
name_generator = utils.naming;

%% Results file filename
if ~LTE_config.trace_simulation_mode
    LTE_config.results_file         = name_generator.results_file(LTE_config);
    LTE_config.network_cache        = name_generator.macroscopic_pathloss_cache(LTE_config);
    LTE_config.UE_cache_file        = name_generator.UE_cache(LTE_config);
end
LTE_config = name_generator.channel_trace_cache(LTE_config);

%% Channel trace generation configuration
switch LTE_config.channel_model.type
    case 'winner+'
        winner_model_path = './Winner Channel Model';
        if LTE_config.recalculate_fast_fading || ~exist([LTE_config.pregenerated_ff_file '.mat'],'file')
            channel_gain_wrappers.winnerChannelFactory.prepare_winner_channel_model_path(winner_model_path); % Add winner channel model path to the Matlab path
        end
        LTE_config.trace_params = channel_gain_wrappers.winnerChannelFactory.get_default_config(LTE_config.frequency,LTE_config.nTX,LTE_config.nRX,LTE_config.UE_speed);
    otherwise
        LTE_config.trace_params = channel_gain_wrappers.pdpChannelFactory.get_default_config(LTE_config.frequency,LTE_config.nTX,LTE_config.nRX,LTE_config.UE_speed,LTE_config.channel_model.correlated_fading,LTE_config.channel_model.type);
end

%% Transmission parameters (used for the throughput calculation)
% We will assume subcarrier a spacing of 15 kHz
switch LTE_config.cyclic_prefix
    case 'normal'
        LTE_config.N_sym = 7;
    case 'extended'
        LTE_config.N_sym = 6;
    otherwise
        error('CP can only be "normal" or "extended"');
end
switch LTE_config.bandwidth
    case 1.4e6
        LTE_config.N_RB = 6;
        LTE_config.fft_points = 128;
        switch LTE_config.cyclic_prefix
            case 'normal'
                LTE_config.CP_length_samples = 9;
            case 'extended'
                LTE_config.CP_length_samples = 32;
        end
    case 3e6
        LTE_config.N_RB = 15;
        LTE_config.fft_points = 256;
        switch LTE_config.cyclic_prefix
            case 'normal'
                LTE_config.CP_length_samples = 18;
            case 'extended'
                LTE_config.CP_length_samples = 64;
        end
    case 5e6
        LTE_config.N_RB = 25;
        LTE_config.fft_points = 512;
        switch LTE_config.cyclic_prefix
            case 'normal'
                LTE_config.CP_length_samples = 36;
            case 'extended'
                LTE_config.CP_length_samples = 128;
        end
    case 10e6
        LTE_config.N_RB = 50;
        LTE_config.fft_points = 1024;
        switch LTE_config.cyclic_prefix
            case 'normal'
                LTE_config.CP_length_samples = 72;
            case 'extended'
                LTE_config.CP_length_samples = 256;
        end
    case 15e6
        LTE_config.N_RB = 75;
        LTE_config.fft_points = 1536;
        switch LTE_config.cyclic_prefix
            case 'normal'
                LTE_config.CP_length_samples = 108;
            case 'extended'
                LTE_config.CP_length_samples = 384;
        end
    case 20e6
        LTE_config.N_RB = 100;
        LTE_config.fft_points = 2048;
        switch LTE_config.cyclic_prefix
            case 'normal'
                LTE_config.CP_length_samples = 144;
            case 'extended'
                LTE_config.CP_length_samples = 512;
        end
    otherwise
        error('Bandwidth not supported');
end
LTE_config.Ntot = LTE_config.N_RB*12;
LTE_config.fs = 15e3*LTE_config.fft_points;
switch LTE_config.nTX   % number of reference symbols
    case 0
        numb = 0;
    case 1
        numb = 4;
    case 2
        numb = 8;
    case 4
        numb = 12;
    case 8
        numb = 8;
    otherwise
        error('Not defined for %d TX antennas',LTE_config.nTX);
end
LTE_config.sym_per_RB_nosync = 12*LTE_config.N_sym - numb; % no synchronization signals 
LTE_config.sym_per_RB_sync = 12*LTE_config.N_sym -(numb+12); % take also synchronization signals into account
LTE_config.scheduler_params.overhead_ref = numb;
LTE_config.scheduler_params.overhead_sync = numb+12;

%% BLER curves
LTE_config.BLER_curves.folder = fullfile(pwd,'data_files','AWGN_BLERs');
LTE_config.BLER_curves.filenames = {
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi1.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi2.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi3.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi4.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi5.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi6.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi7.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi8.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi9.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi10.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi11.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi12.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi13.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi14.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi15.mat')
};

%% MIESM configuration
LTE_config.SINR_averaging.algorithm = 'MIESM';
LTE_config.SINR_averaging.BICM_capacity_tables = './data_filesBICM_capacity_tables_20000_realizations.mat';
LTE_config.SINR_averaging.betas = [3.07,4.41,0.6,1.16,1.06,1.06,0.87,1.01,1.04,1.03,1.11,1.01,1.07,1,1.05];

%% CQI parameters
% CQI 1 is index 1, CQI 2 is index 2, etc...
LTE_config.CQI_params(1).CQI = 1;
LTE_config.CQI_params(1).modulation = 'QPSK';
LTE_config.CQI_params(1).modulation_order = 2;
LTE_config.CQI_params(1).coding_rate_x_1024 = 78;
LTE_config.CQI_params(1).efficiency = 0.1523;

LTE_config.CQI_params(2).CQI = 2;
LTE_config.CQI_params(2).modulation = 'QPSK';
LTE_config.CQI_params(2).modulation_order = 2;
LTE_config.CQI_params(2).coding_rate_x_1024 = 120;
LTE_config.CQI_params(2).efficiency = 0.2344;

LTE_config.CQI_params(3).CQI = 3;
LTE_config.CQI_params(3).modulation = 'QPSK';
LTE_config.CQI_params(3).modulation_order = 2;
LTE_config.CQI_params(3).coding_rate_x_1024 = 193;
LTE_config.CQI_params(3).efficiency = 0.3770;

LTE_config.CQI_params(4).CQI = 4;
LTE_config.CQI_params(4).modulation = 'QPSK';
LTE_config.CQI_params(4).modulation_order = 2;
LTE_config.CQI_params(4).coding_rate_x_1024 = 308;
LTE_config.CQI_params(4).efficiency = 0.6016;

LTE_config.CQI_params(5).CQI = 5;
LTE_config.CQI_params(5).modulation = 'QPSK';
LTE_config.CQI_params(5).modulation_order = 2;
LTE_config.CQI_params(5).coding_rate_x_1024 = 449;
LTE_config.CQI_params(5).efficiency = 0.8770;

LTE_config.CQI_params(6).CQI = 6;
LTE_config.CQI_params(6).modulation = 'QPSK';
LTE_config.CQI_params(6).modulation_order = 2;
LTE_config.CQI_params(6).coding_rate_x_1024 = 602;
LTE_config.CQI_params(6).efficiency = 1.1758;

LTE_config.CQI_params(7).CQI = 7;
LTE_config.CQI_params(7).modulation = '16QAM';
LTE_config.CQI_params(7).modulation_order = 4;
LTE_config.CQI_params(7).coding_rate_x_1024 = 378;
LTE_config.CQI_params(7).efficiency = 1.4766;

LTE_config.CQI_params(8).CQI = 8;
LTE_config.CQI_params(8).modulation = '16QAM';
LTE_config.CQI_params(8).modulation_order = 4;
LTE_config.CQI_params(8).coding_rate_x_1024 = 490;
LTE_config.CQI_params(8).efficiency = 1.9141;

LTE_config.CQI_params(9).CQI = 9;
LTE_config.CQI_params(9).modulation = '16QAM';
LTE_config.CQI_params(9).modulation_order = 4;
LTE_config.CQI_params(9).coding_rate_x_1024 = 616;
LTE_config.CQI_params(9).efficiency = 2.4063;

LTE_config.CQI_params(10).CQI = 10;
LTE_config.CQI_params(10).modulation = '64QAM';
LTE_config.CQI_params(10).modulation_order = 6;
LTE_config.CQI_params(10).coding_rate_x_1024 = 466;
LTE_config.CQI_params(10).efficiency = 2.7305;

LTE_config.CQI_params(11).CQI = 11;
LTE_config.CQI_params(11).modulation = '64QAM';
LTE_config.CQI_params(11).modulation_order = 6;
LTE_config.CQI_params(11).coding_rate_x_1024 = 567;
LTE_config.CQI_params(11).efficiency = 3.3223;

LTE_config.CQI_params(12).CQI = 12;
LTE_config.CQI_params(12).modulation = '64QAM';
LTE_config.CQI_params(12).modulation_order = 6;
LTE_config.CQI_params(12).coding_rate_x_1024 = 666;
LTE_config.CQI_params(12).efficiency = 3.9023;

LTE_config.CQI_params(13).CQI = 13;
LTE_config.CQI_params(13).modulation = '64QAM';
LTE_config.CQI_params(13).modulation_order = 6;
LTE_config.CQI_params(13).coding_rate_x_1024 = 772;
LTE_config.CQI_params(13).efficiency = 4.5234;

LTE_config.CQI_params(14).CQI = 14;
LTE_config.CQI_params(14).modulation = '64QAM';
LTE_config.CQI_params(14).modulation_order = 6;
LTE_config.CQI_params(14).coding_rate_x_1024 = 873;
LTE_config.CQI_params(14).efficiency = 5.1152;

LTE_config.CQI_params(15).CQI = 15;
LTE_config.CQI_params(15).modulation = '64QAM';
LTE_config.CQI_params(15).modulation_order = 6;
LTE_config.CQI_params(15).coding_rate_x_1024 = 948;
LTE_config.CQI_params(15).efficiency = 5.5547;

%% Plot options
LTE_config.plots.BLER_curves                       = 1;
LTE_config.plots.CQI_mapping                       = 2;
LTE_config.plots.antenna_gain_pattern              = 3;
LTE_config.plots.macroscopic_pathloss              = 4;
LTE_config.plots.macroscopic_pathloss_sector       = 5;
% Leave space for some extra plots
LTE_config.plots.shadow_fading_loss                = 22;
% LTE_config.plots.shadow_fading_loss_histogram      = 23;
% LTE_config.plots.initial_UE_positions              = 24;
LTE_config.plots.user_positions                    = 25;
LTE_config.plots.sector_SINR                       = 26;
LTE_config.plots.sector_SINR_no_shadowing          = 27;
LTE_config.plots.sector_SINR_cdf                   = 28;
LTE_config.plots.sector_spectral_densities         = 29;
LTE_config.plots.sector_spectral_densities2        = 30;
LTE_config.plots.sector_spectral_densities_shadow  = 31;
LTE_config.plots.sector_spectral_densities2_shadow = 32;
LTE_config.plots.capesso_vs_pathloss_models        = 33;
LTE_config.plots.FFR_FR_SINR                       = 50;
LTE_config.plots.FFR_PR_SINR                       = 51;
LTE_config.plots.FFR_FR_sector_assignment          = 52;
LTE_config.plots.FFR_PR_sector_assignment          = 53;
LTE_config.plots.capesso_maps_begin                = 100;

%% Traffic models
if isfield(LTE_config,'traffic_models')
    if LTE_config.traffic_models.usetraffic_model && ~(strcmp(LTE_config.scheduler,'constrained') || strcmp(LTE_config.scheduler,'alpha fair') || strcmp(LTE_config.scheduler,'prop fair traffic')|| strcmp(LTE_config.scheduler,'round robin traffic')|| strcmp(LTE_config.scheduler,'CarScheduler'))
        warning('Traffic models are just supported with the constrained scheduler - deactivating traffic models');
        LTE_config.traffic_models.usetraffic_model = false;
    end
    if isfield(LTE_config.traffic_models,'type')
        if strcmp(LTE_config.traffic_models.type,'MLaner')
            lambda_pois = 0.0031*exp(1.085*log10(LTE_config.traffic_models.av_cell_TP));
            warning('Overruling the number of UEs according to the traffic model');
            randn('state',round(cputime*10));
            rand('state',round(cputime*3*10));
            N_BS = sum(6*(1:LTE_config.nr_eNodeB_rings))+1;
            if LTE_config.traffic_models.user_distribution % generate users from poission distribution
                LTE_config.UE_per_eNodeB = poissrnd(lambda_pois,N_BS,3); % user number for each Sector of the eNodeBs
            else % generate users from uniform distribution
                N_UE = [floor(lambda_pois),ceil(lambda_pois)];
                for b_ = 1:N_BS
                    for s_ = 1:3
                        coin_toss = rand(1);
                        if coin_toss <= lambda_pois-floor(lambda_pois)
                            LTE_config.UE_per_eNodeB(b_,s_) = N_UE(2);
                        else
                            LTE_config.UE_per_eNodeB(b_,s_) = N_UE(1);
                        end
                    end
                end
            end
        end
    end
else
    LTE_config.traffic_models.usetraffic_model = false;
end

%% ICI emulation parameters
if ~isfield(LTE_config,'ICI_emulation')
     LTE_config.ICI_emulation = false;
end
LTE_config.ICI_power = zeros(length(LTE_config.BF_samples),LTE_config.Ntot/6);
if LTE_config.ICI_emulation 
%     if LTE_config.ICI_emulation && (~(LTE_config.tx_mode == 1) && ~(LTE_config.tx_mode == 4) && ~(LTE_config.tx_mode  == 9) && ~(LTE_config.tx_mode  == 8))
%         error('ICI emulation (LTE_config.ICI_emulation) not implemented for this transmission mode')
%     end
    
    fd = LTE_config.UE_speed/299792458*LTE_config.frequency; % Doppler frequency 

    T_bar = 1e-3/(2*LTE_config.N_sym);
    bessel_store = zeros(LTE_config.Ntot,1);
    for kk = 1:LTE_config.fft_points
        bessel_store(kk) = besselj(0,2*pi*fd*kk*T_bar/LTE_config.fft_points);
    end
    ICI_power = zeros(1,LTE_config.Ntot);
    for ii = 1:LTE_config.Ntot % ICI power according to "ON CHANNEL ESTIMATION AND DETECTION FOR MULTICARRIER SIGNALS"
        sum_val = 0;
        for pp = 1:LTE_config.Ntot
            if pp ~= ii
                sum_val = sum_val + LTE_config.fft_points;
                temp_sum = 0;
                for kk = 1:LTE_config.fft_points                    
                    temp_sum = temp_sum + (LTE_config.fft_points-kk)*cos(2*pi*kk*(pp-ii)/LTE_config.fft_points)*bessel_store(kk);
                end
                sum_val = sum_val + 2*temp_sum;
            end
        end
        sum_val = 1/LTE_config.fft_points^2*sum_val;
        ICI_power(ii) = sum_val;
    end
%     figure(2)
%     plot(ICI_power)
    sc_samples = repmat([1:6],length(LTE_config.ICI_power),1) + repmat(6*[0:length(LTE_config.ICI_power)-1].',1,6);
    LTE_config.ICI_power = repmat(mean(ICI_power(sc_samples),2).',length(LTE_config.BF_samples),1);
%     sc_samples = 3:6:LTE_config.Ntot;
%     LTE_config.ICI_power = repmat(ICI_power(sc_samples),length(LTE_config.BF_samples),1);
end

