classdef example_capesso < utils.simulatorConfig
    methods (Static)
        function LTE_config = apply_parameters(LTE_config)
            % Load the LTE System Level Simulator config parameters
            % (c) Josep Colom Ikuno, INTHFT, 2008
            % www.nt.tuwien.ac.at
            
            %% Debug options
            LTE_config.debug_level = 1;  % 0=no output
            % 1=basic output
            % 2=extended output
            
            %% Plotting options
            LTE_config.show_network = 0; % 0= show no plots
            % 1= show some plots
            % 2= show ALL plots (moving UEs)
            % 3=plot even all of the pregenerated fast fading
            
            %% General options
            LTE_config.frequency       = 2.14e9;           % Frequency in Hz
            LTE_config.bandwidth       = 20e6;           % Frequency in Hz
            
            LTE_config.nTX           = 2;
            LTE_config.nRX           = 2;
            LTE_config.tx_mode       = 4;
            
            LTE_config.always_on     = false; % Controls whether the eNodeB is always radiating power (default and worse-case scenario) or no power is used when no UEs are attached
            
            %% Random number generation options
            LTE_config.seedRandStream = false;
            LTE_config.RandStreamSeed = 0;      % Only used if the latter is set to 'true'
            
            %% Simulation time
            LTE_config.simulation_time_tti = 100; % Simulation time in TTIs
            
            %% Cache options. Saves the generated eNodeBs, Pathloss map and Shadow fading map to a .mat file
            LTE_config.cache_network = false;
            %LTE_config.network_cache = 'capesso_mka_lte';%'auto';
            LTE_config.network_cache = 'auto';
            LTE_config.delete_pathloss_at_end = false; % Reduces the amount needed to store the traces by deleting the pathloss data and shadow fading data (when present) from the results file
            LTE_config.delete_ff_trace_at_end = true; % Reduces the amount needed to store the traces by deleting the fading parameters trace from the results file
            LTE_config.UE_cache      = false;  % Option to save the user position to a file. This works in the following way:
            %   - cache=true and file exists: read position from file
            %   - cache=true and file does not exist: create UEs and save to cache
            %   - cache=false: do not use cache at all
            LTE_config.UE_cache_file = 'auto';
            
            %% How to generate the network. If the map is loaded, this parameters will be overridden by the loaded map
            % Use network planning tool Capesso by placing 'capesso'
            LTE_config.network_source = 'capesso';
            
            % Configure the network source. Overridden if a pregenerated network pathloss map is used.
            switch LTE_config.network_source
                case 'generated'
                    % Network size
                    LTE_config.inter_eNodeB_distance = 500; % In meters. When the network is generated, this determines the
                    % distance between the eNodeBs.
                    LTE_config.map_resolution = 5; % In meters/pixel. Also the resolution used for initial user creation
                    LTE_config.nr_eNodeB_rings = 0; % Number of eNodeB rings
                    LTE_config.minimum_coupling_loss = 70; % Minimum Coupling Loss: the parameter describing the minimum
                    % loss in signal [dB] between BS and UE or UE and UE in the worst
                    % case and is defined as the minimum distance loss including
                    % antenna gains measured between antenna connectors.
                    % Recommended in TS 36.942 are 70 dB for urban areas, 80 dB for rural.                                              % is
                    
                    % Models to choose
                    % Available are:
                    %  'free space': (more something for testing purposes than to really use it...)
                    %  'cost231'
                    %  'TS36942': Recommended by TS 36.942, subclause 4.5
                    %  'TS25814': Recommended by TS 25.814 (Annex). The same as in HSDPA
                    LTE_config.macroscopic_pathloss_model = 'TS25814';
                    
                    % Additional pathloss model configuration parameters. Will depend on which model is chosen.
                    % Available options are:
                    %  'urban_micro' (COST231)
                    %  'urban_macro' (COST231)
                    %  'suburban_macro' (COST231)
                    %  'urban' (TS36942)
                    %  'rural' (TS36942)
                    LTE_config.macroscopic_pathloss_model_settings.environment = '';
                    
                    % eNodeB settings
                    LTE_config.eNodeB_tx_power = 20; % eNodeB's transmit power, in Watts.
                    % Recommended by TS.36.814 are:
                    % 43 dBm for 1.25, 5 MHz carrier
                    % 46/49 dBm for 10, 20 MHz carrier
                    
                    % Add here cases for other sources (eg. netowrk planning tools)
                case 'capesso'
                    % Define basic parameters needed also in functions different from network generation
                    % All other parameters defined in "LTE_init_generate_capesso_network"
                    
                    %% Example data
                    % warning('Check maps resolution when using data from example maps and Capesso data');
                    LTE_config.map_resolution =  5;  % In meters/pixel. Use for default scenario
                    
                    LTE_config.manually_set_ROI = false; % If needed, false by default
                    LTE_config.rescale_factor =  1; % Resolution rescale factor for ROI
                    LTE_config.macroscopic_pathloss_model = 'Capesso';
                    LTE_config.macroscopic_pathloss_model_settings.environment = '';
                    
                    % Configure Capesso import params
                    LTE_config.capesso_params.planning_tool              = 'capesso';
                    LTE_config.capesso_params.pathloss_data_folder     = './data_files/CapessoExample/pathloss';
                    
                    LTE_config.capesso_params.enable_dtm                 = true;
                    LTE_config.capesso_params.dtm_folder                 = './data_files/CapessoExample/dtm';
                    LTE_config.capesso_params.dtm_file_name              = 'exampleDTM.bil';
                    LTE_config.capesso_params.dtm_hdr_file_name          = 'exampleDTM.hdr';
                    
                    LTE_config.capesso_params.kathrein_antenna_folder    = './data_files/KATHREIN_antenna_files';
                    
                    % CapessoExample
                    LTE_config.capesso_params.cell_atoll_filename        = 'exampleCluster_Cell.txt';
                    LTE_config.capesso_params.site_atoll_filename        = 'exampleCluster_Sites.txt';
                    LTE_config.capesso_params.transmitter_atoll_filename = 'exampleCluster_Transmitter.txt';
                    
                    LTE_config.capesso_params.rx_height                  = 1.5;   % Height of receiver in meter from ground level
                    LTE_config.capesso_params.use_default_tilt_value     = false; % Use common mechanical tilt value for all sites. When not using this option, values will be extracted from Capesso data
                    LTE_config.capesso_params.default_mechanical_tilt    = 0;     % Default mechanical tilt in degree - usage decreases simulation time significantly
                    LTE_config.capesso_params.eNodeB_ROI_increase_factor = 0;   % Increase ROI relative to size
                    
                    % Debugging and plotting
                    LTE_config.capesso_params.plot_antenna_gain_patterns = 1;     % 0 ... false, 1 ... true - don't use boolean values
                    LTE_config.capesso_params.enable_debug_plotting      = false; % This is the extended DEBUG plotting for each of the subfunctions. Probably A LOT of plots!!!
                    LTE_config.capesso_params.number_of_sectors          = 3 ;    % Chosen static for now. In the future, an arbitrary number of sectors/eNodeB may be supported
                otherwise
                    error([LTE_config.network_source ' network source not supported']);
            end
            
            %% Generation of the shadow fading
            LTE_config.shadow_fading_type = 'none'; % Right now only 2D space-correlated shadow fading maps implemented
            
            % Configure the network source
            switch LTE_config.shadow_fading_type
                case 'claussen'
                    LTE_config.shadow_fading_map_resolution = 5; % Recommended value for 8 neighbors
                    LTE_config.shadow_fading_n_neighbors    = 8; % Either 4 or 8
                    LTE_config.shadow_fading_mean           = 0;
                    LTE_config.shadow_fading_sd             = 10;
                    LTE_config.r_eNodeBs                    = 0.5; % inter-site shadow fading correlation
                case 'none'
                    % do not use shadow fading
                otherwise
                    error([LTE_config.shadow_fading_type ' shadow fading type not supported']);
            end
            
            %% Microscale Fading Generation config
            % Microscale fading trace to be used between the eNodeB and its attached UEs.
            LTE_config.channel_model.type = 'winner+'; % 'PedB' 'extPedB' --> the PDP to use
            LTE_config.channel_model.trace_length = 10; % Length of the trace in seconds. Be wary of the size you choose, as it will be loaded in memory.
            LTE_config.channel_model.correlated_fading = true;
            LTE_config.pregenerated_ff_file           = 'auto';
            
            % With this option set to 'true', even if cache is present, the channel trace will be recalculated
            LTE_config.recalculate_fast_fading = false;
            
            %% UE (users) settings
            % note that for reducing trace sizes, the UE_id is stored as a uint16, so
            % up to 65535 users in total are supported. To change that, modify the scheduler class.
            
            LTE_config.UE_distribution                          = 'traffic map';
            LTE_config.traffic_map_config.traffic_map_upscaling = 100;
            LTE_config.traffic_map_config.udtm_folder           = './data_files/CapessoExample/trafficmaps';
            LTE_config.traffic_map_config.udtm_filename         = 'trafficmap_generic';
            
            LTE_config.UE.receiver_noise_figure   = 9;    % Receiver noise figure in dB
            % LTE_config.UE.thermal_noise_density = -131.59; % Thermal noise density in dBm/Hz (-174 dBm/Hz is the typical value) - use this (-131.59) for the SISO LL-SL comparison
            % LTE_config.UE.thermal_noise_density = -134.89; % Thermal noise density in dBm/Hz (-174 dBm/Hz is the typical value) - use this (-134.89) for the MIMO LL-SL comparison
            LTE_config.UE.thermal_noise_density = -174;  % Thermal noise density in dBm/Hz (-174 dBm/Hz is the typical value)
            LTE_config.UE_speed                 = 5/3.6; % Speed at which the UEs move. In meters/second: 5 Km/h = 1.38 m/s
            
            %% Scheduler options
            LTE_config.scheduler        = 'round robin';  % 'round robin', 'best cqi', 'max min', 'max TP', 'resource fair', 'proportional fair' or 'prop fair Sun'
            LTE_config.scheduler_params.fairness = 0.5; % fairness for the variable fairness scheduler
            % LTE_config.scheduler_params.alpha = 1;          % Additional configuration parameters for the scheduler (defau
            % LTE_config.scheduler_params.beta  = 1;
            LTE_config.power_allocation = 'homogeneous;'; % 'right now no power loading is implemented, so just leave it as 'homogeneous'
            
            %% CQI mapper options
            LTE_config.CQI_mapper.CQI2SNR_method = 1; % 1 to use less conservative method to map from CQI back to SNR
            % uses the value in the middle of the SNR interval corresponding to a CQI instead of the lower boarder value
            % this value will just be used in connection with quantized CQI feedback
            
            
            %% Uplink channel options
            LTE_config.feedback_channel_delay = 3; % In TTIs
            LTE_config.unquantized_CQI_feedback = false;
            
            %% Where to save the results
            LTE_config.results_folder         = './results';
            LTE_config.results_file           = 'auto'; % NOTE: 'auto' assigns a filename automatically
            
            %% Values that should not be changed
            LTE_config.antenna_azimuth_offsett = 30;    % This controls the antenna layout that will be generated. 0 degrees generates hexagonal cells,
            % while 30 degrees hexagonal sectors.
        end
    end
end

