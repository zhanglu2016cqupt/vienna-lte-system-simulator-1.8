classdef hex_grid_tilted_traffic < utils.simulatorConfig
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
            LTE_config.frequency       = 2.14e9;         % Frequency in Hz
            LTE_config.bandwidth       = 20e6;           % Frequency in Hz
            
            LTE_config.nTX           = 2;
            LTE_config.nRX           = 2;
            LTE_config.tx_mode       = 4; % CLSM
            
            LTE_config.always_on     = true; % Controls whether the eNodeB is always radiating power (default and worse-case scenario) or no power is used when no UEs are attached
            
            %% Random number generation options
            LTE_config.seedRandStream = false;
            LTE_config.RandStreamSeed = 0;      % Only used if the latter is set to 'true'
            
            %% Simulation time
            LTE_config.simulation_time_tti = 100; % Simulation time in TTIs
            
            %% Cache options. Saves the generated eNodeBs, Pathloss map and Shadow fading map to a .mat file
            LTE_config.cache_network = true;
            LTE_config.network_cache = 'auto';
            LTE_config.delete_ff_trace_at_end = true; % Reduces the amount needed to store the traces by deleting the fading parameters trace from the results file
            LTE_config.delete_pathloss_at_end = true; % Further reduces the amount of space needed to store the traces by deleting the pathloss maps from the results file
            LTE_config.UE_cache      = true;  % Option to save the user position to a file. This works in the following way:
            %   - cache=true and file exists: read position from file
            %   - cache=true and file does not exist: create UEs and save to cache
            %   - cache=false: do not use cache at all
            LTE_config.UE_cache_file = 'auto';
            
            %% How to generate the network. If the map is loaded, this parameters will be overridden by the loaded map
            % Use network planning tool Capesso by placing 'capesso'
            LTE_config.network_source = 'generated';
            
            % Configure the network source. Overridden if a pregenerated network pathloss map is used.
            switch LTE_config.network_source
                case 'generated'
                    % Network size
                    LTE_config.inter_eNodeB_distance = 500; % In meters. When the network is generated, this determines the
                    % distance between the eNodeBs.
                    LTE_config.map_resolution = 10; % In meters/pixel. Also the resolution used for initial user creation
                    LTE_config.nr_eNodeB_rings = 2; % Number of eNodeB rings
                    LTE_config.minimum_coupling_loss = 70; % Minimum Coupling Loss: the parameter describing the minimum
                    % loss in signal [dB] between BS and UE or UE and UE in the worst
                    % case and is defined as the minimum distance loss including
                    % antenna gains measured between antenna connectors.
                    % Recommended in TS 36.942 are 70 dB for urban areas, 80 dB for rural.
                    
                    % Models to choose
                    % Available are:
                    %  'free space'
                    %  'cost231'
                    %  'TS36942': Recommended by TS 36.942, subclause 4.5
                    %  'TS25814': Recommended by TS 25.814 (Annex). The same as in HSDPA
                    LTE_config.macroscopic_pathloss_model = 'TS36942';
                    
                    % Additional pathloss model configuration parameters. Will depend on which model is chosen.
                    % Available options are:
                    %  'urban_micro' (COST231)
                    %  'urban_macro' (COST231)
                    %  'suburban_macro' (COST231)
                    %  'urban' (TS36942)
                    %  'rural' (TS36942)
                    LTE_config.macroscopic_pathloss_model_settings.environment = 'urban';
                    
                    % eNodeB settings
                    LTE_config.eNodeB_tx_power = 40; % eNodeB's transmit power, in Watts.
                    % Recommended by TS.36.814 are:
                    % 43 dBm for 1.25, 5 MHz carrier
                    % 46/49 dBm for 10, 20 MHz carrier
                    
                    % Add here cases for other sources (eg. network planning tools)
                case 'capesso'
                otherwise
                    error([LTE_config.network_source ' network source not supported']);
            end
            
            %% Generation of the shadow fading
            LTE_config.shadow_fading_type = 'claussen'; % Right now only 2D space-correlated shadow fading maps implemented
            
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
            LTE_config.channel_model.type = 'winner+'; % 'winner+' 'PedB' 'extPedB' 'TU' --> the PDP to use
            LTE_config.channel_model.trace_length = 10; % Length of the trace in seconds. Be wary of the size you choose, as it will be loaded in memory.
            LTE_config.channel_model.correlated_fading = true;
            LTE_config.pregenerated_ff_file            = 'auto';
            
            % With this option set to 'true', even if cache is present, the channel trace will be recalculated
            LTE_config.recalculate_fast_fading = false;
            
            %% UE (users) settings
            
            LTE_config.UE.receiver_noise_figure = 9;    % Receiver noise figure in dB
            % LTE_config.UE.thermal_noise_density = -131.59; % Thermal noise density in dBm/Hz (-174 dBm/Hz is the typical value) - use this (-131.59) for the SISO LL-SL comparison
            % LTE_config.UE.thermal_noise_density = -134.89; % Thermal noise density in dBm/Hz (-174 dBm/Hz is the typical value) - use this (-134.89) for the MIMO LL-SL comparison
            LTE_config.UE.thermal_noise_density = -174;  % Thermal noise density in dBm/Hz (-174 dBm/Hz is the typical value)
            
            LTE_config.UE_distribution = 'constant UEs per cell';
            LTE_config.UE_per_eNodeB   = 20;    % number of users per eNodeB sector (calculates it for the center sector and applies this user density to the other sectors)
            LTE_config.UE_speed        = 5/3.6; % Speed at which the UEs move. In meters/second: 5 Km/h = 1.38 m/s
            
            %% eNodeB options
            % LTE_config.antenna.antenna_gain_pattern = 'berger';
            LTE_config.antenna.antenna_gain_pattern = 'kathreinTSAntenna';  % Additional parameters needed
            % LTE_config.antenna.antenna_gain_pattern = 'TS 36.942'; % As defined in TS 36.942. Identical to Berger, but with a 65° 3dB lobe
            
            % LTE_config.antenna.max_antenna_gain = 14; % For a berger antenna
            % LTE_config.antenna.max_antenna_gain = 15; % LTE antenna, rural area (900 MHz)
            LTE_config.antenna.max_antenna_gain = 15; % LTE antenna, urban area (2000 MHz)
            % LTE_config.antenna.max_antenna_gain = 12; % LTE antenna, urban area (900 MHz)
            
            %% Additional parameters needed when using 'kathreinTSAntenna'
            if strcmp(LTE_config.antenna.antenna_gain_pattern, 'kathreinTSAntenna')
                % Default values used for each site
                LTE_config.site_altitude                = 0;     % Altiude of site [m]
                LTE_config.site_height                 = 20;    % Height of site  [m]
                LTE_config.rx_height                   = 1.5;   % Receiver height [m]
                LTE_config.antenna.mechanical_downtilt = 0; % [°]
                LTE_config.antenna.electrical_downtilt = 8; % [°]
                LTE_config.antenna.kathrein_antenna_folder = './data_files/KATHREIN_antenna_files';
                % Kathrein Antenna Type : Extensions possible
                LTE_config.antenna.antenna_type = '742212';
                % LTE_config.antenna.antenna_type = '742215';
                LTE_config.antenna.frequency = 2140;
            end
            
            %% Scheduler options
            LTE_config.scheduler        = 'prop fair traffic';  % 'round robin', 'best cqi', 'max min', 'max TP', 'resource fair', or 'prop fair Sun'
            LTE_config.scheduler_params.fairness = 0.5; % fairness for the variable fairness scheduler
            % LTE_config.scheduler_params.alpha = 1;          % Additional configuration parameters for the scheduler
            % LTE_config.scheduler_params.beta  = 1;
            LTE_config.scheduler_params.av_window = 25; % Number of TTIs used to calculate the average throughput (exponential filtering)
            % See: P. Viswanath, D. Tse and R.
            % Laroia, "Opportunistic Beamforming using Dumb Antennas", IEEE Transactions on Information Theory, vol. 48(6), June, 2002.
            LTE_config.power_allocation = 'homogeneous;'; % 'right now no power loading is implemented, so just leave it as 'homogeneous'
            
            %% CQI mapper options
            LTE_config.CQI_mapper.CQI2SNR_method = 1; % 1 to use less conservative method to map from CQI back to SNR
            % uses the value in the middle of the SNR interval corresponding to a CQI instead of the lower boarder value
            % this value will just be used in connection with quantized CQI feedback
            
            %% Uplink channel options
            LTE_config.feedback_channel_delay   = 3; % In TTIs
            LTE_config.unquantized_CQI_feedback = false;
            
            %% Where to save the results
            LTE_config.results_folder         = './results';
            LTE_config.results_file           = 'auto'; % NOTE: 'auto' assigns a filename automatically
            
            %% Values that should not be changed
            LTE_config.antenna_azimuth_offsett = 30;    % This controls the antenna layout that will be generated. 0 degrees generates hexagonal cells,
            % while 30 degrees hexagonal sectors.
            
            %% Traffic models
            LTE_config.traffic_models.usetraffic_model = true;
            LTE_config.traffic_models.type = '3GPP';
        end
    end
end

