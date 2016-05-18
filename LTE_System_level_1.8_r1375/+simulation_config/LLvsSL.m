classdef LLvsSL < utils.simulatorConfig
    methods (Static)
        function LTE_config = apply_parameters(LTE_config)
            % Load the LTE System Level Simulator config parameters
            % (c) Josep Colom Ikuno, INTHFT, 2008
            % www.nt.tuwien.ac.at
            
            %% Debug options
            LTE_config.debug_level = 0; 
            % 0=no output
            % 1=basic output
            % 2=extended output
            
            %% Plotting options
            LTE_config.show_network = 0; % 0= show no plots
            % 1= show some plots
            % 2= show ALL plots (moving UEs)
            % 3=plot even all of the pregenerated fast fading
            
            %% General options
            LTE_config.frequency       = 2.14e9;          % Frequency in Hz
            LTE_config.bandwidth       = 1.4e6;           % Frequency in Hz
            
            LTE_config.nTX           = 2;
            LTE_config.nRX           = 2;
            LTE_config.tx_mode       = 4;
            
            LTE_config.always_on     = true; % Controls whether the eNodeB is always radiating power (default and worse-case scenario) or no power is used when no UEs are attached
            
            %% Random number generation options
            LTE_config.seedRandStream = true;
            LTE_config.RandStreamSeed = 0;      % Only used if the latter is set to 'true'
            
            %% Simulation time
            LTE_config.simulation_time_tti = 10; % Simulation time in TTIs
            
            LTE_config.latency_time_scale = 10;   % Number of TTIs used to calculate the average throughput (exponential filtering)
            
            %% Cache options. Saves the generated eNodeBs, Pathloss map and Shadow fading map to a .mat file
            LTE_config.cache_network = true;
            LTE_config.network_cache = 'auto';
            LTE_config.delete_ff_trace_at_end = true; % Reduces the amount needed to store the traces by deleting the fading parameters trace from the results file
            LTE_config.UE_cache      = false;
            LTE_config.UE_cache_file = 'auto';
            
            %% How to generate the network. If the map is loaded, this parameters will be overridden by the loaded map
            % Use network planning tool Capesso by placing 'capesso'
            LTE_config.network_source = 'generated';
            LTE_config.network_geometry = 'regular_hexagonal_grid';
            
            % Configure the network source. Overridden if a pregenerated network pathloss map is used.
            % Network size
            LTE_config.inter_eNodeB_distance = 1000; % In meters. When the network is generated, this determines the
            % distance between the eNodeBs.
            LTE_config.map_resolution = 1; % In meters/pixel. Also the resolution used for initial user creation
            LTE_config.nr_eNodeB_rings = 0; % Number of eNodeB rings
            LTE_config.minimum_coupling_loss = 70;                                            % is
            
            % Using the free space model with an alpha of three
            LTE_config.macroscopic_pathloss_model       = 'free space';
            LTE_config.macroscopic_pathloss_model_settings.alpha = 3;
            
            % eNodeB settings
            LTE_config.eNodeB_tx_power = 5; % eNodeB's transmit power, in Watts. Set to obtain a given SNR range
            
            % Add here cases for other sources (eg. netowrk planning tools)
            
            %% Generation of the shadow fading
            LTE_config.shadow_fading_type = 'none'; % Right now only 2D space-correlated shadow fading maps implemented
            
            %% Microscale Fading Generation config
            % Microscale fading trace to be used between the eNodeB and its attached UEs.
            LTE_config.channel_model.type = 'VehA'; % 'winner+' 'PedB' 'extPedB' 'TU' --> the PDP to use
            LTE_config.channel_model.trace_length = 1; % Length of the trace in seconds. Be wary of the size you choose, as it will be loaded in memory.
            LTE_config.channel_model.correlated_fading = false;
            LTE_config.pregenerated_ff_file            = 'auto';
            
            % With this option set to 'true', even if cache is present, the channel trace will be recalculated
            LTE_config.recalculate_fast_fading = false;
            
            %% UE (users) settings
            % note that for reducing trace sizes, the UE_id is stored as a uint16, so
            % up to 65535 users in total are supported. To change that, modify the scheduler class.
            
            LTE_config.UE.receiver_noise_figure = 9;    % Receiver noise figure in dB
            % LTE_config.UE.thermal_noise_density = -131.59; % Thermal noise density in dBm/Hz (-174 dBm/Hz is the typical value) - use this (-131.59) for the SISO LL-SL comparison
            % LTE_config.UE.thermal_noise_density = -134.89; % Thermal noise density in dBm/Hz (-174 dBm/Hz is the typical value) - use this (-134.89) for the MIMO LL-SL comparison
            LTE_config.UE.thermal_noise_density = -160;  % Thermal noise density in dBm/Hz (-174 dBm/Hz is the typical value)
            LTE_config.UE_distribution = 'radial';
            LTE_config.UE_distribution_overlap_radial_UEs = true;
            LTE_config.UE_distribution_radii = 1:1:700; % Place by default one ring 50m away from the eNodeB
            LTE_config.UE_distribution_nUEs  = 1;  % PHow many UEs to place on each radius
            LTE_config.UE_speed        = 5/3.6; % Speed at which the UEs move. In meters/second: 5 Km/h = 1.38 m/s
            
            %% eNodeB options
            % LTE_config.antenna.antenna_gain_pattern = 'berger';
            LTE_config.antenna.antenna_gain_pattern = 'omnidirectional'; % As defined in TS 36.942. Identical to Berger, but with a 65° 3dB lobe
            
            %% Scheduler options
            LTE_config.scheduler        = 'round robin';  % 'round robin', 'best cqi', 'max min', 'max TP', 'resource fair' or 'proportional fair'
            LTE_config.power_allocation = 'homogeneous;'; % 'right now no power loading is implemented, so just leave it as 'homogeneous'
            
            %% CQI mapper options
            LTE_config.CQI_mapper.CQI2SNR_method = 1; % 1 to use less conservative method to map from CQI back to SNR
            % uses the value in the middle of the SNR interval corresponding to a CQI instead of the lower boarder value
            % this value will just be used in connection with quantized CQI feedback
            
            %% Uplink channel options
            LTE_config.feedback_channel_delay = 0; % In TTIs
            LTE_config.unquantized_CQI_feedback = false;
            
            %% Where to save the results
            LTE_config.results_folder         = './reproducibility/included_simulation_results/SLvsLL';
            LTE_config.results_file           = 'auto'; % NOTE: 'auto' assigns a filename automatically
            
            %% Values that should not be changed
            LTE_config.sector_azimuths = 0; % Just one sector in a single site
            
            LTE_config.compact_results_file = 2;
            
            %% Walking model
            LTE_config.UE.walk = 'SLvsLL';
        end
    end
end

