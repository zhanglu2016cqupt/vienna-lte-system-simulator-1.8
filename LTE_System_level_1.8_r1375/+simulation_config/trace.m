classdef trace < utils.simulatorConfig
    methods (Static)
        function LTE_config = apply_parameters(LTE_config)
            % Load the LTE System Level Simulator config parameters
            % (c) Josep Colom Ikuno, INTHFT, 2013
            % www.nt.tuwien.ac.at
            
            %% Debug options
            LTE_config.debug_level = 1;
            
            %% Plotting options
            LTE_config.show_network = 0;
            
            %% General options
            LTE_config.frequency       = 2.14e9; % Frequency in Hz
            LTE_config.bandwidth       = 10e6;   % Frequency in Hz
            
            LTE_config.nTX           = 1;
            LTE_config.nRX           = 1;
            LTE_config.tx_mode       = 1;
            
            LTE_config.eNodeB_tx_power = 40; % eNodeB's transmit power, in Watts.
            
            LTE_config.always_on     = true; % Controls whether the eNodeB is always radiating power (default and worse-case scenario) or no power is used when no UEs are attached
            
            %% Random number generation options
            LTE_config.seedRandStream = false;
            LTE_config.RandStreamSeed = 0;      % Only used if the latter is set to 'true'
            
            %% Simulation time
            LTE_config.simulation_time_tti = 100; % Simulation time in TTIs
            
            LTE_config.trace_simulation_mode = true;
            LTE_config.trace_filename        = './data_files/UE_pathloss_trace/UE_traces'; % UE pathloss trace
            
            %% Microscale Fading Generation config
            % Microscale fading trace to be used between the eNodeB and its attached UEs.
            LTE_config.channel_model.type = 'winner+'; % 'PedB' 'extPedB' --> the PDP to use
            LTE_config.channel_model.trace_length = 15; % Length of the trace in seconds. Be wary of the size you choose, as it will be loaded in memory.
            LTE_config.channel_model.correlated_fading = true;
            LTE_config.pregenerated_ff_file           = 'auto';
            
            % With this option set to 'true', even if cache is present, the channel trace will be recalculated
            LTE_config.recalculate_fast_fading = false;
            
            LTE_config.UE.receiver_noise_figure = 9;   % Receiver noise figure in dB
            LTE_config.UE.thermal_noise_density = -174;  % Thermal noise density in dBm/Hz (-174 dBm/Hz is the typical value)
            LTE_config.UE_speed                 = 5/3.6; % Channel speed speed. In meters/second: 5 Km/h = 1.38 m/s
            
            %% Scheduler options
            LTE_config.scheduler                 = 'round robin';  % 'round robin', 'best cqi', 'max min', 'max TP', 'resource fair', 'proportional fair' or 'prop fair Sun'
            LTE_config.scheduler_params.fairness = 0.5; % fairness for the variable fairness scheduler
            LTE_config.power_allocation          = 'homogeneous;'; % 'right now no power loading is implemented, so just leave it as 'homogeneous'
            
            %% CQI mapper options
            LTE_config.CQI_mapper.CQI2SNR_method = 1; % 1 to use less conservative method to map from CQI back to SNR
                                                      % uses the value in the middle of the SNR interval corresponding
                                                      % to a CQI instead of the lower boarder value
                                                      % this value will just be used in connection with quantized CQI feedback
            
            %% Uplink channel options
            LTE_config.feedback_channel_delay   = 0; % In TTIs
            LTE_config.unquantized_CQI_feedback = false;
            
            %% Where to save the results
            LTE_config.results_folder = './results';
            LTE_config.results_file   = sprintf('%.0_TTI_trace_scenario',LTE_config.simulation_time_tti);
        end
    end
end

