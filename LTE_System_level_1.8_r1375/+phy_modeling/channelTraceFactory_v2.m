classdef channelTraceFactory_v2
    % Channel trace for run-time application of the precoding. Just stores
    % the channel coefficients.
    % (c) Josep Colom Ikuno, INTHFT, 2011
    
    properties
    end
    
    methods(Static)
        
        function pregenerated_fast_fading = generate_channel_trace(config,H_0,H_i)
            %% Create output fast fading trace
            pregenerated_fast_fading                      = phy_modeling.PregeneratedFastFading_v2;
            pregenerated_fast_fading.trace_length_s       = config.trace_length_s;
            pregenerated_fast_fading.trace_length_samples = config.trace_length_s / 1e-3;
            pregenerated_fast_fading.system_bandwidth     = config.system_bandwidth;
            pregenerated_fast_fading.channel_type         = config.channel_type;
            pregenerated_fast_fading.nTX                  = config.nTX;
            pregenerated_fast_fading.nRX                  = config.nRX;
            pregenerated_fast_fading.UE_speed             = config.UE_speed;
            
            pregenerated_fast_fading.t_step               = 1e-3;
            pregenerated_fast_fading.f_step               = 15e3*6;

            pregenerated_fast_fading.H_0                  = permute(H_0,[1 2 5 3 4]);
            pregenerated_fast_fading.H_i                  = permute(H_i,[1 2 5 3 4]);
        end
        
        function secondary_configs = generate_secondary_trace_configs(LTE_config, primary_trace_config)
            % Generate the config structs of the secondary traces (e.g., RRHs)
            
            % Now, we have to calculate all of the TX-RX combinations which
            % are possible in this scenario. Other cases would be if for
            % example the could be eNodeBs with different antenna counts.
            % Right now, only the DAS case exists, but this would need to
            % be updated if more functionality would be implemented in the
            % simulator.
            
            if ~LTE_config.RRHs_enabled
                % No secondary traces
                secondary_configs = [];
            else
                secondary_configs = primary_trace_config;
                RRH_nTX           = LTE_config.RRH.nTX;
                
                % Fill in config struct
                secondary_configs.nTX = RRH_nTX;
                secondary_configs.trace_params.BS_config.nTx = RRH_nTX;
            end
        end
    end
end