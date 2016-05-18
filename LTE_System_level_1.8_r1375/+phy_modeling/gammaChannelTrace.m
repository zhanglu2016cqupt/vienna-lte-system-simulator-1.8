classdef gammaChannelTrace
    % Class encapsulates all methods to generate Gamma-fading trace
    % (c) Martin Taranetz, ITC 2013
    
    properties
    end
    
    methods(Static)
        
        function gammaVariates = generate_Gamma_variates(k,theta,number_of_time_samples,number_of_frequency_samples)
            % Throw Variates from a Gamma Distribution with shape parameter k and
            % scale parameter theta
            for jj = 1:number_of_frequency_samples
                i_ = floor(k);
                f_ = k-i_;
                r_ = 0;
                r_ = -log(prod(rand(i_,number_of_time_samples),1));
                if f_>0
                    for ii=1:number_of_time_samples
                        while true
                            w1 = rand^(1/f_);
                            w2 = rand^(1/(1-f_));
                            if w1+w2<=1 break; end
                        end
                        r_(ii)=r_(ii)-log(rand)*w1/(w1+w2);
                    end
                end
                gammaVariates_ = r_*theta;
                gammaVariates(jj,:) = gammaVariates_;
            end
        end
        
        function trace_to_fill = trace_SISO(config)
            % Generate the fading trace for the SISO LTE mode
            % Gamma Fading Parameters
            gammaFading_k    = 1;   % Shape parameter
            gammaFading_t    = 1;   % Scale parameter
            
            % Re-create config params from input
            system_bandwidth = config.system_bandwidth;
            channel_type     = config.channel_type;
            nTX              = config.nTX;
            nRX              = config.nRX;
            trace_length_s   = config.trace_length_s;
            UE_speed         = config.UE_speed;
            t_step           = config.t_step;
            f_step           = config.f_step;
            
            nTTIs            = trace_length_s/t_step;
            nSC_samples      = config.n_RB*2; % 2 Samples per resource block
            
            % chi doesn't exist, as there is only one stream being transmitted
            zeta = ones([nSC_samples,nTTIs]); % Already permuted
            psi  = phy_modeling.gammaChannelTrace.generate_Gamma_variates(gammaFading_k,gammaFading_t, nTTIs,nSC_samples);
            theta= phy_modeling.gammaChannelTrace.generate_Gamma_variates(gammaFading_k,gammaFading_t, nTTIs,nSC_samples);
            
            %% Fill in the output trace object
            trace_to_fill                  = phy_modeling.txModeTrace;
            trace_to_fill.tx_mode          = 1;
            trace_to_fill.trace_length_s   = trace_length_s;
            trace_to_fill.system_bandwidth = system_bandwidth;
            trace_to_fill.channel_type     = channel_type;
            trace_to_fill.nTX              = nTX;
            trace_to_fill.nRX              = nRX;
            trace_to_fill.UE_speed         = UE_speed;
            
            trace_to_fill.trace.zeta  = zeta;
            trace_to_fill.trace.psi   = psi;
            trace_to_fill.trace.theta = theta;
        end
        
        function pregenerated_fast_fading = generate_channel_trace(config)
            
            config.t_step                   = 1e-3;
            config.f_step                   = 15e3*6;
            
            switch config.tx_mode
                case 1
                    SISO_trace = phy_modeling.gammaChannelTrace.trace_SISO(config);
                otherwise
                    error('Tx mode %d for Gamma fading not supported yet',config.tx_mode);
             end
            
             %% Create output fast fading trace
            pregenerated_fast_fading                      = phy_modeling.PregeneratedFastFading;
            pregenerated_fast_fading.trace_length_s       = config.trace_length_s;
            pregenerated_fast_fading.trace_length_samples = config.trace_length_s / 1e-3;
            pregenerated_fast_fading.system_bandwidth     = config.system_bandwidth;
            pregenerated_fast_fading.channel_type         = config.channel_type;
            pregenerated_fast_fading.nTX                  = config.nTX;
            pregenerated_fast_fading.nRX                  = config.nRX;
            pregenerated_fast_fading.UE_speed             = config.UE_speed;
            
            pregenerated_fast_fading.t_step               = config.t_step;
            pregenerated_fast_fading.f_step               = config.f_step;
            
            switch config.tx_mode
                case 1
                    % SISO trace (mode 1)
                    pregenerated_fast_fading.traces{1} = SISO_trace;
                otherwise
            end
             
        end
    end
end