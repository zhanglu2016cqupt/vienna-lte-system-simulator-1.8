function fading_trace = LTE_init_get_microscale_fading_SL_trace(LTE_config)
% Generate the fading parameters that model the fast (microscale) fading at system level.
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at.
% (c) 2009/2011 by INTHFT
% www.nt.tuwien.ac.at

%% Config

% Possible Tx modes are:
%   1: Single Antenna
%   2: Transmit Diversity
%   3: Open Loop Spatial Multiplexing
%   4: Closed Loop SM
% Number of antenna ports can either be 2 or 4
% Codebook index specifies the codebook index as of TS.36.211 (for closed loop SM)
% nLayers specifies how many layers (symbols) are going to be transmitted.
% Either 1, 2, 3 or 4

standalone = false;

if standalone
    LTE_config.debug_level            = 1;
    % Initial config (Winner II+)
    config.system_bandwidth           = 1.4e6;
    config.channel_type               = 'winner+';
    config.nTX                        = 2;
    config.nRX                        = 2;
    config.trace_length_s             = 1;
    config.UE_speed                   = 400/3.6; % converted to m/s. High speed ~ uncorrelated
    config.parallel_toolbox_installed = true; % change it to false for testing purposes, as if not you will not be able to debug properly
    config.feedback_channel_delay     = 0;
    config.correlated_fading          = true;
    config.f                          = 2.1400e+009;
    config.TTI_length                 = 1.0000e-003;
    config.tx_mode                    = 4; % OLSM
    
    winner_config                     = load('winner_trace_config');
    config.trace_params               = winner_config.config_trace_params;
    config.trace_params.speed         = config.UE_speed;
else
    % Initial config (simulator-linked)
    config.system_bandwidth           = LTE_config.bandwidth;
    config.channel_type               = LTE_config.channel_model.type;
    config.nTX                        = LTE_config.nTX;
    config.nRX                        = LTE_config.nRX;
    config.trace_length_s             = LTE_config.channel_model.trace_length;
    config.UE_speed                   = LTE_config.UE_speed; % converted to m/s
    config.parallel_toolbox_installed = LTE_config.parallel_toolbox_installed; % change it to false for testing purposes, as if not you will not be able to debug properly
    config.feedback_channel_delay     = LTE_config.feedback_channel_delay;
    config.correlated_fading          = LTE_config.channel_model.correlated_fading;
    config.f                          = LTE_config.frequency;
    config.trace_params               = LTE_config.trace_params;
    config.TTI_length                 = LTE_config.TTI_length;
    config.tx_mode                    = LTE_config.tx_mode;
    config.non_parallel_channel_trace = LTE_config.non_parallel_channel_trace;
    config.trace_params.Ns            = LTE_config.N_sym*2; % number of OFDM symbols per subframe
    
    if LTE_config.BF_short
        config.BF_short = true;
        config.trace_params.BF_samples = LTE_config.BF_samples;
    else
        config.BF_short = false;
        config.trace_params.BF_samples = 0;
    end
    
    % Option for wideband precoding
    if LTE_config.tx_mode==4 || LTE_config.tx_mode==6 || LTE_config.tx_mode==5 || LTE_config.tx_mode==9
        config.wideband_precoding = LTE_config.wideband_precoding;
    end
end

% sigma_n2 = 10^((LTE_config.UE.receiver_noise_figure + LTE_config.UE.thermal_noise_density)/10)/1000;    % Receiver noise variance in Watt

% We now have all of the possible precoding combinations stored
precoding_configs = phy_modeling.miscUtils.get_all_precoding_combinations;

switch LTE_config.trace_version
    case 'v1'
        % Trace mode with precalculated precoding
        % Since v2 incorporates a full channel model it is used as
        % reference. If it exists, it is loaded, if not, it is generated.
        % This way, all tx_modes with the same antenna configuration access
        % the same precalculated channel, allowing for reproducibility
        
        if LTE_config.parallel_network
            t = getCurrentTask();
            if ~isempty(t)
                par_ID = num2str(t.ID);
            else
                par_ID = '';
            end
        else
            par_ID = '';
        end
        
        filename_v2 = fullfile(LTE_config.default_channel_trace_folder,...
            sprintf('%s_BS-%dTX%s_UE-%dRX_%dMHz_%.1fKmph_%ds_%s%s%s_%s.mat',...
            LTE_config.channel_model.type,...
            LTE_config.nTX,...
            LTE_config.RRH_string,...
            LTE_config.nRX,...
            LTE_config.bandwidth/1e6,...
            LTE_config.UE_speed*3.6,...
            LTE_config.channel_model.trace_length,...
            'v2',...
            LTE_config.extra_params_hash,...
            LTE_config.extra_params_trace_id,...
            par_ID));
        
        ff_file_exists = exist(filename_v2,'file');    
        if LTE_config.recalculate_fast_fading || (~ff_file_exists && ~LTE_config.recalculate_fast_fading)
            fprintf('Generating master channel trace first...\n')
            pregenerated_ff = generate_trace_v2(LTE_config, config);
            try
                if exist(filename_v2,'file')
                    throw(MException('LTEsim:cacheExists', 'The cache file was concurrently generated during another simulation run'));
                end
                save(filename_v2,'pregenerated_ff','-v7.3');
            catch err
                fprintf('Channel trace could not be saved. If needed, it will be generated again in the next run (%s).\n',err.message);
            end
        else
            fprintf('Loading master channel trace\n');
            load(filename_v2,'pregenerated_ff');
        end
        fileinfo = dir(filename_v2);
        % Old implementation: Calculate at random
        %[H_0_RB, H_i_RB] = phy_modeling.generateChannelMatrix(config,LTE_config.debug_level);
        
        % load channel from master trace
        H_0_RB_temp = pregenerated_ff.H_0;
        H_0_RB = permute(H_0_RB_temp,[1 2 4 5 3]);
        H_i_RB_temp = pregenerated_ff.H_i;
        H_i_RB = permute(H_i_RB_temp,[1 2 4 5 3]);
        fading_trace     = phy_modeling.channelTraceFactory_v1.generate_channel_trace(config,precoding_configs,H_0_RB,H_i_RB);
        fading_trace.generated_from = filename_v2;
        fading_trace.source_info = fileinfo;
    case 'v2'
        fading_trace = generate_trace_v2(LTE_config, config);
    otherwise
        error('Trace mode %s not implemented',LTE_config.trace_version);
end
end

function fading_trace = generate_trace_v2(LTE_config, config)
% Trace mode with runtime precoding. Implemented like this so as not to "break" the old code

% Main channel trace for the target and interfering channels
[H_0_RB, H_i_RB] = phy_modeling.generateChannelMatrix(...
    config,...    
    LTE_config.debug_level,...
    length(LTE_config.BF_samples));
fading_trace     = phy_modeling.channelTraceFactory_v2.generate_channel_trace(...
    config,...
    H_0_RB,H_i_RB);

% Secondary traces, such as the ones from Remote Radio Heads (RRHs)
secondary_configs = phy_modeling.channelTraceFactory_v2.generate_secondary_trace_configs(LTE_config,config);
for t_=1:length(secondary_configs)
    [H_0_RB_s, H_i_RB_s] = phy_modeling.generateChannelMatrix(...
        secondary_configs(t_),...
        LTE_config.debug_level);
    if t_==1
        fading_trace.secondary_traces = phy_modeling.channelTraceFactory_v2.generate_channel_trace(...
            secondary_configs(t_),...
            H_0_RB_s,H_i_RB_s);
    else
        fading_trace.secondary_traces(t_) = phy_modeling.channelTraceFactory_v2.generate_channel_trace(...
            secondary_configs(t_),...
            H_0_RB_s,H_i_RB_s);
    end
end
end

% function fading_trace = generate_trace_BF(LTE_config, config)
% % Trace mode with runtime precoding. Implemented like this so as not to "break" the old code
% 
% if strcmp(config.channel_type,'winner+')
%     error('Short block fading not yet implemented for winner model')
% end
% 
% % Main channel trace for the target and interfering channels
% [H_0_RB, H_i_RB] = phy_modeling.generateChannelMatrix(...
%     config,...
%     LTE_config.debug_level);
% fading_trace     = phy_modeling.channelTraceFactory_v2.generate_channel_trace(...
%     config,...
%     H_0_RB,H_i_RB);
% 
% % Secondary traces, such as the ones from Remote Radio Heads (RRHs)
% secondary_configs = phy_modeling.channelTraceFactory_v2.generate_secondary_trace_configs(LTE_config,config);
% for t_=1:length(secondary_configs)
%     [H_0_RB_s, H_i_RB_s] = phy_modeling.generateChannelMatrix(...
%         secondary_configs(t_),...
%         LTE_config.debug_level);
%     if t_==1
%         fading_trace.secondary_traces = phy_modeling.channelTraceFactory_v2.generate_channel_trace(...
%             secondary_configs(t_),...
%             H_0_RB_s,H_i_RB_s);
%     else
%         fading_trace.secondary_traces(t_) = phy_modeling.channelTraceFactory_v2.generate_channel_trace(...
%             secondary_configs(t_),...
%             H_0_RB_s,H_i_RB_s);
%     end
% end
% end