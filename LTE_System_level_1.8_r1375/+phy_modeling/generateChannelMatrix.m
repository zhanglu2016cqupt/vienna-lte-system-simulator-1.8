function [ H_0_RB, H_i_RB ] = generateChannelMatrix( config, debug_level, N_samples )
% Generates the channel matrix per (half) RB for the given input
% configuration
%
% Josep Colom Ikuno, jcolom@nt.tuwien.ac.at, 2013

switch config.channel_type
    case 'winner+'
        channel_factory_H0 = channel_gain_wrappers.winnerChannelFactory(config.system_bandwidth,config.trace_params,N_samples);
        channel_factory_H1 = channel_gain_wrappers.winnerChannelFactory(config.system_bandwidth,config.trace_params,N_samples);
    otherwise
        channel_factory_H0 = channel_gain_wrappers.pdpChannelFactory(config.system_bandwidth,config.trace_params);
        channel_factory_H1 = channel_gain_wrappers.pdpChannelFactory(config.system_bandwidth,config.trace_params);
end
if debug_level>=1
    fprintf('Generating %dx%d channel trace of length %3.2fs\n',config.nTX,config.nRX,ceil(config.trace_length_s));
end
H_trace0 = channel_factory_H0.generate_FF_trace(config.trace_length_s/config.TTI_length);

% Interfering channel trace
if debug_level>=1
    fprintf('Generating %dx%d interfering channel trace of length %3.2fs\n',config.nTX,config.nRX,ceil(config.trace_length_s));
end
H_trace1 = channel_factory_H1.generate_FF_trace(config.trace_length_s/config.TTI_length);


% Note: each MIMO channel is normalized to a mean power of one
H_0_RB = H_trace0.H_RB_samples;
H_i_RB = H_trace1.H_RB_samples;

end

