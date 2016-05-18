classdef PregeneratedFastFading_v2 < handle
% Object that stores a pregenerated Fast Fading trace
% (c) Josep Colom Ikuno, INTHFT, 2008

   properties
       % Length of the trace in seconds (1 value for every TTI, as we
       % assume Block fading, 1s-> 1000 samples)
       trace_length_s
       trace_length_samples
       % LTE bandwidth used [Hz]
       system_bandwidth
       % Channel type (eg. 'PedB')
       channel_type
       % Number of Tx antennas (eg. 2)
       nTX
       % Number of Rx antennas (eg. 2)
       nRX
       % UE speed, used for the tie correlation [m/s]. Eg: 5/3.6 m/s for 5 Km/h
       UE_speed
       % Step (s) between each channel realization
       t_step
       % Step (Hz) between each sampled channel realization
       f_step
       
       % Here the channel coefficients trace matrix
       H_0
       
       % An interfering channel trace (if used)
       H_i
       
       % Secondary channel traces (e.g., for DAS). Assumed to be of the
       % same length as the primary trace (i.e., they share the time
       % indices)
       secondary_traces
   end

   methods
   end
end 
