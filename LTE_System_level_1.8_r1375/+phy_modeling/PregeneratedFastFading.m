classdef PregeneratedFastFading < handle
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
       
       % Here the several TX mode traces are stored
       traces
       % 1: Single Antenna
       % 2: Transmit Diversity
       % 3: Open Loop Spatial Multiplexing
       % 4: Closed Loop SM
       
       % 
       generated_from % v2 filename from which this trace was calculated
       source_info    % properties of the file, used to check if it has changed
   end

   methods
   end
end 
