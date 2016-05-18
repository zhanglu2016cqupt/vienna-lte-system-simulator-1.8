classdef txModeTrace
    % Defines a trace that models the physical layer for a given TX mode.
    % Since the number of actual parameters stored in the trace, just basic
    % info is compulsory stored in the actual object. The parameters being
    % then stored in a struct.
    % Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
    % (c) 2009 by INTHFT
    % www.nt.tuwien.ac.at
    
    properties
        % For which Tx mode this trace is generated
        % 1: Single Antenna
        % 2: Transmit Diversity
        % 3: Open Loop Spatial Multiplexing
        % 4: Closed Loop SM
        tx_mode
        % Length of the trace in seconds (1 value for every TTI, as we
        % assume Block fading, 1s-> 1000 samples)
        trace_length_s
        % LTE bandwidth [Hz]
        system_bandwidth
        % Channel type (eg. 'PedB')
        channel_type
        % Number of Tx antennas (eg. 2)
        nTX
        % Number of Rx antennas (eg. 2)
        nRX
        % UE speed, used for the time correlation [m/s]. Eg: 5/3.6 m/s for 5 Km/h
        UE_speed
        
        % The struct which contains the actual trace values. Its structure
        % will depend on the actual mode being modeled.
        trace
    end
    
    methods
    end
    
end

