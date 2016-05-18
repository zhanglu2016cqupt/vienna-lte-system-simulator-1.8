classdef ueSpecificTraces < handle
% Needed to import LL simulation results
% Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

properties
    % Fields updated after every TTI
    ACK
    rv_idx
    biterrors_coded
    biterrors_uncoded
    blocksize_coded
    blocksize_uncoded
    % CQI and RI mapping?
    FER_coded
    FER_uncoded
    throughput_coded
    throughput_uncoded
    
    % Aggregates calculated after the simulation is finisched
    BER_coded           % coded BER (each stream)
    BER_uncoded         % uncoded BER (each stream)
    BER_coded_overall   % coded BER (overall)
    BER_uncoded_overall % uncoded BER (overall)
    BLER                % BLER (each stream)
    BLER_overall        % BLER (overall)

    % Used to signal what entries are valid
    used_codewords
end

   methods
       % Class contructor. Data preallocation
       function obj = ueSpecificTraces(N_subframes,SNR_vector_length,maxStreams)
           obj.ACK                = false(N_subframes,SNR_vector_length,maxStreams);          % true/false -> ACK of the received subframes for BLER calculation
           obj.rv_idx             = zeros(N_subframes,SNR_vector_length,maxStreams,'uint8');  % 0-255      -> redundancy version index of the received subframes
           obj.biterrors_coded    = zeros(N_subframes,SNR_vector_length,maxStreams,'uint32'); % 0-4294967295
           obj.biterrors_uncoded  = zeros(N_subframes,SNR_vector_length,maxStreams,'uint32'); % 0-4294967295
           obj.blocksize_coded    = zeros(N_subframes,SNR_vector_length,maxStreams,'uint32'); % 0-4294967295
           obj.blocksize_uncoded  = zeros(N_subframes,SNR_vector_length,maxStreams,'uint32'); % 0-4294967295
           obj.FER_coded          = false(N_subframes,SNR_vector_length,maxStreams);          % This is actually ~ACK
           obj.FER_uncoded        = false(N_subframes,SNR_vector_length,maxStreams);          %
           obj.throughput_coded   = zeros(N_subframes,SNR_vector_length,maxStreams,'uint32'); % 0-4294967295
           obj.throughput_uncoded = zeros(N_subframes,SNR_vector_length,maxStreams,'uint32'); % 0-4294967295
           obj.used_codewords     = false(N_subframes,SNR_vector_length,maxStreams);          % What codewords were used
       end
   end
end 
