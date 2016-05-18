classdef sinrAverager < handle
    % Defines the abstract classes needed by a class that implements a SINR
    % averaging method
    % (c) Josep Colom Ikuno, INTHFT, 2008

   properties
   end

   methods (Abstract)
       [effective_SINR_lin effective_SINR_dB] = average(SINR_vector,MCSs,varargin)
       [effective_SINR_lin effective_SINR_dB] = average_codeword(SINR_vector,MCSs,varargin)
       [effective_SINR_lin effective_SINR_dB] = average_for_RI(SINR_vector,MCSs,varargin)
   end
end 
