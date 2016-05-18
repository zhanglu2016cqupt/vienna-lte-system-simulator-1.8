classdef eesmAverager < utils.sinrAverager
    % Implements an EESM averager. Given a set of predefined \beta values
    % and the Modulation and Coding Scheme (MCS) used, returns an effective
    % SINR value.
    % (c) Josep Colom Ikuno, INTHFT, 2008

   properties
       % Place where we will store the beta value for each possible MCS
       beta_tables
       % MCSs corresponding to the previous beta values
       MCS_values
       % Indicate which MCSs are valid
       isvalid
   end

   methods
       function obj = eesmAverager(beta_values,MCS_values)
           % Input checking
           if length(beta_values) ~= length(MCS_values)
               error('The vector containg the beta values and the vector specifying the corresponging MCSs are of different length');
           end
           if min(MCS_values)<0
               error('The minimum MCS cannot be lower than 0');
           end
           
           % Translate the MCSs so we can put them in a lookup table.
           % Assume that the first MCS is 0
           for i_=1:length(beta_values)
               MCS_value_idx = MCS_values(i_) + 1;
               obj.MCS_values(MCS_value_idx)  = MCS_values(i_);
               obj.beta_tables(MCS_value_idx) = beta_values(i_);
               obj.isvalid(MCS_value_idx)     = true;
           end
       end
       
       % Average the give SINR vector. varargin contains the following:
       %   - MCS -> values in the range 0:15
       % ALL INPUT VALUES IN LINEAR unless an extra parameter is set to "true"
       % Allows for calculations with multiple beta values
       function [effective_SINR_lin effective_SINR_dB] = average(obj,SINR_vector,MCSs)
           
           if isempty(varargin)
               input_in_dB  = false;
           else
               input_in_dB = varargin{1};
           end

           % Input should be linear
           if input_in_dB
               SINR_vector = 10.^(SINR_vector/10);
           end
           
           % Ensure matrix dimensions
           SINR_vector = SINR_vector(:);
           MCS_idx = MCSs+1;
           
           if min(MCS_idx)<1
               error('CQI cannot be lower than 0');
           elseif max(MCS_idx)>length(obj.isvalid)
               error('CQI cannot be higher than %d',length(obj.isvalid)-1);
           elseif sum(~obj.isvalid(MCS_idx))~=0
               error('SINR averaging not defined for all of the CQIs');
           end
           
           betas = obj.beta_tables(MCS_idx);
           % Needed due to some strange handling of vectors from Matlab:
           % SINR_vector(ja) results in a column vector, while betas(jb) a
           % row.
           if length(betas)>1 && length(SINR_vector)>1
               % Change to allow to calculate EESMs for multiple betas in a single pass (needed for the scheduler)
               na = length(SINR_vector);
               nb = length(betas);
               [ja,jb] = meshgrid(1:na,1:nb);
               effective_SINR = -(betas.').*log(mean(exp(-SINR_vector(ja)./betas(jb)),2));
           elseif length(SINR_vector)==1
               effective_SINR = SINR_vector; % there is nothing to average
           else
               effective_SINR = -betas.*log(mean(exp(-SINR_vector./betas)));
           end
           effective_SINR(effective_SINR == Inf) = 10^10; % if the SINR is too good the averager does not work, as the numerical accuracy is too low 
           
           effective_SINR_lin = effective_SINR;
           effective_SINR_dB  = 10.*log10(effective_SINR);
       end
   end
end 
