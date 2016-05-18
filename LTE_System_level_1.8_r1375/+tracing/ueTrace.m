classdef ueTrace < handle
% This class stores, for each UE the traces that we wanto to store. eg, CQI,
% throughput, etc, etc.
% (c) Josep Colom Ikuno, INTHFT, 2008

   properties
       latency_time_scale % This is used to compute the averate throughput using an exponential filter
       av_temp 
       TTI_length_s       % Length of one TTI
       nCodewords         % Number of codewords used
       RI                 % Rank indicator used
       CQI_sent           % CQI thta was fed back
       attached_site      % Index of the attached site
       attached_eNodeB    % Index of the attached eNodeB
       position           % UE position
       assigned_RBs       % Number of assigned RBs
       assigned_power_W   % Assigned power during this TTI in Watts
       rx_power_tb        % Mean power per-RB of desired signal
       rx_power_interferers % Mean power per-RB of interferers % (1,:) ... mean per-RB power, calculated over all RBs in a TTI, (:,2) ... ID of the interfering eNodeB
       ACK                % ACK/NACK
       TB_CQI             % TB CQI
       TB_size            % TB size
       N_used_bits        % Number of used bits
       BLER               % The used BLER for this TB
       avg_throughput     % Average throughput calculated with an exponential averaging window
       TB_SINR_dB         % The AWGN-equivalent TB SINR as measured by the link quality model. Stored in dB
       wideband_SINR      % Geometry
       
       UE_was_disabled    % Whether this UE was being simulated (enabled) or not during the given TTIs
       
       average_throughput_Mbps                % Average UE throughput in Mb/s
       average_spectral_efficiency_bit_per_cu % The spectral efficiency calculated in bits/channel use
       average_RBs_per_TTI                    % How many RBs did the UE get per TTI in average
       average_energy_per_bit                 % Energy needed per bit in J/bit
       
       % Optional traces (can be turned on/off in the config file)
       trace_SINR = true;
       SINR            % Stores the RS SINRs [dB]
       SNR             % Stores the average prequalization SNR
       
       parent_results_object % The parent results object
       
       % Whether less stuff should be stored in the trace (mainly
       % feedback). Helpful for the cases where you have a lot of UEs
       reduced_feedback_logs
   end

   methods
       function obj = ueTrace(...
               simulation_length_TTI,...
               n_RB,...
               nTX,...
               nRX,...
               maxStreams,...
               unquantized_CQI_feedback,...
               trace_SINR,...
               latency_time_scale,...
               TTI_length_s,...
               reduced_feedback_logs)
           
           obj.reduced_feedback_logs = reduced_feedback_logs;
           
           obj.latency_time_scale   = latency_time_scale;
           obj.TTI_length_s         = TTI_length_s;
           
           obj.attached_site        = zeros(1,simulation_length_TTI,'uint16');
           obj.attached_eNodeB      = zeros(1,simulation_length_TTI,'uint16');
           obj.position             = NaN(2,simulation_length_TTI);
           obj.assigned_RBs         = zeros(1,simulation_length_TTI,'uint8');
           obj.ACK                  = false(maxStreams,simulation_length_TTI);
           obj.TB_CQI               = NaN(maxStreams,simulation_length_TTI,'single');
           obj.TB_size              = zeros(maxStreams,simulation_length_TTI,'uint32'); % Max is 80*6*200 (20 MHz, 1 ms for 1 user)
           obj.N_used_bits          = zeros(maxStreams,simulation_length_TTI,'uint32'); % Max is 80*6*200 (20 MHz, 1 ms for 1 user)
           obj.BLER                 = zeros(maxStreams,simulation_length_TTI);
           obj.TB_SINR_dB           = zeros(maxStreams,simulation_length_TTI);
           obj.avg_throughput       = zeros(maxStreams,simulation_length_TTI);
           obj.wideband_SINR        = zeros(1,simulation_length_TTI);
           obj.UE_was_disabled      = false(1,simulation_length_TTI);
           obj.SNR                  = nan(1,simulation_length_TTI);
           obj.rx_power_tb          = zeros(1,simulation_length_TTI);
           obj.rx_power_interferers = cell(1,simulation_length_TTI);
           obj.assigned_power_W     = zeros(1,simulation_length_TTI,'single');
           
           % Feedback which is not stored if the "reduced feedback" flag is on
           if ~reduced_feedback_logs
               obj.nCodewords       = zeros(1,simulation_length_TTI,'uint8');
               obj.RI               = zeros(1,simulation_length_TTI,'uint8');
               
               % Trace of sent CQIs is initialized to -1 (an invalid value).
               if unquantized_CQI_feedback
                   obj.CQI_sent = zeros(maxStreams,n_RB,simulation_length_TTI,'single')-1;
               else
                   obj.CQI_sent = zeros(maxStreams,n_RB,simulation_length_TTI,'int8')-1;
               end
           end
           
           obj.trace_SINR = trace_SINR;
           if trace_SINR
               max_layers = min([nTX nRX]);
               obj.SINR = nan(max_layers,2*n_RB,max_layers,simulation_length_TTI);
           end
       end
       % Trace this specific TTI
       function store(obj,...
               feedback,...        % UE feedback
               nLayers,...         % Number of layers employed in the TB
               nCodewords,...      % Number of codewords employed in the TB
               attached_eNodeB,... % attached eNodeB object
               position,...        % UE position
               tti_idx,...         % TTI index
               assigned_RBs,...    % numbero of assigned RBs
               assigned_power,...  % assigned power
               TB_CQI,...          % CQI assigned to the UE
               TB_size,...         % TB size
               BLER,...            % TB BLER
               TB_SINR_dB,...      % TB post-equalization SINR
               N_used_bits,...     % number of used bits
               wideband_SINR,...   % Ue wideband SINR
               UE_disabled,...     % whether the UE is disabled
               SNR,...             % 5 UE SNR, as calculated in the LL simulator
               rx_power_rb_in_current_tti,... % Mean power per resource block of desired signal
               rx_power_interferers_rb_in_current_tti,... % Mean power per resource blocks of most significant interferers (number might vary).
               extra_traces)       % extra stuff (mainly the subcarrier SINR)
           
           % Optional varargin variables to trace are:
           %  - extra_traces{1} -> SINR
           %  - extra_traces{2} -> SNR
           obj.position(:,tti_idx)                          = position;
           obj.wideband_SINR(tti_idx)                       = wideband_SINR;
           obj.UE_was_disabled(tti_idx)                     = UE_disabled;

           % TB info
           obj.RI(tti_idx)         = nLayers;
           obj.nCodewords(tti_idx) = nCodewords;

           % Feedback which is not stored if the "reduced feedback" flag is on
           if ~obj.reduced_feedback_logs
               feedback_fields = fieldnames(feedback);
               for f_=1:length(feedback_fields)
                   current_field = feedback_fields{f_};
                   switch current_field
                       case 'CQI'
                           fb_CQI = feedback.CQI;
                           obj.CQI_sent(1:size(fb_CQI,1),:,tti_idx) = fb_CQI;
                   end
               end
           end
           
           if isempty(attached_eNodeB)
               obj.attached_site(tti_idx)   = -1;
               obj.attached_eNodeB(tti_idx) = -1;
           else
               obj.attached_site(tti_idx)   = attached_eNodeB.parent_eNodeB.id;
               obj.attached_eNodeB(tti_idx) = attached_eNodeB.eNodeB_id;
           end
           
           throughput     = zeros(size(obj.avg_throughput,1),1);
           if ~isempty(assigned_RBs)
               fb_ACK = feedback.ACK;
               obj.assigned_RBs(tti_idx)                        = assigned_RBs;
               obj.assigned_power_W(tti_idx)                    = assigned_power;
               obj.ACK(1:length(fb_ACK),tti_idx)                = fb_ACK;
               obj.TB_CQI(1:length(TB_CQI),tti_idx)             = TB_CQI;
               obj.TB_size(1:length(TB_size),tti_idx)           = TB_size;
               obj.N_used_bits(1:length(N_used_bits),tti_idx)   = N_used_bits;
               obj.BLER(1:length(BLER),tti_idx)                 = BLER;
               if ~isempty(TB_SINR_dB)
                   obj.TB_SINR_dB(1:length(TB_SINR_dB),tti_idx) = TB_SINR_dB;
               end
               if ~isempty(rx_power_rb_in_current_tti)
                   obj.rx_power_tb(tti_idx)                     = rx_power_rb_in_current_tti; % Calculate mean RB power to save memory
               end
               if ~isempty(rx_power_interferers_rb_in_current_tti)
                   obj.rx_power_interferers{tti_idx}            = rx_power_interferers_rb_in_current_tti; % Cell of variable length, depending on number of regarded interferers
               end
               throughput(1:length(fb_ACK),1) = sum(TB_size(:).*fb_ACK(:) / obj.TTI_length_s);
           end

           obj.av_temp = min(obj.latency_time_scale,tti_idx);
           if tti_idx==1
               obj.avg_throughput(1:length(throughput),tti_idx) = throughput;
           else
               obj.avg_throughput(1:length(throughput),tti_idx) = (1-1/obj.av_temp)*obj.avg_throughput(1:length(throughput),tti_idx-1) + 1/obj.av_temp*throughput(:);
           end
           obj.SNR(tti_idx) = SNR; % In dB
           if obj.trace_SINR
               obj.SINR(:,:,:,tti_idx) = extra_traces{1}; % In dB
           end
       end
       
       % Calculates the spectral efficiency based on the data stored on the ACK, assigned_RBs, and TB_size variables
       function calculate_final_average_spectral_efficiency_bit_per_cu(obj)
           TTIs_to_ignore              = obj.parent_results_object.TTIs_to_ignore_when_calculating_aggregates;
           total_bits                  = sum(double(obj.ACK).*double(obj.TB_size),1);
           TTIs_to_account_for         = true(1,length(obj.assigned_RBs));
           TTIs_to_account_for(1:TTIs_to_ignore)    = false; % Ignore TTIs where no feedback information was available
           TTIs_to_account_for(obj.UE_was_disabled) = false;
           
           sum_total_bits              = sum(total_bits(TTIs_to_account_for));
           channel_uses_total          = sum(double(obj.assigned_RBs(TTIs_to_account_for)))*12*14;
           obj.average_spectral_efficiency_bit_per_cu = sum_total_bits ./ channel_uses_total;
       end
       
       % Calculates the average UE throughput based on the TB_size, ACK, and assigned_RBs variables
       function calculate_final_average_throughput_Mbps(obj)
           TTIs_to_ignore                           = obj.parent_results_object.TTIs_to_ignore_when_calculating_aggregates;
           total_bits                               = sum(double(obj.ACK).*double(obj.TB_size),1);
           TTIs_to_account_for                      = true(1,length(obj.assigned_RBs));
           TTIs_to_account_for(1:TTIs_to_ignore)    = false; % Ignore TTIs where no feedback information was available
           TTIs_to_account_for(obj.UE_was_disabled) = false;

           sum_total_bits              = sum(total_bits(TTIs_to_account_for));
           accounted_TTIs              = sum(TTIs_to_account_for);
           obj.average_throughput_Mbps = sum_total_bits / (accounted_TTIs*obj.TTI_length_s) / 1e6;
       end
       
       % Calculates the average number of RBs used per TTI
       function calculate_final_average_RBs_per_TTI(obj)
           TTIs_to_ignore                           = obj.parent_results_object.TTIs_to_ignore_when_calculating_aggregates;
           TTIs_to_account_for                      = true(1,length(obj.assigned_RBs));
           TTIs_to_account_for(1:TTIs_to_ignore)    = false; % Ignore TTIs where no feedback information was available
           TTIs_to_account_for(obj.UE_was_disabled) = false;

           obj.average_RBs_per_TTI     = mean(obj.assigned_RBs(TTIs_to_account_for));
       end
       
       % Calculate the average energy per bit employed for this UE
       function calculate_final_average_energy_per_bit(obj)
           TTIs_to_ignore                           = obj.parent_results_object.TTIs_to_ignore_when_calculating_aggregates;
           TTIs_to_account_for                      = true(1,length(obj.assigned_RBs));
           TTIs_to_account_for(1:TTIs_to_ignore)    = false; % Ignore TTIs where no feedback information was available
           TTIs_to_account_for(obj.UE_was_disabled) = false;

           total_bits             = sum(double(obj.ACK).*double(obj.TB_size),1);
           sum_total_bits         = sum(total_bits(TTIs_to_account_for));
           sum_total_energy       = sum(obj.assigned_power_W(TTIs_to_account_for))*obj.TTI_length_s;
           
           obj.average_energy_per_bit = sum_total_energy / sum_total_bits;
       end
   end
end 
