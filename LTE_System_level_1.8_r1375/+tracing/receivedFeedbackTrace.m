classdef receivedFeedbackTrace < handle
% This class stores all of the feedbacks received by eNodeBs in the
% network. It is done like this in order to be able to preallocate the space
% (you don't a priori know where each user will be).
% (c) Josep Colom Ikuno, INTHFT, 2008

   properties
       % To which UE this entry belongs
       UE_id
       
       % To which eNodeB-sector this belongs
       eNodeB_id
       sector_id
       
       % When was this received
       tti_idx
       
       % Idx that is now available to write
       current_idx
       
       % Info transmitted by the UE
       CQI          % CQI reports
       RI           % RI report (for SM modes)
       ACK          % TB ACK/NACK reports
       TB_size      % TB size
       UE_scheduled % Whether the UE was scheduled
       
       parent_results_object % The parent results object
   end

   methods
       function obj = receivedFeedbackTrace(simulation_length_TTI,numUEs,n_RB,maxStreams,unquantized_CQI_feedback)
           % Trace of sent CQIs is initialized to -1 (an invalid value)
           obj.UE_id        = zeros(1,simulation_length_TTI*numUEs,'int32')-1;
           obj.eNodeB_id    = zeros(1,simulation_length_TTI*numUEs,'uint16');
           obj.sector_id    = zeros(1,simulation_length_TTI*numUEs,'uint8');
           obj.tti_idx      = zeros(1,simulation_length_TTI*numUEs)-1;
           
           % Info transmitted by the UE
           if unquantized_CQI_feedback
               obj.CQI      = zeros(n_RB,maxStreams,simulation_length_TTI*numUEs,'single')-1;
           else
               obj.CQI      = zeros(n_RB,maxStreams,simulation_length_TTI*numUEs,'int8')-1;
           end
           
           obj.ACK          = false(maxStreams,simulation_length_TTI*numUEs); % UE ACKs
           obj.TB_size      = zeros(maxStreams,simulation_length_TTI*numUEs); % TB size
           obj.RI           = zeros(1,simulation_length_TTI*numUEs,'uint8');  % Rank indicator (for SM modes)
           obj.UE_scheduled = false(1,simulation_length_TTI*numUEs);          % Whether the UE was scheduled in this TTI or not
           
           obj.current_idx  = 1;
       end
       function store(obj,feedback,UE_id,eNodeB_id,sector_id,tti_idx)
           % Store the values  
           the_current_idx                = obj.current_idx;
           obj.UE_id(the_current_idx)     = UE_id;
           obj.eNodeB_id(the_current_idx) = eNodeB_id;
           obj.sector_id(the_current_idx) = sector_id;
           obj.tti_idx(the_current_idx)   = tti_idx;
           
           obj.UE_scheduled(1,the_current_idx)                     = feedback.UE_scheduled;
           obj.TB_size(1:length(feedback.TB_size),the_current_idx) = feedback.TB_size;
           obj.ACK(1:length(feedback.ACK),the_current_idx)         = feedback.ACK;
           
           % Store UE feedback (if it exists)
           feedback_fields = fieldnames(feedback.feedback);
           for f_=1:length(feedback_fields)
               current_field = feedback_fields{f_};
               switch current_field
                   case 'CQI'
                       obj.CQI(:,1:size(feedback.feedback.CQI,1),the_current_idx) = feedback.feedback.CQI';
                   case 'RI'
                       obj.RI(1,the_current_idx) = feedback.feedback.RI;
                   otherwise
                       % Do nothing
               end
           end
           
           % Advance the counter 1 position
           obj.current_idx = obj.current_idx + 1;
       end
   end
end 
