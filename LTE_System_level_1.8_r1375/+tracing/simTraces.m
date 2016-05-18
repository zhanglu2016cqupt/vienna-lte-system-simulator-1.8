classdef simTraces < handle
% This class stores in an ordered way the traces of all of the enodeBs and
% UEs.
% (c) Josep Colom Ikuno, INTHFT, 2008

   properties
       % Traces of what the UE sends and receives
       UE_traces
       % Traces of what the eNodeb sends
       eNodeB_tx_traces
       % Traces of the UEs feedback
       eNodeB_rx_feedback_traces
       % Traces from the schedulers
       scheduler_traces
       
       % Used to store extra UE-related info, whatever this could be
       extra_UE_info
       
       % When calculating aggregates, it will ignore the following number of TTIs
       TTIs_to_ignore_when_calculating_aggregates = 0;
   end
   
   methods
       function calculate_UE_aggregates(obj)
           for s_=1:length(obj.eNodeB_tx_traces)
               obj.eNodeB_tx_traces(s_).calculate_final_average_BLER;
           end
           for u_=1:length(obj.UE_traces)
               obj.UE_traces(u_).calculate_final_average_throughput_Mbps;
               obj.UE_traces(u_).calculate_final_average_spectral_efficiency_bit_per_cu;
               obj.UE_traces(u_).calculate_final_average_RBs_per_TTI;
               obj.UE_traces(u_).calculate_final_average_energy_per_bit;
           end
       end
   end
end 
