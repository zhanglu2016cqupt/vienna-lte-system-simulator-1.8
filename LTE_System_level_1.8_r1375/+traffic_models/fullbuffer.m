classdef fullbuffer < traffic_models.generic_tm
% This class is used for full buffer simulations
% Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2010 by INTHFT 
properties
    type = 'fullbuffer';
end

methods
    function obj = fullbuffer(UE,HARQ_delay)
        obj = obj@traffic_models.generic_tm(UE,HARQ_delay,'fullbuffer');
        obj.is_fullbuffer = true;
    end
    
    function decrease_packets(obj,N_bits)
    end
end

end