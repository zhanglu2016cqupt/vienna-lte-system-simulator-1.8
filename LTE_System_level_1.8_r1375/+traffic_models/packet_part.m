classdef packet_part < handle
% this class is a part of a full data_packet that is already transmitted
% but needs to be acknowledged
% author: Stefan Schwarz

properties
    data_packet_id
    part_size
    id
end

methods
    function obj = packet_part(data_packet_id,size,part_id)
        obj.part_size = size; 
        obj.data_packet_id = data_packet_id;
        obj.id = part_id;
    end
    
    function restore_packet_part(obj)
%         if isvalid(obj.data_packet)
%             obj.data_packet.restore_packet_part(obj.id);
%         end
    end
    
    function id = get_id(obj)
        id = obj.id;
    end
    
    function [packet_done,packet_id] = acknowledge_packet_part(obj,state)
        [packet_done,packet_id] = obj.data_packet.acknowledge_packet_part(obj.id,state);
    end
    
    function back_size = get_size(obj)
        back_size = obj.part_size;
    end
end

end