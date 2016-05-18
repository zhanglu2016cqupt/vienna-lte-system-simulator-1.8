classdef data_packet < handle
% this class represents a data packet
% author: Stefan Schwarz
% contact: sschwarz@nt.tuwien.ac.at
% 2010

properties
    id
    size
    origin_TTI
    part_id = 1;
    packet_parts
    stop_packet = false;
end
    
methods
    function obj = data_packet(size,origin,id)
        obj.size = size;
        obj.origin_TTI = origin;
        obj.id = id;
        temp(1:100) = traffic_models.packet_part(obj.id,0,0);
        obj.packet_parts = temp;
    end
    
    function size_back = get_size(obj)
        size_back = obj.size;
    end
    
    function [TTI_back] = get_origin(obj)
        TTI_back = obj.origin_TTI;
    end
    
    function [back_size,part] = send_data(obj,send_size)
        if ~obj.stop_packet
            part = traffic_models.packet_part(obj.id,min(obj.size,send_size),obj.part_id);
%             obj.packet_parts = [obj.packet_parts,part];
            obj.packet_parts(obj.part_id) = part;
            back_size = max(0,send_size-obj.size);
            obj.size = max(obj.size-send_size,0);
%             obj.part_id = obj.part_id + 1;
            obj.part_id = mod(obj.part_id,100) + 1;
%             obj.size
        else
            obj.size = 0;
            back_size = send_size;
            part = [];
        end
    end
    
    function restore_packet_part(obj,part_id)
        ind = find([obj.packet_parts.id] == part_id);
        if ~obj.stop_packet
            obj.size = obj.size+obj.packet_parts(ind).part_size;
        end
        obj.packet_parts(ind) = traffic_models.packet_part(obj.id,0,0);
%         delete(obj.packet_parts(ind));
%         obj.part_id = obj.part_id - 1;
    end
    
    function [packet_done,packet_id] = acknowledge_packet_part(obj,part_id,state)
%         delete(obj.packet_parts([obj.packet_parts.id] == part_id));
        obj.packet_parts([obj.packet_parts.id] == part_id) = traffic_models.packet_part(obj.id,0,0);
        if ~state
            obj.stop_packet = true;
        end
        if sum(obj.packet_parts.get_size) == 0 && obj.size == 0
            packet_done = true;
        else
            packet_done = false;
        end
%         if (obj.size == 0 && isempty(obj.packet_parts)) || (isempty(obj.packet_parts) && ~state)
%             packet_done = true;
%         else
%             packet_done = false;
%         end
        packet_id = obj.id;
%         obj.part_id = obj.part_id - 1;
    end
    
    function back_id = get_id(obj)
        back_id = obj.id;
    end
end

end
        