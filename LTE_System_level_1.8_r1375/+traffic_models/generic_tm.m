classdef generic_tm < handle
% This is a generic class for all traffic models
% Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2010 by INTHFT

properties
    UE
    id_counter = 1;
    packet_buffer
    delay_buffer = [];
    HARQ_delay
    bit_count = 0; 
    read_start = 1;
    
    % Added for easy and faster identification in the link performance model without having to call strcmp
    is_fullbuffer = false; % Set to false by default
end

methods
    function obj = generic_tm(UE,HARQ_delay,type)
        obj.UE = UE;
        obj.HARQ_delay = HARQ_delay;
        switch type
            case {'ftp','http','car'}
                temp = traffic_models.data_packet(0,0,0);
            otherwise
                temp(1:100) = traffic_models.data_packet(0,0,0);
        end
        obj.packet_buffer = temp;
    end
       
    function generate_packet(obj,dat_size,type,varargin)
        switch type
            case {'voip','gaming','MLaner'}
                packet = traffic_models.data_packet(dat_size,obj.UE.clock.current_TTI,obj.id_counter);
            case {'ftp','http','video'}
                packet_size = obj.eval_cmf(dat_size,varargin{1});
                packet = traffic_models.data_packet(packet_size,obj.UE.clock.current_TTI,obj.id_counter);
            case {'car'}
                packet = traffic_models.data_packet(dat_size,obj.UE.clock.current_TTI,obj.id_counter);                
        end
%         obj.packet_buffer = [obj.packet_buffer,packet];
%         obj.id_counter = obj.id_counter +1;
        obj.packet_buffer(obj.id_counter) = packet;
        obj.id_counter = mod(obj.id_counter,length(obj.packet_buffer)) + 1;
    end
    
    function check_TTI(obj)
    end
    
    function size = get_buffer_length(obj)
        if strcmp(obj.type,'fullbuffer')
            size = Inf;
        else
            size = 0;
            for pp = 1:length(obj.packet_buffer)
                size = size + sum(obj.packet_buffer(pp).get_size);
            end
        end
    end
    
    function remove_packet(obj,packet_id,success)
        temp_packet = obj.packet_buffer(obj.get_packet_ids == packet_id);
        if obj.read_start == packet_id % oldest packet removed, start reading the buffer at the next one
            if obj.get_buffer_length
                while(1)
                    obj.read_start = obj.read_start+1;
                    if obj.read_start > length(obj.packet_buffer) 
                        obj.read_start = 1;
                        break;
                    end
                    if obj.packet_buffer(obj.read_start).id ~= 0 % if next packet is not acknowledged start reading here, else increase read_start
                        break;
                    end
                end
            else
                obj.read_start = obj.id_counter;
            end
        end
%         obj.packet_buffer([obj.packet_buffer.get_id] == packet_id) = [];
        obj.packet_buffer(obj.get_packet_ids == packet_id) = traffic_models.data_packet(0,0,0);
        if success
            obj.delay_buffer = [obj.delay_buffer,obj.UE.clock.current_TTI-temp_packet.get_origin-obj.HARQ_delay];
        else 
            obj.delay_buffer = [obj.delay_buffer,Inf];
        end
%         delete(temp_packet);
%         obj.id_counter = obj.id_counter - 1;
    end
    
    function X = eval_cmf(obj,cmf,x)
        coin_toss = rand;
        larger = cmf/coin_toss>ones(size(cmf));
        [C,I] = max(larger);
        toss = rand;
        if I < 2
            I = 2;
        end
        if toss > 0.5   % instead of just rounding
            X = ceil(x(I-1)); 
        else
            X = floor(x(I-1));
        end
    end
       
    function bits_left = coarse_decrease(obj,N_bits,varargin)
        if ~isempty(varargin)
            obj.bit_count = obj.get_buffer_length;
        end
        obj.bit_count = max(0,obj.bit_count-N_bits);
        bits_left = obj.bit_count;
    end
    
    function ids = get_packet_ids(obj)
        ids = zeros(1,length(obj.packet_buffer));
        for i = 1:length(obj.packet_buffer)
            ids(i) = obj.packet_buffer(i).get_id;
        end
    end
end

end