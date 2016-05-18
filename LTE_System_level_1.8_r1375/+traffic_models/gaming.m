classdef gaming < traffic_models.generic_tm
% This class is used for gaming simulations
% Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2010 by INTHFT 
properties
    type = 'gaming';
    UDP_header = 2; % header overhead
    next_packet; 
    initial_packet_a = 0;
    initial_packet_b = 40;
    packet_size_a = 120;
    packet_size_b = 36;
    packet_time_a = 55;
    packet_time_b = 6;
    delay_constraint = 60;
    arrival_rate = (120+36*-psi(1))*8/(55+6*-psi(1));
%     state = true;
end

methods
    function obj = gaming(UE,HARQ_delay)
        obj = obj@traffic_models.generic_tm(UE,HARQ_delay,'gaming');
        obj.next_packet = randi([obj.initial_packet_a,obj.initial_packet_b],1)+1; % intial packet of the game    
        obj.is_fullbuffer = false;
    end
    
    function check_TTI(obj)
        obj.next_packet = obj.next_packet-1;
        if obj.next_packet == 0
            coin_toss = rand;
            obj.generate_packet(round((obj.packet_size_a-obj.packet_size_b*log(-log(coin_toss))+obj.UDP_header)*8),obj.type);
            coin_toss = rand;
            obj.next_packet = round(obj.packet_time_a-obj.packet_time_b*log(-log(coin_toss)));
        end
        obj.bit_count = obj.get_buffer_length;
    end
    
    function rate = get_arrival_rate(obj)
        rate = obj.arrival_rate;
%         if obj.state 
%             rate = obj.arrival_rate_voice;
%         else
%             rate = obj.arrival_rate_silence;
%         end
    end
    
    function packet_parts = decrease_packets(obj,N_data_bits)
        packet_parts = [];
        for bb = 1:length(obj.packet_buffer)
            read_temp = obj.read_start+bb-1;
            if read_temp > length(obj.packet_buffer)
                read_temp = mod(read_temp,length(obj.packet_buffer));
            end
            [N_data_bits,part] = obj.packet_buffer(read_temp).send_data(N_data_bits);
            packet_parts = [packet_parts,part];
            if N_data_bits <= 0
                break;
            end
        end
        obj.bit_count = obj.get_buffer_length;
    end
    
end
end