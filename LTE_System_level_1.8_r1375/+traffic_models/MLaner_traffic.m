classdef MLaner_traffic < traffic_models.generic_tm
% This class is used to generate traffic according to Markus Laners
% measurements
% Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2011 by INTHFT 
properties
    type = 'MLaner';
    mean_log_rate = 1.3525;
    std_log_rate = 0.1954;
    mean_rate
    sim_time
    number_IPpackets
    generation_times
    counter = 0;
    IP_packet_size = 12000;
    last_packet_size
end

methods
    function obj = MLaner_traffic(UE,HARQ_delay,total_time)
        obj = obj@traffic_models.generic_tm(UE,HARQ_delay,'MLaner');
        obj.is_fullbuffer = false;
        
        obj.sim_time = total_time;
        obj.mean_rate = min(10^7,10^(lognrnd(obj.mean_log_rate,obj.std_log_rate)));
        obj.number_IPpackets = ceil(obj.mean_rate/obj.IP_packet_size*obj.sim_time*10^-3);
        obj.last_packet_size = round(obj.number_IPpackets*obj.IP_packet_size-obj.mean_rate*obj.sim_time*10^-3);
        obj.generation_times = sort(randi(obj.sim_time,obj.number_IPpackets,1),'ascend');
    end
    
    function check_TTI(obj)
        obj.counter = obj.counter+1;
        packet_num = sum(obj.counter == obj.generation_times);
        if packet_num > 0
            if obj.counter == obj.generation_times(end)
                obj.generate_packet(obj.last_packet_size,obj.type);
            else
                obj.generate_packet(obj.IP_packet_size*packet_num,obj.type);
            end
        end
        obj.bit_count = obj.get_buffer_length;
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