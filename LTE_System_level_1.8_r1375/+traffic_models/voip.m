classdef voip < traffic_models.generic_tm
% This class is used for voip simulations according to R1-070674
% Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2010 by INTHFT 

properties
    type = 'voip';
    state = true;
    c = 0.01;  % data according to RAN R1-070674
    d = 0.99;
    iit_voice = 20;
    iit_silence = 160;
    delay_constraint = 50;
    voice_size = 40*8;
    silence_size = 15*8;
    arrival_rate_voice
    arrival_rate_silence
    iit_offset
end

methods
    function obj = voip(UE,HARQ_delay)
        obj = obj@traffic_models.generic_tm(UE,HARQ_delay,'voip');
        obj.arrival_rate_voice = obj.voice_size/obj.iit_voice;
        obj.arrival_rate_silence = obj.silence_size/obj.iit_silence; 
        obj.is_fullbuffer = false;
        obj.iit_offset = randi(20)-1;
    end
    
    function check_TTI(obj)  
        if ~mod(obj.UE.clock.current_TTI-1-obj.iit_offset,20)    % here comes the two state markov model assumed for the speech process (state changes every 20 ms)
            coin_toss = rand;
            if obj.state    % active state
               if coin_toss < obj.d
                   obj.state = true;
               else
                   obj.state = false;
               end
            else        % inactive state
                if coin_toss < obj.c
                    obj.state = true;
                else
                    obj.state = false;
                end
            end
            if obj.state % generate a new speech packet in the active state 
                  obj.generate_packet(obj.voice_size,obj.type);
            else        % in the inactive state
                if ~mod(obj.UE.clock.current_TTI-1-obj.iit_offset,160) % generate a new silence descriptor (SID) packet if the deadline is zero
                      obj.generate_packet(obj.silence_size,obj.type);
                end
            end
            obj.bit_count = obj.get_buffer_length;
        end
    end
        
    function rate = get_arrival_rate(obj)
        if obj.state 
            rate = obj.arrival_rate_voice;
        else
            rate = obj.arrival_rate_silence;
        end
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