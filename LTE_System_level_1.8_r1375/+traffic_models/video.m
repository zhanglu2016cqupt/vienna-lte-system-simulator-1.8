classdef video < traffic_models.generic_tm
% This class is used for video streaming simulations
% Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2010 by INTHFT 
properties
    type = 'video';
    slice_mean = 100;    % mean slice size
    slice_max = 250;    % max slice size in bytes
    slice_min = 53;    % min slice size
    inter_mean = 6;     % mean slice interarrival time (encoder delay)
    inter_max = 12.5;   % max slice interarrival time
    inter_min = 2.5;   % min slice interarrival time
    slice_x
    inter_x
    slices = zeros(8,1); % eight slices per frame
    inters = zeros(7,1);
    slice_cmf
    inter_cmf
    counter = 1;
    slice_counter;
    iat;
    rate_constraint = 64; % 64 kbit/s minimum rate constraint
    delay_constraint = 100; % in ms, mean arrival time for all slices in a frame is 48ms, interarrival between frames is 100ms
    iit_offset
end

methods
    function obj = video(UE,HARQ_delay)
        obj = obj@traffic_models.generic_tm(UE,HARQ_delay,'video');
        obj.is_fullbuffer = false;
        
        % Truncated Pareto distribution for the slice size
        a = 1.2;
        k = obj.slice_min;
        m = obj.slice_max;
        obj.slice_x = linspace(obj.slice_min,obj.slice_max,1000);
        pmf = a*k^a./obj.slice_x.^(a+1)*(obj.slice_x(end)-obj.slice_x(end-1));
%         pmf(end) = pmf(end)+(k/m)^a;
        pmf = pmf/sum(pmf);
        obj.slice_cmf = cumsum(pmf);
        obj.slice_x = obj.slice_x*8;
        obj.iit_offset = randi(100)-1;
%         obj.generate_slices;
        
        % Truncated Pareto distribution for the interarrival time
        a = 1.2;
        k = obj.inter_min;
        m = obj.inter_max;
        obj.inter_x = linspace(obj.inter_min,obj.inter_max,1000);
        pmf = a*k^a./obj.inter_x.^(a+1)*(obj.inter_x(end)-obj.inter_x(end-1));
%         pmf(end) = pmf(end)+(k/m)^a;
        pmf = pmf/sum(pmf);
        obj.inter_cmf = cumsum(pmf);
%         obj.generate_inters;
        
    end
    
    function check_TTI(obj)
        obj.counter = obj.counter + 1;
        if obj.UE.clock.current_TTI > obj.iit_offset
            if ~mod(obj.UE.clock.current_TTI-1-obj.iit_offset,100)
                obj.generate_packet(obj.slice_cmf,obj.type,obj.slice_x);
                obj.iat = obj.eval_cmf(obj.inter_cmf,obj.inter_x);
                obj.counter = 0;
                obj.slice_counter = 1;
            end
            if obj.iat <= obj.counter && obj.slice_counter < 8
                obj.generate_packet(obj.slice_cmf,obj.type,obj.slice_x);
                obj.iat = obj.eval_cmf(obj.inter_cmf,obj.inter_x);
                obj.counter = 0;
                obj.slice_counter = obj.slice_counter + 1;
            end
        end
        obj.bit_count = obj.get_buffer_length;
    end
    
   function rate = get_arrival_rate(obj)
        rate = obj.slice_mean*8*8/0.1; % first 8: byte->bit second 8: 8 slices per frame
    end
    
    function packet_parts = decrease_packets(obj,N_data_bits)
        packet_parts = [];
        for bb = 1:length(obj.packet_buffer)
            read_temp = obj.read_start+bb-1;
            if read_temp > length(obj.packet_buffer)
                read_temp = mod(read_temp,length(obj.packet_buffer));
            end
            if obj.packet_buffer(read_temp).id 
                [N_data_bits,part] = obj.packet_buffer(read_temp).send_data(N_data_bits);
                packet_parts = [packet_parts,part];
                if N_data_bits <= 0
                    break;
                end
            end
        end
        obj.bit_count = obj.get_buffer_length;
    end
         
end

end