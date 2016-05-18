classdef ftp < traffic_models.generic_tm
% This class is used for ftp simulations
% Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2010 by INTHFT 

properties
    type = 'ftp';
    mean_file = 2*10^6; % mean file size before truncation
    max_file = 5*10^6;  % maximum file size (truncation value)
    sigma_file = 0.722*10^6; % file size standard deviation before truncation
    mean_reading_time = 180; % reading time in sec.
    data_cmf
    reading_cmf
    data_x
    reading_x
    waiting_time = 1;
    state = false;
end

methods
    function obj = ftp(UE,HARQ_delay)
        obj = obj@traffic_models.generic_tm(UE,HARQ_delay,'ftp');
        obj.is_fullbuffer = false;
        
        % Log normally distributed file size
        obj.data_x = linspace(obj.max_file/10000,obj.max_file,10000); 
        sigma1 = 0.35^2;
        mu1 = 14.45;
        pmf = 1./sqrt(2*pi*sigma1*obj.data_x.^2).*exp(-(log(obj.data_x)-mu1).^2/(2*sigma1))*(obj.data_x(end)-obj.data_x(end-1)); % according to RAN R1-070674 (mean = 2 Mbyte)
        pmf = pmf/sum(pmf);  % renormalization due to truncation to max_file
        obj.data_cmf = cumsum(pmf);
        obj.data_x = obj.data_x * 8;
       
        % Exponentially distributed reading time
        lambda = 1/obj.mean_reading_time;
        percentile9995 = 1/lambda + 1/lambda*10;
        obj.reading_x = linspace(percentile9995/10^4,percentile9995,10^4);
        pmf = lambda*exp(-lambda*obj.reading_x)*(obj.reading_x(end)-obj.reading_x(end-1));
        pmf(end) = pmf(end)+1-sum(pmf);  % renormalization due to truncation to max_file
        obj.reading_cmf = cumsum(pmf);

    end
      
    function check_TTI(obj)  
        if ~obj.state
            obj.waiting_time = obj.waiting_time-1;
            if obj.waiting_time <= 0  % if the end of the reading state is reached
                obj.generate_packet(obj.data_cmf,obj.type,obj.data_x);
                obj.state = true;
                obj.bit_count = obj.get_buffer_length;
            end
        end
    end
    
    function decrease_packets(obj,N_data_bits)
        diff = sum([obj.packet_buffer.size])-N_data_bits;
        obj.packet_buffer.size = max(0,diff);
        if diff <= 0
            obj.waiting_time = obj.eval_cmf(obj.reading_cmf,obj.reading_x);
            obj.waiting_time = round(obj.waiting_time/10^-3); % transform from sec. to TTIs
            obj.state = false;  % enter the reading state
            obj.packet_buffer = traffic_models.data_packet(0,0,0);
        end
        obj.bit_count = obj.get_buffer_length;
    end
      
end

end