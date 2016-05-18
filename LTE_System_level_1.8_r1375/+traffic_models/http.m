classdef http < traffic_models.generic_tm
% This class is used for http simulations
% Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2010 by INTHFT 
properties
    type = 'http';
    main_mean_file = 10710; % mean file size before truncation
    main_max_file = 2*10^6;  % maximum file size (truncation value)
    main_sigma_file = 25032; % file size standard deviation before truncation
    main_min_file = 100;    
        
    emb_mean_file = 7758; % mean file size before truncation
    emb_max_file = 2*10^6;  % maximum file size (truncation value)
    emb_sigma_file = 126168; % file size standard deviation before truncation
    emb_min_file = 50;
    
    numb_mean = 5.64;
    numb_max = 53;
    
    mean_reading_time = 30; % reading time in sec.
    
    main_cmf
    emb_cmf
    numb_cmf
    reading_cmf
    main_x
    emb_x
    numb_x
    reading_x
    numb = Inf;
    waiting_time = 0;
    state = 0;
end

methods
    function obj = http(UE,HARQ_delay)
        obj = obj@traffic_models.generic_tm(UE,HARQ_delay,'http');
        obj.is_fullbuffer = false;
        
        % Log normally distributed main object size
        obj.main_x = linspace(obj.main_min_file,obj.main_max_file,500000); 
        sigma1 = 1.37^2;
        mu1 = 8.37;
%         sigma1 = log(1+obj.main_sigma_file/obj.main_mean_file^2);
%         mu1 = log(obj.main_mean_file)-0.5*log(1+obj.main_sigma_file/obj.main_mean_file^2);
        pmf = 1./sqrt(2*pi*sigma1*obj.main_x.^2).*exp(-(log(obj.main_x)-mu1).^2/(2*sigma1))*(obj.main_x(end)-obj.main_x(end-1)); % according to RAN R1-070674 (mean = 2 Mbyte)
        pmf = pmf/sum(pmf);  % renormalization due to truncation to max_file
        obj.main_cmf = cumsum(pmf);
        obj.main_x = obj.main_x*8;
        
        % Log normally distributed embedded object size
        obj.emb_x = linspace(obj.emb_min_file,obj.emb_max_file,500000); 
        sigma1 = 2.36^2;
        mu1 = 6.17;
%         sigma1 = log(1+obj.emb_sigma_file/obj.emb_mean_file^2);
%         mu1 = log(obj.emb_mean_file)-0.5*log(1+obj.emb_sigma_file/obj.emb_mean_file^2);
        pmf = 1./sqrt(2*pi*sigma1*obj.emb_x.^2).*exp(-(log(obj.emb_x)-mu1).^2/(2*sigma1))*(obj.emb_x(end)-obj.emb_x(end-1)); % according to RAN R1-070674 (mean = 2 Mbyte)
        pmf = pmf/sum(pmf);  % renormalization due to truncation to max_file
        obj.emb_cmf = cumsum(pmf);
        obj.emb_x = obj.emb_x*8;
        
        % Truncated Pareto distribution for the number of embedded objects
        a = 1.1;
        k = 2;
        m = obj.numb_max+2;
        obj.numb_x = linspace(k,m,1000);
        pmf = a*k^a./obj.numb_x.^(a+1)*(obj.numb_x(end)-obj.numb_x(end-1));
%         pmf(end) = pmf(end)+(k/m)^a;
        pmf = pmf/sum(pmf);
        obj.numb_cmf = cumsum(pmf);
        obj.numb_x = obj.numb_x-k;
        
        % Exponentially distributed reading time
        lambda = 1/obj.mean_reading_time;
        percentile9995 = 1/lambda + 1/lambda*10;
        obj.reading_x = linspace(percentile9995/10^4,percentile9995,10^4);
        pmf = lambda*exp(-lambda*obj.reading_x)*(obj.reading_x(end)-obj.reading_x(end-1));
        pmf = pmf/sum(pmf);  % renormalization due to truncation to max_file
        obj.reading_cmf = cumsum(pmf);

    end
      
    function check_TTI(obj)
        switch obj.state
            case 0
                obj.generate_packet(obj.main_cmf,obj.type,obj.main_x);
                obj.state = 1;
            case 1
                if obj.bit_count <= 0 
                    if isinf(obj.numb)
                        obj.numb = obj.eval_cmf(obj.numb_cmf,obj.numb_x);
                    end
                    if  obj.numb > 0
                        obj.generate_packet(obj.emb_cmf,obj.type,obj.emb_x);
                        obj.numb = obj.numb-1;
                    else
                        obj.numb = Inf;
                        obj.state = 2;
                        obj.waiting_time = obj.eval_cmf(obj.reading_cmf,obj.reading_x);
                        obj.waiting_time = round(obj.waiting_time/10^-3); % transform from sec. to TTIs
                    end
                end
            case 2 
                obj.waiting_time = obj.waiting_time-1;
                if obj.waiting_time <= 0  % if the end of the reading state is reached
                    obj.state = 0;
                end
        end
        obj.bit_count = obj.get_buffer_length;
    end
    
    function decrease_packets(obj,N_data_bits)
        diff = obj.packet_buffer.size-N_data_bits;
        obj.packet_buffer.size = max(0,diff);
        if diff <= 0
            obj.packet_buffer = traffic_models.data_packet(0,0,0);
        end
        obj.bit_count = obj.get_buffer_length;
    end
      
end

end