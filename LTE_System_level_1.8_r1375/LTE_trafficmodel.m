function model = LTE_trafficmodel(trafficmodel_params,UE,HARQ_delay,varargin)
% This function chooses the actual traffic model for the user
% Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2010 by INTHFT

aPrioriPdf = [0.1,0.2,0.2,0.3,0.2]; % a priori pdf according to which a traffic model is picked (defined in RAN R1-070674); just use single digit numbers for this, otherwise set will not work
% [ftp,http,video,voip,gaming]
% aPrioriPdf = [1,0,0,0,0,0];
if trafficmodel_params.usetraffic_model % if the traffic models shall be used
    if isempty(varargin)
        set = [];
        for i = 1:length(aPrioriPdf)
            set = [set,i*ones(1,aPrioriPdf(i)*10)];
        end
        index = randi([1,length(set)]);
    else
        index = 1;
        set = varargin{1};
    end
    switch set(index)   % randomly pack one of the traffic models according to the aPrioriPdf
        case 1
            model = traffic_models.ftp(UE,HARQ_delay);
        case 2
            model = traffic_models.http(UE,HARQ_delay);
        case 3
            model = traffic_models.video(UE,HARQ_delay);
        case 4
            model = traffic_models.voip(UE,HARQ_delay);
        case 5
            model = traffic_models.gaming(UE,HARQ_delay);
        case 6
            model = traffic_models.fullbuffer(UE,HARQ_delay);
        case 7
            if varargin{2}
                model = traffic_models.MLaner_traffic(UE,HARQ_delay,varargin{3});     
            else
                model = traffic_models.fullbuffer(UE,HARQ_delay);
            end
        case 8
            model = traffic_models.car(UE,HARQ_delay);
    end
else
    model = traffic_models.fullbuffer(UE,HARQ_delay);
end
end
