clc;
clear all
close all
sectors{1}= 1;
sectors{2}= 2;
sectors{3}= 2;
sectors{4}= [1,3];
sectors{5}= [1,2,3];
sectors{6}= 2;
sectors{7} = 3;

ccc = 0;
for i1 = 1:52
    i1
    load(['uniform (' num2str(i1) ').mat']);
    for b_i = 1:5
        for s_i = 1:length(sectors{i1})
            UE_ids = [eNodeBs(b_i).sectors(s_i).UEs_attached_last_TTI.id];
            for u_i = 1:length(UE_ids)
                ue_id = UE_ids(u_i);
                ccc = ccc+1;
                TP(ccc) = sum(simulation_traces.UE_traces(ue_id).TB_size(1,:).*uint32(simulation_traces.UE_traces(ue_id).ACK(1,:)))/size(simulation_traces.UE_traces(ue_id).TB_size,2)*10^-3;
            end
        end
    end
end
[F1,X1] = ecdf(TP);

% for i1 = 1
%     i1
%     load(['poisson (' num2str(i1) ').mat']);
%     TP(i1) = sum(simulation_traces.UE_traces(i1).TB_size(1,:).*uint32(simulation_traces.UE_traces(i1).ACK(1,:)))/size(simulation_traces.UE_traces(i1).TB_size,2)*10^-3;
% %     clear -except TP i1 F1 X1
% end
% 
% [F2,X2] = ecdf(TP);