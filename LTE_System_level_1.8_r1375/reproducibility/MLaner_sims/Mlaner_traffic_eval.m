close all
N_sector = length(eNodeBs_sectors);
N_UE = length(UEs);
color = [1,0,0;0,1,0;0,0,1];

for uu = 1:N_UE
    eNB_ind = simulation_traces.UE_traces(uu).attached_eNodeB(1);
    sec_ind = simulation_traces.UE_traces(uu).attached_sector(1);
    if eNB_ind ~= 5
        continue
    end
    uu
    TP = simulation_traces.UE_traces(uu).TB_size(1,:).*uint32(simulation_traces.UE_traces(uu).ACK(1,:));
    figure(sec_ind)
    hold on
    grid on
    plot(TP,'r');
end