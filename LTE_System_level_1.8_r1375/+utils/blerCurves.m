classdef blerCurves < handle
% This class store the BLER curves obtained from Link Level simulations
% (c) Josep Colom Ikuno, INTHFT, 2008

   properties
       BLER_data
       SNR_data
   end

   methods
       % Class constructor. As input receives two [Nx1] cell arrays containing
       % the BLER curves for N CQI values
       function obj = blerCurves(BLER_curves,SNR_data)
           obj.BLER_data = BLER_curves;
           obj.SNR_data = SNR_data;
       end
       
       % Get the BLER given a CQI, and SNR. changed to support a query with
       % multiple CQIs. Input values (SINRs) in LOG scale.
       function BLER = get_BLER(obj,CQI,SNR)
           BLER = zeros(length(CQI),length(SNR));
           
           for CQI_idx = 1:length(CQI)
               current_SNR = SNR;
               minimum_SNR = obj.SNR_data{CQI(CQI_idx)}(1);
               maximum_SNR = obj.SNR_data{CQI(CQI_idx)}(end);
               current_SNR(current_SNR < minimum_SNR) = minimum_SNR;
               current_SNR(current_SNR > maximum_SNR) = maximum_SNR;
               BLER(CQI_idx,:) = interp1(obj.SNR_data{CQI(CQI_idx)},obj.BLER_data{CQI(CQI_idx)},current_SNR);
           end
       end
       
       % Get the BLER curves for the specified CQI values with the
       % specified SINR values. Each SINR value corresponds to one CQI. ie.
       % both vectors must have the same length
       function BLER = get_BLER_CQI(obj,CQI,SNR)
           BLER = zeros(length(CQI),1);
           
           for CQI_idx = 1:length(CQI)
               current_SNR = SNR(CQI_idx);
               minimum_SNR = obj.SNR_data{CQI(CQI_idx)}(1);
               maximum_SNR = obj.SNR_data{CQI(CQI_idx)}(end);
               current_SNR(current_SNR < minimum_SNR) = minimum_SNR;
               current_SNR(current_SNR > maximum_SNR) = maximum_SNR;
               BLER(CQI_idx) = interp1(obj.SNR_data{CQI(CQI_idx)},obj.BLER_data{CQI(CQI_idx)},current_SNR);
           end
       end
       
       % Plot all of the CQI curves
       function plot(obj,fig_number)
           cqi_figure = figure(fig_number);
           axes('Parent',cqi_figure,'YScale','log','YMinorTick','on');
           hold all;
           box on;
           initial_cqi = 1;
           final_cqi = length(obj.BLER_data);
           for cqi=initial_cqi:final_cqi
               semilogy(obj.SNR_data{cqi},obj.BLER_data{cqi},'DisplayName',['CQI ' num2str(cqi)],'LineWidth',2);
               ylim([1e-3,1]);
           end
           grid on;
           xlabel('SNR [dB]');
           ylabel('BLER');
           title(['LTE BLER for CQIs ' num2str(initial_cqi) ' to ' num2str(final_cqi)]);
           legend('Location','SouthEastOutside');
       end
   end
end 
