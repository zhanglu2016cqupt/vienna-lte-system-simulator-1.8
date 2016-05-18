classdef blerCurvesFast
% This class store the BLER curves obtained from Link Level simulations.
% Optimized version with respect ot the original version included with v1.0
% of the LTE System Level Simulator
% (c) Josep Colom Ikuno, INTHFT, 2010

   properties
       BLER_data
       SNR_data
       resolution = 0.05; % We will set the BLER curve to this fixed resolution to speed up BLER lookup
       SINR_range
       multipliers
   end

   methods
       % Class constructor. As input receives two [Nx1] cell arrays containing
       % the BLER curves for N CQI values
       function obj = blerCurvesFast(BLER_curves,SNR_data)
           N_CQIs = length(BLER_curves);
           min_SINR = min(SNR_data{1});
           max_SINR = max(SNR_data{1});
           for i_=1:N_CQIs
               min_SINR = min([min_SINR SNR_data{i_}]);
               max_SINR = max([max_SINR SNR_data{i_}]);
           end
           SINR_range  = min_SINR:obj.resolution:max_SINR;
           BLER_matrix = zeros(N_CQIs,length(SINR_range));
           % Fill in the BLER data
           for i_=1:N_CQIs
               current_SNR  = SNR_data{i_};
               current_BLER = BLER_curves{i_};
               if current_SNR(1)>min_SINR
                   current_SNR  = [min_SINR current_SNR];
                   current_BLER = [current_BLER(1) current_BLER];
               end
               if current_SNR(end)<max_SINR
                   current_SNR  = [current_SNR max_SINR];
                   current_BLER = [current_BLER current_BLER(end)];
               end
               BLER_matrix(i_,:) = interp1(current_SNR,current_BLER,SINR_range);
           end
           obj.BLER_data   = BLER_matrix';
           obj.SNR_data    = SINR_range;
           obj.SINR_range  = [min_SINR max_SINR];
           obj.multipliers = reshape(0:length(SINR_range):length(SINR_range)*(N_CQIs-1),[],1);
       end
       
       % Get the BLER given a CQI, and SNR. changed to support a query with
       % multiple CQIs. Input values (SINRs) in LOG scale.
       function BLER = get_BLER(obj,CQI,SNR)
           SNR_idxs = round((SNR-obj.SINR_range(1))/obj.resolution + 1);
           SNR_length = length(obj.SNR_data);
           SNR_idxs(SNR_idxs<1) = 1;
           SNR_idxs(SNR_idxs>SNR_length) = SNR_length;
           CQI_multipliers = obj.multipliers(CQI);
           if length(SNR_idxs)>1
               SNR_idxs = repmat(SNR_idxs,[length(CQI) 1]);
               CQI_multipliers = repmat(CQI_multipliers,[1 size(SNR_idxs,2)]);
           end
           BLER = obj.BLER_data(SNR_idxs+CQI_multipliers);
       end
       
       % Get the BLER curves for the specified CQI values with the
       % specified SINR values. Each SINR value corresponds to one CQI. ie.
       % both vectors must have the same length
       function BLER = get_BLER_CQI(obj,CQI,SNR)
           SNR_idxs = round((SNR-obj.SINR_range(1))/obj.resolution + 1);
           SNR_length = length(obj.SNR_data);
           SNR_idxs(SNR_idxs<1) = 1;
           SNR_idxs(SNR_idxs>SNR_length) = SNR_length;
           BLER = zeros(length(CQI),1);
           
           BLER_indexes = SNR_idxs(:)+obj.multipliers(CQI);
           BLER = obj.BLER_data(BLER_indexes);
       end
       
       % Plot all of the CQI curves
       function plot(obj,fig_number)
           cqi_figure = figure(fig_number);
           axes('Parent',cqi_figure,'YScale','log','YMinorTick','on');
           hold all;
           box on;
           initial_cqi = 1;
           final_cqi = size(obj.BLER_data,2);
           for cqi=initial_cqi:final_cqi
               semilogy(obj.SNR_data,obj.BLER_data(:,cqi),'DisplayName',['CQI ' num2str(cqi)]);
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