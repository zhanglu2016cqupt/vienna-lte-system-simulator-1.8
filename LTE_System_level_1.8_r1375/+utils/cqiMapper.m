classdef cqiMapper < handle
% This class abstracts the mapping from SNR to CQI and viceversa. This one
% works with a linear approximation.
% (c) Josep Colom Ikuno, INTHFT, 2008

   properties
       max_CQI
       min_CQI
       p % p(1) and p(2) represent the coefficients
       lesscons
   end

   methods
       % Class constructor
       function obj = cqiMapper(p,min_CQI,max_CQI,lesscons)
           obj.p = p;
           obj.max_CQI = max_CQI;
           obj.min_CQI = min_CQI;
           obj.lesscons = lesscons;
       end
       % Convert from SINR to CQI
       function CQIs = SINR_to_CQI(obj,SINRs)
           % Clipped CQI feedback (from obj.min_CQI to obj.max_CQI)
           CQIs = obj.p(1)*SINRs + obj.p(2);
           min_CQI = obj.min_CQI;
           max_CQI = obj.max_CQI;
           CQIs((~isfinite(CQIs))|CQIs<min_CQI) = min_CQI;
           CQIs(CQIs>max_CQI)                   = max_CQI;
       end
       
       function CQIs_back = SINR_to_CQI2(obj,SINRs,CQIs)
            % table = [-7;-5.661;-3.595;-1.572;0.537;2.545;4.593;6.408;8.479;10.34;12.22;14.13;15.85;17.79;19.86];
            table = [-6.934;-5.147;-3.18;-1.254;0.761;2.70;4.697;6.528;8.576;10.37;12.3;14.18;15.89;17.82;19.83];
            for i = 1:length(SINRs)
                temp = zeros(size(CQIs));
                temp(table(CQIs) <= SINRs(i)) = 1;
                CQIs_back(i) = find(temp,1,'last');
                if isempty(CQIs_back(i))
                    CQIs_back(i) = 1;
                end
            end
       end
           
       % Convert from CQI to SINR [dB]
       % Please take into account that NO INPUT CHECKING IS DONE!!
       % Just take care that the input CQI values are correct, if not
       % results could be inconsistent.
       function SINRs = CQI_to_SINR(obj,CQIs)
%            table = [-6.934;-5.147;-3.18;-1.254;0.761;2.70;4.697;6.528;8.576;10.37;12.3;14.18;15.89;17.82;19.83;21];
           if obj.lesscons
               SINR_1 = (CQIs-obj.p(2)) / obj.p(1);
               SINR_2 = (CQIs+1-obj.p(2)) / obj.p(1);
               SINRs  = SINR_1+(SINR_2-SINR_1)/2.5; % output in dB
%                     SINRs
%                 SINR_1 = table(CQIs);
%                 SINR_2 = table(CQIs+1);
%                 SINRs = SINR_1+(SINR_2-SINR_1)/3.5; % output in dB
%                 SINRs
%                 error('j')
           else
               SINRs = (CQIs-obj.p(2)) / obj.p(1); % output in dB
           end
       end
   end
end 
