function indLOS=LOSprobability(wimpar,linkpar,fixpar,iterpar)
%LOSPROBABILITY Random LOS/NLOS condition generation for WIM
%   INDLOS=LOSPROBABILITY(WIMPAR,LINKPAR) is a vector defining links
%   propagation condition. Vector elements have values '0' or '1' and length is
%   number of links. '0' stands for NLOS and '1' stands for LOS link,
%   vector length is number of links. LOS/NLOS condition is drawn randomly
%   for each link according to LOS probabilities defined in [1, Table 4-7].
%
%   Ref. [1]: D1.1.1 V1.0, "WINNER II interim channel models"
%
%   See also WIM, WIMPARSET, LAYOUTPARSET, ANTPARSET.

%   Authors: Pekka Kyösti (EBIT), Mikko Alatossava (CWC/UOULU)

%   Modifications: Parameters from [1, Table 4-7] are used:         15.1.2007 Marko
%                  Divider for C1 changed to 200, used to be 500    12.2.2007 MikkoA

NumLinks = length(iterpar.UserIndeces);
Scenario = iterpar.Scenario;
MsBsDistance = linkpar.MsBsDistance(iterpar.UserIndeces);


% Probability of LOS
switch upper(Scenario)      % See equations 3.20-25 in [1, sect 3.1.6 ]
    
    case {'A1'}
        dBp = 2.5;      % max distance of LOS probability 1 [m]
        pLOS = ones(size(MsBsDistance));
        indBp = MsBsDistance > dBp;
        pLOS(indBp) = 1-0.9*(1-(1.24-0.61*log10(MsBsDistance(indBp))).^3).^(1/3);        
        
    case {'B1'}
        if isnan(linkpar.Dist1) % NaN default -> will be drawn randomly
            StreetWidth = linkpar.StreetWidth(iterpar.UserIndeces);
            Dist1 = 1;
            Dist2 = 1;
            while (Dist1>5000 | Dist1<10)
                Dist2 = (MsBsDistance-StreetWidth/2).*rand(1,length(MsBsDistance)) + StreetWidth/2;
                Dist1 = sqrt(MsBsDistance.^2 - Dist2.^2);
            end
        else
            StreetWidth = linkpar.StreetWidth(iterpar.UserIndeces);
            Dist1 = linkpar.Dist1(iterpar.UserIndeces);
            for linkNum = 1:length(MsBsDistance) % check applicability and change Dist1 if needed
                if (MsBsDistance(linkNum)^2 < (StreetWidth(linkNum)/2)^2 + Dist1(linkNum)^2)
                    Dist1(linkNum) = sqrt(MsBsDistance(linkNum)^2 - (StreetWidth(linkNum)/2)^2);
                end
            end
            
            Dist2 = sqrt(MsBsDistance.^2 - Dist1.^2);
        end
        
        dBp = 15;      % max distance of LOS probability 1 [m]
        pLOS = ones(size(MsBsDistance));
        indBp = MsBsDistance > dBp;
        pLOS(indBp) = 1-(1-(1.56-0.48*log10(sqrt(Dist1(indBp).^2+Dist2(indBp).^2))).^3).^(1/3);

    case {'B3'}
        % Scenario B3 has two LOS probability functions, default is eq 3.22
        
        % P(LOS) for big factory halls, airport and railway stations
        dBp = 10;      % max distance of LOS probability 1 [m]
        pLOS = ones(size(MsBsDistance));
        indBp = MsBsDistance > dBp;
        pLOS(indBp) = exp(-(MsBsDistance(indBp)-10)/45); 
        
        % P(LOS) for big lecture halls or conference halls
%         dBp = 5;      % max distance of LOS probability 1 [m]
%         pLOS = ones(size(MsBsDistance));
%         indBp = MsBsDistance > dBp;
%         pLOS(indBp) = 1-(MsBsDistance(indBp)-5)/150;
%         if sum(MsBsDistance > 40)>0;
%             warning('MsBsDistance exceeds maximum 40 m in B3 P(LOS) computation')
%         end

    case {'C1'}
        pLOS = exp(-MsBsDistance/200);
        
    case {'C2'}
        pLOS = min(18./MsBsDistance,ones(size(MsBsDistance))).*...
               (1-exp(-MsBsDistance/63)) + exp(-MsBsDistance/63);
        
    case {'D1', 'D2a'}
        pLOS = exp(-MsBsDistance/1000);
        
end

% output, 0 for NLOS and 1 for LOS links
indLOS = rand(1,NumLinks)<pLOS;