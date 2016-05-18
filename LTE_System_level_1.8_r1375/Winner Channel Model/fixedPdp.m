%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed delays,  powers and cluster-wice K-factors
%  for different scenarios.
% Needed when wimpar.FixedPdpUsed='yes'
%
% NOTE! Pprime power parameter values from CDL tables are changed to
% (non-dominant) ray powers. The information is same, but now the parameter
% values can be read directly from CDL tables. No change in B5 scenarios.  14.8.2008 PK

function [taus_LoS,Pprime_LoS,Kcluster_LoS, ...
        taus_NLoS,Pprime_NLoS,Kcluster_NLoS] = fixedPdp(wimpar,iterpar)

LoSConnectionLinks = iterpar.LoSConnectionLinks;
NLoSConnectionLinks = iterpar.NLoSConnectionLinks;
Scenario = iterpar.Scenario;


switch Scenario
   
    case {'A1'}
        
        %for LoS links
        taus_LoS = [0 10 25 50 65 75 75 115 115 145 195 350]*1E-9;   % [s]
        Pprime_LoS = 10.^(-[22.9 28.8 26.5 25.1 32.2 36.5 31.3 36.4 42.2 27.2 34.6 36.4]/10);          % lin.
        Kcluster_LoS = [4.7 ;...                % K-factors for CDL clusters [dB]
                       1];                      % cluster number

        %for NLoS links
        taus_NLoS = [0 5 5 5 15 15 15 20 20 35 80 85 110 115 150 175]*1E-9;      % [s]
        Pprime_NLoS = 10.^(-[15.2 19.7 15.1 18.8 16.3 17.7 17.1 21.2 13.0 14.6 23.0 25.1 25.4 24.8 33.4 29.6]/10); 
        Kcluster_NLoS = [-1000000;...                     % K-factors for CDL clusters [dB]
                           1];                            % cluster number
    case {'A2'}

        taus_NLoS = [0 0 5 10 35 35 65 120 125 195 250 305]*1E-9;   % [s]
        Pprime_NLoS = 10.^(-[13.0 21.7 16.7 24.9 29.2 19.9 13.4 23.3 33.7 29.1 34.0 35.9]/10);          % lin.
        Kcluster_NLoS = [-10000000  ;...           % K-factors for CDL clusters [dB]
                       1   ];                      % cluster number
        taus_LoS = NaN;
        Pprime_LoS = NaN;
        Kcluster_LoS = NaN;
        
    case {'B1'}
        
        taus_LoS = [0 30 55 60 105 115 250 460]*1E-9;
        Pprime_LoS = 10.^(-[24.7 20.5 27.8 23.6 26.9 30.8 32.6 44.4]/10);
        Kcluster_LoS = [3.3 ;...                   % K-factors for CDL clusters [dB]
                       1];                      % cluster number

        %this one is not updated accordingly in D111, 15.2.2007. Values from LH
        taus_NLoS = [0 90 100 115 230 240 245 285 390 430 460 505 515 595 600 615]*1E-9;   % [s]
        Pprime_NLoS = 10.^(-[14.0 13.0 13.9 21.1 21.6 24.7 25.0 25.9 32.6 36.9 35.1 38.6 36.4 45.2 44.7 42.9]/10); % lin.
        Kcluster_NLoS = [-100000 ;...            % K-factors for CDL clusters [dB]
                      1 ];                       % cluster number
        
    case {'B2'}
        
        taus_NLoS = [0 35 135 190 350 425 430 450 470 570 605 625 625 630 1600 2800]*1E-9;
        Pprime_NLoS  = 10.^(-[13.0 18.4 15.0 21.2 34.8 38.5 41.7 33.8 43.7 47.9 47.5 44.5 48.3 50.5 18.7 20.7]/10);
        Kcluster_NLoS = [-100000 ;...            % K-factors for CDL clusters [dB]
                       1];                       % cluster number 
        
        taus_LoS = NaN;
        Pprime_LoS = NaN;
        Kcluster_LoS = NaN;

    case {'B3'}
        taus_LoS = [0 0 15 25 40 40 90 130 185 280]*1E-9;
        Pprime_LoS = 10.^(-[24.5 19.6 27.6 25.8 26.8 24.1 25.6 28.2 36.4 40.7]/10);
        Kcluster_LoS = [0.3 ;...                   % K-factors for CDL clusters [dB], mistake in D111 (0.3dB)
                       1   ];                      % cluster number

        taus_NLoS = [0 5 5 10 20 20 30 60 60 65 75 110 190 290 405]*1E-9;                                % [s]
        Pprime_NLoS = 10.^(-[19.6 13.0 24.0 14.3 20.1 15.7 17.3 28.1 19.2 22.1 18.5 24.1 24.8 30.1 37.9]/10); % lin.
        Kcluster_NLoS = [-10000000;...             % K-factors for CDL clusters [dB]
                      1   ];                       % cluster number
        
      case {'B4'}

        taus_NLoS = [0 0 5 10 35 35 65 120 125 195 250 305]*1E-9;   % [s]
        Pprime_NLoS = 10.^(-[13.0 21.7 16.7 24.9 29.2 19.9 13.4 23.3 33.7 29.1 34.0 35.9]/10);          % lin.
        Kcluster_NLoS = [-10000000  ;...            % K-factors for CDL clusters [dB]
                       1   ];                       % cluster number
        taus_LoS = NaN;
        Pprime_LoS = NaN;
        Kcluster_LoS = NaN;
                 

    case {'B5a'}
        taus_LoS = [0 10 20 50 90 95 100 180 205 260]*1E-9;
        Pprime_LoS = 10.^([-0.39 -20.6 -26.8 -24.2 -15.3 -20.5 -28.0 -18.8 -21.6 -19.9 ]/10);
        Kcluster_LoS = [21.8 ;...                   % K-factors for CDL clusters [dB]
                     1 ];                       % cluster number
                 
        taus_NLoS = NaN;
        Pprime_NLoS = NaN;
        Kcluster_NLoS = NaN;
        
    case {'B5b'}    
        if wimpar.range==1
            taus_LoS = [0 5 15 20 40 45 50 70 105 115 125 135 140 240 300 345 430 440 465 625]*1E-9;                                % [s]
            Pprime_LoS = 10.^([-0.37 -15.9 -22.2 -24.9 -26.6 -26.2 -22.3 -22.3 -29.5 -17.7 -29.6 -26.6 -23.4 -30.3 -27.7 -34.8 -38.5 -38.6 -33.7 -35.2 ]/10); % lin.
            Kcluster_LoS = [ 20.0;...                    % K-factors for CDL clusters [dB]
                          1   ];                     % cluster number
        elseif wimpar.range==2
            taus_LoS = [0 5 30 45 75 90 105 140 210 230 250 270 275 475 595 690 855 880 935 1245]*1E-9;                                % [s]
            Pprime_LoS = 10.^([-1.5 -10.2 -16.6 -19.2 -20.9 -20.6 -16.6 -16.6 -23.9 -12.0 -23.9 -21.0 -17.7 -24.6 -22.0 -29.2 -32.9 -32.9 -28.0 -29.6]/10); % lin.
            Kcluster_LoS = [ 13.0;...                    % K-factors for CDL clusters [dB]
                          1   ];                     % cluster number
        elseif wimpar.range==3
            taus_LoS = [0 10 90 135 230 275 310 420 630 635 745 815 830 1430 1790 2075 2570 2635 2800 3740]*1E-9;                                % [s]
            Pprime_LoS = 10.^([-2.6 -8.5 -14.8 -17.5 -19.2 -18.8 -14.9 -14.9 -22.1 -10.3 -22.2 -19.2 -16.0 -22.9 -20.3 -27.4 -31.1 -31.2 -26.3 -27.8]/10); % lin.
            Kcluster_LoS = [ 10.0;...                    % K-factors for CDL clusters [dB]
                          1   ];                     % cluster number
        end

        taus_NLoS = NaN;
        Pprime_NLoS = NaN;
        Kcluster_NLoS = NaN;
        
    case {'B5c'}        % NOTE! PDP is a copy of B1 LOS
        taus_LoS = [0 30 55 60 105 115 250 460]*1E-9;
        Pprime_LoS = 10.^([0 -11.7 -14.8 -14.8 -13.9 -17.8 -19.6 -31.4]/10);
        Kcluster_LoS = [3.3;...                   % K-factors for CDL clusters [dB]
                     1   ];                        % cluster number
                 
        taus_NLoS = NaN;
        Pprime_NLoS = NaN;
        Kcluster_NLoS = NaN;
        
    case {'B5f'}     % NOTE! Same with B5a, except 1st tap -15dB
        taus_NLoS = [0 10 20 50 90 95 100 180 205 260]*1E-9;
        Pprime_NLoS = 10.^([-0.1 -5.3 -11.5 -8.9 0.0 -5.2 -12.7 -3.5 -6.3 -4.6]/10);
        Kcluster_NLoS = [-10000000 ;...              % K-factors for CDL clusters [dB]
                     1 ];                       % cluster number
                 
        taus_LoS = NaN;
        Pprime_LoS = NaN;
        Kcluster_LoS = NaN;
        
    case {'C1'}

        taus_LoS = [0 85 135 135 170 190 275 290 290 410 445 500 620 655 960]*1E-9;
        Pprime_LoS = 10.^(-[33.1 34.7 39.3 38.1 38.4 35.0 42.2 34.3 36.2 45.2 39.5 45.1 41.5 43.5 45.6]/10);
        Kcluster_LoS = [12.9;...                 % K-factors for CDL clusters [dB]
                     1 ];                        % cluster number

        taus_NLoS = [0 25 35 35 45 65 65 75 145 160 195 200 205 770]*1E-9;                                % [s]
        Pprime_NLoS = 10.^(-[13.0 20.5 23.5 16.2 16.1 27.0 19.4 16.1 17.6 21.0 20.2 16.1 22.5 35.4]/10); % lin.
        Kcluster_NLoS = [-10000000;...           % K-factors for CDL clusters [dB]
                            1   ];               % cluster number

    case {'C2'}
        
        taus_LoS = [0 0 30 85 145 150 160 220]*1E-9;
        Pprime_LoS = 10.^(-[30.6 26.2 28.3 29.7 28.2 31.2 28.3 36.1]/10);
        Kcluster_LoS = [7.0;...                   % K-factors for CDL clusters [dB]
                        1 ];                      % cluster number
        
        taus_NLoS = [0 60 75 145 150 190 220 335 370 430 510 685 725 735 800 960 1020 1100 1210 1845]*1E-9;   % [s]
        Pprime_NLoS = 10.^(-[19.5 16.4 15.0 13.0 14.9 16.4 13.4 17.7 20.8 20.8 22.3 25.0 21.5 26.2 24.2 33.8 27.5 24.7 30.2 29.7]/10);                      % lin.
        Kcluster_NLoS = [-1000000;...               % K-factors for CDL clusters [dB]
                           1];                      % cluster number
        
    case {'C3'}

        taus_NLoS = [0 5 35 60 160 180 240 275 330 335 350 520 555 555 990 1160 1390 1825 4800 7100]*1E-9;
        Pprime_NLoS = 10.^(-[16.5 22.0 17.6 22.2 13.0 14.7 15.7 20.0 18.9 19.7 14.3 18.3 17.9 22.4 25.3 25.2 33.8 38.4 22.7 26.0]/10);
        Kcluster_NLoS = [-100000 ;...               % K-factors for CDL clusters [dB]
                          1 ];                      % cluster number  
                   
        taus_LoS = NaN;
        Pprime_LoS = NaN;
        Kcluster_LoS = NaN;

    case {'C4'}

        taus_NLoS = [0 15 95 145 195 215 250 445 525 815 1055 2310]*1E-9;
        Pprime_NLoS = 10.^(-[13.0 19.9 16.6 29.3 21.5 28.9 19.9 27.1 13.8 26.6 30.8 45.2 ]/10);
        Kcluster_NLoS = [-100000 ;...               % K-factors for CDL clusters [dB]
                          1 ];                      % cluster number  
                   
        taus_LoS = NaN;
        Pprime_LoS = NaN;
        Kcluster_LoS = NaN;
        
    case {'D1'}    
        taus_LoS = [0 20 20 25 45 65 65 90 125 180 190]*1E-9;
        Pprime_LoS = 10.^(-[22.8 28.5 29.2 25.3 33.5 31.9 34.2 36.6 39.1 42.4 41.3]/10);
        Kcluster_LoS = [5.7 ;...                   % K-factors for CDL clusters [dB]
                         1 ];                      % cluster number

        taus_NLoS = [0 0 5 10 20 25 55 100 170 420]*1E-9;                                % [s]
        Pprime_NLoS = 10.^(-[13.0 14.8 16.3 14.8 18.3 20.1 22.0 17.2 25.4 39.5]/10); % lin.
        Kcluster_NLoS = [-10000000;...             % K-factors for CDL clusters [dB]
                          1   ];                   % cluster number
    case {'D2a'}
        taus_LoS = [0 45 60 85 100 115 130 210]*1E-9;   % [s]
        Pprime_LoS = 10.^(-[28.8 27.8 30.2 29.5 28.1 28.7 30.8 30.3]/10);                      % lin.
        Kcluster_LoS = [ 6.0;...                         % K-factors for CDL clusters [dB]
                         1 ];                            % cluster number
                   
        taus_NLoS = NaN;
        Pprime_NLoS = NaN;
        Kcluster_NLoS = NaN;

end % switch