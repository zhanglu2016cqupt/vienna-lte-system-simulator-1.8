%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed AoAs for different scenarios
% Needed when wimpar.FixedAnglesUsed='yes'
function [aoas_los,aoa_clusterAS_los, ...
          aoas_nlos,aoa_clusterAS_nlos] = fixedAoas(wimpar,iterpar)

Scenario = iterpar.Scenario;

switch Scenario

    case {'A1'}

        aoas_los = [0 -110 102 -134 121 -134 -118 -134 149 105 129 -134];
        aoa_clusterAS_los = 5;      % Cluster ASA [deg], [1, table 7-1]

        aoas_nlos = [41 -70 39 66 -49 59 -55 -78 0 95 86 95 -96 -94 123 -111];
        aoa_clusterAS_nlos = 5;      % Cluster ASA [deg], [1, table 7-2]

    case {'A2'}
        aoas_nlos = [0 32 -21 37 -43 28 -49 -34 -49 43 49 51];
        aoa_clusterAS_nlos = 5;      % Cluster ASA [deg], [1, table 7-3]
        
        aoas_los = NaN;
        aoa_clusterAS_los = NaN; 
        
    case {'B1'}
        aoas_los = [0 45 63 -69 61 -69 -73 92];
        aoa_clusterAS_los = 18;      % Cluster ASA [deg], [1, table 7-4]

        aoas_nlos = [-20 0 57 -55 57 67 -68 70 -86 -95 -92 -99 94 111 110 -107];
        aoa_clusterAS_nlos = 22;      % Cluster ASA [deg], [1, table 7-5]

    case {'B2'}

        aoas_nlos = [0 -46 -92 57 -92 -100 -106 90 -110 -117 -116 -111 -118 121 15 -25];
        aoa_clusterAS_nlos = 22;      % Cluster ASA [deg], [1, table 7-6]
        
        aoas_los = NaN;
        aoa_clusterAS_los = NaN; 

    case {'B3'}
        aoas_los = [0 -53 -79 -74 76 80 -73 80 -100 -108];
        aoa_clusterAS_los = 5;      % Cluster ASA [deg], [1, table 7-7]

        aoas_nlos = [-73 0 -94 -46 75 -46 -59 107 71 86 67 95 98 117 142];
        aoa_clusterAS_nlos = 13;      % Cluster ASA [deg], [1, table 7-8]

    case {'B4'}

        aoas_nlos = [0 102 -66 -119 139 91 157 -111 157 138 158 165];
        aoa_clusterAS_nlos = 8;      % Cluster ASA [deg], [1, table 7-9]

        aoas_los = NaN;
        aoa_clusterAS_los = NaN; 

    case {'B5a'}
        aoas_los = [0 0.2 1.5 2.0 0 3.6 -0.7 4.0 -2.0 -4.1];
        aoa_clusterAS_los = 0.5;      %KTH cluster AS at  MS [deg], [1, table 7-17]
   
        aoas_nlos = NaN;
        aoa_clusterAS_nlos = NaN; 
        
    case {'B5b'}            
        if wimpar.range==1
            aoas_los = [0 70 -27.5 106.4 94.8 -94.0 48.6 -96.6 41.7 -83.3 176.8 93.7 -6.4 160.3 -50.1 -149.6 161.5 68.7 41.6 142.2];
            aoa_clusterAS_los = 2;      %KTH cluster AS at  MS [deg], [1, table 7-19]
        elseif wimpar.range==2
            aoas_los = [0 70 -27.5 106.4 94.8 -94.0 48.6 -96.6 41.7 -83.3 176.8 93.7 -6.4 160.3 -50.1 -149.6 161.5 68.7 41.6 142.2];
            aoa_clusterAS_los = 2;      %KTH cluster AS at  MS [deg], [1, table 7-20]
        elseif wimpar.range==3
            aoas_los = [0 70 -27.5 106.4 94.8 -94.0 48.6 -96.6 41.7 -83.3 176.8 93.7 -6.4 160.3 -50.1 -149.6 161.5 68.7 41.6 142.2];
            aoa_clusterAS_los = 2;      %KTH cluster AS at  MS [deg], [1, table 7-21]
        end
        
        aoas_nlos = NaN;
        aoa_clusterAS_nlos = NaN; 
        
    case {'B5c'}        % Parameters same with B1 LOS
        aoas_los = [0 45 63 -69 61 -69 -73 92];
        aoa_clusterAS_los = 18; 
        
        aoas_nlos = NaN;
        aoa_clusterAS_nlos = NaN; 
        
    case {'B5f'}    % Copy of B5a
        aoas_nlos = [0 0.2 1.5 2.0 0 3.6 -0.7 4.0 -2.0 -4.1];
        aoa_clusterAS_nlos = 0.5;  
        
        aoas_los = NaN;
        aoa_clusterAS_los = NaN; 
            
    case {'C1'}
        aoas_los = [0 -144 -159 155 156 -146 168 -176 149 -176 -159 -176 -165 -171 177];
        aoa_clusterAS_los = 5;      % Cluster ASA [deg], [1, table 7-10]

        aoas_nlos = [0 -71 -84 46 -66 -97 -66 -46 -56 73 70 -46 -80 123];
        aoa_clusterAS_nlos = 10;      % Cluster ASA [deg], [1, table 7-11]

    case {'C2'}
        aoas_los = [0 -120 129 -135 -129 141 -129 -158];
        aoa_clusterAS_los = 12; 
        
        aoas_nlos = [61 44 -34 0 33 -44 -67 52 -67 -67 -73 -83 -70 87 80 109 91 -82 99 98];
        aoa_clusterAS_nlos = 15;      % Cluster ASA [deg], [1, table 7-12]
     
     case {'C3'}
        aoas_nlos = [-52 -83 -60 -85 0 -36 46 74 68 -72 -62 -64 -62 85 -98 -97 127 140 25 40];
        aoa_clusterAS_nlos = 15;      % Cluster ASA [deg], [1, table 7-13]

        aoas_los = NaN;
        aoa_clusterAS_los = NaN; 

    case {'C4'}
        aoas_nlos = [0 -91 65 -139 101 -138 -91 130 -146 128 -146 196];
        aoa_clusterAS_nlos = 8;      % Cluster ASA [deg], [1, table 7-13]

        aoas_los = NaN;
        aoa_clusterAS_los = NaN; 

    case {'D1'}        
        aoas_los = [0 44 -45 -48 50 -48 51 -54 57 -60 59];
        aoa_clusterAS_los = 3;      % Cluster ASA [deg], [1, table 7-14]

        aoas_nlos = [0 28 38 -55 48 -55 62 42 -73 107];
        aoa_clusterAS_nlos = 3;      % Cluster ASA [deg], [1, table 7-15]
            
    case {'D2a'}
        aoas_los = [0 -80 86 84.4 87.5 -82.2 87.5 86.2];
        aoa_clusterAS_los = 3;      % Cluster ASA [deg], [1, table 7-16]
        
        aoas_nlos = NaN;
        aoa_clusterAS_nlos = NaN; 

end % switch