%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed AoDs for different scenarios
% Needed when wimpar.FixedAnglesUsed='yes'
function [aods_los,aod_clusterAS_los, ...
          aods_nlos,aod_clusterAS_nlos]=fixedAods(wimpar,iterpar)
    
Scenario = iterpar.Scenario;

switch Scenario

    case {'A1'}
        aods_los = [0 -107 -100 131 118 131 116 131 146 102 -126 131];
        aod_clusterAS_los = 5;      % Cluster ASD [deg], [1, table 6-1]


        aods_nlos = [45 77 43 72 54 -65 -60 85 0 -104 95 -104 -105 103 -135 -122];
        aod_clusterAS_nlos = 5;      % Cluster ASD [deg], [1, table 6-2]

        
     case {'A2'}
        
        aods_nlos = [0 102 -66 -119 139 91 157 -111 157 138 158 165];
        aod_clusterAS_nlos = 8;      % Cluster ASD [deg], [1, table 6-3]
        
        aods_los = NaN;
        aod_clusterAS_los = NaN;      
        
    
    case {'B1'}
        
        aods_los = [0 5 8 8 7 8 -9 11];
        aod_clusterAS_los = 3;      % Cluster ASD [deg], [1, table 6-4]

        aods_nlos = [8 0 -24 -24 -24 29 29 30 -37 41 -39 -42 -40 47 47 46];
        aod_clusterAS_nlos = 10;      % Cluster ASD [deg], [1, table 6-5]

        
     case {'B2'}
       
        aods_nlos = [0 20 40 25 40 -44 -46 39 -48 -51 -51 -48 -51 53 -110 75];
        aod_clusterAS_nlos = 10;      % Cluster ASD [deg], [1, table 6-6]
     
        aods_los = NaN;
        aod_clusterAS_los = NaN;    

        
    case {'B3'}
       
        aods_los = [0 -23 -34 -32 33 -35 32 -35 -43 47];
        aod_clusterAS_los = 5;      % Cluster ASD [deg], [1, table 6-7]

        aods_nlos = [-16 0 -21 -10 17 -10 -13 -24 -16 19 -15 -21 22 -26 -32];
        aod_clusterAS_nlos = 6;      % Cluster ASD [deg], [1, table 6-8]

        
     case {'B4'}

        aods_nlos = [29 0 20 -18 18 20 29 24 29 -21 36 46];
        aod_clusterAS_nlos = 5;      % Cluster ASD [deg], [1, table 6-3]
        
        aods_los = NaN;
        aod_clusterAS_los = NaN;      

     
    case {'B5a'} %KTH Wireless feeder, Rooftop-Rooftop
        aods_los =[0 0.9 0.3 -0.3 3.9 -0.8 4.2 -1.0 5.5 7.6];
        aod_clusterAS_los = 0.5;      %KTH ZDSC AS at FS [deg], [1, table 6-19]
        
        aods_nlos = NaN;
        aod_clusterAS_nlos = NaN; 
        
    case {'B5b'}    
        if wimpar.range==1
            aods_los = [0 -71.7 167.4 -143.2 34.6 -11.2 78.2 129.2 -113.2 -13.5 145.2 -172.0 93.7 106.5 -67.0 -95.1 -2.0 66.7 160.1 -21.8];
            aod_clusterAS_los = 2;      %KTH cluster AS at BS [deg], [1, table 6-21]
            
        elseif wimpar.range==2
            aods_los = [0 -71.7 167.4 -143.2 34.6 -11.2 78.2 129.2 -113.2 -13.5 145.2 -172.0 93.7 106.5 -67.0 -95.1 -2.0 66.7 160.1 -21.8];
            aod_clusterAS_los = 2;      %KTH cluster AS at BS [deg], [1, table 6-22]
            
        elseif wimpar.range==3
            aods_los = [0 -71.7 167.4 -143.2 34.6 -11.2 78.2 129.2 -113.2 -13.5 145.2 -172.0 93.7 106.5 -67.0 -95.1 -2.0 66.7 160.1 -21.8];
            aod_clusterAS_los = 2;      %KTH cluster AS at BS [deg], [1, table 6-23]
        end
    
        aods_nlos = NaN;
        aod_clusterAS_nlos = NaN;  
    
    case {'B5c'} %KTH Wireless feeder, Rooftop-Rooftop. Parameters same with B1 LOS
        aods_los = [0 5 8 8 7 8 -9 11];
        aod_clusterAS_los = 3; 
    
        aods_nlos = NaN;
        aod_clusterAS_nlos = NaN; 
        
    case {'B5f'} % Copy of B5a
        aods_nlos = [0 0.9 0.3 -0.3 3.9 -0.8 4.2 -1.0 5.5 7.6];
        aod_clusterAS_nlos = 0.5; 
        
        aods_los = NaN;
        aod_clusterAS_los = NaN;  
        
    case {'C1'}

        aods_los = [0 -29 -32 -31 31 29 -33 35 -30 35 -32 35 33 34 35];
        aod_clusterAS_los = 5;      % Cluster ASD [deg], [1, table 6-9]

        aods_nlos = [0 13 -15 -8 12 -17 12 -8 -10 -13 12 8 14 22];
        aod_clusterAS_nlos = 2;      % Cluster ASD [deg], [1, table 6-10]

    case {'C2'}
        aods_los = [0 -24 26 -27 26 28 26 -32];
        aod_clusterAS_los = 6;      % Cluster ASD [deg], [1, table 6-11]
        
        aods_nlos = [11 -8 -6 0 6 8 -12 -9 -12 -12 13 15 -12 -15 -14 19 -16 15 18 17];
        aod_clusterAS_nlos = 2;      % Cluster ASD [deg], [1, table 6-12]
        
     case {'C3'}

        aods_nlos = [-9 14 -10 -14 0 -6 7 -12 11 -12 -10 -10 -10 14 16 16 21 -23 -135 80];
        aod_clusterAS_nlos = 2;      % Cluster ASD [deg], [1, table 6-13]

        aods_los = NaN;
        aod_clusterAS_los = NaN;
        
    case {'C4'}
        aods_nlos = [0 28 -20 43 -31 43 28 -40 45 -39 45 -61];
        aod_clusterAS_nlos = 5;      % Cluster ASD [deg], [1, table 6-14]

        aods_los = NaN;
        aod_clusterAS_los = NaN; 
    
    case {'D1'}
        aods_los = [0 17 17 18 -19 28 -19 -20 -22 23 -22];
        aod_clusterAS_los = 2;      % Cluster ASD [deg], [1, table 6-15]

        aods_nlos = [0 -8 -10 15 13 15 -17 -12 20 29];
        aod_clusterAS_nlos = 2;      % Cluster ASD [deg], [1, table 6-16]
      
    case {'D2a'}
        aods_los = [0 12.7 -13.6 13.4 -13.9 -13 -13.9 13.7];
        aod_clusterAS_los = 2;      % Cluster ASD [deg], [1, table 6-17]
        
        aods_nlos = NaN;
        aod_clusterAS_nlos = NaN; 

end % switch