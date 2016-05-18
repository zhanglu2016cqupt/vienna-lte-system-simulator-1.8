%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed Doppler frequencies of scatterers
% Needed for stationary feeder scenarios B5
function scatter_freq=fixedScatterFreq(wimpar,iterpar) %KTH

M=wimpar.NumSubPathsPerPath;
Scenario = iterpar.Scenario;
switch lower(Scenario) %KTH
    case {'b5a'}%KTH
        scatter_freq=[41.6,-21.5,-65.2,76.2,10.5,-20.2,1.3,2.2,-15.4,48.9]*wimpar.CenterFrequency/5.25e9*1e-3;%KTH
        scatter_freq=[scatter_freq;zeros(M-1,length(scatter_freq))];
    case {'b5b'}%KTH
        scatter_freq=[744,-5,-2872,434,295,118,2576,400,71,3069,1153,-772,1298,-343,-7,-186,...
                     -2287,26,-1342,-61]*wimpar.CenterFrequency/5.25e9*1e-3;%KTH
        scatter_freq=[scatter_freq;zeros(M-1,length(scatter_freq))];
        
    case {'b5c'}    % NOTE! copied from B5a, removed two last clusters
        scatter_freq=[-127,385,-879,0,0,-735,-274,691]*wimpar.CenterFrequency/5.25e9*1e-3;
        scatter_freq=[scatter_freq;zeros(M-1,length(scatter_freq))];
        scatter_freq(:,4) = [45:0.5:54.5].'*wimpar.CenterFrequency/5.25e9*1e-3;
        scatter_freq(:,5) = [-55:-0.5:-64.5].'*wimpar.CenterFrequency/5.25e9*1e-3;

    case {'b5f'}    % NOTE! copy of B5a 
        scatter_freq=[41.6,-21.5,-65.2,76.2,10.5,-20.2,1.3,2.2,-15.4,48.9]*wimpar.CenterFrequency/5.25e9*1e-3;%KTH
        scatter_freq=[scatter_freq;zeros(M-1,length(scatter_freq))];
end; %switch