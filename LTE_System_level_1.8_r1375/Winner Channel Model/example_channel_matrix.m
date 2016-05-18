
% Example of channel matrix generation for 10 MS-BS links (semi-random layout)
% 
% Authors: Milan Narandži? (TUI), Martin Käske (TUI)

%% Generation of 3D-AA has to be performed in pre-processing phase
%Arrays=example_syntetic_arrays;

%% Network layout
MsAAIdx = [1 1 2 3];
BsAAIdxCell = {[1 3]; [2]; [1 1 2]};
layoutpar=layoutparset(MsAAIdx, BsAAIdxCell,10,Arrays);
layoutpar.PropagConditionVector=zeros(1,10);  % (NLOS=0/LOS=1)
layoutpar.ScenarioVector=3*ones(1,10);  % B1 scenario
layoutpar.Stations(2).Rot=[0 0 pi]; %Rotate sector array -180
layoutpar.Stations(5).Rot=[0 0 -pi/2]; % 90
layoutpar.Stations(6).Rot=[0 0 2*pi/3]; % -120

% Visualization
NTlayout(layoutpar);
%% Model and simulation settings
wimpar=wimparset;
wimpar.NumTimeSamples=1000;     % 100 time samples per link

%% Generation of channel matrix
[H1,delays,out]=wim(wimpar, layoutpar);

%% External initialization of structural parameters in consequtive calls
[H2,delays,out]=wim(wimpar,layoutpar,out);