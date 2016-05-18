function linkpar=layout2link(layoutpar)
%LAYOUT2LINK Layout to link parameter conversion for WIM.
%   LINKPAR=LAYOUT2LINK(LAYOUTPAR) returns extended set of link parameters
%   in the case of layout parameters are defined. It converts layout
%   parameters to Ms/Bs distances, LOS directions etc. LAYOUT2LINK is used with 
%   WIM the following way: [..] = wim(wimparset,layout2link(layoutpar),antpar).
%
%   The output parameters are:
%
%   Stations     - directly from LAYOUTPAR
%   NofSect      - directly from LAYOUTPAR 
%   Pairing      - directly from LAYOUTPAR
%   ThetaBs      - calculated from LAYOUTPAR
%   ThetaMs      - calculated from LAYOUTPAR
%   MsBsDistance - calculated from LAYOUTPAR
%   MsVelocity   - calculated from LAYOUTPAR 
%   MsDirection  - calculated from LAYOUTPAR
%   StreetWidth  - defined in LAYOUTPAR
%   Dist1        - defined in LAYOUTPAR
%
%   See [1, Fig. 6.1 and 6.2].
%
%   Ref. [1]: D1.1.1 V1.0, "WINNER II interim channel models"
%
%   See also WIM, LAYOUTPARSET, WIMPARSET, ANTPARSET.

%   Authors: Pekka Kyösti (EBIT), Martin Käske (TUI)
%               
NofSect = layoutpar.NofSect;
Pairing = layoutpar.Pairing;
Stations = layoutpar.Stations;
 
% linkpar struct with layout parameters included
linkpar=struct( 'Stations', Stations,...
                'NofSect', [],...
                'Pairing', [],...
                'ScenarioVector',[],...
                'PropagConditionVector',[],...
                'ThetaBs', [],...
                'ThetaMs', [],...
                'MsHeight', [],...
                'BsHeight', [],...
                'MsBsDistance', [],...
                'MsVelocity', [],...
                'MsDirection', [],...
                'StreetWidth', [],...
                'NumFloors', [],...
                'NumPenetratedFloors', [],...
                'Dist1', []);
                

% MS-BS distance
linkpar.MsBsDistance = StationDistXY(Stations(Pairing(1,:)), Stations(Pairing(2,:)));
BsHeight = [Stations(Pairing(1,:)).Pos]; BsHeight=BsHeight(3,:);
MsHeight = [Stations(Pairing(2,:)).Pos]; MsHeight=MsHeight(3,:);
% LOS direction from BS array to MS array
[ThetaBs, ThetaMs]=StationDirectionXY(Stations(Pairing(1,:)),Stations(Pairing(2,:)));
ThetaBs=prin_value(ThetaBs*180/pi);
ThetaMs=prin_value(ThetaMs*180/pi);

% the rest of the link parameters
[linkpar.MsDirection,linkpar.MsVelocity]=StationVelocityXY(Stations(Pairing(2,:)));
linkpar.MsDirection=linkpar.MsDirection*180/pi;
linkpar.ScenarioVector = layoutpar.ScenarioVector;
linkpar.PropagConditionVector = layoutpar.PropagConditionVector;
linkpar.StreetWidth = layoutpar.StreetWidth;
linkpar.NumFloors = layoutpar.NumFloors;
linkpar.NumPenetratedFloors = layoutpar.NumPenetratedFloors;
linkpar.Dist1       = layoutpar.Dist1;

% convert actual layoutpar to linkpar
linkpar.ThetaBs = ThetaBs;
linkpar.ThetaMs = ThetaMs;
linkpar.BsHeight = BsHeight;
linkpar.MsHeight = MsHeight;
linkpar.NofSect = NofSect;
linkpar.Pairing = Pairing;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function that maps inputs from (-inf,inf) to (-180,180)
function y=prin_value(x)
y=mod(x,360);
y=y-360*floor(y/180);
