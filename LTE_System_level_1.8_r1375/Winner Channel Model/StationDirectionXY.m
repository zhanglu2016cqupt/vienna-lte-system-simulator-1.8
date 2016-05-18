% [thetaS1,thetaS2]=StationDirectionXY(S1,S2)
% returns direction angles
% thetaS1 is the direction of S2 as "seen" from S1
%
% Author: Martin Käske (TUI)
function [thetaS1,thetaS2]=StationDirectionXY(S1,S2)
PosS1=[S1.Pos]; PosS2=[S2.Pos];
thetaS1=-atan2(PosS2(2,:)-PosS1(2,:),PosS2(1,:)-PosS1(1,:))+pi/2;
if(nargout > 1)
    thetaS2=pi+thetaS1;
end;