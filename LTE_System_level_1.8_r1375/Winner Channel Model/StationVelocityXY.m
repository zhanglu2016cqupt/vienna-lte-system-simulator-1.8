% [dir, mag] = StationVelocityXY(S)
% returns magnitude and direction 
%
% Author: Martin Käske (TUI)
function [dir, mag]=StationVelocityXY(S)
v=[S.Velocity];
mag=sqrt(sum(v.^2));
if(nargout > 1)
    dir=-atan2(v(2,:),v(1,:))+pi/2;
end;
