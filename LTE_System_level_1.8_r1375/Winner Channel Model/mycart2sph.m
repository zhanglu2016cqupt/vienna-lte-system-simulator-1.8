% MYCART2SPH() transforms from carthesian
% to spherical coordinates, the "right" way
%
% Author: Martin Käske (TUI)
function [phi, theta, r]=mycart2sph(x,y,z)
    [phi, theta, r]=cart2sph(x,y,z);
    theta=pi/2-theta;