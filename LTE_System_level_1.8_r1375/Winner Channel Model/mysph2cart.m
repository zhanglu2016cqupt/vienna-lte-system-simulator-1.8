% MYSPH2CART transforms spherical to cartesian coordinates the "right" way
% that means theta is the angle from the positive y-axis and phi the angle
% from the positive x-axis
%
% Author: Martin Käske (TUI)
function [varargout]=mysph2cart(phi, theta, r)
    x=r.*sin(theta).*cos(phi);
    y=r.*sin(theta).*sin(phi);
    z=r.*cos(theta);
    if(nargout <= 1)
        varargout{1}=[x;y;z];
    elseif(nargout==3)
        varargout{1}=x;
        varargout{2}=y;
        varargout{3}=z;        
    else
        error(['incorrect number(' num2str(nargout) ') of output arguments']);
    end;