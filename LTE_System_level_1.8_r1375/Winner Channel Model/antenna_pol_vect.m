% [i_phi, i_theta] = antenna_pol_vect(phi, theta)
%
% antenna_pol returns the polarization vectors of an antenna
% for azimuth/elevation phi/theta
%
% Author: Martin Käske (TUI)
function [i_phi, i_theta]=antenna_pol_vect(phi,theta)
i_theta=[ cos(theta).*cos(phi);
    cos(theta).*sin(phi);
    -sin(theta)];
% don't care for length of vector only direction, same with i_phi
i_theta=i_theta./repmat(sqrt(sum(i_theta.^2)),3,1);

i_phi=[-sin(theta).*sin(phi);
    sin(theta).*cos(phi);
    zeros(1,length(phi))];
% at the poles i_phi is zero, but we know that it can be created by the
% cross product of i_theta and the DoA/DoD
d=sqrt(sum(i_phi.^2));
idx=find(d==0);
if(~isempty(idx))
    i_phi(:,idx)=cross(mysph2cart(phi(idx),theta(idx),ones(1,length(idx))), i_theta(:,idx));
    d(idx)=sqrt(sum(i_phi(:,idx).^2)); %recalculate norm if vector
end;
i_phi=i_phi./repmat(d,3,1);