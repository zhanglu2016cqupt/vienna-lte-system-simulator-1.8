% g=AntennaResponse(arrays, phi, theta)
% arguments: array, phi, theta, lambda
function g=AntennaResponse(arrays, phi,varargin)
sz_phi=size(phi);

if(nargin <= 2)
    theta=ones(sz_phi).*pi/2;
else
    theta=varargin{1};
end;

if(sz_phi~=size(theta))
    error('phi and theta must have the same size');
end;
if(length(arrays) ~= sz_phi(1))
    error('first dimension of phi/theta must match length of arrays');
end;

g=cell(1,length(arrays));
for s=1:length(arrays)
    array=arrays(s);
    % make phi/theta for this station col vector, reshape later
    phi_s=phi(s,:);
    theta_s=theta(s,:);
    % pol of the wave in GCS
    [i_phi, i_theta]=antenna_pol_vect(phi_s,theta_s);

    % DoA in GCS
    %[wave_dir_gc(1,1), wave_dir_gc(2,1), wave_dir_gc(3,1)]=mysph2cart(phi, theta, 1);
    wave_dir_gc=mysph2cart(phi_s, theta_s, ones(1,length(phi_s)));
    % DoA in ACS
    wave_dir_ac=unrotate_vector(wave_dir_gc, array.Rot);

    NumElem=length(array.Element);
    g{s}=complex(zeros(NumElem,3, length(phi_s)));
    fp=complex(zeros(length(phi_s)*NumElem,2));
    fp_gc=complex(zeros(length(phi_s)*NumElem,2));
    i_phi_ac=zeros(3,length(phi_s)*NumElem);
    i_theta_ac=zeros(3,length(phi_s)*NumElem);
    if(isfield(array,'Aperture')) %preprocessing was used
        %disp('with preprocessing');
        [phi_ac, theta_ac, r]=mycart2sph(wave_dir_ac(1,:),wave_dir_ac(2,:),wave_dir_ac(3,:));
        if(~isempty(array.Aperture))
            fptmp=interpbp(array.Aperture, phi_ac.', 0).';
            fp=reshape(fptmp(:), NumElem*length(phi_ac),2);
        else
            fp(:,:)=ones(NumElem*length(phi_ac),2);
        end;
        % fp_ac contains now FP for each element and DoA (pols are VH)
        % fp_ac(1,1:2) is FP for first element/first DoA
        % fp_ac(2,1:2) is FP for first element/second DoA ...

        % pol vectors for each DoA in ACS
        [i_phi_ac,i_theta_ac]=antenna_pol_vect(phi_ac, theta_ac);
        for i=1:NumElem
            g{s}(i,3,:)=calc_dist(array.Element(i).Pos, wave_dir_ac);
        end
        % i_phi_ac is the same for all elements
        i_phi_ac=repmat(i_phi_ac,1,NumElem);
        i_theta_ac=repmat(i_theta_ac,1,NumElem);
    else % no preprocessing
        %disp('without preprocessing');
        elements=array.Element;
        for i=1:NumElem
            %precalculate idx since its used often
            idx=(i-1)*length(phi_s)+1:i*length(phi_s);
            % DoA in ECS
            wave_dir_ec=unrotate_vector(wave_dir_ac, elements(i).Rot);
            [phi_ec, theta_ec, r]=mycart2sph(wave_dir_ec(1,:), wave_dir_ec(2,:), wave_dir_ec(3,:));
            % 1. phi 2. theta
            if(~isfield(elements(i),'Aperture') || isempty(elements(i).Aperture))
                [fp(idx,1), fp(idx,2)]=InterpHelper(array.CommonAperture, phi_ec, theta_ec);
            else
                [fp(idx,1), fp(idx,2)]=InterpHelper(elements(i).Aperture, phi_ec, theta_ec);
            end;
            % pol in ECS
            [i_phi_ec,i_theta_ec]=antenna_pol_vect(phi_ec, theta_ec);
            % pol in ACS
            i_phi_ac(:,idx)=rotate_vector(i_phi_ec, elements(i).Rot);
            i_theta_ac(:,idx)=rotate_vector(i_theta_ec, elements(i).Rot);
            g{s}(i,3,:)=calc_dist(elements(i).Pos, wave_dir_ac);
        end;
        %now we've reached a state that looks like preprocessing was used
    end;
    % pol in GCS
    i_phi_gc=rotate_vector(i_phi_ac, array.Rot);
    i_theta_gc=rotate_vector(i_theta_ac, array.Rot);
    % make i_phi/i_theta compatible (wrt. to dimensions) with
    % i_phi_ac/i_theta_ac
    i_phi=repmat(i_phi,1,NumElem);
    i_theta=repmat(i_theta,1,NumElem);
    % project rotated pol-vectors of antenna onto pol-vectors of wave
    % calculation works like this (<,> denotes scalar product):
    % fp_gc_phi  =<i_phi_gc,i_phi>  *fp_ac_phi+<i_theta_gc,i_phi>  *fp_ac_theta
    % fp_gc_theta=<i_phi_gc,i_theta>*fp_ac_phi+<i_theta_gc,i_theta>*fp_ac_theta
    % since i_phi/i_theta are matrices the scalar product cannot be
    % computed by i_phi'*i_phi_ac instead sum() and element-wise product is
    % used
    %i_phi_proj=[sum(i_phi.*i_phi_gc).' sum(i_phi.*i_theta_gc).'];
    %i_theta_proj=[sum(i_theta.*i_phi_gc).' sum(i_theta.*i_theta_gc).'];
    i_phi_proj=[sum(i_phi.*i_theta_gc).' sum(i_phi.*i_phi_gc).'];
    i_theta_proj=[sum(i_theta.*i_theta_gc).' sum(i_theta.*i_phi_gc).'];
    % fp_gc contains V/H for all elements and all DoAs at once
    fp_gc(:,2)=sum(i_phi_proj.*fp,2);
    fp_gc(:,1)=sum(i_theta_proj.*fp,2);

    % reshape
    %g(i,:,:)=[fp_gc(2);fp_gc(1);phase_term];
    g{s}(:,1,:)=reshape(fp_gc(:,1), length(phi_s), NumElem).';
    g{s}(:,2,:)=reshape(fp_gc(:,2), length(phi_s), NumElem).';

    % reshape to fit size of input phi/theta
    g{s}=reshape(g{s}(:,:,:),[NumElem 3 sz_phi(2:end)]);
end;

% little helper function to change layout of interpbp's output
% g_theta is vertical pol
% g_phi is horizontal pol
function [g_theta, g_phi]=InterpHelper(aperture, phi, theta)
%     [g_theta, g_phi]=my_fieldpattern(phi, theta, pi/7);
%    g_theta=1;
%    g_phi=0;
g=interpbp(aperture, phi.',theta).';% NumElem*2
g=reshape(g(:), aperture.elements*length(phi),2);
g_theta=g(:,1);
g_phi=g(:,2);

% calc_dist() calculates the projection
% of a to b where b is normalized to have length==1
function d=calc_dist(a, b)
% b has to have length==1
b=b./repmat(sqrt(sum(b.^2)),3,1);
d=a.'*b;