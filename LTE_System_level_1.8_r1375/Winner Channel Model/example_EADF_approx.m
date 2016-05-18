% Example of EADF representation of radiation pattern
%
% Authors: Milan Narandži? (TUI), Martin Käske (TUI)

%% 
NAz=120; %3 degree sampling interval
Az=linspace(-180,180-1/NAz,NAz);

%% Radiation pattern
%slanted dipole
pattern(1,:,1,:)=dipole(Az,12); % slanted by 12 degree
Arrays=AntennaArray('ULA',2,0.01,'FP-ECS',pattern,'Azimuth',Az); %ULA-2 1cm spacing

% isotropic anenna with Horizontal polarization, XPD=Inf
% pattern=complex(cat(2,ones(1,1,1,length(Az)),zeros(1,1,1,length(Az)))); %[NumElem Pols(2) NumEle NumAz]
% Arrays=AntennaArray('FP-ECS',pattern,'Azimuth',Az);

%% Comparison
ar=AntennaResponse(Arrays,Az);

for el=1:1 %first element in array
    for pol=1:2
        subplot(1,2,2*(el-1)+pol)
        fp=squeeze(pattern(1,pol,1,:));
        h(1)=polar(Az'*pi/180,fp);
        set(h(1),'LineWidth',2)
        hold on


        h(2)=polar(Az'*pi/180,abs(squeeze(ar{1}(el,pol,:))),'rx');
        set(h(2),'LineWidth',2)
        legend(h,{'original','EADF approx'})
        xlabel({['antenna element: ',int2str(el)];['polarization: ',int2str(pol)]}) 
        hold off
    end %for pol
end %for el