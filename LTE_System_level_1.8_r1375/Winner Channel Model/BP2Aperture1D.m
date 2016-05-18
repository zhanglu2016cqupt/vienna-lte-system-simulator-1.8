% [aperture]=BP2Aperture1D(BP,Az)
% BP [ELNUM POL EL AZ]
% POL denotes vertical polarization first
% Az is a vector containing the azimuth angles (in degrees) the beam pattern is sampled
% at
% ATTENTION: the definition of azimuth used here is the one where 0 degree
% is the positive x-axis and angles increase in mathematical positive
% direction (+90 degree positive y-axis, -90 degree negative y-axis)
%
% Author: Markus Landmann (TUI), Martin Käske (TUI)
function [aperture]=BP2Aperture1D(BP,Az)

NAz=length(Az);
Options.Typ='BP2Aperture';
Options.pol=size(BP,2);
Options.save=0;
Options.positioning_correction=0;
Options.thrs=-320;
Options.darst=0;
%search for optimal NAz over whole range of azimuth samples
Options.NAzh=(2:2:NAz);
Options.Sampleflag=1; %BP is already periodic
Options.NoFreeSpaceCorrection=1;
Options.phase_optimisation=0;
RawData.Az=Az;
RawData.Ele=90*ones(1,NAz); %actually it doesn't matter
% YM are [ELNUM*POL length(Az)] matrices
% first ELNUM rows denote one pol
% second denote the other pol
% if BP contains only one element squeeze would remove this dimension also,
% this works around it
if(size(BP,1)>1)    
    RawData.YM1=[squeeze(BP(:,1,1,:));squeeze(BP(:,2,1,:))];
else
    RawData.YM1=[squeeze(BP(:,1,1,:)).';squeeze(BP(:,2,1,:)).'];
end;
RawData.YM2=RawData.YM1;
RawData.YM3=RawData.YM1;
%disp(['EADF is being calculated ...']);
[arraydef Nopt] = Aperture_Calc(RawData,Options);
if(~isempty(Nopt))    
    Options.NAzh=Nopt;
    arraydef = Aperture_Calc(RawData,Options);
    disp(['Optimal aperture size is Nopt=' num2str(Nopt)]);
end;
aperture=arraydef.apertur;
