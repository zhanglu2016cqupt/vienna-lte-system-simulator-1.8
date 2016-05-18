% [aperture]=BP2Aperture(BP,Az)
% BP [ELNUM POL EL AZ]
% POL denotes vertical polarization first
% Az is a vector containing the azimuth angles (in degrees) the beam pattern is sampled
% at
%
% Author: Markus Landmann (TUI), Martin Käske (TUI)
function [aperture]=BP2Aperture(BP)
clear RawData
clear Options
Options.Typ='BP2Aperture';
Options.pol=1;
Options.save=0;
Options.positioning_correction=0;
Options.thrs=-320;
Options.darst=0;
%Options.NAzh=(2:2:60);
%Options.NAzh=42;
%Options.NEleh=(2:2:60);
Options.Sampleflag=0;
Options.NoFreeSpaceCorrection=1;
Options.phase_optimisation=0;
NAz=260;
NEle=200;
RawData.Az=(0:NAz-1)*360/NAz-180;
RawData.Ele=(0:NEle)*180/NEle;
%RawData.Ele=90*ones(1,NAz);
Az=RawData.Az;
Ele=RawData.Ele;
% RawData.Ele=90*ones(1,NAz);
% RawData.YM1=ones(2,NAz);
%RawData.YM1=[cos(RawData.Az/180*pi);sin(RawData.Az/180*pi)].*[exp(j*2*pi*0.5*cos(RawData.Az/180*pi));exp(j*2*pi*0.5*cos(RawData.Az/180*pi))];
%RawData.YM1=[cos(RawData.Az/180*pi);sin(RawData.Az/180*pi)].*[exp(j*2*pi*0.5*cos(RawData.Az/180*pi));exp(j*2*pi*0.5*cos(RawData.Az/180*pi))];
%BP1=cos(RawData.Az/180*pi).'*sin(RawData.Ele/180*pi).*exp(j*2*pi*0.5*cos(RawData.Az/180*pi).'*ones(size(RawData.Ele)));
BP1=cos(10*RawData.Az/180*pi).'*sin(10*RawData.Ele/180*pi);
RawData.YM1=[BP1(:).';BP1(:).'];
RawData.YM2=RawData.YM1;
RawData.YM3=RawData.YM1;
RawData.Az=repmat(RawData.Az.',NAz,1).';
RawData.Ele=(RawData.Ele.'*ones(1,NAz)).';
RawData.Ele=RawData.Ele(:).';
% load('X:\makae\src\PHK2_UCA16_2D_20-Jul-2006_p_25-Jul-2006.mat');
%Options.NAzh = [120*ones(size(40:2:120)) (4:2:30)*4];  % only even numbers + number/4 is also even
%Options.NAzh=[120];
Options.NAzh=NAz;
%Options.NEleh = [40:2:120 120*ones(size((4:2:30)*4))]; % only even numbers
Options.NEleh=NEle;
Options.indele=1:NEle;
Options.indaz=1:NAz;
[arraydef] = Aperture_Calc(RawData,Options);
aperture=arraydef.apertur;