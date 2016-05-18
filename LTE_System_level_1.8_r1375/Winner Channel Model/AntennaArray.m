% AA=AntennaArray(varargin)
% see arrayparset.m for examples
% 
% Author: Milan Narandzic (TUI), Martin Käske (TUI)
function AA=AntennaArray(varargin)
% Default array geometry is single antenna (does not require any parameters)
% Default element field-pattern is isotropic, verticaly polarized with
% XPD=Inf;
% if Azimuth is supplied, it should be in degrees
% ATTENTION: the definition of azimuth used here is the one where 0 degree
% is the positive x-axis and angles increase in mathematical positive
% direction (+90 degree positive y-axis, -90 degree negative y-axis)
%         1     2     3     4      5        6         7
Labels={'Pos','Rot','UCA','ULA','FP-ACS','FP-ECS','Azimuth'};
AA = struct('Name', [],...
            'Pos', [0;0;0],...
            'Rot', [0;0;0],...
            'Element', []...
            );

%Analyze inputs
lb=[];
for i=1:length(varargin)
    if ischar(varargin{i});
        tmp=find(strcmp(varargin{i},Labels));
        if ~isempty(tmp)
            lb=[lb; i tmp];
        else
            warning(['ignoring unknown parameter: ' varargin{i}]);
        end
    end
end
lb=[lb;length(varargin)+1 0]; %dummy to make get_arg() work

%first define Array Geometry

arg_i=get_arg(lb,1); %Pos overrides ULA or UCA setings
if ~isempty(arg_i)
    Pos=varargin{arg_i(1)};  %only one argument is expected - take the first one
    N=size(Pos,1);
    for n=1:N
        AA.Element(n).Pos=Pos(n,:)';
    end
    AA.Name=['Custom-' num2str(N)];
    arg_i=get_arg(lb,2); %Rot overrides ULA or UCA
    if ~isempty(arg_i)
        Rot=varargin{arg_i(1)};  %only one argument is expected - take the first one
        rn=size(Rot,1);
        if rn~=length(AA.Element)
            error('AntennaArray: Mismatch in ''Pos'' and ''Rot'' dimensions');
        end
        for n=1:rn
            AA.Element(n).Rot=Rot(n,:)';
        end
    else
        for n=1:length(AA.Element)
            AA.Element(n).Rot=[0;0;0];
        end
    end
else
    arg_i=get_arg(lb,3); %UCA
    if ~isempty(arg_i)
        narg=length(arg_i);
        N=varargin{arg_i(1)};
        if narg>1
            R=varargin{arg_i(2)};
        else
            R=1;
        end;
        disp('Creating UCA');
        for n=1:N
            phi=2*pi/N*(n-1);
            [x,y,z]=sph2cart(phi,0,R);
            AA.Element(n).Pos=[x; y; z];
            AA.Element(n).Rot=[0;0;phi];
        end;
        AA.Name=['UCA-' num2str(N)];
    else
        arg_i=get_arg(lb,4); %ULA
        if ~isempty(arg_i)
            narg=length(arg_i);
            N=varargin{arg_i(1)};
            %disp('creating ULA, NOTE: for testing purpose the elements are not centered around origin of CS');
            disp('Creating ULA');
            if narg>1
                delta_d=varargin{arg_i(2)};
            else
                delta_d=1/N;
            end;
            for n=1:N
                %place elements along x-axis that the center of the array is at [0;0;0]
                AA.Element(n).Pos=[(n-1)*delta_d-(N-1)*delta_d/2; 0; 0];
                %AA.Element(n).Pos=[(n-1)*delta_d; 0; 0];
                AA.Element(n).Rot=[0;0;0];
            end;
            AA.Name=['ULA-' num2str(N)];
        else % use default single antenna
          AA.Element(1).Pos=[0; 0; 0];
          AA.Element(1).Rot=[0; 0; 0];
          AA.Name='Single';
        end
    end
end

        
% Express Field Patterns as EADF

N=length(AA.Element); % Number of elements is defined previously with Antenna Geometry
Azimuth=[];
ECS=0;
 
% Recover azimuth angles
arg_i=get_arg(lb,7); %Azimuth
if ~isempty(arg_i)
    Azimuth=varargin{arg_i(1)};
end

arg_i=get_arg(lb,5); %FP-ACS
if ~isempty(arg_i)
    FP=adjust_FP(varargin{arg_i(1)},N);
    if isempty(Azimuth) && ~isempty(FP)
        No_Azimuth_Angles=size(FP,4); %Azimuth in 4th DIM
        Azimuth=linspace(-180,180-1/No_Azimuth_Angles,No_Azimuth_Angles); % Azimuth angles are starting from -180 deg
    end
    disp(['Calculating EADF in ACS ...']);
    AA.Aperture=BP2Aperture1D(FP,Azimuth); %take first N
else
    arg_i=get_arg(lb,6); %FP-ECS
    if ~isempty(arg_i)
        ECS=1;
        FP=varargin{arg_i(1)}; %only one argument is expected - take the first one
   
        if isempty(Azimuth) && ~isempty(FP)
            No_Azimuth_Angles=size(FP,4); %Azimuth in 4th DIM
            Azimuth=linspace(-180,180-1/No_Azimuth_Angles,No_Azimuth_Angles); % Azimuth angles are starting from -180 deg
        end
        
        nFP=size(FP,1);
        disp(['Calculating EADF in ECS ...']);
        if nFP==1
            AA.CommonAperture=BP2Aperture1D(FP,Azimuth);
        elseif nFP<N
            error('AntennaArray: Field patterns (in ECS) are not provided for all array elements');
        else
            for n=1:N
                AA.Element(n).Aperture=BP2Aperture1D(FP(n,:,:,:),Azimuth);
            end
        end
    else % default: izotropic FP and make CommonAperture;
        disp(['Creating default FP: Vertical-Isotropic, XPD=Inf ...']);
        FP(1,1,1,:) = ones(360,1); %Verical polarization
        FP(1,2,1,:) = zeros(360,1); %XPD=Inf;
        Azimuth=linspace(-180,180-1/360,360);
        AA.Aperture=BP2Aperture1D(adjust_FP(FP,N),Azimuth);
    end
end

if ECS
    disp(['FP rotation ECS->ACS: recalculating EADF in ACS ...']);
    AA=ArrayPreprocess(AA);
end


% sub-functions
% returns indices of 'codes's arguments
function arg_i=get_arg(lb,code)
i=find(lb(:,2)==code); 
if ~isempty(i)
    %lb(i+1,1) will always work since last row of lb is a dummy entry
    %containing length(varargin)+1
    arg_i=(lb(i,1)+1):(lb(i+1,1)-1); % Allocate input arguments
else
    arg_i=[];
end

function FP=adjust_FP(FP,N)
if ~isempty(FP)
    nFP=size(FP,1);
    if nFP==1 && N~=1
        FP=repmat(FP,N,1);
    elseif nFP<N
        error('AntennaArray: Field patterns are not provided for all array elements');
    elseif nFP>N
        FP=FP(1:N,:,:,:);
    end
end