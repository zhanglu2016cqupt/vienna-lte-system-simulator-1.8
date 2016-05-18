function layoutpar=layoutparset(varargin)
%LAYOUTPARSET Link parameter configuration for WIM
%   LAYOUTPAR=LAYOUTPARSET(MsIdx,BsIdx,SectPerBs,K,Arrays) is a struct consisting 
%   of randomly generated network layout parameters. BS and MS positions 
%   are set and a pairing matrix with K links between is generated.
%   Input parameters:
%   Arrays - vector of array definitions, the actual MS and BS
%   arrays are selected from this vector
%   MsIdx - vector of indices into ArrayGeometries
%   BsIdx - vector of indices into ArrayGeometries
%           
%   LAYOUTPAR=LAYOUTPARSET(MsIdx,BsIdx,K,Arrays,RMAX) uses layout range
%   RMAX for generation of MS and BS positions on cartesian co-ordinate
%   system (default: 500 meters).
%
%
%   LAYOUTPAR=LAYOUTPARSET(...,SEED) sets the random seed used in layout
%   parameter generation. 
%
%   The parameters and their defaults are:
%
%   Stations        - vector of struct describing both Ms and Bs, created
%                     by using Arrays
%                     first NofBs elements contains Bs, rest Ms
%   NofSect         - vector of number of sectors in each of the BSs, default=ones(1,NofBs)
%   Pairing         - matrix defining which links are modelled, 2xK
%   ScenarioVector        - maps scenario names to links (see ScenarioMapping.m)
%   PropagConditionVector - maps propagation condition (NLOS=0/LOS=1) to links
%   NumFloors   - For scenarios A2/B4 this determines the floor number of BS/MS 
%   NumPenetratedFloors   - Number of floor between BS/MS for the A1 path loss (default is zero)
%   Dist1       - Distance from BS to "the street crossing" (last LOS point), default NaN -> will be drawn randomly
%                    
%   StreetWidth     - 25 meters
%   Dist2           - NaN default -> will be drawn randomly
%
%   See [1, Fig. 6.1 and 6.2].
%   
%   Some notes about the parameters:
%
%   - Co-ordinates of Bs and Ms should be given in meters with resolution
%     of 1 meter. One meter resolution is assumed in auto-correlation
%     generation of LScorrelation.m function.
%   - Pairing is a matrix with dimensions 2 x K,
%     the first row denotes the Bs index into Stations the second row the Ms index
%     i.e. Stations(Pairing(:,1)) returns the two stations comprising the first link
%   - StreetWidth, this is utilized only with path loss model in [1, sec 5.4.1.2]
%   - Dist2 is defined in [1, Figure 5.16] and generated randomly if empty
%
%   Ref. [1]: D1.1.1 V1.0, "WINNER II interim channel models"
%
%   See also WIM, LAYOUT2LINK, WIMPARSET, ANTPARSET.

%   Authors: Pekka Kyösti (EBIT), Martin Käske (TUI), Milan Narandži? (TUI)
%
%   Updates to Phase II model:  Added new (linkpar) parameters ScenarioVector,
%              PropagConditionVector, NumFloors             (29.5.07  PekKy)


% defaults
NofMs=1;       % number of BMs
NofBs=1;       % number of BSs
K=1;           % number of links
rmax=500;      % layout range [m]
SectPerBs=1;
%BSrmin=10;     % minimum distance of BSs [m]

% inputs
ni=length(varargin);
if ni>0, if (~isempty(varargin{1})), MsAAIdx=varargin{1}; NofMs=length(MsAAIdx); end, end
if ni>1, 
    if (~isempty(varargin{2})), 
        BsAAIdxCell=varargin{2}; 
        NofSect = cellfun(@length,BsAAIdxCell)'; 
        BsAAIdx=cell2mat(cellfun(@(x) x(:),BsAAIdxCell,'UniformOutput',false))';
        NofBs=length(BsAAIdx);  
    end, 
end
if ni>2, if (~isempty(varargin{3})), K=varargin{3}; end, end
if ni>3, if (~isempty(varargin{4})), Arrays=varargin{4}; end, end
if ni>4, if (~isempty(varargin{5})), rmax=varargin{5}; end, end
if ni>5, if (~isempty(varargin{6})), seed=varargin{6}; rand('state',floor(seed)); end, end
if ni>6, error('Too many input arguments!'), end

% check input SectPerBs
% if length(SectPerBs)==1 NofSect=repmat(SectPerBs,1,NofBs); % SectPerBs [scalar]: NofSect equal in any BS
% elseif length(SectPerBs)~=NofBs                            % SectPerBs missing in some BSs
%     SectPerBs = 1; NofSect=repmat(SectPerBs,1,NofBs);        % set default values
%     warning('Number of sectors required for each Bs')
%     disp(['SectPerBs set to default ' mat2str(SectPerBs)])
% else                                                       % SectPerBs [1*NofBs], NofSect differ among BSs
%     NofSect=SectPerBs;
% end

if(max([MsAAIdx(:);BsAAIdx(:) ])>length(Arrays))
    error('MsAAIdx/BsAAIdx out of supported Antenna Array bounds !');
end;

% create stations
for i=1:length(BsAAIdx)
    tmpStation=Arrays(BsAAIdx(i));
    tmpStation.Name=['BS' num2str(i) ' ' tmpStation.Name];
    tmpStation.Pos=[round(rand(2,1)*rmax); 32];
    %tmpStation.Rot=[0;0;(2*rand-1)*pi];
    tmpStation.Velocity=[0;0;0];
    Stations(i)=tmpStation;
end;    
for i=(1:length(MsAAIdx))+length(BsAAIdx)
    tmpStation=Arrays(MsAAIdx(i-length(BsAAIdx)));
    tmpStation.Name=['MS' num2str(i-NofBs) ' ' tmpStation.Name];
    tmpStation.Pos=[round(rand(2,1)*rmax); 1.5];
    %tmpStation.Rot=[0;0;(2*rand-1)*pi];
    tmpStation.Velocity=rand(3,1)-0.5;
    tmpStation.Velocity=tmpStation.Velocity./sqrt(sum(abs(tmpStation.Velocity).^2))*10; % 10 m/s
    Stations(i)=tmpStation;
end;

% outputs
layoutpar=struct('Stations', Stations,...
                 'NofSect', NofSect,...
                 'Pairing', fillpairing(length(MsAAIdx),length(BsAAIdx),K)+[zeros(1,K);length(BsAAIdx)*ones(1,K)],... % Ms are located after Bs in Stations
                 'StreetWidth',20*ones(1,K),...
                 'Dist1',repmat(NaN,1,K),...
         		 'ScenarioVector',1*ones(1,K),... % A1, A2, B1, B2, B3, B4, B5a, B5c, B5f, C1, C2, C3, D1, D2a
                 'PropagConditionVector',round(rand(1,K)),...
                 'NumFloors', 1*ones(1,K),...   % The ground floor is number 1
                 'NumPenetratedFloors',0*ones(1,K)); % Number of floor for the A1 path loss (Default is zero)

function A=fillpairing(NofMs,NofBs,K)
% FILLPAIRING
%   A=FILLPAIRING(NOFMS,NOFBS,K) generates pairing matrix
%   first row denotes the Bs, second row the Ms, the actual link are the
%   columns of A (e.g. for the k'th link: A(:,k)==[BsAAIdx;MsAAIdx])
%   as of now only links from Bs to Ms are supported
%   support for P2P links (Ms-Ms), will be added in future

%   Authors: Pekka Kyösti (EBIT), Daniela Laselva (EBIT), 
%            Martin Käske (TUI)

% avoid duplicate links (same Ms/Bs)
if K>NofBs*NofMs
    K=NofBs*NofMs;
    warning('Number of modelled links limited by the layout')
    disp(['Number of links set to ' mat2str(K)])
end
A=zeros(2,K);
tmp=randperm(NofBs*NofMs);
% link id is created like BsAAIdx*MsAAIdx
% the following two lines extract BsAAIdx and MsAAIdx from link id
A(1,:)=floor((tmp(1:K)-1)/NofMs)+1; % BsAAIdx
A(2,:)=mod(tmp(1:K)-1,NofMs)+1; % MsAAIdx


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function that maps inputs from (-inf,inf) to (-180,180)
function y=prin_value(x)
y=mod(x,360);
y=y-360*floor(y/180);
   