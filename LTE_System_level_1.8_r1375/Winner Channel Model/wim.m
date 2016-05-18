function [H, delays, full_output]=wim(wimpar,layoutpar,initvalues)
%WIM  WINNER Phase II Channel Model (WIM2)
%   H=WIM(WIMPAR,LINKPAR,ANTPAR) is a 5D-array of channel coefficients. For
%   explanation of the input parameter structs, see WIMPARSET, LAYOUT2LINK, LAYOUTPARSET,
%   and ANTPARSET. H is a cell array of size K (number of links.
%   SIZE(H{i})=[U S N T], where U is the number of MS (RX)
%   elements, S is the number of BS (TX) elements, N is the number of paths,
%   T is the number of time samples
%
%   [H DELAYS]=WIM(...) outputs also a [KxN] matrix of path delays (in seconds).
%
%   [H DELAYS BULKPAR]=WIM(...) outputs also the struct BULKPAR, whose fields
%   are as follows:
%
%   With NLOS propagation condition:
%
%   delays          - path delays in seconds [KxN]
%   path_powers     - relative path powers [KxN]
%   aods            - angles of departure in degrees over (-180,180) [KxNxM]
%   aoas            - angles of arrival in degrees over (-180,180) [KxNxM]
%   subpath_phases  - final phases for subpaths in degrees over (0,360) [KxNxM]
%   path_losses     - path losses in linear scale [Kx1]
%   MsBsDistance    - distances between MSs and BSs in meters [1xK]
%   shadow_fading   - shadow fading losses in linear scale [Kx1]
%   delta_t         - time sampling intervals for all links [Kx1]
%
%   In addition, when LOS condition (in addition to the above):
%
%   K_factors       - K factors for all links [Kx1]
%   Phi_LOS         - final phases for LOS paths in degrees over (-180,180) [Kx1]
%
%   [H ...]=WIM(...,INIT_VALUES) uses initial values given in the struct
%   INIT_VALUES, instead of random parameter generation. INIT_VALUES has
%   the same format as BULKPAR, except that SUBPATH_PHASES are now the
%   initial phases. Also, time sampling intervals (delta_t) are not used
%   (they are recalculated for every call of WIM).
%
%   Ref. [1]: D1.1.2 V1.2, "WINNER II channel models"
%        [2]: 3GPP TR 25.996 v6.1.0 (2003-09)
%
%   See also WIMPARSET, LAYOUT2LINK, LAYOUTPARSET and ANTPARSET

%   Authors: Jari Salo (HUT), Giovanni Del Galdo (TUI), Pekka Kyösti (EBIT),
%   Daniela Laselva (EBIT), Marko Milojevic (TUI), Christian Schneider (TUI)
%   Lassi Hentilä (EBIT), Mikko Alatossava (CWC/UOULU), Martin Käske (TUI)


% Note: all units are in degrees, meters, Hertz (1/s) and meters/second (m/s)




ni=nargin;
if (ni<2 || ni>3)
    error('WIM requires two or three input arguments !')
end

linkpar=layout2link(layoutpar);

% Read fixed scenario dependent parameters from a table
fixpar = ScenParTables(linkpar.StreetWidth(1)); %same street width for all links

% WIM parameters, common to all links
SampleDensity=wimpar.SampleDensity;
NumTimeSamples=wimpar.NumTimeSamples;
%N=wimpar.NumPaths;
M=wimpar.NumSubPathsPerPath;
CenterFrequency=wimpar.CenterFrequency;
DelaySamplingInterval=wimpar.DelaySamplingInterval;
PathLossModel=wimpar.PathLossModel;
RandomSeed=wimpar.RandomSeed;
UniformTimeSampling=wimpar.UniformTimeSampling;
PathLossModelUsed=wimpar.PathLossModelUsed;
ShadowingModelUsed=wimpar.ShadowingModelUsed;
FixedPdpUsed = wimpar.FixedPdpUsed;
FixedAnglesUsed = wimpar.FixedAnglesUsed;
PolarisedArrays = wimpar.PolarisedArrays;

% link parameters
ScenarioVector = linkpar.ScenarioVector;
PropagConditionVector = linkpar.PropagConditionVector;
MsBsDistance=linkpar.MsBsDistance;
ThetaBs=linkpar.ThetaBs;
ThetaMs=linkpar.ThetaMs;
MsVelocity=linkpar.MsVelocity;
MsDirection=linkpar.MsDirection;
StreetWidth=linkpar.StreetWidth;
NumFloors = linkpar.NumFloors;
Dist1=linkpar.Dist1;
if(isfield(linkpar,'Stations'))
  Stations=linkpar.Stations;
end;
if(isfield(linkpar,'Pairing'))
  Pairing=linkpar.Pairing;
end;

% extract the number of links
NumLinks=length(MsBsDistance);

% Check that the struct linkpar has the same number of parameters in
% each of its fields. This is also the number of links/users.
if (    NumLinks ~= length(ThetaBs)     ||...
        NumLinks ~= length(ThetaMs)     ||...
        NumLinks ~= length(MsVelocity)  ||...
        NumLinks ~= length(MsDirection)  ||...
        NumLinks ~= length(StreetWidth)  ||...
        NumLinks ~= length(NumFloors)  ||...
        NumLinks ~= length(Dist1))
    
    error('All fields in input struct LINKPAR must be of same size!')
end

% Set random seeds if given
if (isempty(RandomSeed)==0)
    rand('state',RandomSeed);
    randn('state',RandomSeed);
else
    rand('state',sum(100*clock));
    randn('state',sum(101*clock));
end



% GENERATION OF RANDOM "BULK" PARAMETERS FOR ALL LINKS
switch (ni)

    case (2)    % do the basic thing

        % generate bulk parameters for all links
        %bulkpar=generate_bulk_par_polarised(wimpar,linkpar,antpar,fixpar);
        bulkpar=generate_bulk_par(wimpar,linkpar,fixpar);

        % get number of clusters from bulk parameters (located here because
        %  for the case of FixedPdpUsed
        N = size(bulkpar.delays,2);

        % for interpolation
        aods=bulkpar.aods;
        aoas=bulkpar.aoas;


    case (3)    % do not generate random link parameters, use initial values

        % take bulk parameters from input struct
        bulkpar=initvalues;
        
        % This IF is added to remove intra cluster delay spred effects from
        % initial values (spread takes effect in wim_core.m)
        if strcmp(wimpar.IntraClusterDsUsed,'yes')
            for k=1:NumLinks
                % Remove intra cluster delay values from initial values
                tmp = initvalues.delays(k,:);
                tmp([initvalues.IndexOfDividedClust(k,1)+[1:2],initvalues.IndexOfDividedClust(k,2)+2+[1:2]])= [];
                tmpDelay(k,:) = tmp;
                % Remove intra cluster power values from initial values
                tmp = initvalues.path_powers(k,:);
                tmp([initvalues.IndexOfDividedClust(k,1)+[1:2],initvalues.IndexOfDividedClust(k,2)+2+[1:2]])= [];
                tmpPower(k,:) = tmp;
            end
            bulkpar.delays = tmpDelay;
            bulkpar.path_powers = tmpPower;
        end

        % get number of clusters from bulk parameters (located here because
        %  for the case of FixedPdpUsed
        N = size(bulkpar.delays,2);

        % for interpolation
        aods=bulkpar.aods;
        aoas=bulkpar.aoas;


end



% ANTENNA FIELD PATTERN INTERPOLATION
% Interpolation is computationally intensive, so avoid it if possible.
% Since elevation will not be supported, dismiss the elevation dimension (for now)
% NOTE: aods/aoas should be given in degrees.
BsGainIsScalar=0;
MsGainIsScalar=0;

% {NumLinks}[NumElems 3(Pol+Phase) N M]
BsGainPatternInterpolated=AntennaResponse(Stations(Pairing(1,:)), pi/2-aods*pi/180);
MsGainPatternInterpolated=AntennaResponse(Stations(Pairing(2,:)), pi/2-aoas*pi/180);

%% Do antenna field pattern interpolation for the LOS path
% Polarised arrays case added 19.12.2005, PekKy
BsGain_Theta_BS=AntennaResponse(Stations(Pairing(1,:)), pi/2-ThetaBs(:)*pi/180);
MsGain_Theta_MS=AntennaResponse(Stations(Pairing(2,:)), pi/2-ThetaMs(:)*pi/180);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Channel Matrix Generation    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H=cell(NumLinks,1);
if ~strcmpi(wimpar.TimeEvolution,'yes') % NO time evolution

    %%%% Separate processing for LOS and NLOS scenarios %%%%
    %
    if sum(bulkpar.propag_condition)>0      % LOS links

        PCind = find(bulkpar.propag_condition);   % Propagation condition index, LOS

        % CHANNEL MATRIX GENERATION
        [Htmp delta_t FinalPhases FinalPhases_LOS] = wim_core( wimpar,...
            linkpar,...
            bulkpar,...
            BsGainPatternInterpolated,...
            BsGain_Theta_BS,...             % gain of LOS path
            MsGainPatternInterpolated,...
            MsGain_Theta_MS,...             % gain of LOS path
            0,...                           % offset time (not used typically)
            BsGainIsScalar,...
            MsGainIsScalar,...
            PCind);

        H(PCind) = Htmp;

        % final phases
        bulkpar.subpath_phases(PCind,:,:,:)=FinalPhases;
        % time sampling grid
        bulkpar.delta_t(PCind)=delta_t;
    end

    if (length(bulkpar.propag_condition)-sum(bulkpar.propag_condition))>0   % NLOS links

        PCind = find(bulkpar.propag_condition==0);   % Propagation condition index, NLOS

        % CHANNEL MATRIX GENERATION
        [Htmp delta_t FinalPhases] = wim_core( wimpar,...
            linkpar,...
            bulkpar,...
            BsGainPatternInterpolated,...
            BsGain_Theta_BS,...             % gain of LOS path
            MsGainPatternInterpolated,...
            MsGain_Theta_MS,...             % gain of LOS path
            0,...                           % offset time (not used typically)
            BsGainIsScalar,...
            MsGainIsScalar,...
            PCind);

        H(PCind) = Htmp;

        % final phases
        bulkpar.subpath_phases(PCind,:,:,:)=FinalPhases;
        % time sampling grid
        bulkpar.delta_t(PCind)=delta_t;
    end

else    % YES time evolution

    %% Time evolution is not finalized yet...
    error('In current version of WIM2 time evolution is not supported!')

    % Check parameters
    if sum(ScenarioVector(1)~=ScenarioVector)>0
        error('With time evolution option all Scenarios must be same!')
    end
    if sum(PropagConditionVector(1)~=PropagConditionVector)>0
        error('With time evolution option all PropagConditions must be same (LOS/NLOS)!')
    end
    if length(ScenarioVector)<2
        error('With time evolution option number of links must be >1 !')
    end
    if strcmpi(wimpar.FixedPdpUsed,'yes')
        error('With time evolution option FixedPdp is not feasible!')
    end
    if (linkpar.ScenarioVector(1)>=7 & linkpar.ScenarioVector(1)<=9)
        error('With time evolution option B5 scenarios are not feasible!')
    end
    if nargin>3
        warning('Initial values as an input parameter not supported with time evolution option!')
    end

    NofC = size(bulkpar.delays,2);      % Number of clusters
    NofDr = length(ScenarioVector);     % Number of drops (links)
    [tmp,Pind] = sort(bulkpar.path_powers');    % Indices to power sorted clusters
    Pind = flipud(Pind);
    Pind = Pind + (ones(NofC,1)*[1:NofC:NofDr*NofC]-1);

    NofSeg = NofC*NofDr-NofC;   % Number of evolution segments

    for k=1:NofSeg
        [c,r] = ind2sub(size(bulkpar.delays'),Pind(k:k+NofC));  % Note! c,r order is not a mistake
        IndSeg = sub2ind(size(bulkpar.delays),r,c);  % index to clusters of current evolution segment

        bulkpar_evo.delays(k,:)           = bulkpar.delays(IndSeg);
        bulkpar_evo.path_powers(k,:)      = bulkpar.path_powers(IndSeg);
        tmp = permute(bulkpar.aods,[3 1 2]); tmp = tmp(:,:);
        bulkpar_evo.aods(k,:,:)           = permute(tmp(:,IndSeg),[2 1]);
        tmp = permute(bulkpar.aoas,[3 1 2]); tmp = tmp(:,:);
        bulkpar_evo.aoas(k,:,:)           = permute(tmp(:,IndSeg),[2 1]);
        bulkpar_evo.propag_condition      = repmat(bulkpar.propag_condition(1),NofSeg,1);

        % Note! this is not correct
        tmp = permute(bulkpar.subpath_phases,[3 1 2]); tmp = tmp(:,:);
        bulkpar_evo.subpath_phases(k,:,:) = permute(tmp(:,IndSeg),[2 1]);
    end

    % interpolate path loss, K-factor and shadowing to evolution segments
    for k = 1:NofDr-1
        bulkpar_evo.path_losses((k-1)*NofC+1 : k*NofC,1) =...
            linspace(bulkpar.path_losses(k),bulkpar.path_losses(k+1),NofC)';
        bulkpar_evo.Kcluster((k-1)*NofC+1 : k*NofC,1) =...
            linspace(bulkpar.Kcluster(k),bulkpar.Kcluster(k+1),NofC)';
        bulkpar_evo.shadow_fading((k-1)*NofC+1 : k*NofC,1) =...
            linspace(bulkpar.shadow_fading(k),bulkpar.shadow_fading(k+1),NofC)';
        bulkpar_evo.MsBsDistance((k-1)*NofC+1 : k*NofC,1) =...
            linspace(bulkpar.MsBsDistance(k),bulkpar.MsBsDistance(k+1),NofC)';
    end

    % Note, not correct!
    bulkpar_evo.Phi_LOS               = bulkpar.Phi_LOS(1);
    
    %% Call for channel matrix generation
    for k=1:NofSeg
        PCind = ones(1,NofSeg);   % Propagation condition index, all ones
        
        % CHANNEL MATRIX GENERATION
        [Htmp delta_t FinalPhases] = wim_core( wimpar,...
            linkpar,...
            bulkpar,...
            BsGainPatternInterpolated,...
            BsGain_Theta_BS,...             % gain of LOS path
            MsGainPatternInterpolated,...
            MsGain_Theta_MS,...             % gain of LOS path
            0,...                           % offset time (not used typically)
            BsGainIsScalar,...
            MsGainIsScalar,...
            PCind);
        
        H(PCind) = Htmp;
        
        % final phases
        bulkpar.subpath_phases(PCind,:,:,:)=FinalPhases;
        % time sampling grid
        bulkpar.delta_t(PCind)=delta_t;
    end

end     % end IF time evolution



% If path loss and shadowing are to be multiplied into the output
if(strcmpi(PathLossModelUsed,'yes'))
    H=cellfun(@(x,y)(x*y),H,num2cell(sqrt(bulkpar.path_losses)),'UniformOutput',false);
end;
if(strcmpi(ShadowingModelUsed,'yes'))
    H=cellfun(@(x,y)(x*y),H,num2cell(sqrt(bulkpar.shadow_fading)),'UniformOutput',false);
end;


% GENERATE OUTPUT
no=nargout;

if strcmpi(wimpar.IntraClusterDsUsed,'yes')
    bulks = bulkpar;
    bulkpar.delays = repmat(NaN,size(bulks.delays,1),size(bulks.delays,2)+4);
    bulkpar.path_powers = repmat(NaN,size(bulks.delays,1),size(bulks.delays,2)+4);
    for link = 1:NumLinks
        B5ind = find((linkpar.ScenarioVector>=7 & linkpar.ScenarioVector<=9));
        if link==B5ind
            bulkpar.delays(link,1:length(bulks.delays(link,:))) = bulks.delays(link,:);
            bulkpar.path_powers(link,1:length(bulks.delays(link,:))) = bulks.path_powers(link,:);
        else
            P = bulks.path_powers(link,:); P(isnan(P)) = -Inf;
            SortedPower = fliplr(sort(P,2)); P(isinf(P)) = NaN;
            %SubClustInd = P > SortedPower(3); % Index of the cluster to be divided
            % Find index to the clusters to be divided
            [tmp tmpind] = sort(P); SubClustInd = tmpind(end-1:end);
            SubClustInd = zeros(1,size(P,2)); SubClustInd(tmpind(end-1:end)) = 1;
            SubClustDelays = [0  5e-009  10e-009]';
            taus = repmat(bulks.delays(link,:),length(SubClustDelays),1) + repmat(SubClustDelays,1,length(P));
            taus(2:3,~SubClustInd) = NaN;
            taus = reshape(taus,1,length(SubClustDelays)*size(P,2));
            taus(isnan(taus)) = [];
            powers = repmat(bulks.path_powers(link,:),length(SubClustDelays),1);
            SubClustP = [10/20 6/20 4/20]';
            powers(:,find(SubClustInd==1)) = powers(:,find(SubClustInd==1)).*repmat(SubClustP,1,2);
            powers(2:3,~SubClustInd) = NaN;
            powers = reshape(powers,1,length(SubClustP)*size(P,2));
            powers(isnan(powers)) = [];
            IndexOfDividedClust(link,:) = find(SubClustInd==1);


            if (no>1)
                delays(link,1:length(taus)) = taus;
                if (no>2)
                    bulkpar.delays(link,1:length(taus)) = taus;
                    bulkpar.path_powers(link,1:length(taus)) = powers;
                    bulkpar.IndexOfDividedClust = IndexOfDividedClust;
                    bulkpar.aods = aods;
                    bulkpar.aoas = aoas;
                    full_output = bulkpar;
                end
            end
            clear taus powers
        end
    end

else
    if (no>1)
        delays = bulkpar.delays;
        if (no>2)
            if (sum(bulkpar.propag_condition)>0)     % At least one LOS links included
                bulkpar.Phi_LOS=FinalPhases_LOS;
                full_output=bulkpar;
            else    % Only NLOS links included
                full_output=bulkpar;
            end
        end
    end
end











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function that maps inputs from (-inf,inf) to (-180,180)
function y=prin_value(x)
y=mod(x,360);
y=y-360*floor(y/180);