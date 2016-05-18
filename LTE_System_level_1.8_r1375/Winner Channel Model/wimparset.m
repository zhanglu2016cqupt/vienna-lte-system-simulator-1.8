function wimpar=wimparset(varargin)
%WIMPARSET Model parameter configuration for WIM
%   WIMPAR=WIMPARSET sets default parameters for the input struct WIMPAR 
%   (see WIM). 
%
%   WIMPARSET parameters [ {default} ]:
%
%   Scenario                - WINNER scenario [ {A1} | A2 | B1 | B2 | B3 | B4 | B5a | B5b | B5c | B5f| C1 | C2 | C4 | D1 | D2a]
%   range                   - if Scenario='B5b', the path-loss ranges 1, 2 and 3 are defined, see [1]
%   PropagCondition         - Line of sight condition [ LOS | {NLOS} ]
%   SampleDensity           - number of time samples per half wavelength [ {2} ]
%   NumTimeSamples          - number of time samples [ {100} ]
%   UniformTimeSampling     - Use same time sampling grid for all links [ yes | {no} ] 
%   NumSubPathsPerPath      - number of subpaths per path [ {10} ] (cannot be changed)
%   FixedPdpUsed            - nonrandom path delays and powers [ yes | {no}]
%   FixedAnglesUsed         - nonrandom AoD/AoAs [ yes | {no} ]
%   PolarisedArrays         - usage of dual polarised arrays [ yes | {no} ]
%   TimeEvolution           - usage of time evolution  [ yes | {no} ]
%   CenterFrequency         - carrier frequency in Herz [ {5.25e9} ]
%   DelaySamplingInterval   - delay sampling grid [ {5e-9} ]
%   PathLossModelUsed       - usage of path loss model [ yes | {no} ]
%   ShadowingModelUsed      - usage of shadow fading model [ yes | {no} ]
%   PathLossModel           - path loss model function name [ {pathloss} ]
%   RandomSeed              - sets random seed [ {[empty]} ]
%   UseManualPropCondition  - whether to use manual propagation condition (los/nlos) setting or not. 
%                             If not, the propagation condition is drawn from probabilities.  
%
%   Notes about parameters:
%   - For successful Doppler analysis, one should select SampleDensity > 1.
%     The time sample interval is calculated from CenterFrequency and
%     MsVelocity (see LAYOUT2LINK) according to wavelength/(MsVelocity*SampleDensity).
%     The calculated time sample interval for each link is included in the optional 
%     output argument of WIM. 
%   - If UniformTimeSampling is 'yes' all links will be sampled at
%     simultaneous time instants. In this case, the time sample interval is
%     the same for all links it is calculated by replacing MsVelocity with
%     MAX(MsVelocity), where the maximum is over all links. 
%   - Number of rays is fixed to 20. This is because the AoD/AoAs for
%     subpaths in WIM have fixed angle spread. see [1, Table 4-1].
%   - If FixedPdpUsed='yes', the delays and powers of paths are taken from
%     a table [1, Table 6-1...6-26].  
%   - If FixedAnglesUsed='yes', the AoD/AoAs are taken from a table 
%     [1, Table 6-1...6-26]. Random pairing of AoDs and AoAs is not used.
%   - If PolarisedArrays='yes', single channel coefficient of impulse
%     response turns to 2x2 coefficient matrix, with elements [VV VH;HV HH].
%     Where V stands for vertical polarisation and H for horizontal.
%   - If TimeEvolution='yes', the transition between adjacent channel segments
%     is enabled. Transition from segment to segment is carried out by
%     replacing clusters of the "old" segment by the clusters of the "new"
%     segment, one by one. The route between adjacent channel segments is
%     divided to number of sub-intervals equal to maximum number of
%     clusters within the channel segments. During each sub-interval the
%     power of one old cluster ramps down and one new cluster ramps up.
%     Power ramps are linear. Clusters from the old and new segments are
%     coupled based on their power. If number of clusters is different in
%     the channel segments the weakest clusters are ramped up or down
%     without a pair from other cluster. See [1, Sec 3.4].
%   - CenterFrequency affects path loss and time sampling interval.
%   - DelaySamplingInterval determines the sampling grid in delay domain.
%     All path delays are rounded to the nearest grid point. It can also 
%     be set to zero. 
%   - When PathLossModelUsed is 'no' the path losses are still computed for
%     each link but they are not multiplied into the channel matrices. If
%     ShadowingModelUsed is also 'no', each channel matrix element has unit
%     mean power (summed over delay domain). In other words,
%     MEAN(MEAN(ABS(SUM(H,3)).^2,4),5) is a matrix of (approximately) ones 
%     when isotropic unit-gain antennas are used. Exception: with
%     'polarized' option (and default antennas) the mean power is two.
%   - Path loss model is implemented in a separate function, whose name is
%     defined in PathLossModel. For syntax, see PATHLOSS. 
%   - Even fixing the random seed may not result in fully repeatable
%     simulations due to differences in e.g. MATLAB versions. 
%
%   Ref. [1]: D1.1.2 V1.2, "WINNER II channel models"
%        [2]: 3GPP TR 25.996 v6.1.0 (2003-09)
%
%   See also WIM, LAYOUT2LINK, LAYOUTPARSET, ANTPARSET.

%   Authors: Jari Salo (HUT), Pekka Kyösti (EBIT), Daniela Laselva (EBIT), 
%   Giovanni Del Galdo (TUI), Marko Milojevic (TUI), Christian Schneider (TUI)
%   Lassi Hentilä (EBIT), Mikko Alatossva (CWC/UOULU)

%   Modifications:
%   PathLossOption added                                    4.5.2006 HentLas
%   subScenario removed and combined to Scenario,           15.5.2006 PekKy
%   PropagCondition removed                                 12.2.2007 MikkoA
%   Scenario moved to linkparset -function                  12.2.2007 MikkoA
%   TimeEvolution parameter added, Ansi-C default changed   15.2.2007 PekKy
%   Ansi-C support removed, references updated              15.8.2008 HentLas

if length(varargin)>0
    error('No such functionality yet. Try ''wimpar=wimparset'' instead.')
end

% Set the default values
wimpar=struct(  'range',1,...  
                'end_time',1,...                        % Observation end time for B5 - time points are taken as:  wimpar.TimeVector=linspace(0,wimpar.end_time,T);
                'SampleDensity', 2,...                  % in samples/half-wavelength
                'NumTimeSamples',100,...         
                'UniformTimeSampling','no',... 
                'IntraClusterDsUsed','yes',...          % Two strongest clusters are divided into three subclusters
                'NumSubPathsPerPath',20,...             % only value supported is 20.
                'FixedPdpUsed','no',...                 % Use fixed delays and path powers
                'FixedAnglesUsed','no',...              % Use fixed AoD/AoAs
                'PolarisedArrays','yes',...             % Obsolite - always polarised arrays are used!
                'TimeEvolution','no',...                % use of time evolution option
                'CenterFrequency',5.25e9,...            % in Herz
                'DelaySamplingInterval',5e-9,...        
                'PathLossModelUsed','no',...            
                'ShadowingModelUsed','no',...           
                'PathLossModel','pathloss',...
                'PathLossOption','CR_light',...         % 'CR_light' or 'CR_heavy' or 'RR_light' or 'RR_heavy', CR = Corridor-Room, RR = Room-Room nlos  
                'RandomSeed',[],...                     % if empty, seed is not set. 
                'UseManualPropCondition','yes');
