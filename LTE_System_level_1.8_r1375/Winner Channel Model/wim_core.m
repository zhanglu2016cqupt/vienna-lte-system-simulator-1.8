%WIM_CORE Channel coefficient computation for a geometric channel model
%   [H DELTA_T FINAL_PHASES FINAL_PHASES_LOS]=WIM_CORE(WIMPAR,LINKPAR,BULKPAR,BSGAIN,BSGAIN_LOS,MSGAIN,MSGAIN_LOS,OFFSET_TIME, BSGAINISSCALAR, MSGAINISSCALAR, PCind)
%   This is the wim_core aka the big for loop. It implements the formula in [2, Eq. 4.14, 4.17, 4.19].
%
%   Outputs:
%
%   H               - cell array of size K
%   H{i}            - [UxSxNxT] array of channel coefficients
%   DELTA_T         - time sampling intervals (in seconds) for all links
%   FINAL_PHASES    - final phases of all subpaths in degrees over (-180,180)
%   FINAL_PHASES_LOS- final phases for LOS paths in degrees over (-180,180)
%
%   Inputs:
%
%   WIMPAR          - input struct, see WIMPARSET
%   LINKPAR         - input struct, see LAYOUT2LINK and LAYOUTPARSET
%   ANTPAR          - input struct, see ANTPARSET
%   BULKPAR         - input BULKPAR, see GENERATE_BULK_PAR
%   BSGAIN          - {K}[SxNxM] cell array of interpolated antenna field
%                     patterns (complex)
%   BSGAIN_LOS      - {K}[S] cell array of interpolated antenna field patterns
%                     (complex) for LOS paths. Only used with the LOS
%                     option; it is set to scalar otherwise.
%   MSGAIN          - {K}[UxNxM] cell array of interpolated antenna field
%                     patterns (complex)
%   MSGAIN_LOS      - {K}[U] cell array of interpolated antenna field patterns
%                     (complex) for LOS paths. Only used with the LOS
%                     option; it is set to scalar otherwise.
%   OFFSET_TIME     - time offset added to the initial phase (set to zero by default)
%   BSGAINISSCALAR  - this is 1 if BsGain is uniform over azimuth, 0 otherwise.
%   MSGAINISSCALAR  - this is 1 if MsGain is uniform over azimuth, 0 otherwise.
%
%   With 'polarized' option:
%
%   BSGAIN          - [KxSx2xNxM] array of interpolated antenna field
%                     patterns (complex), where the third dimension are
%                     the patterns for [V H] polarizations.
%   MSGAIN          - [KxUx2xNxM] array of interpolated antenna field
%                     patterns (complex), where the third dimension are
%                     the patterns for [V H] polarizations.
%
%
%   Ref. [1]: 3GPP TR 25.996 v6.1.0 (2003-09)
%        [2]: D1.1.2 V1.2, "WINNER II channel models"
%
%   Authors: Giovanni Del Galdo (TUI), Marko Milojevic (TUI), Jussi Salmi (HUT),
%   Christian Schneider (TUI), Jari Salo (HUT), Pekka Kyösti (EBIT),
%   Daniela Laselva (EBIT), Lassi Hentilä (EBIT), Martin Käske (TUI)

% Bug fixes:
%   Erroneus variable name s changed to u on line 573. Caused error with
%    asymmetric MIMO configuration (e.g.2x4) & CDL models or B5 scenario. (22.8.2006 PekKy)
%   Polarised arrays, non-ANSI-C version, power scaling changed harmonised
%    with ANSI-C version on rows 419-420, 462-468. (22.8.2006 PekKy)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                   --------                     %%
function [H, delta_t, output_SubPathPhases, output_Phi_LOS] = wim_core (wimpar,linkpar,bulkpar,BsGain,BsGain_Theta_BS,MsGain,MsGain_Theta_MS,offset_time, BsGainIsScalar, MsGainIsScalar,PCind)
%%                   --------                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% offset_time [samples] = defines the starting point (in samples) of the
%                         time axis
% Examples: you want to calculate 1000 time samples calling the wim_core
%           twice (everytime for 500 timesamples)
%           H1 = wim_core (..., 0)
%           H2 = wim_core (..., 500)
%



DEBUG_MODE_FLAG   = 0;
PROFILE_MODE_FLAG = 0;
DISPLAY_MODE_FLAG = 0;

N    = size(bulkpar.delays,2);        % number of paths
T    = wimpar.NumTimeSamples;         % number of time samples
K    = length(PCind);                 % number of links
M    = wimpar.NumSubPathsPerPath;     % number of subpaths

% AnsiC_core = wimpar.AnsiC_core;


% intra-cluster delays spread added based on the D1.1.1
if strcmpi(wimpar.IntraClusterDsUsed,'yes')
    AnsiC_core = 'no';      % NOTE! IntraClusterDsUsed forces AnsiC_core = 'no', change this after scm_mex_core.c is modified

    NumRaysPerSubCluster = [10,6,4];
    RayOrder = [1,2,3,4,5,6,7,8,19,20,9,10,11,12,17,18,13,14,15,16];
    P = bulkpar.path_powers; P(isnan(P))=-Inf;
    SortedPower = fliplr(sort(P,2)); P(isinf(P))=NaN;
    for xx = 1:size(P,1)
        SubClusterInd(xx,:) = P(xx,:) > SortedPower(xx,3); % Index of the cluster to be divided
    end

else    % special case, one midpath only
    LM = 1:M;
    LN=M;
    L = length(LN);
    P = bulkpar.path_powers;
    SubClusterInd = zeros(size(P));
end

% initialize H and determine number of array elements for each link
H = cell(1,K);
U=zeros(1,K); % number of Rx(Ms) elements
S=zeros(1,K); % number of Tx(Bs) elements
for k=1:K
    U(k)=size(MsGain{PCind(k)},1); S(k)=size(BsGain{PCind(k)},1);
    H{k}=zeros(U(k),S(k),N+4,T); %"+4" is due to 4 extra sub-clusters
end;

% Set internal parameters
speed_of_light=2.99792458e8;
wavelength=speed_of_light/wimpar.CenterFrequency;

% dummy
output_Phi_LOS       = zeros(K,1);

% let's make the time axis - for that we need to check UniformTimeSampling
% and the MSs' velocities
% Note: SampleDensity is samples per half wavelength.
if strcmp(wimpar.UniformTimeSampling,'yes')

    max_vel = max(linkpar.MsVelocity);
    delta_t = repmat((wavelength / max_vel)/2/wimpar.SampleDensity,K,1);

else % 'UniformTimeSampling' is 'no'

    delta_t = (wavelength ./ linkpar.MsVelocity(PCind).')./2/wimpar.SampleDensity ;

end
t = repmat(delta_t,1,T).*repmat([0:T-1]+offset_time,K,1); % matrix containing the time axes for all links [KxT]

% Time axis generation for fixed feeder links (B5 scenarios)
tmp = zeros(1,length(linkpar.ScenarioVector));
tmp(PCind) = 1;
B5ind = find((linkpar.ScenarioVector>=7 & linkpar.ScenarioVector<=9).*tmp);

if length(B5ind)>0
    H = zeros(U,S,N+4,T,K);  %"+4" is due to 4 extra sub-clusters (even though not used with B5)
    SubClusterInd(B5ind,:) = zeros(length(B5ind),size(bulkpar.delays,2));
    AnsiC_core = 'no';      % NOTE! B5 forces AnsiC_core = 'no', change this after scm_mex_core.c is modified
    for k=1:length(B5ind)
        B5ind2(k) = find((B5ind(k)==PCind));
    end
    % not final
    wimpar.TimeVector=linspace(0,wimpar.end_time,T);
    % not final

    t(B5ind2,:) = repmat(wimpar.TimeVector,length(B5ind),1);%KTH
    linkpar.MsVelocity(B5ind)=zeros(length(B5ind),1);
    delta_t(B5ind2)=repmat(wimpar.TimeVector(2)-wimpar.TimeVector(1),length(B5ind),1); %% Dummy value
end

k_CONST = 2*pi/wavelength;      % wave number


if ~strcmpi(wimpar.PolarisedArrays,'yes')

    if DISPLAY_MODE_FLAG
        disp('entering main loop...');
    end

    output_SubPathPhases = zeros(K,N,M);

    for kk = 1:K % cycles links
        k = PCind(kk);
        for u = 1:U(kk) % cycles (MS) antennas

            for s = 1:S(kk) % cycles Tx (BS) atennas

                LN_index = 0; %
                for n = 1:N % cycles paths

                    LM_index = 0; %

                    if SubClusterInd(kk,n)==1; L=3; LM = RayOrder; LN = NumRaysPerSubCluster;
                        path_powers = [P(kk,n)*10/20 P(kk,n)*6/20 P(kk,n)*4/20];
                    else L=1; LN=M; LM=1:M;
                        path_powers = P(kk,n);
                    end

                    for km = 1:L % cycles subclusters

                        LN_index = LN_index+1; %Running index of the clusters (including sub-clusters)

                        temp = zeros(M,T);   % oversized, just to keep it always the same size

                        for m=1:LN(km) % cycles rays

                            LM_index = LM_index+1; % Running index of the rays

                            % Calculate Doppler frequency nu of scatterer m
                            if sum(k==B5ind)    % IF current link is B5
                                nu = bulkpar.scatterer_freq(k,n,m);
                            else    % IF not B5
                                nu = (linkpar.MsVelocity(k) * cos((bulkpar.aoas(k,n,LM(LM_index)) - linkpar.MsDirection(k))*pi/180))/wavelength;
                            end
                            temp(m,:) =  BsGain{k}(s,1,n,LM(LM_index)) *...
                                exp(j*(...
                                k_CONST * BsGain{k}(s,3,n,LM(LM_index)) +...
                                (bulkpar.subpath_phases(k,n,LM(LM_index))*pi/180)+...
                                k_CONST * MsGain{k}(u,3,n,LM(LM_index)) ...
                                )) *...
                                MsGain{k}(u,1,n,LM(LM_index)) * exp(1j*2*pi*nu * t(kk,:) );


                        end % rays

                        H{kk}(u,s,LN_index,:) =  sqrt(path_powers(km) / LN(km)) * sum(temp,1);

                    end % subclusters

                end % paths

            end % Tx antennas
        end % Rx antennas
    end % links


    for kk = 1:K % cycles links  % Of course the for loop could be avoided.
        k = PCind(kk);
        for n = 1:N % cycles paths
            for m=1:M % cycles supaths

                % SIMPLE LOOP
                output_SubPathPhases(kk,n,m) =  k_CONST * linkpar.MsVelocity(k) * cos((bulkpar.aoas(k,n,m) - linkpar.MsDirection(k))*pi/180) * (delta_t(kk)*T);

            end % subpaths
        end % paths
    end % links

    output_SubPathPhases = prin_value((output_SubPathPhases*180/pi + bulkpar.subpath_phases(PCind,:,:)));



else % it's polarized!
    if DISPLAY_MODE_FLAG
        disp('entering polarized option...');
    end

    output_SubPathPhases = zeros(K,4,N,M);

    % Set polarisation matrix powers according to XPRs
    % Assume Pvv+Pvh=Phh+Phv=1. In this case
    %         Pvv = bulkpar.xprV./(1+bulkpar.xprV);
    %         Pvh = 1-bulkpar.xprV./(1+bulkpar.xprV);
    %         Phh = bulkpar.xprH./(1+bulkpar.xprH);
    %         Phv = 1-bulkpar.xprH./(1+bulkpar.xprH);
    % Pxy replaced by R_ni (8.5.2006 PekKy)
    r_n1 = 1 ./ bulkpar.xpr;
    r_n2 = r_n1; % D1.1.2 definition
    %r_n2 = 1 ./ bulkpar.xprH;


    % BsGain has must have size: [K x S x 2 x N x M]
    % the first dimension in the polarization must be vertical
    %
    % MsGain has must have size: [K x U x 2 x N x M]
    % the first dimension in the polarization must be vertical
    %
    % subpath_phases has size: [K x 4 x N x M]
    % bulkpar.xpd has size: [K,2,N]
    temp = zeros(M,T);
    for kk = 1:K % cycles links
        k = PCind(kk);
        for u = 1:U(kk) % cycles (MS) antennas

            for s = 1:S(kk) % cycles Tx (BS) atennas

                LN_index = 0; %
                for n = 1:N % cycles paths

                    LM_index = 0; %

                    if SubClusterInd(kk,n)==1; L=3; LM = RayOrder; LN = NumRaysPerSubCluster;
                        path_powers = [P(kk,n)*10/20 P(kk,n)*6/20 P(kk,n)*4/20];
                    else L=1; LN=M; LM=1:M;
                        path_powers = P(kk,n);
                    end

                    for km = 1:L % cycles midpaths

                        LN_index = LN_index+1; %Running index of the clusters (including sub-clusters)

                        temp = zeros(M,T);   % oversized, just to keep it always the same size

                        for m=1:LN(km) % cycles subpaths

                            LM_index = LM_index+1;

                            %                                 % Assume Pvv+Pvh=Phh+Phv=1. In this case        % Commented 22.8.2006, PekKy
                            %                                 % Pvv=xprV/(1+xprV) and Pvh=1-xprV/(1+xprV),
                            %                                 % Phh=xprH/(1+xprH) and Phv=1-xprH/(1+xprH)
                            %                                 temp(m,:) =  [BsGain(k,s,1,n,LM(LM_index)) BsGain(k,s,2,n,LM(LM_index))] * ...
                            %                                     [ sqrt(Pvv(k,n,m)) * exp(j*bulkpar.subpath_phases(k,1,n,LM(LM_index))*pi/180),   sqrt(Pvh(k,n,m)) * exp(j*bulkpar.subpath_phases(k,2,n,LM(LM_index))*pi/180);...
                            %                                       sqrt(Phv(k,n,m)) * exp(j*bulkpar.subpath_phases(k,3,n,LM(LM_index))*pi/180),   sqrt(Phh(k,n,m)) * exp(j*bulkpar.subpath_phases(k,4,n,LM(LM_index))*pi/180)] *...
                            %                                     [MsGain(k,u,1,n,LM(LM_index)); MsGain(k,u,2,n,LM(LM_index))] * ...
                            %                                     exp(j * (k_CONST * ds * sin((bulkpar.aods(k,n,LM(LM_index)))*pi/180))) * ...
                            %                                     exp(j * (k_CONST * du * sin((bulkpar.aoas(k,n,LM(LM_index)))*pi/180))) * ...
                            %                                     exp(j * k_CONST * linkpar.MsVelocity(k) * cos((bulkpar.aoas(k,n,LM(LM_index)) - linkpar.MsDirection(k))*pi/180) * t(k,:));

                            % Calculate Doppler frequency nu of scatterer m
                            if sum(k==B5ind)    % IF current link is B5
                                nu = bulkpar.scatterer_freq(k,n,m);
                            else    % IF not B5
                                nu = (linkpar.MsVelocity(k) * cos((bulkpar.aoas(k,n,LM(LM_index)) - linkpar.MsDirection(k))*pi/180))/wavelength;
                            end

                            temp(m,:) =  [BsGain{k}(s,1,n,LM(LM_index)) BsGain{k}(s,2,n,LM(LM_index))] * ...
                                [ exp(j*bulkpar.subpath_phases(k,1,n,LM(LM_index))*pi/180)  ,  sqrt(r_n1(k,n,m)) * exp(j*bulkpar.subpath_phases(k,2,n,LM(LM_index))*pi/180);...
                                sqrt(r_n2(k,n,m)) * exp(j*bulkpar.subpath_phases(k,3,n,LM(LM_index))*pi/180)  ,  exp(j*bulkpar.subpath_phases(k,4,n,LM(LM_index))*pi/180)] *...
                                [MsGain{k}(u,1,n,LM(LM_index)); MsGain{k}(u,2,n,LM(LM_index))] * ...
                                exp(j * (k_CONST * BsGain{k}(s,3,n,LM(LM_index)))) * ...
                                exp(j * (k_CONST * MsGain{k}(u,3,n,LM(LM_index)))) *...
                                exp(1j*2*pi*nu * t(kk,:) );


                        end % rays

                        H{kk}(u,s,LN_index,:) =  sqrt(path_powers(km) / LN(km)) * sum(temp,1);

                    end % subclusters
                end % paths
            end % Tx antennas
        end % Rx antennas
    end % links

    %         for k = 1:K % cycles links  % Of course the for loop could be avoided.
    %             for n = 1:N % cycles paths
    %                 for m=1:M % cycles supaths
    %
    %
    %                     output_SubPathPhases(k,:,n,m) =  (k_CONST * linkpar.MsVelocity(k) * cos((bulkpar.aoas(k,n,m) - linkpar.MsDirection(k))*pi/180) * (t(k,end)+delta_t(k))) * ones(1,4);
    %
    %                 end % subpaths
    %             end % paths
    %         end % links
    %
    %         output_SubPathPhases = prin_value((output_SubPathPhases*180/pi + bulkpar.subpath_phases));


    for kk = 1:K % cycles links  % Of course the for loop could be avoided.
        k = PCind(kk);
        for n = 1:N % cycles paths
            for m=1:M % cycles supaths

                % SIMPLE LOOP
                output_SubPathPhases(kk,:,n,m) =  (k_CONST * linkpar.MsVelocity(k) * cos((bulkpar.aoas(k,n,m) - linkpar.MsDirection(k))*pi/180) * (delta_t(kk)*T)) * ones(1,4);

            end % subpaths
        end % paths
    end % links

    output_SubPathPhases = prin_value((output_SubPathPhases*180/pi + bulkpar.subpath_phases(PCind,:,:,:)));



end % is it polarized?




%%% LOS OPTION %%%

% index to LOS but not B5 links
LosNonB5ind = find(bulkpar.propag_condition'.*(linkpar.ScenarioVector<7 | linkpar.ScenarioVector>9));

% If LOS, but not B5 link and not 'FixedPDP'
if (bulkpar.propag_condition(PCind(1))==1) & length(LosNonB5ind)>0

    % Take the values of K factors
    K_factors       = bulkpar.Kcluster;%;K_factors;

    ThetaBs      = linkpar.ThetaBs; ThetaBs=ThetaBs(:).';
    ThetaMs      = linkpar.ThetaMs; ThetaMs=ThetaMs(:).';

    output_Phi_LOS       = zeros(length(linkpar.ScenarioVector),1);

    for kk = 1:length(LosNonB5ind) % cycles links
        k = LosNonB5ind(kk);
        k_ind = find(k==PCind);     % index to current LOS/NLOS links which are LOS but not B5
        output_Phi_LOS(k,1) = k_CONST * linkpar.MsVelocity(k) * cos((ThetaMs(k) - linkpar.MsDirection(k))*pi/180) * (t(k_ind,end)+delta_t(k_ind));
        for u = 1:U(kk) % cycles (MS) antennas

            for s = 1:S(kk) % cycles (BS) antennas

                temp =  BsGain_Theta_BS{k}(s,1) * exp(j * k_CONST * BsGain_Theta_BS{k}(s,3)).* ...
                    MsGain_Theta_MS{k}(u,1) * exp(j * (k_CONST * MsGain_Theta_MS{k}(u,3) + bulkpar.Phi_LOS(k,1) * pi/180 )) * ...
                    exp(j * k_CONST * linkpar.MsVelocity(k) * cos((ThetaMs(k) - linkpar.MsDirection(k))*pi/180) * t(find(k==PCind),:));


                H{k_ind}(u,s,1,:)= (sqrt(1/(K_factors(k)+1)) * squeeze(H{k_ind}(u,s,1,:)) + sqrt(K_factors(k)/(K_factors(k)+1)) * temp.').';
            end % Rx antennas

        end % Tx antennas

        H{k_ind}(:,:,2:end,:)= sqrt(1/(K_factors(k)+1)) *  H{k_ind}(:,:,2:end,:);
    end % links

    output_Phi_LOS = prin_value((output_Phi_LOS*180/pi + bulkpar.Phi_LOS));
end % if 'LOS' propagation condition


%% B5 links %%%
%   direct rays are added and cluster powers adjusted
%   according to cluster-wise K-factors given in [2, tables 7.17-24]
%
% index to links which are LOS and B5
LosB5ind = find(bulkpar.propag_condition'.*(linkpar.ScenarioVector>=7 & linkpar.ScenarioVector<=9));

if (bulkpar.propag_condition(PCind(1))==1) & length(LosB5ind)>0

    Kcluster = bulkpar.Kcluster;             % read cluster-wise K-factors

    %    output_Phi_LOS       = zeros(K,1);

    for kk = 1:length(LosB5ind) % cycles links
        k = LosB5ind(kk);
        k_ind = find(k==PCind);     % index to current LOS/NLOS links which are LOS and B5
        for u = 1:U % cycles (MS) antennas
            du = antpar.MsElementPosition(u)*wavelength;
            for s = 1:S % cycles (BS) atennas

                n = 1;     % index to cluster with a direct ray, in WIM2 always 1st cluster
                ds = antpar.BsElementPosition(s)*wavelength;

                aod_direct = bulkpar.aods(k,n,1)-bulkpar.aods(k,n,2);     % AoD for the direct ray (middle)
                aoa_direct = bulkpar.aoas(k,n,1)-bulkpar.aoas(k,n,2);     % AoA for the direct ray (middle)

                % antenna gain of direct ray is approximated by linear interpolation
                BsGain_direct = mean(BsGain(k,s,n,1:2));
                MsGain_direct = mean(MsGain(k,u,n,1:2));        % 22.8.2006 PekKy, index s corrected to u

                nu = 0;     % LOS ray has always 0 Hz Doppler in B5 scenarios

                temp =  BsGain_direct * exp(j * k_CONST * ds * sin( aod_direct*pi/180))* ...
                    MsGain_direct * exp(j * (k_CONST * du * sin( aoa_direct*pi/180  ) + bulkpar.Phi_LOS(k,1) * pi/180 )) * ...
                    exp(1j*2*pi*nu * t(k_ind,:) );

                H(u,s,n,:,k_ind) = (sqrt(1/(Kcluster(k,1)+1)) * squeeze(H(u,s,n,:,k_ind)) + sqrt(Kcluster(k,1)/(Kcluster(k,1)+1)) * temp.').';

                output_Phi_LOS(k_ind,n) = 0;

            end % Rx antennas
        end % Tx antennas

    end % links

    output_Phi_LOS(k_ind,:) = prin_value((output_Phi_LOS(k_ind,:)*180/pi + bulkpar.Phi_LOS(k_ind,:)));

end     % end of B5 part

%%%%%%%%%%%%%%%%
%%%%%%%%%%%
%%%%%%%%
%%%%%        That's all folks !!!
%%
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function that maps inputs from (-inf,inf) to (-180,180)
function y=prin_value(x)
y=mod(x,360);
y=y-360*floor(y/180);