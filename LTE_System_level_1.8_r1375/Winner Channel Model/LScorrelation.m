%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generation of correlated large scale (LS) parameters DS,ASA,ASD
%  and SF for all links
%
%  If layout parameters (co-ordinates etc.) are given, auto-correlation
%  between channel segments is generated. If layout parameters are not
%  defined, only cross-correlation between LS parameters is generated.
%
%  Notes! Auto-correlation is generated for each BS separately. If two MSs
%  are linked to same BS, correlation is distance dependent. If one MS is
%  linked to two sectors of single BS, correlation is full.
%
% The procedure and notations are from [1]
% Ref: [1]
% Authors: Pekka Kyösti (EBIT), Martin Käske (TUI)

% Updated by Pekka Kyösti 8.6.2007: Auto-correlation generation changed to
% 2D- filtering of grid of samples (based on Per Zetterberg's original
% proposal). Cross correlation generation cleaned to simple linear
% transformation (multiplication by R^0.5). Distribution transformations
% simplified, all distributions are log-Normal in WINNER II.

function sigmas = LScorrelation(wimpar,linkpar,fixpar,iterpar)

% scenario parameters
Scenario = iterpar.Scenario;
PropagCondition = iterpar.PropagCondition;
LoSConnectionLinks = iterpar.LoSConnectionLinks;
NLoSConnectionLinks = iterpar.NLoSConnectionLinks;
NumLoSConnectionLinks = length(LoSConnectionLinks);
NumNLoSConnectionLinks = length(NLoSConnectionLinks);
N = iterpar.N;
N_max = max(N(1),N(2));
NumLinks = length(iterpar.UserIndeces);

for loop_index = 1:2 %los and nlos links are handled separately
    if loop_index == 1 & LoSConnectionLinks
        Condition = 'LoS';
        NumLinks = NumLoSConnectionLinks;
    elseif loop_index == 2 & NLoSConnectionLinks
        Condition = 'NLoS';
        NumLinks = NumNLoSConnectionLinks;
    else
        continue;
    end      

    evalstr = sprintf('DS_lambda = iterpar.%s.DS_lambda;',Condition); eval(evalstr);
    evalstr = sprintf('AS_D_lambda = iterpar.%s.AS_D_lambda;',Condition); eval(evalstr);
    evalstr = sprintf('AS_A_lambda = iterpar.%s.AS_A_lambda;',Condition); eval(evalstr);
    evalstr = sprintf('SF_lambda = iterpar.%s.SF_lambda;',Condition); eval(evalstr);
    evalstr = sprintf('KF_lambda = iterpar.%s.KF_lambda;',Condition); eval(evalstr);

    if 1  % layout parameters are always given
        % layout parameters
        NofSect = linkpar.NofSect;
        Stations = linkpar.Stations;
        Pairing = linkpar.Pairing(:,iterpar.UserIndeces);
        if loop_index == 1  % LOS condition
            Pairing = Pairing(:,LoSConnectionLinks);  % Take only LOS links
        else    % NLOS condition
            Pairing = Pairing(:,NLoSConnectionLinks);  % Take only NLOS links
        end 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Auto-correlation generation %%%
        % determine auto-correlation separately for each BS
        for jj=1:max(Pairing(1,:))
            indLink = find(Pairing(1,:)==jj);   % indices to links of BS k
            tmpMs = Pairing(2,indLink);     % indices to MSs linked to BS k
            if length(tmpMs)>1  % if only 1 MS linked to BS, no auto-correlation
                
                % Generate grid of iid Gaussian random numbers ~N(0,1) with
                %  100 extra samples to all directions (+-x and +-y coordinates)
                xtra = 100;
                MsPos=[Stations(tmpMs).Pos];
                gridn = randn(max(MsPos(2,:))-min(MsPos(2,:))+2*xtra+1,...
                              max(MsPos(1,:))-min(MsPos(1,:))+2*xtra+1 ,5);
                % Index to MS locations on the grid  
                gind = sub2ind(size(gridn),MsPos(2,:)-min(MsPos(2,:))+xtra+1,MsPos(1,:)-min(MsPos(1,:))+xtra+1);
                
                % Define auto-correlation filter for each of the 4 LS parameters
                delta = [DS_lambda AS_D_lambda AS_A_lambda SF_lambda KF_lambda];
                d = 0:100;
                h = exp(-1*repmat(d',1,5)./repmat(delta,length(d),1)); % d>0
                h = h./repmat(sum(h),length(d),1);  
                
                % Filter Gaussian grid in 2D to get exponential auto-correlation
                for ii=1:5
                    if(delta(ii)~=0)
                        tmp = filter(h(:,ii),1,gridn(:,:,ii),[],1); 
                        grida = filter(h(:,ii),1,tmp,[],2);
                        % Pick correlated MS locations from the grid 
                        ksi(ii,indLink) = grida(gind); clear tmp grida
                    else
                        ksi(ii,indLink) = gridn(gind);
                    end;
                end
                
            elseif ~isempty(indLink)    % only one MS linked, no auto-correlation
                
                ksi(:,indLink) = randn(5,1);
                
            end % end if >1 MSs

        end % end for BS sites

    else  % if layout parameters not given

        % generate non-correlated 4xK Gaussian random numbers
        ksi = randn(5,NumLinks);

    end % if (layout parameters given)

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Cross-correlation generation %%%
    % Extract cross correlation parameters from input
    evalstr = sprintf('a = iterpar.%s.asD_ds;',Condition); eval(evalstr); % departure AS vs delay spread
    evalstr = sprintf('b = iterpar.%s.asA_ds;',Condition); eval(evalstr); % arrival AS vs delay spread
    evalstr = sprintf('c = iterpar.%s.asA_sf;',Condition); eval(evalstr); % arrival AS vs shadowing std
    evalstr = sprintf('d = iterpar.%s.asD_sf;',Condition); eval(evalstr); % departure AS vs shadoving std
    evalstr = sprintf('e = iterpar.%s.ds_sf;',Condition); eval(evalstr);  % delay spread vs shadoving std
    evalstr = sprintf('f = iterpar.%s.asD_asA;',Condition); eval(evalstr); % departure AS vs arrival AS
    evalstr = sprintf('g = iterpar.%s.asD_kf;',Condition); eval(evalstr); % departure AS vs k-factor
    evalstr = sprintf('h = iterpar.%s.asA_kf;',Condition); eval(evalstr); % arrival AS vs k-factor
    evalstr = sprintf('k = iterpar.%s.ds_kf;',Condition); eval(evalstr); % delay spread vs k-factor
    evalstr = sprintf('l = iterpar.%s.sf_kf;',Condition); eval(evalstr); % shadowing std vs k-factor

    
    
    % Cross-correlation matrix 
    % Order of rows and columns is ds asD asA sf
    A = [ 1  a  b  e  k ;...
          a  1  f  d  g ;...
          b  f  1  c  h ;...
          e  d  c  1  l ;...
          k  g  h  l  1 ];

    % get A^(0.5)    (equal to sqrtm(A))
    R_sqrt = sqrtm(A);
    ksi = R_sqrt*ksi;   % generate cross-correlation by linear transformation
    
    evalstr = sprintf('a = iterpar.%s.DS_mu;',Condition); eval(evalstr); 
    evalstr = sprintf('b = iterpar.%s.DS_sigma;',Condition); eval(evalstr); 
    evalstr = sprintf('c = iterpar.%s.AS_D_mu;',Condition); eval(evalstr); 
    evalstr = sprintf('d = iterpar.%s.AS_D_sigma;',Condition); eval(evalstr);
    evalstr = sprintf('e = iterpar.%s.AS_A_mu;',Condition); eval(evalstr); 
    evalstr = sprintf('f = iterpar.%s.AS_A_sigma;',Condition); eval(evalstr); 
    evalstr = sprintf('g = iterpar.%s.SF_sigma;',Condition); eval(evalstr);
    evalstr = sprintf('h = iterpar.%s.KF_mu;',Condition); eval(evalstr);
    evalstr = sprintf('k = iterpar.%s.KF_sigma;',Condition); eval(evalstr);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Transform Normal distributed random numbers to scenario specific distributions
    
    sigma_ds  = 10.^(b*ksi(1,:).' + a);      % Log-Normal 
    sigma_asD = 10.^(d*ksi(2,:).' + c);      % Log-Normal 
    sigma_asA = 10.^(f*ksi(3,:).' + e);      % Log-Normal 
    sigma_sf  = 10.^(0.1*g*ksi(4,:).');      % Log-Normal dB
    sigma_kf  = 10.^(0.1*(k*ksi(5,:).'+ h));   % Log-Normal dB

    
    % output
    evalstr= sprintf('sigmas_%s = [sigma_asD sigma_asA sigma_ds sigma_sf sigma_kf];',Condition);
    eval(evalstr);
    clear ksi indBs
end  % END loop_index (LOS/NLOS)

sigmas = NaN*ones(length(iterpar.UserIndeces),5);
if LoSConnectionLinks
    sigmas(LoSConnectionLinks,:) = sigmas_LoS;
end
if NLoSConnectionLinks
    sigmas(NLoSConnectionLinks,:) = sigmas_NLoS;
end
