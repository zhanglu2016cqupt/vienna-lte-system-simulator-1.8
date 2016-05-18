close all; clc; clear all;

% Martin Taranetz, Jan. 2013
 
% Use codesnippets at the end of the code to extract channel taps.
% 1. ) Evaluate statistics of composite shadow and microscale fading
%      Take only first column of composite_fading_taps since shadow fading is
%      the same for all RBs.
% 2. 1. ) Evaluate statistics of aggregated interference

load('Interference Statistics_Rayleigh.mat')
fadingData     = sort(aggregated_interference_taps); 
fadingData     = fadingData(1:200:end);
% fadingData     = fadingData(1:end-floor(length(microscale_fading_taps)*0.01));
% fadingData_log = log(fadingData);

%% Replace dfittool by manual steps
% Scaled histogram in linear scale

fadingData=theta;
fadingData=fadingData(:,1);

dmin = min(fadingData);
dmax = max(fadingData);
binCount = 1e4;
binWidth = (dmax - dmin)/binCount;
nSamples = length(fadingData);
[n xout] = hist(fadingData, binCount);
figure('name', 'Histogram'); hold on;
xlim([0 5]);
plot(xout,n/(binWidth*nSamples));

% ***** EXPONENTIAL
parmhat_e = expfit(fadingData);
Y_exp   = pdf('exp', xout, parmhat_e(1));
plot(xout, Y_exp,'r');
% ***** GAMMA
parmhat_g = gamfit(fadingData);
Y_gam   = pdf('gam', xout, parmhat_g(1), parmhat_g(2));
plot(xout, Y_gam, 'm');
% ***** LOGNORMAL
parmhat_l = lognfit(fadingData);
Y_logn  = pdf('logn', xout, parmhat_l(1), parmhat_l(2));
plot(xout, Y_logn,'g');
hold off;

% Scaled histogram in logarithmic scale
% dmin_log = min(fadingData_log);
% dmax_log = max(fadingData_log);
% BinWidth_log = (dmax_log - dmin_log)/BinCount;
% [n_log xout_log] = hist(fadingData_log, BinCount);
% figure('name', 'Histogram in logarithmic scale');
% plot(xout_log, n_log/(BinWidth_log*NSamples));


% Distribution Fitting
% X_lin = dmin:binWidth:dmax; % Axis for pdf in linear scale
% X_log = dmin_log:0.001:dmax_log; % Axis for pdf in log scale

% Fit to Linear Data
% figure('name', 'Fit distribution in linear scale');
% **** LOGNORMAL
% parmhat_logn = lognfit(fadingData)
% Y_logn = pdf('logn', X_lin, parmhat_logn(1), parmhat_logn(2));
% figure('name', 'Lognormal fit in linear scale');
% plot(X_lin,Y_logn);
% **** GAMMA
% parmhat  = gamfit(FadingData);
% Y_pdf = pdf('gam', exp(X_log), parmhat(1), parmhat(2));
% figure('name', 'Distribution fitted on linear data - log plot');
% plot(X_log, Y_pdf.*exp(X_log));

% Fit to Logarithmic Data
% [mu_norm sigma_norm] = normfit(fadingData_log)
% Y_normal = pdf('norm', X_log, mu_norm, sigma_norm);
% figure('name', 'Fit distribution in log scale');
% plot(X_log, Y_normal);

% Fitting gamma distribution to Lognormal distribution
% mu and sigma are the parameters of the Gaussian distribution
% corresponding to the lognormal distribution
% a,b are the parameters of the gamma distribution
% a = (exp(sigma_norm^2)-1)^(-1)
% b = exp(mu_norm + sigma_norm^2/2)*(exp(sigma_norm^2)-1)

%% Stable distribution
% p_stbl = stblfit(FadingData)
% Y_stbl  = stblpdf(X_lin, p_stbl(1), p_stbl(2), p_stbl(3), p_stbl(4));
% Y_stbl_log_fit  = stblpdf(exp(X_log), p_stbl(1), p_stbl(2), p_stbl(3), p_stbl(4));
% Y_stbl_log_calc = stblpdf(exp(X_log),  0.33, 1, 7.74e-10, -1e-11);
% figure('name', 'Stable Distribution in log scale'); hold on;
% plot(X_log, Y_stbl);
%plot(X_log, Y_stbl_log_fit.*exp(X_log));
% plot(X_log, Y_stbl_log_calc.*exp(X_log),'g'); hold off;

%% Add this to UE.m in the "if there_are_interferers" branch
% Temporarily safe interference power on per RB block basis; For tracing and statistical evaluation
% of interference.

% if (obj.clock.current_TTI==1)&(~obj.deactivate_UE)
%     % Assume unit power, i.e. transmit power per RB = 1
%     % Take only taps from first RB as representative
%     % Reason: Shadow fading is constant for all RBs.
%     microscale_fading_taps_temp  = microscale_interfering_thetas(:,1);
%     shadow_fading_taps_temp      = interfering_shadow_fading_loss_linear;
%     composite_fading_taps_temp   = microscale_fading_taps_temp.*shadow_fading_taps_temp;
%     aggregated_interference_temp = sum(microscale_interfering_thetas./temp_macro_mat'./temp_shadow_mat',1);
%     if (exist('Interference Statistics.mat','file'))
%         load('Interference Statistics.mat', 'microscale_fading_taps', 'shadow_fading_taps', 'composite_fading_taps', 'aggregated_interference_taps');
%         microscale_fading_taps       = [microscale_fading_taps;       microscale_fading_taps_temp];
%         shadow_fading_taps           = [shadow_fading_taps;           shadow_fading_taps_temp];
%         composite_fading_taps        = [composite_fading_taps;        composite_fading_taps_temp];
%         aggregated_interference_taps = [aggregated_interference_taps; aggregated_interference_temp(1)];
%         save('Interference Statistics.mat', 'microscale_fading_taps', 'shadow_fading_taps','composite_fading_taps','aggregated_interference_taps');
%     else
%         microscale_fading_taps       = microscale_fading_taps_temp;
%         shadow_fading_taps           = shadow_fading_taps_temp;
%         composite_fading_taps_temp   = composite_fading_taps_temp;
%         aggregated_interference_taps = aggregated_interference_temp(1);
%         save('Interference Statistics.mat', 'microscale_fading_taps', 'shadow_fading_taps','composite_fading_taps','aggregated_interference_taps');
%     end
% end
