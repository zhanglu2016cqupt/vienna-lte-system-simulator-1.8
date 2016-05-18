classdef channelTraceFactory_v1
    % This class encapsulates all of the methods needed to generate the channel traces (old method)
    % (c) Josep Colom Ikuno, INTHFT, 2011
    
    properties
    end
    
    methods(Static)
        
        function trace_to_fill = trace_SISO( config,H_trace_normalized,H_trace_interf_normalized,debug_output )
            % Generate the fading trace for the SISO LTE mode
            % Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at.
            % (c) 2009 by INTHFT
            % www.nt.tuwien.ac.at
            
            % Re-create config params from input
            system_bandwidth = config.system_bandwidth;
            channel_type     = config.channel_type;
            nTX              = config.nTX;
            nRX              = config.nRX;
            trace_length_s   = config.trace_length_s;
            UE_speed         = config.UE_speed;
            
            size_H = size(H_trace_normalized);
            
            nTTIs       = size_H(4);
            nSC_samples = size_H(5);
            Symbol_samples = size_H(3);

%             H_trace_normalized = permute(H_trace_normalized,[1,2,4,3,5]);
%             H_trace_interf_normalized = permute(H_trace_interf_normalized,[1,2,4,3,5]);
            
            % chi doesn't exist, as there is only one stream being transmitted
            zeta = ones([Symbol_samples,nSC_samples,nTTIs]); % Already permuted
            switch nRX
                case 1
                    % Directly calculate the inverse
                    inv_H = 1./H_trace_normalized;
                    psi   = reshape(abs(inv_H).^2,[Symbol_samples,nTTIs,nSC_samples]);
                    theta = reshape(abs(inv_H .* H_trace_interf_normalized).^2,[Symbol_samples,nTTIs,nSC_samples]);
                otherwise
                    % Calculate MRC for SIMO case
                    h_H_h   = sum(conj(H_trace_normalized).*H_trace_normalized,1);
                    h_i_H_h = sum(conj(H_trace_normalized).*H_trace_interf_normalized,1);
                    psi     = reshape(sum(abs(conj(H_trace_normalized)./h_H_h(ones(1,nRX),:,:,:,:)).^2,1),[Symbol_samples,nTTIs,nSC_samples]);
                    theta   = reshape(abs(h_i_H_h ./ h_H_h).^2,[Symbol_samples,nTTIs,nSC_samples]);
            end
            
            % Rearrange dimensions in the new form (v1.2)

            psi   = permute(psi,[1 3 2]);
            theta = permute(theta,[1 3 2]);
            
            %% Fill in the output trace object
            trace_to_fill                  = phy_modeling.txModeTrace;
            trace_to_fill.tx_mode          = 1;
            trace_to_fill.trace_length_s   = trace_length_s;
            trace_to_fill.system_bandwidth = system_bandwidth;
            trace_to_fill.channel_type     = channel_type;
            trace_to_fill.nTX              = nTX;
            trace_to_fill.nRX              = nRX;
            trace_to_fill.UE_speed         = UE_speed;
            
            trace_to_fill.trace.zeta  = zeta;
            trace_to_fill.trace.psi   = psi;
            trace_to_fill.trace.theta = theta;
            
            %% Some plotting
            
            if debug_output
                figure;
                hold on;
                plot(squeeze(zeta(1,:)),'r','Displayname','\zeta (RX power)');
                plot(squeeze(psi(1,:)),'b','Displayname','\psi (noise enhancement)');
                plot(squeeze(theta(1,:)),'m','Displayname','\theta (inter-cell interference)');
                set(gca,'Yscale','log');
                title('SISO, subcarrier 1');
                grid on;
                legend('show','Location','best');
            end
        end
        
        function trace_to_fill = trace_TxD_2x2( config,H_trace_normalized,H_trace_interf_normalized,precoding_matrix,debug_output )
            % Generate the fading trace for the 2x2 TxD LTE mode
            % Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at.
            % (c) 2009 by INTHFT
            % www.nt.tuwien.ac.at
            
            fprintf('TxD: 2x2\n');
            
            % Re-create config params from input
            system_bandwidth = config.system_bandwidth;
            channel_type     = config.channel_type;
            nTX              = config.nTX;
            nRX              = config.nRX;
            trace_length_s   = config.trace_length_s;
            UE_speed         = config.UE_speed;
            
            % 1/sqrt(2)*precoding_matrix.Z
            % We will use the equivalent channel expression instead of the matrix
            % representation in the standard.
            
            % Preallocate
            TTI_length = trace_length_s*1000;
            nLayers    = precoding_matrix.nLayers;
            SC_samples = size(H_trace_normalized,5);
            Symbol_samples = size(H_trace_normalized,3);
            
            H_tilde1_pinv = zeros(nRX,2*nTX,Symbol_samples,SC_samples,TTI_length);
            A             = zeros(nRX,nRX,Symbol_samples,SC_samples,TTI_length);
            C             = zeros(nRX,nRX,Symbol_samples,SC_samples,TTI_length);
            
            blockLength      = 500;
            nBlocks_start    = 1:blockLength:TTI_length;
            nBlocks_end      = nBlocks_start+blockLength-1;
            nBlocks_end(end) = TTI_length; 
            nBlocks          = length(nBlocks_start);
            output           = cell(1,nBlocks);
            parfor blockIdx=1:nBlocks
                block_H        = H_trace_normalized(:,:,:,nBlocks_start(blockIdx):nBlocks_end(blockIdx),:);
                block_H_interf = H_trace_interf_normalized(:,:,:,nBlocks_start(blockIdx):nBlocks_end(blockIdx),:);
                
                size_H = size(H_tilde1_pinv);
                size_A = size(A);
                size_C = size(C);
                H_tilde1_pinv_block = zeros([size_H(1) size_H(2) size_H(3) size_H(4) size(block_H,4)]);
                A_block = zeros([size_A(1) size_A(2) size_A(3) size_A(4) size(block_H,4)]);
                C_block = zeros([size_C(1) size_C(2) size_C(3) size_C(4) size(block_H,4)]);
                for TTI_ = 1:size(block_H,4)
                    TTI_H       = block_H(:,:,:,TTI_,:);
                    TTI_H_inter = block_H_interf(:,:,:,TTI_,:);
                    
                    [   H_tilde1_pinv_block(:,:,:,:,TTI_),...
                        A_block(:,:,:,:,TTI_),...
                        C_block(:,:,:,:,TTI_) ] = phy_modeling.channelTraceFactory_v1.calculate_TTI_params_TxD_2x2(TTI_H,TTI_H_inter,nTX,nRX,SC_samples,Symbol_samples);
                end
                output{blockIdx}.H_tilde1_pinv_block = H_tilde1_pinv_block;
                output{blockIdx}.A_block             = A_block;
                output{blockIdx}.C_block             = C_block;
            end
            
            for blockIdx=1:nBlocks
                H_tilde1_pinv(:,:,:,:,nBlocks_start(blockIdx):nBlocks_end(blockIdx)) = output{blockIdx}.H_tilde1_pinv_block;
                A(:,:,:,:,nBlocks_start(blockIdx):nBlocks_end(blockIdx))             = output{blockIdx}.A_block;
                C(:,:,:,:,nBlocks_start(blockIdx):nBlocks_end(blockIdx))             = output{blockIdx}.C_block;
            end
            
            %% Extract now the fading parameters
            
            zeta  = zeros(nLayers,Symbol_samples,SC_samples,TTI_length);  % Scales the received signal (1 for perfect channel knowledge)
            for layer_idx = 1:nLayers
                zeta(layer_idx,:,:,:) = squeeze(abs(A(layer_idx,layer_idx,:,:,:)).^2);
            end
            chi   = reshape(sum(abs(A).^2,2),nLayers,Symbol_samples,SC_samples,TTI_length) - zeta;      % Represents inter-layer interference (0 for perfect channel knowledge)
            psi   = reshape(sum(abs(H_tilde1_pinv).^2,2),nLayers,Symbol_samples,SC_samples,TTI_length); % Scales the noise
            theta = reshape(sum(abs(C).^2,2),nLayers,Symbol_samples,SC_samples,TTI_length);             % Scales the interference
            
            %% Fill in the output trace object
            trace_to_fill                  = phy_modeling.txModeTrace;
            trace_to_fill.tx_mode          = 3;
            trace_to_fill.trace_length_s   = trace_length_s;
            trace_to_fill.system_bandwidth = system_bandwidth;
            trace_to_fill.channel_type     = channel_type;
            trace_to_fill.nTX              = nTX;
            trace_to_fill.nRX              = nRX;
            trace_to_fill.UE_speed         = UE_speed;
            
            trace_to_fill.trace.zeta  = zeta;
            trace_to_fill.trace.chi   = chi;
            trace_to_fill.trace.psi   = psi;
            trace_to_fill.trace.theta = theta;
            
            %% Some plotting
            
            if debug_output
                for i_=1:nLayers
                    figure;
                    hold on;
                    plot(norms_H0,'k','Displayname','||H_0||_{F}^2 (Channel norm)');
                    plot(norms_H1,'k:','Displayname','||H_1||_{F}^2 (Interf channel norm)');
                    plot(squeeze(zeta(i_,1,:)),'r','Displayname','\zeta (RX power)');
                    plot(squeeze(chi(i_,1,:)),'g','Displayname','\chi (inter-layer interf)');
                    plot(squeeze(psi(i_,1,:)),'b','Displayname','\psi (noise enhancement)');
                    plot(squeeze(theta(i_,1,:)),'m','Displayname','\theta (inter-cell interference)');
                    set(gca,'Yscale','log');
                    title(sprintf('TxD, Layer %d/%d, subcarrier 1',i_,nLayers));
                    grid on;
                    legend('show','Location','best');
                end
            end
        end
        
        function [H_tilde1_pinv A C] = calculate_TTI_params_TxD_2x2(H_trace_normalized,H_trace_interf_normalized,nTX,nRX,SC_samples,Symbol_samples)
            % H_trace_normalized(:,:,1,N_SCs)
            % H_trace_interf_normalized(:,:,1,N_SCs)
            
            H_tilde1_pinv = zeros(nRX,2*nTX,Symbol_samples,SC_samples);
            A             = zeros(nRX,nRX,Symbol_samples,SC_samples);
            C             = zeros(nRX,nRX,Symbol_samples,SC_samples);
            
            for SC_sample = 1:SC_samples
                for SY_sample = 1:Symbol_samples
                    H_tilde1 = 1/sqrt(2)*[
                        H_trace_normalized(1,:,SY_sample,1,SC_sample).'  H_trace_normalized(2,:,SY_sample,1,SC_sample).'
                        H_trace_normalized(2,:,SY_sample,1,SC_sample)'  -H_trace_normalized(1,:,SY_sample,1,SC_sample)'
                        ];
                    H_tilde2 = 1/sqrt(2)*[
                        H_trace_interf_normalized(1,:,SY_sample,1,SC_sample).'  H_trace_interf_normalized(2,:,SY_sample,1,SC_sample).'
                        H_trace_interf_normalized(2,:,SY_sample,1,SC_sample)'  -H_trace_interf_normalized(1,:,SY_sample,1,SC_sample)'
                        ];
                    H_tilde1_pinv(:,:,SY_sample,SC_sample) = pinv(H_tilde1);
                    A(:,:,SY_sample,SC_sample) = H_tilde1_pinv(:,:,SY_sample,SC_sample)*H_tilde1;
                    C(:,:,SY_sample,SC_sample) = H_tilde1_pinv(:,:,SY_sample,SC_sample)*H_tilde2;
                end
            end
        end
        
        function [ trace_to_fill ] = trace_OLSM( config,H_trace_normalized,H_trace_interf_normalized,precoding_matrix,debug_output,nLayers )
            % Generate the fading trace for the 2x2 OLSM LTE mode
            % Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at.
            % (c) 2009 by INTHFT
            % www.nt.tuwien.ac.at
            
            % Re-create config params from input
            system_bandwidth = config.system_bandwidth;
            channel_type     = config.channel_type;
            nTX              = config.nTX;
            nRX              = config.nRX;
            trace_length_s   = config.trace_length_s;
            UE_speed         = config.UE_speed;
            
            fprintf('OLSM %dx%d: %d Layers\n',nTX,nRX,precoding_matrix.nLayers);
            
            switch nTX
                case 2
                    WDU = precoding_matrix.W * precoding_matrix.D * precoding_matrix.U;
                otherwise
                    num_W_matrices = size(precoding_matrix.W,3);
                    WDU = zeros(nTX,nLayers,num_W_matrices);
                    for w_=1:num_W_matrices
                        WDU(:,:,w_) = precoding_matrix.W(:,:,w_) * precoding_matrix.D * precoding_matrix.U;
                    end
            end
            
            % Preallocate
            TTI_length = trace_length_s*1000;
            nLayers    = precoding_matrix.nLayers;
            SC_samples = size(H_trace_normalized,5);
            Symbol_samples = size(H_trace_normalized,3);
            
            H0F      = zeros(nRX,nLayers,Symbol_samples,SC_samples,TTI_length);
            H1F      = zeros(nRX,nLayers,Symbol_samples,SC_samples,TTI_length);
            H0F_pinv = zeros(nLayers,nRX,Symbol_samples,SC_samples,TTI_length);
            A        = zeros(nLayers,nLayers,Symbol_samples,SC_samples,TTI_length);
            C        = zeros(nLayers,nLayers,Symbol_samples,SC_samples,TTI_length);
            
            for TTI_ = 1:TTI_length
                TTI_H       = H_trace_normalized(:,:,:,TTI_,:);
                TTI_H_inter = H_trace_interf_normalized(:,:,:,TTI_,:);
                [   H0F(:,:,:,:,TTI_)...
                    H1F(:,:,:,:,TTI_)...
                    H0F_pinv(:,:,:,:,TTI_)...
                    A(:,:,:,:,TTI_)...
                    C(:,:,:,:,TTI_) ] = phy_modeling.channelTraceFactory_v1.calculate_TTI_params_OLSM(TTI_H,TTI_H_inter,WDU,nRX,nLayers,SC_samples,Symbol_samples);
            end
            
            %% Extract now the fading parameters
            
            zeta  = zeros(nLayers,Symbol_samples,SC_samples,TTI_length); % Scales the received signal (1 for perfect channel knowledge)
            for layer_idx = 1:nLayers
                zeta(layer_idx,:,:,:) = squeeze(abs(A(layer_idx,layer_idx,:,:,:)).^2);
            end
            chi   = reshape(sum(abs(A).^2,2),nLayers,Symbol_samples,SC_samples,TTI_length) - zeta;                      % Represents inter-layer interference (0 for perfect channel knowledge)
            psi   = reshape(sum(abs(H0F_pinv).^2,2),size(H0F_pinv,1),size(H0F_pinv,3),size(H0F_pinv,4),size(H0F_pinv,5)); % Scales the noise
            theta = reshape(sum(abs(C).^2,2),size(C,1),size(C,3),size(C,4),size(C,5));
            
            % For the 2-layer case, the D matrix is cyclically changed between
            % [1 0;0-1] and [1 0;0 1]
            % Thus, the noise and interference enhancement is averaged between the 2 layers.
            % Since the rest or the values are (for now) 1 or 0, they are not averaged
            psi   = mean(psi,1);
            theta = mean(theta,1);
            
            % Trick to still have matching dimensions
            psi   = repmat(psi,[nLayers 1 1 1]);
            theta = repmat(theta,[nLayers 1 1 1]);
            
            %% Fill in the output trace object
            trace_to_fill                  = phy_modeling.txModeTrace;
            trace_to_fill.tx_mode          = 3;
            trace_to_fill.trace_length_s   = trace_length_s;
            trace_to_fill.system_bandwidth = system_bandwidth;
            trace_to_fill.channel_type     = channel_type;
            trace_to_fill.nTX              = nTX;
            trace_to_fill.nRX              = nRX;
            trace_to_fill.UE_speed         = UE_speed;
            
            %% Slight change to reduce memory consumption
            trace_to_fill.trace.zeta  = [];
            trace_to_fill.trace.chi   = [];
            trace_to_fill.trace.psi   = psi;
            trace_to_fill.trace.theta = theta;
        end
        
        function [H0F H1F H0F_pinv A C] = calculate_TTI_params_OLSM(H_trace_normalized,H_trace_interf_normalized,WDU,nRX,nLayers,SC_samples,Symbol_samples)
            
            H0F      = zeros(nRX,nLayers,Symbol_samples,SC_samples);
            H1F      = zeros(nRX,nLayers,Symbol_samples,SC_samples);
            H0F_pinv = zeros(nLayers,nRX,Symbol_samples,SC_samples);
            A        = zeros(nLayers,nLayers,Symbol_samples,SC_samples);
            C        = zeros(nLayers,size(H1F,2),Symbol_samples,SC_samples);
            
            for SC_sample = 1:SC_samples
                for SY_sample = 1:Symbol_samples
                    current_H0 = H_trace_normalized(:,:,SY_sample,1,SC_sample);
                    current_H1 = H_trace_interf_normalized(:,:,SY_sample,1,SC_sample);
                    %for W_idx = 1:num_W_matrices
                    W_idx = mod(SC_sample-1,size(WDU,3))+1;
                    % Get the pinv(HWDU) matrix for each TTI, SC sample and precoding matrix
                    H0F(:,:,SY_sample,SC_sample)      = current_H0 * WDU(:,:,W_idx);
                    H1F(:,:,SY_sample,SC_sample)      = current_H1 * WDU(:,:,W_idx);
                    H0F_pinv(:,:,SY_sample,SC_sample) = pinv(H0F(:,:,SY_sample,SC_sample));
                    A(:,:,SY_sample,SC_sample)        = H0F_pinv(:,:,SY_sample,SC_sample) * H0F(:,:,SY_sample,SC_sample);
                    C(:,:,SY_sample,SC_sample)        = H0F_pinv(:,:,SY_sample,SC_sample) * H1F(:,:,SY_sample,SC_sample);
                end
            end
        end
        
        function [ trace_to_fill ] = trace_CLSM( config,H_trace_normalized,H_t,H_trace_interf_normalized,precoding_matrix,debug_output)
            % Generate the fading trace for the 2x2 CLSM LTE mode using the PMI
            % feedback proposed in "Mutual Information based calculation of the
            % precoding matrix indicator for 3GPP UMTS/LTE", ITG WSA 2010
            % Author: Stefan Schwarz, sschwarz@nt.tuwien.ac.at.
            % (c) 2010 by INTHFT
            % www.nt.tuwien.ac.at
            % Re-create config params from input
            system_bandwidth   = config.system_bandwidth;
            channel_type       = config.channel_type;
            nTX                = config.nTX;
            nRX                = config.nRX;
            trace_length_s     = config.trace_length_s;
            wideband_precoding = config.wideband_precoding;
            
            fprintf('CLSM %dx%d: %d Layers, parallelized over ',nTX,nRX,precoding_matrix.nLayers);
            
            UE_speed = config.UE_speed;
            nLayers  = precoding_matrix.nLayers;
            
            % Preallocate
            TTI_length = trace_length_s*1000;
            SC_samples = size(H_trace_normalized,5);
            Symbol_samples = size(H_trace_normalized,3);
            
            H0F       = zeros(nRX,nLayers,Symbol_samples,SC_samples,TTI_length);
            H1F       = zeros(nRX,nLayers,Symbol_samples,SC_samples,TTI_length);
            H0F_pinv  = zeros(nLayers,nRX,Symbol_samples,SC_samples,TTI_length);
            A         = zeros(nLayers,nLayers,Symbol_samples,SC_samples,TTI_length);
            C         = zeros(nLayers,nLayers,Symbol_samples,SC_samples,TTI_length);
            PMI_trace = zeros(SC_samples/2,TTI_length);
            PMI       = zeros(SC_samples/2,TTI_length+config.feedback_channel_delay);
            
            TTI_values_loop   = 1:TTI_length+config.feedback_channel_delay;
            if config.parallel_toolbox_installed
                number_of_labs = max (matlabpool('size'),1);
            else
                number_of_labs = 1;
            end
            
            TTIs_per_lab      = floor(length(TTI_values_loop)/number_of_labs)*ones(1,number_of_labs);
            TTIs_per_lab(end) = TTIs_per_lab(end)+(length(TTI_values_loop)-sum(TTIs_per_lab));
            lab_end_idxs      = cumsum(TTIs_per_lab);
            lab_start_idxs    = [1 (lab_end_idxs(1:(end-1))+1)];
            
            PMI_tmp    = struct('PMI',[]);
            PMI_tmp    = PMI_tmp(ones(1,number_of_labs));
            TTI_params = struct('H0F',[],'H1F',[],'H0F_pinv',[],'A',[],'C',[],'PMI_trace',[]);
            TTI_params = TTI_params(ones(1,number_of_labs));
            
            fprintf('%g labs\n',number_of_labs);
            
            % NOTE: assumed is number_of_labs<<<TTI_length, which should be the case unless you have thousands of cores.
            
            % Calculate PMI by means of a parfor loop (the most time-consuming part) and post-process
            fprintf('  Feedback calculation\n');
%             H_trace_normalized(:,:,2)
            parfor slice_num = 1:number_of_labs
                TTI_vect               = lab_start_idxs(slice_num):lab_end_idxs(slice_num);
                PMI_tmp(slice_num).PMI = zeros(SC_samples/2,length(TTI_vect));
                
                for TTI_idx = 1:length(TTI_vect)
                    TTI_ = TTI_vect(TTI_idx);
                    if TTI_ <= TTI_length
                        TTI_H_t       = reshape(H_t(:,:,:,TTI_,:),[size(H_t,1) size(H_t,2) size(H_t,3) size(H_t,5)]);
                    else
                        TTI_H_t       = reshape(H_t(:,:,:,mod(TTI_,TTI_length),:),[size(H_t,1) size(H_t,2) size(H_t,3) size(H_t,5)]);
                    end
                    [PMI_tmp(slice_num).PMI(:,TTI_idx)] = feedback_calculation.Rel8_CLSM_PMI_layer_specific(TTI_H_t,nLayers,precoding_matrix.W,nTX,'ZF',wideband_precoding,[]);
                end
            end
            for slice_num = 1:number_of_labs
                slice_TTIs = lab_start_idxs(slice_num):lab_end_idxs(slice_num);
                PMI(:,slice_TTIs) = PMI_tmp(slice_num).PMI;
            end
%             PMI(:,2)
            % Calculate the rest of the TTI parameters with a second parfor
            fprintf('  TTI parameters\n');
            parfor slice_num = 1:number_of_labs
                TTI_vect                        = lab_start_idxs(slice_num):lab_end_idxs(slice_num);
                TTI_params(slice_num).H0F       = zeros(nRX,nLayers,Symbol_samples,SC_samples,length(TTI_vect));
                TTI_params(slice_num).H1F       = zeros(nRX,nLayers,Symbol_samples,SC_samples,length(TTI_vect));
                TTI_params(slice_num).H0F_pinv  = zeros(nLayers,nRX,Symbol_samples,SC_samples,length(TTI_vect));
                TTI_params(slice_num).A         = zeros(nLayers,nLayers,Symbol_samples,SC_samples,length(TTI_vect));
                TTI_params(slice_num).C         = zeros(nLayers,nLayers,Symbol_samples,SC_samples,length(TTI_vect));
                TTI_params(slice_num).PMI_trace = zeros(SC_samples/2,length(TTI_vect));
                
                for TTI_idx = 1:length(TTI_vect)
                    TTI_ = TTI_vect(TTI_idx);
                    if TTI_ <= TTI_length
                        TTI_H         = H_trace_normalized(:,:,:,TTI_,:);
                        TTI_H_inter   = H_trace_interf_normalized(:,:,:,TTI_,:);
                    else
                        TTI_H         = H_trace_normalized(:,:,:,mod(TTI_,TTI_length),:);
                        TTI_H_inter   = H_trace_interf_normalized(:,:,:,mod(TTI_,TTI_length),:);
                    end
                    if TTI_ >  config.feedback_channel_delay
                        TTI__ = TTI_ - config.feedback_channel_delay;
                        PMI_act = PMI(:,TTI__);
                        
                        [   TTI_params(slice_num).H0F(:,:,:,:,TTI_idx)...
                            TTI_params(slice_num).H1F(:,:,:,:,TTI_idx)...
                            TTI_params(slice_num).H0F_pinv(:,:,:,:,TTI_idx)...
                            TTI_params(slice_num).A(:,:,:,:,TTI_idx)...
                            TTI_params(slice_num).C(:,:,:,:,TTI_idx)...
                            TTI_params(slice_num).PMI_trace(:,TTI_idx)] = phy_modeling.channelTraceFactory_v1.calculate_TTI_params_CLSM(TTI_H,TTI_H_inter,precoding_matrix.W,nRX,nLayers,SC_samples,Symbol_samples,PMI_act);
                    end
                end
            end
            for slice_num = 1:number_of_labs
                slice_TTIs = lab_start_idxs(slice_num):lab_end_idxs(slice_num);
                if slice_num==number_of_labs
                    slice_TTIs = slice_TTIs(1:(end-config.feedback_channel_delay));
                    tail_cut   = config.feedback_channel_delay;
                else
                    tail_cut   = 0;
                end
                H0F(:,:,:,:,slice_TTIs)       = TTI_params(slice_num).H0F(:,:,:,:,1:(end-tail_cut));
                H1F(:,:,:,:,slice_TTIs)       = TTI_params(slice_num).H1F(:,:,:,:,1:(end-tail_cut));
                H0F_pinv(:,:,:,:,slice_TTIs)  = TTI_params(slice_num).H0F_pinv(:,:,:,:,1:(end-tail_cut));
                A(:,:,:,:,slice_TTIs)         = TTI_params(slice_num).A(:,:,:,:,1:(end-tail_cut));
                C(:,:,:,:,slice_TTIs)         = TTI_params(slice_num).C(:,:,:,:,1:(end-tail_cut));
                PMI_trace(:,slice_TTIs)     = TTI_params(slice_num).PMI_trace(:,1:(end-tail_cut));
            end
            clear PMI_tmp
            clear TTI_params

            % Extract now the fading parameters
            zeta  = zeros(nLayers,Symbol_samples,SC_samples,TTI_length);  % Scales the received signal (1 for perfect channel knowledge)
            for layer_idx = 1:nLayers
                zeta(layer_idx,:,:,:) = squeeze(abs(A(layer_idx,layer_idx,:,:,:)).^2);
            end
            chi   = reshape(sum(abs(A).^2,2),nLayers,Symbol_samples,SC_samples,TTI_length) - zeta;                      % Represents inter-layer interference (0 for perfect channel knowledge)
            psi   = reshape(sum(abs(H0F_pinv).^2,2),size(H0F_pinv,1),size(H0F_pinv,3),size(H0F_pinv,4),size(H0F_pinv,5)); % Scales the noise
            theta = reshape(sum(abs(C).^2,2),size(C,1),size(C,3),size(C,4),size(C,5));                             % Scales the interference
            
            %% Fill in the output trace object
            trace_to_fill                  = phy_modeling.txModeTrace;
            trace_to_fill.tx_mode          = 4;
            trace_to_fill.trace_length_s   = trace_length_s;
            trace_to_fill.system_bandwidth = system_bandwidth;
            trace_to_fill.channel_type     = channel_type;
            trace_to_fill.nTX              = nTX;
            trace_to_fill.nRX              = nRX;
            trace_to_fill.UE_speed         = UE_speed;
            
            %% Slight change to reduce memory consumption
            trace_to_fill.trace.zeta  = [];
            trace_to_fill.trace.chi   = [];
            trace_to_fill.trace.psi   = psi;
            trace_to_fill.trace.theta = theta;
            trace_to_fill.trace.PMI   = PMI_trace;
        end
        function [H0F H1F H0F_pinv A C PMI] = calculate_TTI_params_CLSM(H_trace_normalized,H_trace_interf_normalized,W,nRX,nLayers,SC_samples,Symbol_samples,PMI)
            % H_trace_normalized(:,:,1,N_SCs)
            % H_trace_interf_normalized(:,:,1,N_SCs)
            
            H0F      = zeros(nRX,nLayers,Symbol_samples,SC_samples);
            H1F      = zeros(nRX,nLayers,Symbol_samples,SC_samples);
            H0F_pinv = zeros(nLayers,nRX,Symbol_samples,SC_samples);
            A        = zeros(nLayers,nLayers,Symbol_samples,SC_samples);
            C        = zeros(nLayers,size(H1F,2),Symbol_samples,SC_samples);
            
            for SC_sample = 1:SC_samples
                for SY_sample = 1:Symbol_samples
                    current_H0 = H_trace_normalized(:,:,SY_sample,1,SC_sample);
                    current_H1 = H_trace_interf_normalized(:,:,SY_sample,1,SC_sample);
                    W_idx_H0 = PMI(ceil(SC_sample/2)); % H0 precoding matrix
                    W_idx_H1 = 1;                      % dummy interfering precoding matrix

                    % Get the pinv(HW) matrix for each TTI, SC sample and precoding matrix
                    H0F_      = current_H0 * W(:,:,W_idx_H0);
                    H1F_      = current_H1 * W(:,:,W_idx_H1);
                    H0F_pinv_ = pinv(H0F_);

                    H0F(:,:,SY_sample,SC_sample)      = H0F_;
                    H1F(:,:,SY_sample,SC_sample)      = H1F_;

                    H0F_pinv(:,:,SY_sample,SC_sample) = H0F_pinv_;
                    A(:,:,SY_sample,SC_sample)        = H0F_pinv_ * H0F_;
                    C(:,:,SY_sample,SC_sample)        = H0F_pinv_ * H1F_;
                end
            end
            
        end
        
        function pregenerated_fast_fading = generate_channel_trace(config,precoding_configs,H_trace_normalized,H_trace_interf_normalized)
            
            % The channel trace is stored in H_trace.H_RB_samples, containing a
            % (nRX,nTX,subframe_num,sample_num) matrix. We will calculate the SINR
            % trace for each of these samples (every 6 subcarriers).
            
            show_trace_generation_debug = false;
            
            max_layers_SM = min(config.nTX,config.nRX);
            if config.parallel_toolbox_installed && ~matlabpool('size') && ~config.non_parallel_channel_trace && config.tx_mode~=1
                try
                    matlabpool open;
                catch 
                    fprintf('Failed to open matlabpool. Maybe already open?\n');
                end
            end
            
            switch config.tx_mode
                case 1
                    % SISO trace
                    SISO_trace = phy_modeling.channelTraceFactory_v1.trace_SISO(config,H_trace_normalized,H_trace_interf_normalized,show_trace_generation_debug);
                case 2
                    % TxD (up to 2x2)
                    clear precoding_matrix
                    precoding_matrix = precoding_configs{config.nTX,config.nTX,config.tx_mode};
                    
                    %% Call trace generation function
                    TxD_2x2_trace = phy_modeling.channelTraceFactory_v1.trace_TxD_2x2(config,H_trace_normalized,H_trace_interf_normalized,precoding_matrix,show_trace_generation_debug);
                case 3
                    % OLSM
                    OLSM_trace = cell(max_layers_SM,1);
                    switch config.nTX
                        case {2,4}
                            sliced_precoding_matrices = [precoding_configs{1:config.nTX,config.nTX,config.tx_mode}];
                        otherwise
                            error('2 or 4 TX antennas supported');
                    end
                    
                    %% Call trace generation function
                    parfor layer_i = 1:min(config.nTX,config.nRX)
                        OLSM_trace{layer_i} = phy_modeling.channelTraceFactory_v1.trace_OLSM(config,H_trace_normalized,H_trace_interf_normalized,sliced_precoding_matrices(layer_i),show_trace_generation_debug,layer_i);
                    end
                case {4, 5, 6,9}
                    % CLSM
                    if config.tx_mode ==6 || config.tx_mode == 5
                        max_layers_SM = 1;
                    end
                    
                    if config.tx_mode == 9 && config.nTX ~= 8
                       error('for TxMode 9 only 8 Tx antennas are possible'); 
                    end
                    
                    CLSM_trace = cell(max_layers_SM,1);                    
                    

                    switch config.nTX
                        case {2,4}
                            sliced_precoding_matrices = [precoding_configs{1:config.nTX,config.nTX,config.tx_mode}]; % Precoding matrices for CLSM, 2 or 4 antenna ports
                        case 8
                            sliced_precoding_matrices = [precoding_configs{1:config.nTX,config.nTX,config.tx_mode}]; % Precoding matrices for Tx Mode 9
                        otherwise
                            error('only 2, 4 or 8 TX antennas supported');
                    end
                    
                    % Calculate the average H from each RB (needed for the PMI calculation)
                    H_size   = size(H_trace_normalized); % channel matrix sizes
                    N_sc     = 2;                        % number of subcarriers per resource block (just 2 samples are picked per RB)
                    N_rb     = H_size(5)/N_sc;           % number of resource blocks
                    N_s      = H_size(3);                % number of symbols per TTI
                    H_t = zeros([H_size(1) H_size(2) H_size(3) H_size(4) N_rb]); % average channel of the RB                    
                    for i_=1:H_size(4)
                        H_TTI = reshape(H_trace_normalized(:,:,:,i_,:),[H_size(1) H_size(2) H_size(3) H_size(5)]);
                        H_t(:,:,:,i_,:) = reshape(mean(reshape(H_TTI,[H_size(1)*H_size(2) N_s N_sc N_rb]),3),[H_size(1) H_size(2) N_s N_rb]);
                    end
                     if config.tx_mode ==6 ||config.tx_mode == 5
                         CLSM_trace{1} = phy_modeling.channelTraceFactory_v1.trace_CLSM(config,H_trace_normalized,H_t,H_trace_interf_normalized,sliced_precoding_matrices(1),show_trace_generation_debug);
                     
                     else
%                      plot([squeeze(abs(H_trace_normalized(1,1,:,1,1))).',squeeze(abs(H_trace_normalized(1,1,:,2,1))).',squeeze(abs(H_trace_normalized(1,1,:,3,1))).',squeeze(abs(H_trace_normalized(1,1,:,4,1))).',squeeze(abs(H_trace_normalized(1,1,:,5,1))).',squeeze(abs(H_trace_normalized(1,1,:,6,1))).'])
                    %% Call trace generation function
                    for layer_i = 1:min(config.nTX,config.nRX)
                        CLSM_trace{layer_i} = phy_modeling.channelTraceFactory_v1.trace_CLSM(config,H_trace_normalized,H_t,H_trace_interf_normalized,sliced_precoding_matrices(layer_i),show_trace_generation_debug);
                    end
                     end
                otherwise
                    error('Tx mode %d not supported',config.tx_mode);
            end
%             if config.parallel_toolbox_installed && ~config.non_parallel_channel_trace && config.tx_mode~=1
%                 matlabpool close;
%             end
            
            %% Create output fast fading trace
            pregenerated_fast_fading                      = phy_modeling.PregeneratedFastFading;
            pregenerated_fast_fading.trace_length_s       = config.trace_length_s;
            pregenerated_fast_fading.trace_length_samples = config.trace_length_s / 1e-3;
            pregenerated_fast_fading.system_bandwidth     = config.system_bandwidth;
            pregenerated_fast_fading.channel_type         = config.channel_type;
            pregenerated_fast_fading.nTX                  = config.nTX;
            pregenerated_fast_fading.nRX                  = config.nRX;
            pregenerated_fast_fading.UE_speed             = config.UE_speed;
            pregenerated_fast_fading.source_info          = [];
            pregenerated_fast_fading.generated_from       = '';
            pregenerated_fast_fading.t_step               = 1e-3;
            pregenerated_fast_fading.f_step               = 15e3*6;
            
            switch config.tx_mode
                case 1
                    % SISO trace (mode 1)
                    pregenerated_fast_fading.traces{1} = SISO_trace;
                case 2
                    % TxD trace (mode 2)
                    pregenerated_fast_fading.traces{2} = TxD_2x2_trace;
                case 3
                    % OLSM trace (mode 3)
                    pregenerated_fast_fading.traces{3} = OLSM_trace;
                case 4
                    % CLSM trace (mode 4)
                    pregenerated_fast_fading.traces{4} = CLSM_trace;
                case 5
                    % MU-MIMO needs precoders at scheduling time -> no
                    % precalculations
                    pregenerated_fast_fading.traces{5} = CLSM_trace;
                case 6
                    % Rank 1 CLSM trace (mode 6)
                    pregenerated_fast_fading.traces{6} = CLSM_trace;
                case 9
                    % Up to 8 layer transmission (mode 9)
                    pregenerated_fast_fading.traces{9} = CLSM_trace;
            end
        end
    end
end