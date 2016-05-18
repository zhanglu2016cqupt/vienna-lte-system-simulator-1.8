function feedback = CQI_and_rank(...
    config,...
    tx_mode,...
    SINR_dB,...
    nRB,...
    feedback_PMIs,...
    DL_signaling,...
    SINR_averager,...
    CQI_mapper,...
    unquantized_CQI_feedback,...
    adaptive_RI)
% Calculate the feedback values based on the input. This function
% is called from the link quality model and is separated for
% convenience and readability. The results of the feedback
% calculation are stored in the following variables:
% - obj.feedback.CQI:              CQI feedback
% - obj.feedback.RI:               Rank Indicator feedback (when applicable)
% - obj.link_quality_model_output: SINR values
%
% As input parameters you have one SINR per RB
% (SINRs_to_map_to_CQI) or all of the SINRs the SL simulator
% traces, which are currently two per RB (SINR_dB)
%
% Take a subset of the SINRs for feedback calculation
% For SM we send 2 CQIs, one for each of the codewords (which in the 2x2
% case are also the layers). For TxD, both layers have the same SINR
% The CQI is calculated as a linear averaging of the SINRs in
% dB. This is done because like this the Tx has an "overall
% idea" of the state of the RB, not just a sample of it.
%
% (c) Stefan Schwarz, Josep Colom Ikuno ITC, 2013
% www.nt.tuwien.ac.at

% switch tx_mode
%     case 1 % SISO
%         SINRs_to_map_to_CQI = (SINR_dB(1:2:end)+SINR_dB(2:2:end))/2;  % NOTE: linear averaging of logarithmic SNR - what the heck???        
%     case 2 % TxD
%         % Both layers have the same SINR
%         SINRs_to_map_to_CQI = (SINR_dB(1,1:2:end)+SINR_dB(1,2:2:end))/2;
%     case {3,4,5,6,9} % OLSM, CLSM
%         SINRs_to_map_to_CQI = (SINR_dB(:,1:2:end,:)+SINR_dB(:,2:2:end,:))/2;
%     otherwise
%         error('TX mode not yet supported');
% end

size_vec = size(SINR_dB);
size_vec = [size_vec,1]; % to make sure that the length of size_vec is at least 4 (not the case for SISO otherwise)

% use BICM SINR averaging instead of the strange linear average above
SINRs_to_map_to_CQI = NaN(size_vec(1),size_vec(3)/2,size_vec(4));
SINR_temp = cat(2,SINR_dB(:,:,1:2:end,:),SINR_dB(:,:,2:2:end,:));
% CQI_SNR_table = [-6.934;-5.147;-3.18;-1.254;0.761;2.70;4.697;6.528;8.576;10.37;12.3;14.18;15.89;17.82;19.83;21];
% if size_vec(2) > 1  % this is to mimic the per-slot CQI feedback of the link level simulator
%     for rb_i = 1:nRB % RBs
%         for r_i = 1:size_vec(1) % rank_possibilities
%             for l_i = 1:r_i % layers
%                 SINR_slot1 = squeeze(SINR_temp(l_i,1:ceil(size_vec(2)/2),rb_i,r_i));
%                 SINR_slot2 = squeeze(SINR_temp(l_i,ceil(size_vec(2)/2)+1:end,rb_i,r_i));
%                 SINR_av1 = SINR_averager.average(SINR_slot1,15,true);
%                 SINR_av2 = SINR_averager.average(SINR_slot2,15,true);
%                 SINR_av1 = CQI_SNR_table(max(floor(CQI_mapper.SINR_to_CQI(SINR_av1)),1));
%                 SINR_av2 = CQI_SNR_table(max(floor(CQI_mapper.SINR_to_CQI(SINR_av2)),1));
%                 SINRs_to_map_to_CQI(l_i,rb_i,r_i) = SINR_averager.average([SINR_av1,SINR_av2],15,true);
%             end
%         end
%     end
% else
% size(SINR_temp)
    for rb_i = 1:nRB % RBs
        for r_i = 1:size_vec(4) % rank_possibilities
            for l_i = 1:r_i % layers  
                SINRs_to_map_to_CQI(l_i,rb_i,r_i) = SINR_averager.average(squeeze(SINR_temp(l_i,:,rb_i,r_i)),15,true); % MCS 15 has beta ~= 1 --> averaging using the BICM capacity without extra calibration
            end
        end
    end
% end

max_rank = size(SINRs_to_map_to_CQI,3);
if (tx_mode==3) || (tx_mode==4) || (tx_mode==9) %|| (tx_mode==5)||(tx_mode==6) % Rank decision for SM
    MCSs_all = 1:15;
    Is_MCSs  = SINR_averager.SINR_to_I(SINRs_to_map_to_CQI,MCSs_all);
    
    if max_rank==1
        Is_MCSs = reshape(Is_MCSs,[size(Is_MCSs,1) size(Is_MCSs,2) 1 size(Is_MCSs,3)]);
    end
    
    % Compute the per-layer mutual-information sum
    Is_dims                         = size(Is_MCSs);
    Is_MCSs_no_nans                 = zeros(size(Is_MCSs));
    Is_finite_idxs                  = isfinite(Is_MCSs);
    Is_MCSs_no_nans(Is_finite_idxs) = Is_MCSs(Is_finite_idxs);
    Is_sum_MCSs_per_layer           = reshape(sum(Is_MCSs_no_nans,1),Is_dims(2:end));
    
    % Optional: take only the N best values, as measured by the last MCS (if not one would need to calcualte it for every rank-MCS pair, making it too computationally costly!)
    Is_mean_MCSs_per_rank = zeros(length(MCSs_all),max_rank);
    last_MCS_Is_sum       = Is_sum_MCSs_per_layer(:,:,MCSs_all(end));
    
    %% Calculate mean MI value for each rank and MCS pair based on the best N MI values for each rank
    for r_=1:max_rank
        multiplier_matrix = kron(1:max_rank,ones(length(MCSs_all),1));
        if adaptive_RI==1 && ~isempty(DL_signaling.adaptive_RI) && ~isempty(DL_signaling.adaptive_RI.avg_MI) && ~isempty(DL_signaling.adaptive_RI.min_MI)
            % Take only the RBs with similar spectral efficiency as previously scheduled (measured in BICM capacity)
            spectral_eff_threshold = DL_signaling.adaptive_RI.min_MI;
            
            number_of_RBs_to_take   = DL_signaling.num_assigned_RBs;
            if number_of_RBs_to_take==0
                number_of_RBs_to_take = nRB;
            end
            [sort_values,sort_idxs]    = sort(last_MCS_Is_sum(:,r_));
            bigger_than_threshold      = sort_values>=spectral_eff_threshold;
            bigger_than_threshold_idxs = sort_idxs(bigger_than_threshold);
            begin_idx                  = max(length(bigger_than_threshold_idxs)-number_of_RBs_to_take+1,1);
            end_idx                    = length(bigger_than_threshold_idxs);
            RBs_to_average_idx_rank    = bigger_than_threshold_idxs(begin_idx:end_idx);
            
            RBs_to_average = Is_sum_MCSs_per_layer(RBs_to_average_idx_rank,r_,:);
            
            if isempty(RBs_to_average)
                RBs_to_average = Is_sum_MCSs_per_layer(:,r_,:); % All values
            end
        elseif adaptive_RI==2 && ~isempty(DL_signaling.adaptive_RI) && ~isempty(DL_signaling.adaptive_RI.RBs_for_feedback)
            RBs_to_average = Is_sum_MCSs_per_layer(DL_signaling.adaptive_RI.RBs_for_feedback,r_,:); % The indicated values
        else
            RBs_to_average = Is_sum_MCSs_per_layer(:,r_,:);
        end
        
        Is_mean_MCSs_per_rank(:,r_) = reshape(mean(RBs_to_average,1),[Is_dims(end) 1]) ./ multiplier_matrix(:,r_);
    end
    
    %% Rank Indicator: Decide based on the number of transmitted data bits for a rank value
    SINR_av_dB_for_RI     = SINR_averager.average_for_RI(Is_mean_MCSs_per_rank,1:15);
    CQI_temp_all          = floor(CQI_mapper.SINR_to_CQI(SINR_av_dB_for_RI));
    all_CQIs              = reshape(1:15,[],1);
    all_CQIs              = all_CQIs(:,ones(1,max_rank));
    temp_var              = CQI_temp_all-all_CQIs;
    temp_var(temp_var<0)  = Inf;
    [~, CQI_layer_all]    = min(temp_var);
    out_of_range          = CQI_layer_all<1;
    CQI_layer_all(out_of_range) = 1;
    bits_layer_config = (1:max_rank).*(8*round(1/8*[config.CQI_params(CQI_layer_all).modulation_order] .* [config.CQI_params(CQI_layer_all).coding_rate_x_1024]/1024 * config.sym_per_RB_nosync * config.N_RB*2)-24);
    bits_layer_config(out_of_range) = 0;
    [~,optimum_rank] = max(bits_layer_config); % Choose the RI for which the number of bits is maximized
    
    %% Calculate CQI feedback on a per-codeword basis
    
    % CQI reporting Layer mappings according to TS 36.211
    switch optimum_rank
        case 1
            SINRs_to_CQI_CWs = SINRs_to_map_to_CQI(1,:,1);
        case 2
            SINRs_to_CQI_CWs = SINRs_to_map_to_CQI(1:2,:,2);
        case 3
            % Manually set to two Codewords. Layer-to-codeword mapping according to TS 36.211 and done with the last CQI
            codeword2_SINRs_dB_avg = SINR_averager.average_codeword(Is_MCSs(2:3,:,optimum_rank,MCSs_all(end)),MCSs_all(end));
            SINRs_to_CQI_CWs       = [SINRs_to_map_to_CQI(1,:,3); codeword2_SINRs_dB_avg];
        case 4
            % Manually set to two Codewords. Layer-to-codeword mapping according to TS 36.211
            codeword1_SINRs_dB_avg = SINR_averager.average_codeword(Is_MCSs(1:2,:,optimum_rank,MCSs_all(end)),MCSs_all(end));
            codeword2_SINRs_dB_avg = SINR_averager.average_codeword(Is_MCSs(3:4,:,optimum_rank,MCSs_all(end)),MCSs_all(end));
            SINRs_to_CQI_CWs       = [codeword1_SINRs_dB_avg; codeword2_SINRs_dB_avg];
        case 5
            % Manually set to two Codewords. Layer-to-codeword mapping according to TS 36.211
            codeword1_SINRs_dB_avg = SINR_averager.average_codeword(Is_MCSs(1:2,:,optimum_rank,MCSs_all(end)),MCSs_all(end));
            codeword2_SINRs_dB_avg = SINR_averager.average_codeword(Is_MCSs(3:5,:,optimum_rank,MCSs_all(end)),MCSs_all(end));
            SINRs_to_CQI_CWs       = [codeword1_SINRs_dB_avg; codeword2_SINRs_dB_avg];
        case 6
            % Manually set to two Codewords. Layer-to-codeword mapping according to TS 36.211
            codeword1_SINRs_dB_avg = SINR_averager.average_codeword(Is_MCSs(1:3,:,optimum_rank,MCSs_all(end)),MCSs_all(end));
            codeword2_SINRs_dB_avg = SINR_averager.average_codeword(Is_MCSs(4:6,:,optimum_rank,MCSs_all(end)),MCSs_all(end));
            SINRs_to_CQI_CWs       = [codeword1_SINRs_dB_avg; codeword2_SINRs_dB_avg];
        case 7
            % Manually set to two Codewords. Layer-to-codeword mapping according to TS 36.211
            codeword1_SINRs_dB_avg = SINR_averager.average_codeword(Is_MCSs(1:3,:,optimum_rank,MCSs_all(end)),MCSs_all(end));
            codeword2_SINRs_dB_avg = SINR_averager.average_codeword(Is_MCSs(4:7,:,optimum_rank,MCSs_all(end)),MCSs_all(end));
            SINRs_to_CQI_CWs       = [codeword1_SINRs_dB_avg; codeword2_SINRs_dB_avg];
        case 8
            % Manually set to two Codewords. Layer-to-codeword mapping according to TS 36.211
            codeword1_SINRs_dB_avg = SINR_averager.average_codeword(Is_MCSs(1:4,:,optimum_rank,MCSs_all(end)),MCSs_all(end));
            codeword2_SINRs_dB_avg = SINR_averager.average_codeword(Is_MCSs(5:8,:,optimum_rank,MCSs_all(end)),MCSs_all(end));
            SINRs_to_CQI_CWs       = [codeword1_SINRs_dB_avg; codeword2_SINRs_dB_avg];
    end
    
    feedback.RI  = optimum_rank;
else
    feedback.RI  = 1;
    SINRs_to_CQI_CWs = SINRs_to_map_to_CQI; % I have to check whether this also holds for TxD because of the matrix dimensions
end

if (tx_mode==4 || tx_mode==9) && ~isempty(feedback_PMIs)
    feedback.PMI = feedback_PMIs.max(:,optimum_rank);
elseif (tx_mode==5 || tx_mode==6) && ~isempty(feedback_PMIs)
    feedback.PMI = feedback_PMIs.max(:,1);
    feedback.nonstandard = feedback_PMIs.nonstandard;
else
    feedback.PMI = nan(nRB,1);
end

% Send as feedback the CQI for each RB.
% Flooring the CQI provides much better results than
% rounding it, as by rounding it to a higher CQI you will
% very easily jump the BLER to 1. The other way around it
% will jump to 0.

if unquantized_CQI_feedback
    CQIs = CQI_mapper.SINR_to_CQI(SINRs_to_CQI_CWs);
else
    CQIs = floor(CQI_mapper.SINR_to_CQI(SINRs_to_CQI_CWs));
end

feedback.CQI     = CQIs;
feedback.tx_mode = tx_mode;
% CQIs
end