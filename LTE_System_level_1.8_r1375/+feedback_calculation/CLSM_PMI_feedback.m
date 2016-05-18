function [ PMI_feedback,TB_SINR_feedback_dB ] = CLSM_PMI_feedback(...
    precoding_codebook,...
    H_0, TX_power_RB_0, noise_W_half_RB,...
    H_i, TX_W_RB_i, W_i, ICI)
% Calculate CLSM feedback
% (c) Josep Colom Ikuno, INTHFT, 2013

% Calculate optimum precoder (SU-MIMO, no coordination)
PMI_feedback.max = feedback_calculation.Rel8_CLSM_PMI(precoding_codebook,H_0,noise_W_half_RB + mean(mean(ICI)));
% Calculate now the TB SINR for each layer possibility and that's the feedback
max_rank = min([size(H_0,1) size(H_0,2) max([precoding_codebook.nLayers])]);
TB_SINR_feedback_lin = NaN([max_rank,size(H_0,4),size(H_0,3),max_rank]);

for r_=1:max_rank
    PMI_0 = PMI_feedback.max(:,r_);
    PMI_0 = reshape([PMI_0, PMI_0].',[],1);
    W_0   = precoding_codebook(r_).W(:,:,PMI_0);
    
    TB_SINR_feedback_lin(1:r_,:,:,r_) = phy_modeling.post_equalization_SINR(...
        sum(TX_power_RB_0, 1),...
        H_0, W_0,...
        noise_W_half_RB,...
        TX_W_RB_i,...
        H_i, W_i, ICI);
end

TB_SINR_feedback_dB = 10*log10(TB_SINR_feedback_lin);

end

