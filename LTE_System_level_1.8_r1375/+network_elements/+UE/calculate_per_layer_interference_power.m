function [ interf_enhancement ] = calculate_per_layer_interference_power(P, H_i, precoders_i ) %#codegen
% Unfortunately, we need these loops here. From each interferer, each RB
% could be allocated to a different UE, thus with a different layer and
% precoder

nSC  = size(H_i,3);
nSY = size(H_i,4);
nInt = size(H_i,5);
nLayers = size(P,1);
nUEs = size(precoders_i, 3); %number of different UEs that can have different precoders at interferer
interf_enhancement = zeros(nLayers,nSC,nSY,nInt, nUEs);

for SC_idx = 1:nSC
    for SY_idx = 1:nSY  
        for I_idx = 1:nInt
            for UE_idx = 1:nUEs
                % PH_interf = P(:,:,SC_idx) * H_i(:,:,SC_idx,I_idx) * precoding_codebook(rank_i(SC_idx,I_idx)).W(:,:,PMI_i(SC_idx,I_idx));
                interf_enhancement(:,SC_idx,SY_idx,I_idx, UE_idx) = sum(abs(P(:,:,SC_idx,SY_idx) * H_i(:,:,SC_idx,SY_idx,I_idx) * precoders_i(SC_idx,I_idx, UE_idx).W).^2,2); % diag(PH_interf*PH_interf');
            end
        end
    end
end
end
