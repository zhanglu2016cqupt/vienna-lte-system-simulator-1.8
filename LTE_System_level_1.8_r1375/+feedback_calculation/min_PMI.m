function min_PMI = min_PMI(precoding,H_0, max_PMI)
% Calculate the optimum precoder according to the passed Rel'8 (or other) codebook
% (c) Josep Colom Ikuno, jcolom@nt.tuwien.ac.at ITC, 2013
% www.nt.tuwien.ac.atnTX        
nTX        = size(H_0,2);
nRX        = size(H_0,1);
max_layers = min([nTX nRX precoding(end).nLayers]);

% Calculate the average H from each RB (needed for the PMI calculation)
H_0_RB = (H_0(:,:,1:2:end)+H_0(:,:,2:2:end))/2;
N_rb   = size(H_0_RB,3);
min_PMI    = zeros(N_rb,max_layers); % PMI for each layer choice

% Call PMI-calculation function (CLSM)
for layer_idx = 1:max_layers
    % As it is built in the struct creation, it is sure that precoding.W is always complex
    min_PMI(:,layer_idx) = feedback_calculation.min_PMI_layer_specific(H_0_RB,layer_idx,precoding(layer_idx).W,nTX,'ZF',false, max_PMI);
end
end
