function [PMI] = Rel8_CLSM_PMI_layer_specific(H_t_precalculated,nLayers,W_all,nAtPort,receiver,wideband_precoder, max_PMI)
% Calculates the PMI based on the given codebook for the LTE feedback
% Author:  Stefan Schwarz, Josep Colom Ikuno (block diagonalization)
% Contact: stefan.schwarz@nt.tuwien.ac.at, jcolom@nt.tuwien.ac.at

%% Some initializations
switch nLayers % precoder indices as in the standard (zero-indexed)
    case 1
        if nAtPort == 2
           W_idxs = 0:3;
        elseif nAtPort == 4 
           W_idxs = 0:15;
        else
           W_idxs = 0:255;
        end
    case 2
        if nAtPort == 2
           W_idxs = 0:1;
        elseif nAtPort == 4 
           W_idxs = 0:15;
        else
           W_idxs = 0:255;
        end
    case 3
        if nAtPort == 4
            W_idxs = 0:15;
        else
            W_idxs = 0:63;
        end
    case 4
        if nAtPort == 4
            W_idxs = 0:15;
        else
            W_idxs = 0:31;
        end
    case {5 6 7}
        W_idxs = 0:3;
    case 8
        W_idxs = 0;
        
    otherwise
        error('Number of layers not supported');
end

% Take the subset of precoders we will use
W_all = W_all(:,:,W_idxs+1);
if isreal(W_all)
    W_all = complex(W_all);
end

% Generate block matrix version of the precoders
N_W         = size(W_all,3);
N_TX        = size(W_all,1);
W_all_block = complex(zeros(N_TX*N_W,nLayers*N_W));

N_TX_vect     = (1:N_TX);
N_layers_vect = (1:nLayers);

% Generate the block diagonal matrix containing the precoders in case we
% use the sparse-optimized function
for w_=1:size(W_all,3)
    W_cols = (w_-1)*nLayers + N_layers_vect;
    W_rows = (w_-1)*N_TX     + N_TX_vect;
    W_all_block(W_rows,W_cols) = W_all(:,:,w_);
end
W_all_block = sparse(W_all_block);

%try
    %I = feedback_calculation.codebook_capacity.non_sparse_mex(H_t_precalculated,W_all,receiver);
%catch err
    I = feedback_calculation.codebook_capacity.sparse_min(H_t_precalculated,W_all_block,N_W,nLayers,receiver, max_PMI);
%end

% Take into account that the PMI is 1-indexed, not 0-indexed, as in the standard!
if ~wideband_precoder
    [~,PMI] = min(I,[],2);
else
    [~,PMI] = min(sum(I,1));
end

PMI = PMI(:); % make column vector
