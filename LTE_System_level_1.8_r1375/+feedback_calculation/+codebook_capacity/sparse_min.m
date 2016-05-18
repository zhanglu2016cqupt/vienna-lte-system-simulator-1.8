function I = sparse_min( H_t_all, W_all_block, N_W, nLayers, receiver, max_PMI )
% Calculate the capacity of a given channel matrix with each one of the
% given precoders in block diagonal matrix form. This function works with
% sparse matrices so as to increase speed. There is also a non-sparse mode
% to be (mainly) intended for C code generation

N_rb = size(H_t_all,3); % number of resource blocks
I    = zeros(N_rb,N_W);

N_TX        = size(H_t_all,2);
N_RX        = size(H_t_all,1);

N_TX_vect     = (1:N_TX);
N_RX_vect     = (1:N_RX);

% Calculate feedback for each RB
for RB_i = 1:N_rb
    
    H_t = H_t_all(:,:,RB_i);
    
    % Generate block diagonal matrix with the repeated current channel matrix
    H_t_block = complex(zeros(N_RX*N_W,N_TX*N_W));
%     for w_=1:N_W
%         H_rows = (w_-1)*N_RX     + N_RX_vect;
%         H_cols = (w_-1)*N_TX     + N_TX_vect;
%         H_t_block(H_rows,H_cols) = H_t;
%     end
    H_t_block = kron(eye(N_W), H_t);
    H_t_block = sparse(H_t_block);
    % Equalize with all of the possible precoders in one step, with sparsity and the mldivide operator (faster)
    P_all = H_t_block*W_all_block;
    
    % Right now only ZF receiver supported in SL, however
    if strcmp(receiver,'ZF')
        F_all = (P_all'*P_all) \ P_all'; % ZF receiver
    elseif strcmp(receiver,'MMSE')
        temp  = P_all'*P_all;
        F_all = (temp+0.01*eye(size(temp))) \ P_all'; % MMSE receiver assuming 0.01 noise variance
    else
        error('Not implemented');
    end
    
    K_all      = F_all*P_all;
    diag_K_all = diag(K_all);
    
    signal      = reshape(abs(diag_K_all).^2,nLayers,[]);
    inter_layer = reshape(sum(abs(K_all-diag(diag_K_all)).^2,2),nLayers,[]);
    noise       = reshape(0.01.*sum(abs(F_all).^2,2),nLayers,[]);
    SNR_all     = signal ./ (inter_layer + noise);
    I(RB_i,:)   = sum(log2(1+SNR_all),1); % rate of one resource block for different precoders and this specific rank (nLayers)
end
end

