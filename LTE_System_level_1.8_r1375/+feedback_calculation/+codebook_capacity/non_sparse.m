function I = non_sparse( H_t_all, W_all, receiver ) %#codegen
% Calculate the capacity of a given channel matrix with each one of the
% given precoders in block diagonal matrix form. This function works with
% sparse matrices so as to increase speed. There is also a non-sparse mode
% to be (mainly) intended for C code generation

N_rb = size(H_t_all,3); % number of resource blocks
N_W  = size(W_all,3);
I    = zeros(N_rb,N_W);

% Calculate feedback for each RB
for RB_i = 1:N_rb
    H_t = H_t_all(:,:,RB_i);
    for w_=1:N_W
        % Effective channel matrix
        P = H_t*W_all(:,:,w_);
        
        % Right now only ZF receiver supported in SL, however
        if strcmp(receiver,'ZF')
            F = (P'*P) \ P'; % ZF receiver
        elseif strcmp(receiver,'MMSE')
            temp  = P'*P;
            F = (temp+0.01*eye(size(temp))) \ P'; % MMSE receiver assuming 0.01 noise variance
        else
            error('Not implemented');
        end
        
        K           = F*P;
        signal      = abs(diag(K)).^2;
        
        % Avoid using the diag function, as codegen had some serious problems with it, as of 2013a.
        for k_idx=1:min(size(K,1),size(K,2))
            K(k_idx,k_idx) = 0;
        end
        
        inter_layer = sum(abs(K).^2,2);
        noise       = 0.01.*sum(abs(F).^2,2);
        SNR         = signal ./ (inter_layer + noise);
        I(RB_i,w_)  = sum(log2(1+SNR),1); % rate
    end
end

end

