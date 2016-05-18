% Calculates the post-equalization SINR
function [ TB_SINR_feedback_lin, TX_W_half_RB_resized, TX_W_half_RB_i ] = post_equalization_SINR(...
    TX_power_RB,...         % Desired tx power (row vector)
    H_0, W_0,...            % Channel and precoder of signal
    noise_W_half_RB,...     % Noise per half RB
    TX_W_RB_i,...           % Tx power of interferers (each interferer equals one column)
    H_i, W_i,ICI_power)               % Channel and precoder of interferers

there_are_interferers = ~(isempty(TX_W_RB_i)||isempty(H_i)||isempty(W_i));

% Avoids complexity problems with the C version
if isreal(W_0)
    W_0 = complex(W_0);
end

% Calculate receiver filter
% ZF-Filter:                P    ( = (HW'*HW) \ HW' )
% Gain of desired signal:   zeta (all ones by design because of ZF)
% Self interference:        xi   (all zeros by design since different RBs are orthogonal)
% Noise enhancement:        psi  (= sum(abs(P).^2,2))
try
    [P,HW_0,zeta,xi,psi] = network_elements.UE.calculate_effective_channel_and_receiver_mex(H_0,W_0);
catch err
    [P,HW_0,zeta,xi,psi] = network_elements.UE.calculate_effective_channel_and_receiver(H_0,W_0);
end

% Interferers (if present)
% theta:                   similar to psi, but with sum(abs(P*H_iW_i).^2, 2)
if there_are_interferers

   try
        theta = network_elements.UE.calculate_per_layer_interference_power_mex( P, H_i, W_i );
   catch err
        theta = network_elements.UE.calculate_per_layer_interference_power( P, H_i, W_i );
   end
    
    % Power allocation per half-RB
    TX_W_half_RB_i = zeros(size(theta));
    TX_W_RB_i = permute(repmat(TX_W_RB_i,[1,1,4]),[4,1,3,2]);
    for i_=1:size(theta,1)
        %TX_W_half_RB_i(i_, :, :) = kron(TX_W_RB_i, [1; 1])/(2); 
        %double rows and half power (RB -> half_RB)        
        TX_W_half_RB_i(i_,1:2:end,:, :) = TX_W_RB_i/(2);
        TX_W_half_RB_i(i_,2:2:end,:, :) = TX_W_RB_i/(2);
    end
else
    theta = 0;
    TX_W_half_RB_i = 0;
end


if iscolumn(TX_power_RB)
    TX_power_RB = TX_power_RB';
end
TX_W_half_RB         = kron(sum(TX_power_RB, 1), [1 1])/(2); % double columns and half power (RB->half_RB)
TX_W_half_RB_resized = repmat(TX_W_half_RB(ones(1,size(W_0,2)),:),[1,1,size(zeta,3)]);

ICI_power = permute(repmat(ICI_power,[1,1,size(psi,1)]),[3,2,1]);
% size(ICI_power)
% size(psi)
TB_SINR_feedback_lin = zeta.*TX_W_half_RB_resized ./ (xi.*TX_W_half_RB_resized + psi.*noise_W_half_RB + psi.*ICI_power + sum(theta.*TX_W_half_RB_i,4));
TB_SINR_feedback_lin = permute(TB_SINR_feedback_lin,[1,3,2,4]); 
end
