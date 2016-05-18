function [P,HW,zeta,xi,psi] = calculate_effective_channel_and_receiver(H,W) %#codegen
HW   = complex(zeros([size(H,1) size(W,2) size(H,3) size(H,4)]));
P    = complex(zeros([size(HW,2) size(HW,1) size(HW,3) size(H,4)]));
zeta = ones([size(P,1) size(P,3) size(H,4)]);  % ZF receiver
xi   = zeros([size(P,1) size(P,3) size(H,4)]); % ZF receiver
psi  = zeros([size(P,1) size(P,3) size(H,4)]); % Noise enhancement

for i_=1:size(H,3)
    for j_ =1:size(H,4)
        HW_ = H(:,:,i_,j_)*W(:,:,i_);
        HW(:,:,i_,j_) = HW_;
        P_         = (HW_'*HW_) \ HW_';
        P(:,:,i_,j_)  = P_;
        psi(:,i_,j_) = sum(abs(P_).^2,2);
    end
end
end