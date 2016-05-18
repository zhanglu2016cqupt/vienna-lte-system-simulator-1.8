% a_old is the vector that is to be rotated
% b contains the rotation angles b=[RotX;RotY,RotZ]
% rotation is done in the following order
% 1. RotX 2. RotY 3. RotZ
%
% Author: Martin Käske (TUI)
function [a_new,varargout]=rotate_vector(a_old,b)
    % T=Tz*Ty*Tx
    % since there are a lot of zeros in the actual rotation matrices
    % it is useful to use the "compact" form, i.e. T=Tx*Ty*Tz
    % complexity is reduced from 2*27+9 multiplications
    % to 16+9 multiplications
	T=[cos(b(2))*cos(b(3)) (-cos(b(1))*sin(b(3))+sin(b(1))*sin(b(2))*cos(b(3))) (sin(b(1))*sin(b(3))+cos(b(1))*sin(b(2))*cos(b(3)));
       cos(b(2))*sin(b(3)) (cos(b(1))*cos(b(3))+prod(sin(b))) (-sin(b(1))*cos(b(3))+cos(b(1))*sin(b(2))*sin(b(3)));
       -sin(b(2))           sin(b(1))*cos(b(2))               cos(b(1))*cos(b(2))];
    
    a_new=T*a_old;
    if(nargout>1)
        varargout{1}=T;
    end;