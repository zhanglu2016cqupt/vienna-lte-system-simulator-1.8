% unrotate_vector() assumes that a_old was
% rotated using a_old=rotate_vector(a_new,b)
% and reverses the rotation
function a_new=unrotate_vector(a_old,b)
	b=-1*b; % rotate in opposite direction
    % T=Tx*Ty*Tz
    % since there are a lot of zeros in the actual rotation matrices
    % it is useful to use the "compact" form, i.e. T=Tx*Ty*Tz
    % complexity is reduced from 2*27+9 multiplications
    % to 16+9 multiplications
    T=[cos(b(2))*cos(b(3)) -cos(b(2))*sin(b(3)) sin(b(2));
       (sin(b(1))*sin(b(2))*cos(b(3))+cos(b(1))*sin(b(3))) (-prod(sin(b))+cos(b(1))*cos(b(3))) -sin(b(1))*cos(b(2));
       (-cos(b(1))*sin(b(2))*cos(b(3))+sin(b(1))*sin(b(3))) (cos(b(1))*sin(b(2))*sin(b(3))+sin(b(1))*cos(b(3))) cos(b(1))*cos(b(2))]; 
	% this is how a_old was created a_old=Tz*Ty*Tx*a_new
	% so we have to undo this 
    a_new=T*a_old;

