function v=vel_vector(intensity,direction)
%direction is defined according to the new model WIM-3D

[X,Y] = pol2cart(pi/2-direction,intensity);
v=[X;Y;0];


