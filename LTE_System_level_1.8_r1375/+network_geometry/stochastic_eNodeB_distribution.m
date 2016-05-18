function [ eNodeB_positions ] = LTE_init_create_stochastic_eNodeB_distribution(LTE_config)
% Generates eNodeB positions, which are homogenously distributed over the ROI.
% (c) Martin Taranetz, ITC, 2012

average_eNodeB_distance = LTE_config.average_eNodeB_distance;
network_size            = LTE_config.network_size;

% Equivalent area of hexagonal cell:
A_hex = 3*sqrt(3)/2 * (average_eNodeB_distance/sqrt(3))^2;
% Area of ROI:
% Square with sidelength = 2*network_size*average_eNodeB_distance
% The factor 2 regards both, negative and positive dimension
A_ROI = (2 * average_eNodeB_distance * network_size)^2;
% Determine number of eNodeBs to generate
number_of_eNodeBs_to_generate = max(floor(A_ROI / A_hex),1); % Make sure that at least one eNodeB is generated
% Generate the positions of the eNodeBs (vector of [x,y] values)
normed_random_eNodeB_positions = rand(number_of_eNodeBs_to_generate, 2)-0.5;
% Position of eNodeBs in m (The factor two accounts for the random
% value range [-0.5 0.5]
random_eNodeB_positions        = normed_random_eNodeB_positions*2*network_size*average_eNodeB_distance;

eNodeB_positions = random_eNodeB_positions;