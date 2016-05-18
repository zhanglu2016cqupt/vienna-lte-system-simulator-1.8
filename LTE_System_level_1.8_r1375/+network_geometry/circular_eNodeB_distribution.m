function [ eNodeB_positions ] = LTE_init_create_circular_eNodeB_distribution(LTE_config)
% Generate circular eNodeB distribution
% (c) Martin Taranetz, ITC, 2013

nBSs    = LTE_config.eNodeB_circle_nBSs;
R       = LTE_config.eNodeB_circle_radius;

angle_positions  = [cos(2*pi*(1:nBSs)'/nBSs) sin(2*pi*(1:nBSs)'/nBSs)];

eNodeB_positions = R*angle_positions;
