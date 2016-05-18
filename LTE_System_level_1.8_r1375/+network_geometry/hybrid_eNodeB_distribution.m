function [ eNodeB_positions ] = LTE_init_create_hybrid_eNodeB_distribution(LTE_config)
% Generate eNodeB distribution with stochastic and deterministic parts.
% (c) Martin Taranetz, ITC, 2012

% Radius of the typical cell; assumed to be half of the typical inter
% eNodeB distance
Rc           = LTE_config.average_eNodeB_distance/2;           
% Multiplication factor which determines the size of the ROI corresponding
% to the size of the typical cell
network_size = LTE_config.network_size;
% Area of the region of interest
A_roi        = (network_size * LTE_config.average_eNodeB_distance)^2;
% From [Modeling Heterogeneous Network Interference Using Poisson Point Processes, Heath, 2013.]
BS_density   = 1/(16*Rc^2);
% Determine the number of eNodeBs to generate
number_of_eNodeBs_to_generate = max(round(A_roi * BS_density),1); % Make sure that at least one eNodeB is generated

% Generate the positions of the eNodeBs (vector of [x,y] values)
normed_random_eNodeB_positions = rand(number_of_eNodeBs_to_generate, 2)-0.5;
% Position of eNodeBs in m (The factor two accounts for the random
% value range [-0.5 0.5]
random_eNodeB_positions        = normed_random_eNodeB_positions*2*network_size*Rc;

% Exclusion of eNodeBs within circle of Radius (Rc + Rg) around center
% Guard Radius
Rg = Rc;
% distance of eNodeBs to center of scenario
random_eNodeB_distances        = sqrt(random_eNodeB_positions(:,1).^2 + random_eNodeB_positions(:,1).^2);
% Index of eNodeBs which are outside the exclusion region
eNodeBs_to_keep                = random_eNodeB_distances >=(Rg + Rc);
eNodeB_positions_               = random_eNodeB_positions(eNodeBs_to_keep,:);
% Add strongest interferer at edge exclusion region
% Its position on the circle of the exclusion region is random
if isfield(LTE_config,'generate_strongest_interferer')
   if LTE_config.generate_strongest_interferer 
       strongest_interferer_angle    = rand()*2*pi;
       strongest_interferer_position = (Rc+Rg)*[cos(strongest_interferer_angle) sin(strongest_interferer_angle)];
       eNodeB_positions_        = [strongest_interferer_position; eNodeB_positions_];
   end
end

% Add center eNodeB
eNodeB_positions               = [0,0; eNodeB_positions_];
