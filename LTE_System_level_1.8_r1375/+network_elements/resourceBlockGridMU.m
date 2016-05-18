classdef resourceBlockGridMU < network_elements.resourceBlockGrid
    % This is an extension of the resource block grid to the multi user
    % case. Instead of one UE per RB, multiple UEs per RB can be scheduled.
    % In this case, the Limit is set to the nTX antennas, as this is the
    % max number of linearly independent precoders that can be found. The
    % grid itself is extended to a matrix.
    properties
        max_ues
    end
    
    methods
        function obj = resourceBlockGridMU(n_RB,sym_per_RB_nosync,sym_per_RB_sync, nTX)
           obj = obj@network_elements.resourceBlockGrid(n_RB,sym_per_RB_nosync,sym_per_RB_sync);
           max_codewords = 2;
           obj.max_ues = nTX;
           max_ues = nTX;
           % "Dynamic" information
           obj.user_allocation  = zeros(n_RB,max_ues,'uint16'); % We will assume that the streams cannot be scheduled to different UEs.
           obj.power_allocation = zeros(n_RB,max_ues); % TTI-based power allocation. No slot-based power allocation
           obj.power_allocation_signaling = zeros(n_RB,max_ues); % Simplification: equal for each RB
           obj.size_bits         = zeros(1,max_codewords);
           
           precoder_struct       = struct('W',[]);
           obj.precoder          = precoder_struct(ones(max_ues,n_RB));
           
           % "Static" information
           obj.n_RB              = n_RB;
           obj.sym_per_RB_nosync = sym_per_RB_nosync;
           obj.sym_per_RB_sync   = sym_per_RB_sync;
           obj.max_ues;
       end
       
    end
    
end

