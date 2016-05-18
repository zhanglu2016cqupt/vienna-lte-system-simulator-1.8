classdef enodebTrace < handle
    % This class stores, for each eNodeB's sector the traces that we wanto to store. eg, CQI assignments,
    % throughput, etc, etc.
    % (c) Josep Colom Ikuno, INTHFT, 2008-2012
    
    properties
        clock                  % Tells the object in what TTI we are now
        user_allocation        % whose slot is each
        sent_data              % Sent data this TTI (bits)
        acknowledged_data      % Data that was acknowledged (bits)
        scheduled_RBs          % Number of RBs scheduled this TTI
        RB_grid_size           % Size of the RB grid in RBs
        
        expected_ACKs          % Total number of TBs sent
        received_ACKs          % The ones that were ACKs
        
        % Aggregates
        average_BLER           % Final average cell BLER during the simulation
        
        parent_results_object % The parent results object
    end
    
    methods
        function obj = enodebTrace(the_RB_grid,maxStreams,simulation_length_TTI,clock)
            %time_slots_power_allocation = size(RB_grid_object.power_allocation,1);
            n_RB                  = the_RB_grid.n_RB;
            if isprop(the_RB_grid, 'max_ues')
                obj.user_allocation   = zeros(n_RB,the_RB_grid.max_ues, simulation_length_TTI,'uint16');
                obj.RB_grid_size      = n_RB*the_RB_grid.max_ues;
                obj.scheduled_RBs     = zeros(simulation_length_TTI,1, 'uint16');
            else
                
                obj.user_allocation   = zeros(n_RB, simulation_length_TTI,'uint16');
                obj.RB_grid_size      = n_RB;
                obj.scheduled_RBs     = zeros(1,simulation_length_TTI,'uint16');
            end
            obj.sent_data         = zeros(maxStreams,simulation_length_TTI,'uint32');
            obj.acknowledged_data = zeros(maxStreams,simulation_length_TTI,'uint32');
            obj.expected_ACKs     = zeros(maxStreams,simulation_length_TTI,'uint16');
            obj.received_ACKs     = zeros(maxStreams,simulation_length_TTI,'uint16');
            obj.clock             = clock;
        end
        % Stores traces after the scheduling was done (Rb grid size)
        function store_after_scheduling(obj,the_RB_grid)
            TTI_idx = obj.clock.current_TTI;
            if iscolumn( the_RB_grid.user_allocation)
                
                obj.user_allocation(:,TTI_idx)      = the_RB_grid.user_allocation;
                obj.scheduled_RBs(TTI_idx)          = sum(the_RB_grid.user_allocation~=0);
            else
                obj.user_allocation(:, :,TTI_idx)      = the_RB_grid.user_allocation;
                obj.scheduled_RBs(TTI_idx)          = sum(sum((the_RB_grid.user_allocation~=0), 1));
            end
            nCodewords                          = length(the_RB_grid.size_bits);
            obj.sent_data(1:nCodewords,TTI_idx) = the_RB_grid.size_bits;
        end
        % Add the trace from a feedback report
        function store_ACK_report(obj,feedback)
            nCodewords   = feedback.nCodewords;
            UE_scheduled = feedback.UE_scheduled;
            ACK          = feedback.ACK;
            TB_size      = feedback.TB_size;
            trace_TTI    = feedback.TTI_idx;
            % In case the ACK report is incorrectly filled. The most important check is whether the user was scheduled or not.
            checked_TB_size = uint32(UE_scheduled .* ACK(:) .* TB_size(:));
            obj.acknowledged_data(1:nCodewords,trace_TTI) = obj.acknowledged_data(1:nCodewords,trace_TTI) + sum(checked_TB_size);
            max_streams = size(obj.expected_ACKs,1);
            if UE_scheduled
                filler_ACK_zeros = zeros(max_streams-nCodewords,1,'uint16');
                current_expected_ACKs = [ones(length(ACK),1,'uint16'); filler_ACK_zeros];
                current_ACKs          = [uint16(ACK(:)); filler_ACK_zeros];
                obj.expected_ACKs(:,trace_TTI) = obj.expected_ACKs(:,trace_TTI) + current_expected_ACKs(:);
                obj.received_ACKs(:,trace_TTI) = obj.received_ACKs(:,trace_TTI) + current_ACKs(:);
            end
        end
        % Calculate the cell average BLER
        function calculate_final_average_BLER(obj)
            TTIs_to_ignore                                  = obj.parent_results_object.TTIs_to_ignore_when_calculating_aggregates;
            TTIs_to_account_for                             = true(1,length(obj.scheduled_RBs));
            TTIs_to_account_for(1:TTIs_to_ignore)           = false; % Ignore TTIs where no feedback information was available or dummy data was sent
            TTIs_to_account_for((end-TTIs_to_ignore+1):end) = false;
            expected_ACKs     = obj.expected_ACKs(:,TTIs_to_account_for);
            received_ACKs     = obj.received_ACKs(:,TTIs_to_account_for);
            average_BLER.BLER = 1-(sum(received_ACKs(:))/sum(expected_ACKs(:)));
        end
    end
end
