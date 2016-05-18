classdef walkingModel < handle
    % Defines a walking model
    % (c) Josep Colom Ikuno, INTHFT, 2008
    
    methods (Abstract)
        % Based on the current position, outputs the next TTI's position
        new_pos = move(obj,current_pos)
        % Based on the current position, outputs the last TTI's position
        old_pos = move_back(obj,current_pos)
        % Print some info
        print(obj)
    end
end
