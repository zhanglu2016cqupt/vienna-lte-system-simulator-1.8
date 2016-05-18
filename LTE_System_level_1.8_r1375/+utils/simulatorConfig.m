classdef simulatorConfig
    % This class defines a configuration class, which takes as input the
    % configuration struct of the simulator and adds certain parameters to
    % it.
    % (c) Josep Colom Ikuno, INTHFT, 2012
    
    properties
    end
    
    methods (Static,Abstract)
        LTE_config = apply_parameters(LTE_config)
    end
end

