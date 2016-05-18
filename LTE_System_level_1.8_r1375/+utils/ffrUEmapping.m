classdef ffrUEmapping < handle
    % Maps the UEs to either the FR or the PR part of the BW.
    % Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at.
    % (c) 2011 by INTHFT
    % www.nt.tuwien.ac.at
    
    properties
        % Technically it could support overlapping assignments, but because
        % of the implementation of the constructor, this is not possible
        FR_assignment
        PR_assignment
        N_FR
        N_PR
    end
    
    methods
        function obj = ffrUEmapping(UE_mapping)
            obj.FR_assignment = (UE_mapping==1);
            obj.PR_assignment = (UE_mapping==2);
            obj.N_FR          = sum(obj.FR_assignment);
            obj.N_PR          = sum(obj.PR_assignment);
        end
        
        function [FR PR] = get_assignment(obj,UE_id)
            FR = obj.FR_assignment(UE_id);
            PR = obj.PR_assignment(UE_id);
        end
    end
end

