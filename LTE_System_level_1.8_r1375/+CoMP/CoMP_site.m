classdef CoMP_site < handle
    % (c) Thomas Blazek, 2014, ITC, thomas.blazek@nt.tuwien.ac.at
    
    properties
        id
        attached_eNodeBs
        cooperation_rules
        scheduler
    end
    
    methods
        function obj = CoMP_site(LTE_config, eNodeBs, id)
            obj.id = id;
            obj.attached_eNodeBs = eNodeBs;
            for i_ = 1: length(eNodeBs)
                eNodeBs(i_).attached_CoMP_site = obj;
            end
            %obj.cooperation_rules = cell(length(eNodeBs), 1);
        end
        
        % Example for a feedback implementation
        function feedback_wideband_SIR(obj, UE, eNodeB, SIR, i_eNodeBs)
            SIR = SIR(eNodeB ~= i_eNodeBs);
            rules.dont_schedule =[];
            
            check_matrix = SIR < 5;
            rules.UE_id = UE.id;
            rules.dont_schedule = i_eNodeBs ( check_matrix);
            coop_rules = obj.cooperation_rules;
            if ~isempty(coop_rules) && sum([coop_rules.UE_id] == UE.id)>0
                obj.cooperation_rules([coop_rules.UE_id] == UE.id).dont_schedule = rules.dont_schedule;
            else
                obj.cooperation_rules = [obj.cooperation_rules; rules];
            end
            
        end
        
        function forbidden_UEs =  get_forbidden_UEs(obj, eNodeB)
            if~isempty(obj.cooperation_rules)
                for k_ = 1:length(obj.cooperation_rules)
                    if ~isempty(obj.cooperation_rules(k_).dont_schedule)
                        forbidden_UEs(k_) = logical(sum(obj.cooperation_rules(k_).dont_schedule==eNodeB));
                        
                    else
                        forbidden_UEs(k_) = 0;
                    end
                end
                forbidden_UEs = vertcat(obj.cooperation_rules(logical(forbidden_UEs)).UE_id);
            else
                forbidden_UEs = 0;
            end
        end
    end
    
end

