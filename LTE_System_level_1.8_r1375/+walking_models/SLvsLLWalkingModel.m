classdef SLvsLLWalkingModel < walking_models.walkingModel
    % Defines a walking model that makes the UE walk in a straight line
    % (c) Josep Colom Ikuno, INTHFT, 2008
    
    properties
        direction
        av_numb
        av_count
        distances
        max_dist
        pos_count = 0;
    end
    methods
        % Class constructor
        function obj = SLvsLLWalkingModel
            % If angle is specified, use it, if not assign a random one
            obj.direction = 60;
            obj.av_numb = 100;
            obj.max_dist = 750;
            obj.av_count = obj.av_numb;
%             obj.distances = (linspace((15)^(1/3),(obj.max_dist)^(1/3),41)).^3;
%             obj.distances = [15,17,18,22,25,30,35,45,50,55,60,70,80,90,100,110,120,135,150,165,180,200,240,280,320,360,400,450,500,550,600,650,700,750];
            PL = 70:3:100;
            obj.distances = 10.^(PL/20)*299792458/(4*pi*2e9);
        end
        % Based on the current position, outputs the next TTI's position
        function new_pos = move(obj,current_pos)            
            if obj.av_count == obj.av_numb
                obj.pos_count = obj.pos_count+1;
                new_pos = obj.distances(obj.pos_count)* [ cosd(obj.direction) sind(obj.direction) ];
                obj.av_count = 0;
            else
                new_pos = current_pos;
            end
            obj.av_count = obj.av_count +1;
        end
        % Based on the current position, outputs the last TTI's position
        function old_pos = move_back(obj,current_pos)
                obj.pos_count = obj.pos_count-1;
                old_pos = obj.distances(obj.pos_count)* [ cosd(obj.direction) sind(obj.direction) ];
        end
        % Print some info
        function print(obj)
            fprintf('  direction: %d°, %d m/s*TTI\n',obj.direction,obj.speed);
        end
    end
end
