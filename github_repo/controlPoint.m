classdef controlPoint < handle
    
    properties
        id          %string
        cp          %plot line
        cp_prec     %plot line
        cp_succ     %plot line
        left_leg    %plot line
        right_leg   %plot line
    end
    
    methods
        
        %Costruttore
        % 
        function newCP = controlPoint(name, plotCP, previousCP, idFigure)
            newCP.id = name;
            if(isa(previousCP,'controlPoint'))
                newCP.cp_prec = previousCP;
                previousCP.cp_succ = newCP;
                
                x1 = get(previousCP.cp,'XData');
                y1 = get(previousCP.cp,'YData');
                x2 = get(plotCP,'XData');
                y2 = get(plotCP,'YData');
                
                %creo il nuovo control polygon leg
                figure(idFigure);
                hold on;
                leg = plot([x1 x2],[y1 y2],'-b');
                
                newCP.left_leg = leg;
                previousCP.right_leg = leg;
                
            end
            newCP.cp = plotCP;
        end
        
        function [xCoo, yCoo] = getCurrentCoo(this)
            
            xCoo = get(this.cp, 'XData');
            yCoo = get(this.cp, 'YData');
            
        end
        
    end
    
end

