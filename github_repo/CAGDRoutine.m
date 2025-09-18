classdef CAGDRoutine
    
    methods (Static)
        
        %------------------------------------------------------------------
        
        function [ct] = DeCasteljou(P, t)
            
            n = size(P,2);
            for i=1:n
                for j=1:n-i
                    P(:,j) = (1-t)*P(:,j) + t*P(:,j+1);
                end
            end
            ct = P(:,1);
            
        end
        
        %------------------------------------------------------------------
        %plotFlag   - [1|0] plot or no the curve in the figure(idFigure)
        function [X, plotCurve] = plotCurve_equiT(plotFlag, P, nT, idFigure, curveColor, lineWidth)
            
            X = [];
            t = 0;
            tIncr = 1/nT;
            
            for i=1:nT
                X(:,i) = DeCasteljou(P,t);
                t = i*tIncr;
            end
            
            plotCurve = [];
            if(plotFlag)
                figure(idFigure);
                hold on;
                plotCurve = plot(X(1,:),X(2,:),curveColor,'LineWidth',lineWidth);
                hold off;
            end
            
        end
        
        %------------------------------------------------------------------
    end
    
end

