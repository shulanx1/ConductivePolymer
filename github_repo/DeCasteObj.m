%Traduzione da function -> Obj
% Le global variable diventano le variabile dell'oggetto
%
% argin degli eventi ->
%[objDellaFunzione(this), objEvento, windowEventData, axisHandler]
classdef DeCasteObj < handle
    
    properties
        %------------------- TYPE       - DESCRIPTION
        mainFigure = [];    %figureRef  - container figure
        curvFigure = [];
        dragging = [];      %plotRef    - moving point graph reference (by plot)
        orPos = [];         %[2x1]      - coordinate vector of original position
        tag_CPobj_map = ...
            containers.Map()%map        - dictionary string tag (graphic obj) to control point object related
        CP_obj = [];        %ObjList    - collection of control point objects
        X = [];             %[2 x #t-val]  - collection of point computed by DeCasteljou alg
        plotCurve = [];     %plotRef    - resulting curve, plot ref obj
        
        convexHull = [];    %plotRef    - graphic reference to current convex hull
        convexHullCkbx = [] %uiCheckBox - checkbox reference for display convex hull
    end
    
    methods
        %---------------- Builder -----------------------------------------
        
        function new = DeCasteObj(varargin)
            clf;
            close all;
            
            % default case if no control points are given
            if(isempty(varargin) == 1)
                %SAMPLE POINTS - DEFAULT CASE
                CP = [-3.5 -3 -2 -1 0 1 3 3.5;
                    0    2   1 -1 1 2 2.5 1];
            else
                CP = varargin{1};
            end
            
            %setting figure
            f = figure(1);
            set(f,'units','normalized')
            set(gcf,'color','w');
            set(f, 'Position', [0 0 0.7 1]);
            set(f,'Name','Bezier Curve','NumberTitle','off');
            ax = axes('xlim',[-4 4],'ylim',[-3 3]);
            new.mainFigure = f;
            
            %add ui element
            addControlPointBtn = uicontrol(f,'Style','pushbutton','String','Add Random Control Point',...
                'Position',[10 10 150 40],'Callback',@new.addCPHandler);
            
            degreeElevationBtn = uicontrol(f,'Style','pushbutton','String','Apply Degree Elevation +1',...
                'Position',[180 10 150 40],'Callback',@new.degreeElevationhandler);
            
            cpCkbx = uicontrol(f,'Style','checkbox',...
                'String','Display Control Points',...
                'BackgroundColor','white',...
                'Value',1,'Position',[340 23 130 20],'Callback',@new.cpHideShowHandler);
            
            cpLegCkbx = uicontrol(f,'Style','checkbox',...
                'String','Display Legs',...
                'BackgroundColor','white',...
                'Value',1,'Position',[340 3 130 20],'Callback',@new.legHandler);
            
            new.convexHullCkbx = uicontrol(f,'Style','checkbox',...
                'String','Display Convex Hull',...
                'BackgroundColor','white',...
                'Value',1,'Position',[340 43 130 20],'Callback',@new.convexHullHandler);
            
            
            %draw point on figure and add event handlers
            hold on;
            for i = 1:size(CP,2)
                pi = plot(CP(1,i), CP(2,i), 'or',...
                    'MarkerSize',10,'MarkerFaceColor','r');     %plot current graphic control point
                set(pi,'ButtonDownFcn',@new.dragObjectMOD);     %add handler to graphic control point
                id = ['p' num2str(i)];                          %add unique control point id as tag to graphic control point
                tagm(pi, id);
                
                %instantiation of control point object using ID and pi-graphic cp
                if(i == 1)
                    cpi = controlPoint(id, pi, -1, 1);
                else
                    cpi = controlPoint(id, pi, new.CP_obj(i-1), 1);
                end
                
                new.tag_CPobj_map(id) = cpi;                    %association of unique id-tag-graphic to object control point
                new.CP_obj = [new.CP_obj cpi];                  %add control point object to obejct list
            end
            hold off;
            
            %plot associate curve
            [new.X, new.plotCurve] = CAGDRoutine.plotCurve_equiT(1,CP,400,1,'k',3);
            
            %plot associate convex hull
            figure(1);
            hold on;
            k = convhull(CP(1,:),CP(2,:));
            new.convexHull = plot(CP(1,k),CP(2,k),'r-');
            hold off;
            
            %Link handler to figure
            set(f,'WindowButtonUpFcn', {@new.dropObjectMOD});
            set(f,'WindowButtonMotionFcn',{@new.moveObjectMOD,ax});
            
            axis equal
            
            %grain settings
            grid on
            set(gca,'Xtick',-1000:1:1000) %grid grain
            set(gca,'Ytick',-1000:1:1000) %grid grain
            
            % set(gca,'Xtick',XGridLowBound:XGridStep:XGridUppBound);
            % set(gca,'Ytick',YGridLowBound:YGridStep:YGridUppBound);
            
            set(gca,'XMinorTick','on')
            grid minor
            
            % axis([-4 4 -3 3]);
            axis([min(CP(1,:))-3 max(CP(1,:))+3 min(CP(2,:))-1 max(CP(2,:))+1 ]);
            
        end
        
        %---------------- Curve Functions ---------------------------------
        
        function addPoint(this, xCoo, yCoo)
            
            xMin = min(this.X(1,:));
            xMax = max(this.X(1,:));
            yMin = min(this.X(2,:));
            yMax = max(this.X(2,:));
            
            if(strcmp(xCoo,'null') || strcmp(yCoo,'null') )
                x = this.randomCooRange(xMin, xMax);
                y = this.randomCooRange(yMin, yMax);
            else
                x = xCoo;
                y = yCoo;
            end
            numNewCP = size(this.CP_obj,2)+1;
            
            hold on;
            %-----------
            pi = plot(x, y, 'or',...
                'MarkerSize',10,'MarkerFaceColor','r');     %plot current graphic control point
            set(pi,'ButtonDownFcn',@this.dragObjectMOD);    %add handler to graphic control point
            id = ['p' num2str(numNewCP)];                   %add unique control point id as tag to graphic control point
            tagm(pi, id);
            
            %instantiation of control point object using ID and pi-graphic cp
            cpi = controlPoint(id, pi, this.CP_obj(numNewCP-1), 1);
            
            this.tag_CPobj_map(id) = cpi;                    %association of unique id-tag-graphic to object control point
            this.CP_obj = [this.CP_obj cpi];                  %add control point object to obejct list
            %-----------
            hold off;
        end
        
        
        
        %--------------- Support Function ---------------------------------
        
        function [r] = randomCooRange(this,a,b)
            r = (b-a).*rand() + a;
        end
        
        %!! > Draw curve based on controlPoint object in this.CP_obj
        function reDrawCurve(this, numTValues, lineWidth, lineColor, lineStyle)
            %retrieve all XYCoos of control points by objects ref
            newCP = [];
            for cp = this.CP_obj
                cpID = cp.id;
                [x, y] = cp.getCurrentCoo();
                %disp([cpID ' ' num2str(x) ' ' num2str(y)]); %display XYCoo
                
                newCP = [newCP [x;y]];
            end
            
            %retrive new X points without plot curve but using the exist one
            [this.X, nothing] = CAGDRoutine.plotCurve_equiT(0, newCP, numTValues, 1, '--r', 0.3);
            %update coordinate for existing curve
            set(this.plotCurve,'XData',this.X(1,:));
            set(this.plotCurve,'YData',this.X(2,:));
            %set visual option
            set(this.plotCurve,'LineWidth', lineWidth);
            set(this.plotCurve,'Color', lineColor);
            set(this.plotCurve,'LineStyle', lineStyle);
            
            set(0,'currentFigure',this.mainFigure);
            axis([min(newCP(1,:))-3 max(newCP(1,:))+3 min(newCP(2,:))-1 max(newCP(2,:))+1 ]);
            
        end
        
        function updateVisualLegs(this,vis_left_leg, vis_right_leg,x,y)
            
            if(isa(vis_left_leg,'matlab.graphics.chart.primitive.Line'))
                xData = get(vis_left_leg,'XData');
                xData(2) = x;
                set(vis_left_leg,'XData', xData);
                
                yData = get(vis_left_leg,'YData');
                yData(2) = y;
                set(vis_left_leg,'YData', yData);
            end
            
            if(isa(vis_right_leg,'matlab.graphics.chart.primitive.Line'))
                xData = get(vis_right_leg,'XData');
                xData(1) = x;
                set(vis_right_leg,'XData', xData);
                
                yData = get(vis_right_leg,'YData');
                yData(1) = y;
                set(vis_right_leg,'YData', yData);
            end
            
        end
        
        function XYCoo = retriveXYCoo(this)
            % Output
            % XYCoo     - [2 x #ControlPoints] return the xy coordinate of
            %             current list of control points.
            XYCoo = [];
            for cp = this.CP_obj
                cpID = cp.id;
                [x, y] = cp.getCurrentCoo();
                disp([cpID ' ' num2str(x) ' ' num2str(y)]);
                
                XYCoo = [XYCoo [x;y]];
            end
        end
        
        function hideConvexHull(this)
            
            hold on;
            set(this.convexHull,'visible','off');
            hold off;
            
            set(this.convexHullCkbx,'Value',0);
        end
        
        %---------------- Handler Function --------------------------------
        %!! > varargin{} = [this callingObj windowDatas Axis];
        
        function degreeElevationhandler(varargin)
            
            this = varargin{1};
            this.hideConvexHull();
            
            XYCoo = this.retriveXYCoo();
            newCP = DeCasteObj.degreeElevation(XYCoo);
            
            %update visual coordinate and add last control point
            for i = 1:size(newCP,2)-1
                visualCP = this.CP_obj(i).cp;
                set(visualCP,'XData',newCP(1,i));
                set(visualCP,'YData',newCP(2,i));
                
                %update visual left leg
                vis_left_leg = this.CP_obj(i).left_leg;
                vis_right_leg = this.CP_obj(i).right_leg;
                this.updateVisualLegs(vis_left_leg,vis_right_leg,newCP(1,i),newCP(2,i))
                
                
            end
            
            %add last new control point
            this.addPoint( newCP(1,size(newCP,2)), newCP(2,size(newCP,2)) );
            %draw new curve
            %             this.reDrawCurve(400, 3, 'k', '-');
            
        end
        
        function addCPHandler(varargin)
            this = varargin{1};
            this.hideConvexHull();
            
            this.addPoint('null','null');
            this.reDrawCurve(400, 3, 'k', '-');
        end
        
        function dragObjectMOD(varargin)
            currentObj = varargin{1};
            currentObj.hideConvexHull();
            
            currentObj.dragging = varargin{2};
            currentObj.orPos = get(gcf,'CurrentPoint');
        end
        
        function moveObjectMOD(varargin)
            currentObj = varargin{1};
            %             currentObj.hideConvexHull();
            
            if ~isempty(currentObj.dragging)
                
                % Serves as the buttondownfcn for the axes.
                objFigure = varargin{2};  % Get the structure.
                cooNORM = get(objFigure, 'currentpoint'); % Get the position of the mouse.
                
                objAxis = varargin{4};
                cooAXIS = get(objAxis, 'currentpoint'); % Get the position of the mouse.
                
                %set new position for current control point
                currentPoint = currentObj.dragging;%varargin{4};
                set(currentPoint,'XData',cooAXIS(1,1));
                set(currentPoint,'YData',cooAXIS(1,2));
                
                %-------------set new position for curreng leg-----------------
                %tramite la map prendo il riferimento all'oggetto control point
                tagcp = get(currentPoint,'Tag');
                objCP = currentObj.tag_CPobj_map(tagcp);
                left_leg = objCP.left_leg;
                right_leg = objCP.right_leg;
                currentObj.updateVisualLegs(left_leg,right_leg,cooAXIS(1,1),cooAXIS(1,2));
                
                currentObj.reDrawCurve(50, 0.3, [0.5 0.5 0.5], '--');
                
                %--------------------------------------------------------------
                
                %----------- compute curvature dinamically ----------------
                %curvature
                h = 1/400;
                xt = currentObj.X(1,:);
                yt = currentObj.X(2,:);
                n = size(xt,2);
                %calcolo le derivate prima
                xt1 = (xt(2:n) - xt(1:n-1))/h;
                yt1 = (yt(2:n) - yt(1:n-1))/h;
                
                %calcolo derivata seconda
                xt2 = ( xt(3:n) - 2*xt(2:n-1) + xt(1:n-2) ) / (h^2);
                yt2 = ((yt(3:n) - 2*yt(2:n-1) + yt(1:n-2))) / (h^2);
                
                xt1 = xt1(2:size(xt1,2)); %non esiste per i = 0 la d''
                yt1 = yt1(2:size(yt1,2)); %non esiste per i = 0 la d''
                
                size(xt1,2)
                size(xt2,2)
                
                kt = ( abs((xt1.*yt2) - (yt1.*xt2) )) ./ ( (xt1.^2 + yt1.^2).^(3/2) );
                curv = figure(2);
                set(curv,'Name','Bezier Curve Curvature','NumberTitle','off');
                plot(linspace(0,1,48),kt,'-b')
%                 axis equal
%                 axis([0 1 min(kt)-0.1 max(kt)+0.1])
                axis([0 1 -10 10]);
                grid on;
                set(gca,'Xtick',-1000:0.1:1000) %grid grain
                set(curv,'units','normalized');
                set(curv,'color','w');
                set(curv, 'Position', [0.7 0.3 0.3 0.3]);
                
                %----------------------------------------------------------
            end
        end
        
        function dropObjectMOD(varargin)
            currentObj = varargin{1};
            
            if ~isempty(currentObj.dragging)
                
                %remove reference for dragging control point
                currentObj.dragging = [];
                
                currentObj.reDrawCurve(400, 3, 'k', '-');
                %--------------------------------------------------------------
            end
        end
        
        function convexHullHandler(varargin)
            this = varargin{1};
            currCPCoo = this.retriveXYCoo();
            
            if( strcmp(get(this.convexHull,'Visible'),'off'))
                
                chIncidicies = convhull(currCPCoo(1,:),currCPCoo(2,:));
                set(this.convexHull,'XData',currCPCoo(1,chIncidicies));
                set(this.convexHull,'YData',currCPCoo(2,chIncidicies));
                set(this.convexHull,'Visible','on');
                
            else
                this.hideConvexHull();
            end
        end
        
        function legHandler(varargin)
            checkBox = varargin{2};
            visible = [];
            if(get(checkBox,'value') == 1)
                visible = 'on';
            else
                visible = 'off';
            end
            
            for cp = varargin{1}.CP_obj
                if(isa(cp.left_leg,'matlab.graphics.chart.primitive.Line') == 1)
                    set(cp.left_leg,'visible',visible);
                end
            end
            
        end
        
        function cpHideShowHandler(varargin)
            checkBox = varargin{2};
            visible = [];
            if(get(checkBox,'value') == 1)
                visible = 'on';
            else
                visible = 'off';
            end
            
            for currCP = varargin{1}.CP_obj
                set(currCP.cp,'visible',visible);
            end
            
        end
    end
    
    methods (Static)
        
        %A new set of control point coordinates is generated
        function newCP = degreeElevation(CP)
            n = size(CP,2);
            newCP = [];
            newCP = [newCP [CP(1,1); CP(2,1)]];
            
            for i=2:n
                newCP(:,i) = ( (i-1)/(n+2) )*CP(:,i-1) + (1 - ( (i-1)/(n+2)) )*CP(:,i);
            end
            
            newCP = [newCP [CP(1,n); CP(2,n)]];
            
        end
    end
    
end

