
function deCasteljouDragAndDrop

clc;
clear;
clf;

%global variables
dragging = [];
orPos = [];
tag_CPobj_map = containers.Map()
CP_obj = []; %collection of control point objects

CP = [-3.5 -3 -2 -1 0 1 3 3.5;
      0    2   1 -1 1 2 2.5 1];


f = figure(1);
set(f,'units','normalized')
set(gcf,'color','w');
set(f, 'Position', [0 0 100 100])

ax = axes('xlim',[-4 4],'ylim',[-3 3]);

hold on;
%
for i = 1:size(CP,2)
    
    %plot current graphic control point
    pi = plot(CP(1,i), CP(2,i), 'or','MarkerSize',10,'MarkerFaceColor','r');
    
    %add handler to graphic control point
    set(pi,'ButtonDownFcn',@dragObjectMOD);
    
    %add unique control point id as tag to graphic control point
    id = ['p' num2str(i)];
    tagm(pi, id);
    
    %instantiation of control point object using ID and pi-graphic cp
    if(i == 1)
        cpi = controlPoint(id, pi, -1, 1);
    else
        cpi = controlPoint(id, pi, CP_obj(i-1), 1);
    end
    
    %association of unique id-tag-graphic to object control point
    tag_CPobj_map(id) = cpi;
    
    %add control point object to obejct list
    CP_obj = [CP_obj cpi];
end

%plot associate curve
[X, plotCurve] = CAGDRoutine.plotCurve_equiT(1,CP,400,1,'k',3);

hold off;

%------------------ manual link handlers to plot object

% hold on
% 
% p1 = plot(-3,-1,'or');
% set(p1,'ButtonDownFcn',@dragObjectMOD);
% tagm(p1,'p1');
% cp1 = controlPoint('p1', p1,-1,1);
% tag_CPobj_map('p1') = cp1;
% 
% 
% p2 = plot(0,-2,'or');
% set(p2,'ButtonDownFcn',@dragObjectMOD);
% tagm(p2,'p2');
% cp2 = controlPoint('p2', p2, cp1, 1);
% tag_CPobj_map('p2') = cp2;
% 
% 
% p3 = plot(2,1,'or');
% set(p3,'ButtonDownFcn',@dragObjectMOD);
% tagm(p3,'p3');
% cp3 = controlPoint('p3', p3,cp2,1);
% tag_CPobj_map('p3') = cp3;
% 
% hold off;

%--------------------------------------------


%Link handler to figure
set(f,'WindowButtonUpFcn', @dropObjectMOD);
set(f,'WindowButtonMotionFcn',{@moveObjectMOD,ax});

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



    function dragObjectMOD(varargin)
        dragging = varargin{1};
        orPos = get(gcf,'CurrentPoint');
    end

    function moveObjectMOD(varargin)
        
        if ~isempty(dragging)
            
            % Serves as the buttondownfcn for the axes.
            objFigure = varargin{1};  % Get the structure.
            cooNORM = get(objFigure, 'currentpoint'); % Get the position of the mouse.
            
            objAxis = varargin{3};
            cooAXIS = get(objAxis, 'currentpoint'); % Get the position of the mouse.
            
            %set new position for current control point
            currentPoint = dragging;%varargin{4};
            set(currentPoint,'XData',cooAXIS(1,1));
            set(currentPoint,'YData',cooAXIS(1,2));
            
            %-------------set new position for curreng leg-----------------
            %tramite la map prendo il riferimento all'oggetto control point
            tagcp = get(currentPoint,'Tag');
            objCP = tag_CPobj_map(tagcp);
            left_leg = objCP.left_leg;
            right_leg = objCP.right_leg;
            if(isa(left_leg,'matlab.graphics.chart.primitive.Line'))
                
                xData = get(left_leg,'XData');
                xData(2) = cooAXIS(1,1);
                set(left_leg,'XData', xData);
                
                yData = get(left_leg,'YData');
                yData(2) = cooAXIS(1,2);
                set(left_leg,'YData', yData);
                
            end
            if(isa(right_leg,'matlab.graphics.chart.primitive.Line'))
                
                xData = get(right_leg,'XData');
                xData(1) = cooAXIS(1,1);
                set(right_leg,'XData', xData);
                
                yData = get(right_leg,'YData');
                yData(1) = cooAXIS(1,2);
                set(right_leg,'YData', yData);
            end
            
            %--------------------------------------------------------------
            
            %------------- draw new curve light sketch --------------------
            %retrieve all XYCoos of control points by objects ref
            newCP = [];
            for cp = CP_obj
                cpID = cp.id;
                [x, y] = cp.getCurrentCoo();
                disp([cpID ' ' num2str(x) ' ' num2str(y)]);
                
                newCP = [newCP [x;y]];
            end
            
            %retrive new X points without plot curve but using the exist one
            [X, nothing] = CAGDRoutine.plotCurve_equiT(0, newCP, 50, 1, '--r', 0.3);
            %update coordinate for existing curve
            set(plotCurve,'XData',X(1,:));
            set(plotCurve,'YData',X(2,:));
            %set visual option
            set(plotCurve,'LineWidth', 0.3);
            set(plotCurve,'Color', [0.5 0.5 0.5]);
            set(plotCurve,'LineStyle', '--');
            
            %--------------------------------------------------------------
            
            axis([min(newCP(1,:))-3 max(newCP(1,:))+3 min(newCP(2,:))-1 max(newCP(2,:))+1 ]);
            %axis([-4 4 -3 3]);


            
        end
    end

    function dropObjectMOD(hObject,eventdata)
        if ~isempty(dragging)
            %             newPos = get(gcf,'CurrentPoint');
            %             posDiff = newPos - orPos;
            %             set(dragging,'Position',get(dragging,'Position') + [posDiff(1:2) 0 0]);
            
            %remove reference for dragging control point
            dragging = [];
            
            %------------- plot stable curve on controlPoint drop ---------
            %retrieve all XYCoos of control points by objects ref
            newCP = [];
            for cp = CP_obj
                cpID = cp.id;
                [x, y] = cp.getCurrentCoo();
                disp([cpID ' ' num2str(x) ' ' num2str(y)]);
                
                newCP = [newCP [x;y]];
            end
            
            %retrive new X points without plot curve but using the exist one
            [X, nothing] = CAGDRoutine.plotCurve_equiT(0, newCP, 400, 1, 'k', 3);
            %update coordinate for existing curve
            set(plotCurve,'XData',X(1,:));
            set(plotCurve,'YData',X(2,:));
            %set visual option
            set(plotCurve,'LineWidth', 3);
            set(plotCurve,'Color', 'k');
            set(plotCurve,'LineStyle', '-');
            
            
            %--------------------------------------------------------------
        end
    end

%-------------------------------------------------------------------------


    function dragObject(hObject,eventdata)
        dragging = hObject;
        orPos = get(gcf,'CurrentPoint');
    end

    function moveObject(hObject,eventdata)
        if ~isempty(dragging)
            newPos = get(gcf,'CurrentPoint');
            posDiff = newPos - orPos;
            orPos = newPos;
            set(dragging,'Position',get(dragging,'Position') + [posDiff(1:2) 0 0]);
        end
    end

    function dropObject(hObject,eventdata)
        if ~isempty(dragging)
            newPos = get(gcf,'CurrentPoint');
            posDiff = newPos - orPos;
            set(dragging,'Position',get(dragging,'Position') + [posDiff(1:2) 0 0]);
            dragging = [];
        end
    end


end

