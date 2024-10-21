function barplot_pairwise(y, colors)
% y: 2D matrix (each dataset in each column)
% colors: color to be used, N*3 matrix, each color in each row

% example
% y = rand(100,2);  % one pair of y values
% barplot_pairwise(y);

% y = {rand(100,2), rand(100,2)+1};  % two pair of y values
% barplot_pairwise(y);




if nargin < 2
    colors = [[0,0,0];[119,177,204];[61,139,191];[6,50,99]];
    colors = colors/256;
end
if ~iscell(y)
    size_x = (0.4 + size(y, 2)*0.15)/2; % figure size in inches
    size_y = 1; 
    f = figure();
    f.Position(3) = f.Position(4)*size_x;
    f.Renderer = 'painters';
    mean_y = mean(y,1);
    sem_y = std(y,[],1)/sqrt(size(y,2));
    hold on
    for i = 1:size(y,2)
        color_idx = mod(i, size(colors, 1));
        if color_idx == 0
            color_idx = size(colors, 1);
        end
        bar(i,mean_y(i),'FaceColor','None','EdgeColor',colors(color_idx,:),'LineWidth',1)
        er = errorbar(i,mean_y(i),sem_y(i));    
        er.Color = colors(color_idx,:);                            
        er.LineStyle = 'none';  
        er.LineWidth = 1;
        er.CapSize = 10;
    end
    for i = 1:size(y, 1)
        plot(1:size(y, 2), y(i,:), 'Color', [0.5,0.5,0.5])
    end
    box on
else
    size_x = (0.4 + length(y)*2*0.15)/2; % figure size in inches
    size_y = 1; 
    f = figure();
    f.Position(3) = f.Position(4)*size_x;
    f.Renderer = 'painters';
    data = y;
    for n = 1:length(data)
        y = data{n};
        mean_y = mean(y,1);
        sem_y = std(y,[],1)/sqrt(size(y,2));
        hold on
        for i = 1:size(y,2)
            color_idx = mod(i, size(colors, 1));
            if color_idx == 0
                color_idx = size(colors, 1);
            end
            bar(i+(n-1)*2,mean_y(i),'FaceColor','None','EdgeColor',colors(color_idx,:),'LineWidth',1)
            er = errorbar(i,mean_y(i),sem_y(i));    
            er.Color = colors(color_idx,:);                            
            er.LineStyle = 'none';  
            er.LineWidth = 1;
            er.CapSize = 10;
        end
        for i = 1:size(y, 1)
            plot(1:size(y, 2)+(n-1)*2, y(i,:), 'Color', [0.5,0.5,0.5])
        end
    end
    box on
end

end