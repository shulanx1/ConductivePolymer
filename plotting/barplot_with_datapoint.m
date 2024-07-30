function barplot_with_datapoint(y, colors)
% y: 2D matrix (each dataset in each column)
% colors: color to be used, N*3 matrix, each color in each row

% example
% y = rand(100,2);  % one pair of y values
% bar_pairwise(y);

% y = {rand(100,1), rand(50, 1)+1}; % multiple pairs of y values
% bar_pairwise(y);


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
    data = cell(1, size(y,2));
    mean_y = zeros(1, size(y,2));
    sem_y = zeros(1, size(y,2));
    for i = 1:size(y,2)
        data{i} = y(:,i);
        data{i}(isoutlier(y(:,i))) = [];
        mean_y(i) = mean(data{i});
        sem_y(i) = std(data{i})/sqrt(length(data{i}));
    end
    bar([1:length(data)],mean_y,'FaceColor','None','EdgeColor','k','LineWidth',1)
    hold on
    er = errorbar([1:length(data)],mean_y,sem_y);    
    er.Color = 'k';                            
    er.LineStyle = 'none';  
    for i = 1:length(data)
    color_idx = mod(i, size(colors, 1));
    if color_idx == 0
        color_idx = size(colors, 1);
    end
    scatter(i+(rand(1, length(data{i}))-0.5)/2,data{i},20,'MarkerFaceColor', colors(color_idx,:),'MarkerEdgeColor', 'None')
    end
else
    mean_y = zeros(1, length(y));
    sem_y = zeros(1, length(y));
    data = y;
    for i = 1:size(y,2)
        data{i}(isoutlier(data{i})) = [];
        mean_y(i) = mean(data{i});
        sem_y(i) = std(data{i})/sqrt(length(data{i}));
    end
    f = figure();
    size_x = (0.4 + length(y)*0.15)/2; % figure size in inches
    size_y = 1; 
    f.Position(3) = f.Position(4)*size_x;
    f.Renderer = 'painters';
    bar([1:length(data)],mean_y,'FaceColor','None','EdgeColor','k','LineWidth',1)
    hold on
    er = errorbar([1:length(data)],mean_y,sem_y);    
    er.Color = 'k';                            
    er.LineStyle = 'none';  
    for i = 1:length(data)
    color_idx = mod(i, size(colors, 1));
    if color_idx == 0
        color_idx = size(colors, 1);
    end
    scatter(i+(rand(1, length(data{i}))-0.5)/2,data{i},20,'MarkerFaceColor', colors(color_idx,:),'MarkerEdgeColor', 'None')
    end

end

end