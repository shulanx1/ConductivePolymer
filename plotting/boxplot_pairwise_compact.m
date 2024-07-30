function boxplot_pairwise_compact(y, colors)
% y: 2D matrix (each dataset in each column)
% colors: color to be used, N*3 matrix, each color in each row

% example
% y = rand(100,2);  % one pair of y values
% boxplot_pairwise(y);

% y = {rand(100,2), rand(50, 2)+1}; % multiple pairs of y values
% boxplot_pairwise(y);


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
    boxplot(y, 'Color', colors, 'PlotStyle','compact')
    hold on
    for i = 1:size(y, 1)
        plot(1:size(y, 2), y(i,:), 'Color', [0.5,0.5,0.5])
    end
else
    data = [];
    grp = [];
    count = 1;
    for i = 1:length(y)
        for j = 1:size(y{i},2)
            data = [data;y{i}(:,j)];
            grp = [grp, count*ones(1, size(y{i},1))];
            count = count + 1;
        end
    end
    size_x = (0.4 + length(y)*2*0.15)/2; % figure size in inches
    size_y = 1; 
    f = figure();
    f.Position(3) = f.Position(4)*size_x;
    f.Renderer = 'painters';
    boxplot(data,grp, 'Color', colors, 'PlotStyle','compact')
    hold on
    count = 1;
    for i = 1:length(y)
        for j = 1:size(y{i},1)
            plot(count:count+size(y{i}, 2)-1, y{i}(j,:), 'Color', [0.5,0.5,0.5])
        end
        count = count + size(y{i}, 2);
    end
end

end