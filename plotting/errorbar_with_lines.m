function errorbar_with_lines(x, y, colors)
% x: cell or 1D array, datapoints on the x axis
% Y: cell or 2D array (dataset arranged as column), datapoints on the y
% axis
% colors: color to be used, N*3 matrix, each color in each row

% example:
% x = [1:20]';
% y = cell(1,2);
% for i = 1:20
%     y{1} = [y{1}, x+5*rand(20,1)];
%     y{2} = [y{2}, 2*x+5*rand(20,1)];
% end
% errorbar_with_lines(x, y);

if nargin < 3
    colors = [[0,0,0];[119,177,204];[61,139,191];[6,50,99]];
    colors = colors/256;
end

if ~iscell(y)
    data{1} = y;
    N = 1;
else
    data = y;
    N = length(y);
end

X = cell(1, N);
if ~iscell(x)
    for i = 1:N
        X{i} = x;
    end
end

figure
hold on
for i = 1:N
    color_idx = mod(i, size(colors, 1));
    if color_idx == 0
        color_idx = size(colors, 1);
    end
    sample_size = zeros(size(data{i},1),1);
    for j = 1:size(data{i},1)
        sample_size(j) = length(find(~isnan(data{i}(j,:))));
    end
    errorbar(X{i}, nanmean(data{i}, 2),  nanstd(data{i}, [], 2)./sqrt(sample_size),'o','Color', colors(color_idx,:),  'MarkerSize',10,  'MarkerEdgeColor', colors(color_idx,:),'MarkerFaceColor', 'None')
    hold on
    plot(X{i}, nanmean(data{i}, 2),'Color', colors(color_idx,:), 'Linewidth', 2)
end
end