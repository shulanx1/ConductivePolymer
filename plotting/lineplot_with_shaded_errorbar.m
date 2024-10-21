function lineplot_with_shaded_errorbar(x, y, colors)
% x: cell or 1D array, datapoints on the x axis
% Y: cell or 2D array (dataset arranged as column), datapoints on the y
% axis
% colors: color to be used, N*3 matrix, each color in each row

% example:
% x = [1:10]';
% y = cell(1,2);
% for i = 1:20
%     y{1} = [y{1}, x+5*rand(20,1)];
%     y{2} = [y{2}, 2*x+5*rand(20,1)];
% end
% lineplot_with_shaded_errorbar(x, y);
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
for i = 1:N
    color_idx = mod(i, size(colors, 1));
    if color_idx == 0
        color_idx = size(colors, 1);
    end
    plot(X{i}, nanmean(data{i},2), 'Color', colors(color_idx,:),'Linewidth', 2)
    sample_size = zeros(size(data{i},1),1);
    for j = 1:size(data{i},1)
        sample_size(j) = length(find(~isnan(data{i}(j,:))));
    end
    hold on
    plot(X{i}, nanmean(data{i},2) + nanstd(data{i},[],2)./sqrt(sample_size), 'Color', colors(color_idx,:),'Linewidth', 0.5)
    plot(X{i}, nanmean(data{i},2) - nanstd(data{i},[],2)./sqrt(sample_size), 'Color', colors(color_idx,:),'Linewidth', 0.5)
end
end