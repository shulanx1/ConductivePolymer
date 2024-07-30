function plot_traces_with_gradient(x, y, color)
% x: 1D array of time
% y: ND array, data arranged in columns
% color: the darkest color that gradient needs to be plot
if nargin < 3, color = [0,0,0]; end
len = size(y, 2) + 2;
red = color;
pink = [1,1,1];
colors_p = [linspace(pink(1),red(1),len)', linspace(pink(2),red(2),len)', linspace(pink(3),red(3),len)'];
figure
for i = 1:size(y, 2)
    plot(x, y(:, i), 'Color', colors_p(i+2, :))
    hold on
end
end