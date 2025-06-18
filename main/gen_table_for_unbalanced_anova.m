
function [measurements, gp1, gp2] = gen_table_for_unbalanced_anova(x, y)
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
% datatable = gen_table_for_ranova(x, y)


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


measurements = cell2mat(data);

gp1 = [];
gp2= [];
for i = 1:length(data)
    for j = 1:size(data{i},2)
        gp1 = [gp1; X{i}'];
    end
    gp2 = [gp2; i*ones(length(X{i})*size(data{i},2),1)];
end
measurements = reshape(measurements,[],1);
end