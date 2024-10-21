function [datatable, within_design] = gen_table_for_ranova(x, y)
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

for i = 1:length(data)
    data{i} = data{i}';
end
measurements = cell2mat(data);
datatable = array2table(measurements);
% measurements = [];
% groups = [];
% cells = [];
% Xs = [];
% for i = 1:length(data)
%     measurements = [measurements;reshape(data{i},[],1)];
%     Xs = [Xs; repmat(X{i}',size(data{i},2),1)];
%     for j = 1:size(data{i},2)
%         cells = [cells;j*ones(size(data{i},2), 1)];
%     end
%     groups = [groups; i*ones(size(reshape(data{i},[],1)))];
%     
% end
% datatable = table(groups, cells, Xs, measurements,'VariableNames',{'treatment','cell','voltage','measurements'});
groups_within = [];
Xs_within= [];
for i = 1:length(data)
    groups_within = [groups_within; i*ones(length(X{1}),1)];
    Xs_within = [Xs_within; X{i}'];
end
within_design = table(groups_within, Xs_within,'VariableNames',{'treatment', 'voltage'});
within_design.treatment = categorical(within_design.treatment);
within_design.voltage = categorical(within_design.voltage);
end