clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

datapath = 'E:\data\polymer\bloodvessel';
files = dir(fullfile(datapath, '*_Segment_Statistics.csv'));
M = readtable('E:\data\polymer\pool_bloodvessel.xlsx', 'Sheet','Sheet1');

base_path = 'E:\data';
save_path = 'E:\data\polymer\bloodvessel\stat';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end
vesc_len = cell(size(M,1), size(M,2));
vesc_diam = cell(size(M,1), size(M,2));
vesc_volume = cell(size(M,1), size(M,2));
csv_file = cell(size(M,1), size(M,2));
mouse_id = cell(size(M,1),size(M,2));
ROI_id = cell(size(M,1), size(M,2));

for i = 1:size(M, 2)
    for j = 1:size(M,1)
        if ~isempty(M{j,i}{1})
            ROI_name = M{j,i}{1};
            for ii = 1:length(files)
                if contains(files(ii).name, ROI_name)
                    csv_file{j,i} = fullfile(files(ii).folder, files(ii).name);
                    break
                end
            end
            if ~isempty(csv_file{j,i})
                T = readtable(csv_file{j,i}, 'Delimiter', ';');
                vesc_len{j,i} = [T.length]/0.9;
                vesc_diam{j,i} = [T.diameter]/0.9;
                vesc_volume{j,i} = [T.volume]/0.9;
                splitParts = strsplit(files(ii).name, '_');
                if strcmp(splitParts(2),'m1')
                    mouse_id{j,i} = ones(size(vesc_len{j,i}));

                else
                    mouse_id{j,i} = 2*ones(size(vesc_len{j,i}));
                end
                ROI_id{j,i} = j*ones(size(vesc_len{j,i}));
            end
        end
    end
    
end

treat_id = cell(1, size(M,2));
vessel_id = cell(1, size(M,2));
for i = 1:size(vesc_diam,2)
    vesc_volume_pool{i} = cell2mat(vesc_volume(:,i));
    vesc_len_pool{i} = cell2mat(vesc_len(:,i));
    vesc_diam_pool{i} = cell2mat(vesc_diam(:,i));   
    mouse_id_pool{i} = cell2mat(mouse_id(:,i));
    ROI_id_pool{i} = cell2mat(ROI_id(:,i));
    idx_rmv = unique([find(vesc_len_pool{i}<=160);find(vesc_volume_pool{i}<=500);find(vesc_volume_pool{i}>=10000)]);
    vesc_volume_pool{i}(idx_rmv) = [];
    vesc_len_pool{i}(idx_rmv) = [];
    vesc_diam_pool{i}(idx_rmv) = [];
    mouse_id_pool{i}(idx_rmv) = [];
    ROI_id_pool{i}(idx_rmv) = [];
    vesc_volume_pool{i} = vesc_volume_pool{i}/1e3;
    treat_id{i} = i*ones(size(vesc_volume_pool{i}));
    vessel_id{i} = [1:length(vesc_volume_pool{i})]';
    
end
colors = [[239,153,75];[231,193,215];[220,165,195];[137,62,129];[175,178,207]]/255;
violinplot_with_datapoint(vesc_diam_pool,colors)
[p,~,stats] = anovan(cell2mat(vesc_diam_pool'),{cell2mat(treat_id'),cell2mat(mouse_id_pool'),cell2mat(ROI_id_pool'),cell2mat(vessel_id')},'varnames',{'treatment','mouse','ROI','vessel'},'nested',[[0,0,0,0];[0,0,0,0];[0,1,0,0];[0,0,1,0]],'random',[2,3,4])

% for i = 1:size(vesc_diam,2)
%     vesc_volume_pool{i} = cellfun(@(x,y,z) median(x(unique([find(y<=100);find(z<=500);find(z>=10000)]))), vesc_volume(:,i),vesc_len(:,i),vesc_volume(:,i)); 
%     vesc_len_pool{i} = cellfun(@(x,y,z) median(x(unique([find(y<=100);find(z<=500);find(z>=10000)]))), vesc_len(:,i),vesc_len(:,i),vesc_volume(:,i)); 
%     vesc_diam_pool{i} = cellfun(@(x,y,z) median(x(unique([find(y<=100);find(z<=500);find(z>=10000)]))), vesc_diam(:,i),vesc_len(:,i),vesc_volume(:,i));    
%     
% end
% violinplot_with_datapoint(vesc_diam_pool)
save(fullfile(save_path, 'vessel_morph.mat'))
%%
img_path = '\\10.165.57.13\Sutter_backup\Shulan\2p imaging\polymer_bloodvessel';
M = readtable('E:\data\polymer\pool_bloodvessel.xlsx', 'Sheet','Sheet2');
files = dir(fullfile(img_path, '*.tif'));
hist_stp = 10;
thr_high = 0.95;
ratio = cell(size(M,1), size(M,2));
for i = 1:size(M, 2)
    for j = 1:size(M,1)
        if ~isempty(M{j,i}{1})
            ROI_name = M{j,i}{1};
            for ii = 1:length(files)
                if contains(files(ii).name, ROI_name)
                    tif_file{j,i} = fullfile(files(ii).folder, files(ii).name);
                    break
                end
            end
            if ~isempty(tif_file{j,i})
                info = imfinfo(tif_file{j,i}); % Get information about the TIFF file
                numPages = length(info);
                ratio{j,i} = zeros(numPages,1);
                for ii = 1:numPages
                    V = imread(tif_file{j,i},ii);
                    [a, x] = histcounts(reshape(V, [],1),0:stp:255);
                    bar_x = x(1:end-1)+stp/2;
                    [sort_a,I] = sort(a,'descend');
                    sort_x = bar_x(I);
                    if (sort_x(1) == bar_x(1))||(sort_x(1)==bar_x(end))
                        ratio{j,i}(ii) = mean(double(sort_x(2:3)))/double(quantile(reshape(V, [],1),thr_high));
                    else
                    	ratio{j,i}(ii) = mean(double(sort_x(1:2)))/double(quantile(reshape(V, [],1),thr_high));
                    end
                end
            end
        end
    end
    
end
% ratio_pool = cell(1, size(M,2));
% for i = 1:size(M,2)
%     a = ratio(:,i)';
%     idx_keep = cellfun(@(x) ~isempty(x), a);
%     a = a(idx_keep);
%     ratio_pool{i} = padcat(a{:});
% end
% colors = [[239,153,75];[231,193,215];[220,165,195];[137,62,129];[175,178,207];[30,30,30]]/255;
% lineplot_with_shaded_errorbar(0:5:70*5, ratio_pool, colors);
% 

ratio_pool = cell(1, size(M,2));
for i = 1:size(M,2)
    a = ratio(:,i)';
    idx_keep = cellfun(@(x) ~isempty(x), a);
    ratio_pool{i} = cellfun(@(x) median(x),a(idx_keep));
end
violinplot_with_datapoint(ratio_pool, colors)

[measurements, gp1, gp2] = gen_table_for_unbalanced_anova(1, ratio_pool);
[p,~,stats] = anovan(measurements,{gp1,gp2});
[results,~,~,gnames] = multcompare(stats,"Dimension",[2]);
save(fullfile(save_path, 'BBB_permeability.mat'), 'ratio', 'ratio_pool', 'M', 'results')