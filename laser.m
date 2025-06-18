clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))
base_path = 'E:\data';
M = readtable('E:\data\polymer\pool_temp_calibration.xlsx', 'Sheet','polymer');
base_path = 'E:\data';
save_path = 'E:\data\polymer\temp_calibration';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end

load('E:\data\polymer\temp_calibration\temp_calibration_log.mat', 'mdl')

dT_polymer = nan(size(M,1), 5);
power_polymer = zeros(size(M,1), 5);

for i = 1:size(M,1)
    base_folder = fullfile(base_path, M.folder{i});
    [dT_laser, ~, ~] = laser_temperature_calibration(base_folder, M.date{i}, M.idx(i), mdl, 34, 0, 1, save_path);
    dT_polymer(i,1:length(dT_laser)) = dT_laser;
    power_polymer(i,:) = M.power(i)*ones(1,5);
end

dT_tissue = nan(size(M,1), 5);
power_tissue = zeros(size(M,1), 5);

M1 = readtable('E:\data\polymer\pool_temp_calibration.xlsx', 'Sheet','tissue');
for i = 1:size(M1,1)
    base_folder = fullfile(base_path, M1.folder{i});
    [dT_laser, ~, ~] = laser_temperature_calibration(base_folder, M1.date{i}, M1.idx(i), mdl, 34, 0, 1, save_path);
    dT_tissue(i,1:length(dT_laser)) = dT_laser;
    power_tissue(i,:) = M1.power(i)*ones(1,5);
end

power_pool = [0:10:50,70,90, 150,200,250,300,350,400];
power_mW = [1.9,3.1,4.7,6.6,8.7,11.2,16.5,23.7,51.9,81.5,114,149,184, 220];

dTemp_polymer = cell(1, length(power_pool));
dTemp_tissue = cell(1, length(power_pool));
for i = 1:length(power_pool)
    idxs = find(power_polymer(:,1)==power_pool(i));
    if ~isempty(idxs)
        dTemp_polymer{i} = reshape(dT_polymer(idxs,:),[],1);
    else
        dTemp_polymer{i} = nan(5,1);
    end
end
for i = 1:length(power_pool)
    idxs = find(power_tissue(:,1)==power_pool(i));
    if ~isempty(idxs)
        dTemp_tissue{i} = reshape(dT_tissue(idxs,:),[],1);
    else
        dTemp_tissue{i} = nan(5,1);
    end
end
dT1 = padcat(dTemp_tissue{:});
dT2 = padcat(dTemp_polymer{:});
[data_table, within_design] = gen_table_for_ranova(power_mW(1:8), {dT2(:,1:8),dT1(:,1:8)});
rm = fitrm(data_table,'measurements1-measurements16 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);
lineplot_with_shaded_errorbar(power_mW, {padcat(dTemp_tissue{:})',padcat(dTemp_polymer{:})'})
save(fullfile(save_path, 'laser_temp_change.mat'), 'M','M1','power_polymer', 'dTemp_polymer', 'power_tissue', 'dTemp_tissue','power_pool', 'power_mW','dTemp_polymer', 'dTemp_tissue');
%%
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))
base_path = 'E:\data';
M = readtable('E:\data\polymer\pool_temp_calibration.xlsx', 'Sheet','polymer_invivo');
base_path = 'E:\data';
save_path = 'E:\data\polymer\temp_calibration';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end

load('E:\data\polymer\temp_calibration\temp_calibration_log.mat', 'mdl')

dT_polymer = nan(size(M,1), 5);
power_polymer = zeros(size(M,1), 5);

for i = 1:size(M,1)
    base_folder = fullfile(base_path, M.folder{i});
    [dT_laser, ~, ~] = laser_temperature_calibration(base_folder, M.date{i}, M.idx(i), mdl, 34, 0, 1, save_path);
    dT_polymer(i,1:length(dT_laser)) = dT_laser;
    power_polymer(i,:) = M.power(i)*ones(1,5);
end

dT_tissue = nan(size(M,1), 5);
power_tissue = zeros(size(M,1), 5);

M1 = readtable('E:\data\polymer\pool_temp_calibration.xlsx', 'Sheet','tissue_invivo');
for i = 1:size(M1,1)
    base_folder = fullfile(base_path, M1.folder{i});
    [dT_laser, ~, ~] = laser_temperature_calibration(base_folder, M1.date{i}, M1.idx(i), mdl, 34, 0, 1, save_path);
    dT_tissue(i,1:length(dT_laser)) = dT_laser;
    power_tissue(i,:) = M1.power(i)*ones(1,5);
end

power_pool = [10,40,60:10:140,200:50:500 ];
power_mW = [3.8,4.4,6.4,7.8,9.6,11.6,14,16.7,18.8,22.1,25.5,65,92.3,124.1,155.6,189.1,218.5,247.6];

dTemp_polymer = cell(1, length(power_pool));
dTemp_tissue = cell(1, length(power_pool));
for i = 1:length(power_pool)
    idxs = find(power_polymer(:,1)==power_pool(i));
    if ~isempty(idxs)
        dTemp_polymer{i} = reshape(dT_polymer(idxs,:),[],1);
    else
        dTemp_polymer{i} = nan(5,1);
    end
end
for i = 1:length(power_pool)
    idxs = find(power_tissue(:,1)==power_pool(i));
    if ~isempty(idxs)
        dTemp_tissue{i} = reshape(dT_tissue(idxs,:),[],1);
    else
        dTemp_tissue{i} = nan(5,1);
    end
end
dT1 = padcat(dTemp_tissue{:});
dT2 = padcat(dTemp_polymer{:});
[data_table, within_design] = gen_table_for_ranova(power_mW(1:11), {dT2(:,1:11),dT1(:,1:11)});
rm = fitrm(data_table,'measurements1-measurements22 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);
lineplot_with_shaded_errorbar(power_mW, {padcat(dTemp_tissue{:})',padcat(dTemp_polymer{:})'})
save(fullfile(save_path, 'laser_temp_change.mat'), 'M','M1','power_polymer', 'dTemp_polymer', 'power_tissue', 'dTemp_tissue','power_pool', 'power_mW','dTemp_polymer', 'dTemp_tissue');