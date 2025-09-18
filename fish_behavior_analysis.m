clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

datapath = 'E:\data\polymer\fish_behavior_videos';
files = dir(fullfile(datapath, '*.csv'));

base_path = 'E:\data';
save_path = 'E:\data\polymer\fish_behavior_videos\stat';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end

control = {};
i_control = 1;
monomer = {};
i_monomer = 1;
for i = 1:length(files)
    if contains(files(i).name, 'control')
        temp = csvread(fullfile(files(i).folder, files(i).name),3 ,1);
        control{i_control}.csv_file = fullfile(files(i).folder, files(i).name);
        control{i_control}.raw_coor = temp(:,[1,2,4,5,7,8]);
        i_control = i_control + 1;
    elseif contains(files(i).name, 'monomer')
        temp = csvread(fullfile(files(i).folder, files(i).name),3 ,1);
        monomer{i_monomer}.csv_file = fullfile(files(i).folder, files(i).name);
        monomer{i_monomer}.raw_coor = temp(:,[1,2,4,5,7,8]);
        i_monomer = i_monomer + 1;
    end        
end

%%
max_curve = cell(1,2);
turn_time = cell(1,2);
head_direction = cell(1,2);
tail_direction = cell(1,2);
curvature = cell(1,2);
for i = 1:length(control)
    control{i} = analyze_fish_video(control{i});
    max_curve{1}(i) = control{i}.max_curve*100;
    turn_time{1}(i) = control{i}.turn_time*1e3;
    head_direction{1}(:,i) = real(rad2deg(control{i}.directions(:,1)));
    tail_direction{1}(:,i) = real(rad2deg(control{i}.directions(:,2)));
    curvature{1}(:,i) = control{i}.curvature*100;
end

for i = 1:length(monomer)
    monomer{i} = analyze_fish_video(monomer{i});
    max_curve{2}(i) = monomer{i}.max_curve*100;
    turn_time{2}(i) = monomer{i}.turn_time*1e3;
    head_direction{2}(:,i) = real(rad2deg( monomer{i}.directions(:,1)));
    tail_direction{2}(:,i) = real(rad2deg( monomer{i}.directions(:,2)));
    curvature{2}(:,i) =  monomer{i}.curvature*100;
end
violinplot_with_datapoint(max_curve)
violinplot_with_datapoint(turn_time)

lineplot_with_shaded_errorbar([0:size(head_direction{1},1)-1]*(1/750)*1e3, head_direction);box on
lineplot_with_shaded_errorbar([0:size(head_direction{1},1)-1]*(1/750)*1e3, tail_direction);box on
lineplot_with_shaded_errorbar([0:size(head_direction{1},1)-1]*(1/750)*1e3, curvature);box on

save(fullfile(save_path,'fish_behavior_pool.mat'), 'control', 'monomer', 'max_curve','turn_time', 'head_direction','tail_direction','curvature')