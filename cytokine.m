%% cytokine
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

above_LOD_brain = {'IL-23', 'TNF', 'IL-1b', 'MCP1', 'IL-1a'};
above_LOD_skin = {'IL-23', 'TNF', 'IL-1a', 'IL-1b', 'IL-6', 'IL-27', 'IFN-b', 'MCP1', 'IL-10'};

filename = 'E:\data\polymer\cytokine\Jayant_Brain_Skin_Ifnlamamtory_Cytokines.xlsx';
save_path = 'E:\data\polymer\cytokine\stat';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end
brain = cell(1, length(above_LOD_brain));
skin = cell(1, length(above_LOD_brain));
brain_LOD = zeros(1, length(above_LOD_skin));
skin_LOD = zeros(1, length(above_LOD_skin));
colors = [[0,0,0];[119,177,204];[220,165,195];[138,62,130]];
colors = colors/256;
brain_p = zeros(1, length(above_LOD_brain));
skin_p = zeros(1, length(above_LOD_skin));
  
for i = 1:length(above_LOD_brain)
    brain{i} = xlsread(filename, above_LOD_brain{i}, 'A15:C17'); 
    temp = readtable(filename, 'Sheet',above_LOD_brain{i});
    LOD_temp = temp.Var2(23);
    numStr = regexp(string(LOD_temp{1}), '\d+\.?\d*', 'match');
    brain_LOD(i) = str2double(numStr);
    brain_p(i) = anova1(brain{i});
    
    barplot_with_datapoint_STD(brain{i}, colors, brain_LOD(i))
end

for i = 1:length(above_LOD_skin)
    skin{i} = xlsread(filename, above_LOD_skin{i}, 'H15:K17'); 
    temp = readtable(filename, 'Sheet',above_LOD_skin{i});
    LOD_temp = temp.Var2(23);
    numStr = regexp(string(LOD_temp{1}), '\d+\.?\d*', 'match');
    skin_LOD(i) = str2double(numStr);
    skin_p(i) = anova1(skin{i});
    
    barplot_with_datapoint_STD(skin{i}, colors, skin_LOD(i))
end

save(fullfile(save_path,'cytokine_stat.mat'),'brain','skin','above_LOD_brain','above_LOD_skin','brain_p','skin_p', 'brain_LOD', 'skin_LOD')