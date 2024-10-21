%% dendrite_puff_excitability
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

M = readtable('E:\data\polymer\pool_dstim.xlsx', 'Sheet','excitability');
% M = readtable('E:\data\dendritic patch\pool\pass_filter.xlsx', 'Sheet','morph');
base_path = 'E:\data';
save_path = 'E:\data\polymer\sub&supra_dend';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end
Rm = [];
Ihold = [];
Rm_fit = [];
Cm = [];
f_pre = [];
f_post = [];
AP_amp = [];
AP_width = [];
I = [];
for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder(i));
    [params_pre_temp, ~, f_pre_temp, I_temp, AP_amp_pre_temp,AP_width_pre_temp] = sub_and_supra_classify(data_path, M.date{i}, M.cell(i), M.pre(i), 0, 1, save_path);
    [params_post_temp, ~, f_post_temp, I_temp, AP_amp_post_temp,AP_width_post_temp] = sub_and_supra_classify(data_path, M.date{i}, M.cell(i), M.post(i),0,1,save_path);
    Rm = [Rm;[params_pre_temp(1), params_post_temp(1)]];
    Cm = [Cm;[params_pre_temp(5), params_post_temp(5)]];
    AP_amp = [AP_amp; [median(AP_amp_pre_temp), median(AP_amp_post_temp)]];
    AP_width = [AP_width; [median(AP_width_pre_temp), median(AP_width_post_temp)]];
    Rm_fit = [Rm_fit;[params_pre_temp(6), params_post_temp(6)]];
    Ihold = [Ihold;[params_pre_temp(7), params_post_temp(7)]];
    f_pre = [f_pre,f_pre_temp'];
    f_post = [f_post,f_post_temp'];
    I = [I, I_temp'];
end


colors = [[0,0,0];[119,176,203]/255];
boxplot_pairwise(Rm, colors)
% boxplot_pairwise(Rm_fit.*Cm/1e3)
boxplot_pairwise(Cm, colors)
barplot_pairwise(AP_amp, colors), ylim([80,140])
barplot_pairwise(-Ihold, colors)
lineplot_with_shaded_errorbar(I_temp, {f_pre(:,find(M.withttx==0)),f_post(:,find(M.withttx==0))}, colors), xlim([0,350])

[data_table, within_design] = gen_table_for_ranova(I_temp, {f_pre(:,find(M.withttx==0)),f_post(:,find(M.withttx==0))});

rm = fitrm(data_table,'measurements1-measurements24 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);
p = zeros(1, length(I_temp));
for i = 1:length(I_temp)
    [p(i),~] = signrank(f_pre(i,find(M.withttx==0))',f_post(i,find(M.withttx==0))');
end
save(fullfile(save_path,'excitability_pool.mat'),'Rm','Cm','AP_amp','AP_width', 'Ihold','I_temp','f_pre', 'f_post', 'anova_table', 'p', 'M')
%% d_stim_AP
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

M = readtable('E:\data\polymer\pool_dstim.xlsx', 'Sheet','stim_AP_lowcond');
base_path = 'E:\data';
base_path = 'E:\data';
save_path = 'E:\data\polymer\dstim';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end
power = [];
dist = [];
attenuation = [];
FR = zeros(size(M,1), 3);
soma_dist = [];

for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder(i));
    [pre_temp,stim_temp, post_temp, dur_stim_temp, no_stim_temp] = stim_AP(data_path, M.date{i}, M.cell(i), M.idx_f(i),0,1,save_path);
    attenuation = [attenuation;nanmean(dur_stim_temp')./nanmean(no_stim_temp')];
    FR(i,:) = [pre_temp,stim_temp, post_temp];
    power = [power, M.power(i)];
    dist = [dist,M.dist_to_polymer(i)];
    soma_dist = [soma_dist,M.dist_to_soma(i)];
end

attenuation_p = attenuation(dist==0,:);
power_p = power(dist==0);
[power_p,I] = sort(power_p);
attenuation_p = attenuation_p(I,:);

M1 = readtable('E:\data\polymer\pool_dstim.xlsx', 'Sheet','stim_AP_control');
save_path = 'E:\data\polymer\dstim_control';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end
power_control = [];
attenuation_control = [];
FR_control = zeros(size(M,1), 3);


for i = 1:size(M1, 1)
    data_path = fullfile(base_path,M1.folder(i));
    [pre_temp,stim_temp, post_temp, dur_stim_temp, no_stim_temp] = stim_AP(data_path, M1.date{i}, M1.cell(i), M1.idx_f(i),0,1,save_path);
    attenuation_control = [attenuation_control;nanmean(dur_stim_temp')./nanmean(no_stim_temp')];
    FR_control(i,:) = [pre_temp,stim_temp, post_temp];
    power_control = [power_control, M1.power(i)];
end

attenuation_p_control = attenuation_control;
power_p_control = power_control;
[power_p_control,I] = sort(power_p_control);
attenuation_p_control = attenuation_p_control(I,:);

colors = [[0,0,0];[119,176,203]/255];
edges = 0:2:ceil(max(power_p)/2)*2;
idxs = discretize(power_p,edges);
power_pool = edges(1:end-1)+(edges(2)-edges(1))/2;

attenuation_pool = cell(1, length(power_pool));
for i = 1:length(power_pool)
    attenuation_pool{i} = attenuation_p(find(idxs==i));
end

attenuation_pool_control = cell(1, length(power_pool));
for i = 1:length(power_pool)
    attenuation_pool_control{i} = attenuation_p_control(find(idxs==i));
end

a = padcat(attenuation_pool{:})';
b = padcat(attenuation_pool_control{:})';
errorbar_with_lines(power_pool([1:9,12]), {b([1:9,12],:), a([1:9,12],:)}, colors)

idx = find((power>=6)&(power<=7)&(dist==0));
FR_pool = FR([idx,25],:);
idx1 = find((power_control>=6)&(power_control<=7));
FR_pool_control = FR_control([idx1,1],:);
errorbar_with_lines([1:3], {FR_pool'}, colors(2,:))
hold on
for i = 1:size(FR_pool,1)
    plot([1:3], FR_pool(i,:),'Color',[0.5,0.5,0.5])
end
xlim([0.5,3.5])
errorbar_with_lines([1:3], {FR_pool_control'}, colors(1,:))
hold on
for i = 1:size(FR_pool_control,1)
    plot([1:3], FR_pool_control(i,:),'Color',[0.5,0.5,0.5])
end
xlim([0.5,3.5])


idx2 = find(power_control==149);
idx2(2) = 23;

FR_pool_control_hightemp = FR_control(idx2,:);
errorbar_with_lines([1:3], {FR_pool_control_hightemp'}, colors(1,:))
hold on
for i = 1:size(FR_pool_control_hightemp,1)
    plot([1:3], FR_pool_control_hightemp(i,:),'Color',[0.5,0.5,0.5])
end
xlim([0.5,3.5])

save_path = 'E:\data\polymer\dstim';
save(fullfile(save_path, 'dstim.mat'))
%% d stim subthreshold
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

M = readtable('E:\data\polymer\pool_dstim.xlsx', 'Sheet','stim_rest');
base_path = 'E:\data';
base_path = 'E:\data';
save_path = 'E:\data\polymer\dstim';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end
RMP = zeros(size(M,1),3);

for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder(i));
    [pre_temp, stim_temp, post_temp] = stim_sub(data_path, M.date{i}, M.cell(i), M.file_name{i}, M.idx_f(i),0,1,save_path);
    RMP(i,:) = [pre_temp(end),stim_temp(end), post_temp(end)];
end

M1 = readtable('E:\data\polymer\pool_dstim.xlsx', 'Sheet','stim_rest_control');
save_path = 'E:\data\polymer\dstim_control';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end
RMP_control = zeros(size(M,1),3);

for i = 1:size(M1, 1)
    data_path = fullfile(base_path,M1.folder(i));
    [pre_temp, stim_temp, post_temp] = stim_sub(data_path, M1.date{i}, M1.cell(i), M1.file_name{i}, M1.idx_f(i),0,1,save_path);
    RMP_control(i,:) = [pre_temp(end),stim_temp(end), post_temp(end)];
end
dV = RMP(:,2)-RMP(:,1);
dV_control = RMP_control(:,2)-RMP_control(:,1);
barplot_with_datapoint({dV_control,dV})

save_path = 'E:\data\polymer\dstim';
save(fullfile(save_path, 'dstim_hyperpolarization.mat'))
%%
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

M = readtable('E:\data\polymer\pool_dstim.xlsx', 'Sheet','VC');
base_path = 'E:\data';
save_path = 'E:\data\polymer\dstim';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end

Vstep = nan(size(M,1), 10);
dI = nan(size(M,1), 10);

for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder(i));
    [pre_temp,stim_temp, post_temp, Vstep_temp] = stim_VC(data_path, M.date{i}, M.cell(i), M.idx_f(i),0,1,save_path);
    Vstep(i,1:length(Vstep_temp)) = Vstep_temp;
    a = stim_temp(:,end)-pre_temp(:,end);
    dI(i,1:length(Vstep_temp))  = a';
end
dI_temp = dI(4,:);
dI(4,:) = nan;
dI(4,3:8) = dI_temp(1:6);
colors = [[0,0,0];[119,176,203]/255];
ft = fittype('a*x+b', 'independent', 'x');
curve = fit(Vstep(1,:)', nanmean(dI, 1)', ft, 'Start', [0, 0]);

errorbar_with_fitcurve(Vstep(1,:), {dI'},{curve}, colors(2,:))
box on
xlim([-95,-40])
x1 = xlim();
xplot = x1(1):x1(2);
plot(xplot, zeros(size(xplot)), 'k--')
save(fullfile(save_path, 'dstim_VC.mat'))
%%
% %% dendrite_puff_excitability_with_light_modulation
% clear all
% close all
% 
% addpath(genpath(fullfile(pwd,'main')))
% addpath(genpath(fullfile(pwd,'plotting')))
% 
% M = readtable('Z:\data\shulan\polymer\pool_dstim.xlsx', 'Sheet','excitability_w_light');
% % M = readtable('E:\data\dendritic patch\pool\pass_filter.xlsx', 'Sheet','morph');
% base_path = 'Z:\data\shulan';
% Rm = [];
% Ihold = [];
% Rm_fit = [];
% Cm = [];
% f_pre = [];
% f_post = [];
% AP_amp = [];
% I = [];
% for i = 1:size(M, 1)
%     data_path = fullfile(base_path,M.folder(i));
%     [params_pre_temp, ~, f_pre_temp, I_temp, AP_amp_pre_temp] = sub_and_supra_classify(data_path, M.date{i}, M.cell(i), M.pre(i));
%     [params_post_temp, ~, f_post_temp, I_temp, AP_amp_post_temp] = sub_and_supra_classify(data_path, M.date{i}, M.cell(i), M.post(i));
%     Rm = [Rm;[params_pre_temp(1), params_post_temp(1)]];
%     Cm = [Cm;[params_pre_temp(5), params_post_temp(5)]];
%     AP_amp = [AP_amp; [median(AP_amp_pre_temp), median(AP_amp_post_temp)]];
%     Rm_fit = [Rm_fit;[params_pre_temp(6), params_post_temp(6)]];
%     Ihold = [Ihold;[params_pre_temp(7), params_post_temp(7)]];
%     f_pre = [f_pre,f_pre_temp'];
%     f_post = [f_post,f_post_temp'];
%     I = [I, I_temp'];
% end
% 
% 
% colors = [[0,0,0];[119,176,203]/255];
% boxplot_pairwise(Rm, colors)
% % boxplot_pairwise(Rm_fit.*Cm/1e3)
% boxplot_pairwise(Cm, colors)
% barplot_pairwise(AP_amp, colors), ylim([80,140])
% barplot_pairwise(-Ihold, colors)
% lineplot_with_shaded_errorbar(I_temp, {f_pre(:,find(M.withttx==0)),f_post(:,find(M.withttx==0))}, colors), xlim([0,350])