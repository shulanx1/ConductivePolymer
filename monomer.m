%% initiated monomer puff bAP_Ca
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

M = readtable('E:\data\polymer\pool_monomer.xlsx', 'Sheet','bAP');
base_path = 'E:\data';
save_path = 'E:\data\polymer\monomer_puff_bAP';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end
baseline_f = zeros(size(M,1),4);
evoke_df = zeros(size(M,1),4);

for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder{i});
    load(fullfile(data_path, M.date{i}, sprintf('cell%d',M.cell(i)),'Ca_output', sprintf('data_bAP_incremental_%d.mat',M.pre(i))));
    baseline_f(i,[1,2]) = mean(raw_f{1}(:,1:stim_idx_Ca(1)),2)./mean(raw_f{2}(:,1:stim_idx_Ca(1)),2);
    evoke_df(i,[1,2]) = Ca_amp(:,3)';
    load(fullfile(data_path, M.date{i}, sprintf('cell%d',M.cell(i)),'Ca_output', sprintf('data_bAP_incremental_%d.mat',M.post(i))));
    baseline_f(i,[3,4]) = mean(raw_f{1}(:,1:stim_idx_Ca(1)),2)./mean(raw_f{2}(:,1:stim_idx_Ca(1)),2);
    evoke_df(i,[3,4]) = Ca_amp(:,3)';
    
end

colors = [[0,0,0];[119,176,203]/255]; % color for nPBDF
boxplot_pairwise(evoke_df(:,[1,3]), colors)
boxplot_pairwise(baseline_f(:,[1,3]), colors)
save(fullfile(save_path, 'monomer_puff_bAP_pool.mat'),'baseline_f','evoke_df')

%% excitability
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

M = readtable('E:\data\polymer\pool_monomer.xlsx', 'Sheet','excitability');
base_path = 'E:\data';
save_path = 'E:\data\polymer\monomer_puff_excitability';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end
decay = [];
AP_ratio = [];
Rm = [];
RMP = [];
Ihold = [];
Cm = [];
Rm_fit = [];
f_pre = [];
f_post = [];
Rin1 = [];
I = [];
for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder(i));
    [params_pre_temp,decay_R_pre,AP_ratio_pre, f_pre_temp, I_temp, Rin_pooled1_pre,~] = sub_and_supra_dend(data_path, M.date{i}, M.cell(i), 'sub and supra_Ch2', M.pre(i), 0, 1, save_path);
    [params_post_temp,decay_R_post, AP_ratio_post,f_post_temp, I_temp, Rin_pooled1_post,~] = sub_and_supra_dend(data_path, M.date{i}, M.cell(i), 'sub and supra_Ch2', M.post(i), 0, 1, save_path);
    decay = [decay;[decay_R_pre, decay_R_post]];
    AP_ratio = [AP_ratio;[AP_ratio_pre, AP_ratio_post]];
    Rm = [Rm;[params_pre_temp(1), params_post_temp(1)]];
    Cm = [Cm;[params_pre_temp(5), params_post_temp(5)]];
    Rm_fit = [Rm_fit;[params_pre_temp(6), params_post_temp(6)]];
    RMP = [RMP;[params_pre_temp(8), params_post_temp(8)]];
    Ihold = [Ihold;[params_pre_temp(7), params_post_temp(7)]];
    Rin1 = [Rin1;[Rin_pooled1_pre, Rin_pooled1_post]];
    f_pre = [f_pre,f_pre_temp'];
    f_post = [f_post,f_post_temp'];
    I = [I, I_temp'];
end

colors = [[0,0,0];[119,176,203]/255]; % color for Na rich polymer
% boxplot_pairwise(decay, colors)
boxplot_pairwise(Rm, colors)
boxplot_pairwise(Cm, colors)
boxplot_pairwise(Cm.*Rm_fit/1e3, colors)
lineplot_with_shaded_errorbar(I_temp, {f_pre,f_post}, colors), xlim([0,350])
[data_table, within_design] = gen_table_for_ranova(I_temp((I_temp>0)&(I_temp<=350)), {f_pre((I_temp>0)&(I_temp<=350),:),f_post((I_temp>0)&(I_temp<=350),:)});

rm = fitrm(data_table,'measurements1-measurements14 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);
save(fullfile(save_path,'monomer_puff_excitability_pool.mat'),'decay','AP_ratio','Rm', 'Cm', 'M','RMP', 'Ihold', 'Rm_fit','Rin1', 'I', 'f_pre', 'f_post')
%% excitability_ch1
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

M = readtable('E:\data\polymer\pool_monomer.xlsx', 'Sheet','excitability_Ch1');
base_path = 'E:\data';
save_path = 'E:\data\polymer\monomer_puff_excitability_Ch1';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end
decay = [];
AP_ratio = [];
Rm = [];
Rm_fit = [];
RMP = [];
Ihold = [];
Cm = [];
f_pre = [];
f_post = [];
Rin2 = [];
I = [];
for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder(i));
    [params_pre_temp,decay_R_pre,AP_ratio_pre, f_pre_temp, I_temp, Rin_pooled1_pre,Rin_pooled2_pre] = sub_and_supra_dend_ch1(data_path, M.date{i}, M.cell(i), 'sub and supra_Ch1', M.pre(i), 0, 1, save_path);
    [params_post_temp,decay_R_post, AP_ratio_post,f_post_temp, I_temp, Rin_pooled1_post,Rin_pooled2_post] = sub_and_supra_dend_ch1(data_path, M.date{i}, M.cell(i), 'sub and supra_Ch1', M.post(i), 0, 1, save_path);
    decay = [decay;[decay_R_pre, decay_R_post]];
    AP_ratio = [AP_ratio;[AP_ratio_pre, AP_ratio_post]];
    Rm = [Rm;[params_pre_temp(1), params_post_temp(1)]];
    Cm = [Cm;[params_pre_temp(5), params_post_temp(5)]];
    Rm_fit = [Rm_fit;[params_pre_temp(6), params_post_temp(6)]];
    RMP = [RMP;[params_pre_temp(8), params_post_temp(8)]];
    Ihold = [Ihold;[params_pre_temp(7), params_post_temp(7)]];
    Rin2 = [Rin2;[Rin_pooled2_pre, Rin_pooled2_post]];
    f_pre = [f_pre,f_pre_temp'];
    f_post = [f_post,f_post_temp'];
    I = [I, I_temp'];
end

colors = [[0,0,0];[119,176,203]/255]; % color for Na rich polymer
boxplot_pairwise(decay, colors)
boxplot_pairwise(Rm, colors)
boxplot_pairwise(Cm, colors)
boxplot_pairwise(Cm.*Rm_fit/1e3, colors)
lineplot_with_shaded_errorbar(I_temp, {f_pre,f_post}, colors), xlim([100,450])
[data_table, within_design] = gen_table_for_ranova(I_temp((I_temp>100)&(I_temp<=450)), {f_pre((I_temp>100)&(I_temp<=450),:),f_post((I_temp>100)&(I_temp<=450),:)});

rm = fitrm(data_table,'measurements1-measurements14 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);
save(fullfile(save_path,'monomer_puff_excitability_pool.mat'),'decay','AP_ratio','Rm', 'Cm', 'M','RMP', 'Ihold', 'Rm_fit','Rin2', 'I', 'f_pre', 'f_post')
%% zap_ch2
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

M = readtable('E:\data\polymer\pool_monomer.xlsx', 'Sheet','zap_Ch2');
base_path = 'E:\data';
save_path = 'E:\data\polymer\monomer_puff_zap_Ch2';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end
Z2_pre = [];
Zprop_pre = [];
Z2_post = [];
Zprop_post = [];
f = [];

for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder{i});
    [f_temp, Z1_pre_temp, Z2_pre_temp] = zap(data_path, M.date{i}, M.cell(i), M.pre(i), 0, 1, save_path);
    [f_temp, Z1_post_temp, Z2_post_temp] = zap(data_path, M.date{i}, M.cell(i), M.post(i), 0, 1, save_path);
    Z2_pre(:,i) = Z2_pre_temp;
    Z2_post(:,i) = Z2_post_temp;
    Zprop_pre(:,i) = Z1_pre_temp./Z2_pre_temp;
    Zprop_post(:,i) = Z1_post_temp./Z2_post_temp;    
    f = [f, f_temp];
end

colors = [[0,0,0];[119,176,203]/255]; % color for nPBDF
df = f_temp(2)-f_temp(1);
lineplot_with_shaded_errorbar(f_temp(ceil(1/df):floor(14/df)), {movmean(Z2_pre(ceil(1/df):floor(14/df),:),5,1),movmean(Z2_post(ceil(1/df):floor(14/df),:),5,1)}, colors)
lineplot_with_shaded_errorbar(f_temp(ceil(1/df):floor(14/df)), {movmean(Zprop_pre(ceil(1/df):floor(14/df),:),5,1),movmean(Zprop_post(ceil(1/df):floor(14/df),:),5,1)}, colors)

a = f_temp(ceil(1/df):floor(14/df));
sample_idx = 1:30:length(a);

data_pre = movmean(Z2_pre(ceil(1/df):floor(14/df),:),5,1);
data_post = movmean(Z2_post(ceil(1/df):floor(14/df),:),5,1);
[data_table, within_design] = gen_table_for_ranova(a(sample_idx), {data_pre(sample_idx,:), data_post(sample_idx,:)});
rm = fitrm(data_table,'measurements1-measurements16 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);

data_pre = movmean(Zprop_pre(ceil(1/df):floor(14/df),:),5,1);
data_post = movmean(Zprop_post(ceil(1/df):floor(14/df),:),5,1);
[data_table, within_design] = gen_table_for_ranova(a(sample_idx), {data_pre(sample_idx,:), data_post(sample_idx,:)});
rm = fitrm(data_table,'measurements1-measurements16 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);

save(fullfile(save_path,'monomer_puff_zapch2_pool.mat'),'Z2_pre','Z2_post','Zprop_pre', 'Zprop_post', 'f','df', 'M')
%% zap_ch1
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

M = readtable('E:\data\polymer\pool_monomer.xlsx', 'Sheet','zap_Ch1');
base_path = 'E:\data';
save_path = 'E:\data\polymer\monomer_puff_zap_Ch1';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end
Z1_pre = [];
Zprop_pre = [];
Z1_post = [];
Zprop_post = [];
f = [];

for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder{i});
    [f_temp, Z1_pre_temp, Z2_pre_temp] = zap_Ch1(data_path, M.date{i}, M.cell(i), M.pre(i), 0, 1, save_path);
    [f_temp, Z1_post_temp, Z2_post_temp] = zap_Ch1(data_path, M.date{i}, M.cell(i), M.post(i), 0, 1, save_path);
    Z1_pre(:,i) = Z1_pre_temp;
    Z1_post(:,i) = Z1_post_temp;
    Zprop_pre(:,i) = Z2_pre_temp./Z1_pre_temp;
    Zprop_post(:,i) = Z2_post_temp./Z1_post_temp;    
    f = [f, f_temp];
end

colors = [[0,0,0];[119,176,203]/255]; % color for nPBDF
df = f_temp(2)-f_temp(1);
lineplot_with_shaded_errorbar(f_temp(ceil(1/df):floor(14/df)), {movmean(Z1_pre(ceil(1/df):floor(14/df),2:4),5,1),movmean(Z1_post(ceil(1/df):floor(14/df),2:4),5,1)}, colors)
lineplot_with_shaded_errorbar(f_temp(ceil(1/df):floor(14/df)), {movmean(Zprop_pre(ceil(1/df):floor(14/df),2:4),5,1),movmean(Zprop_post(ceil(1/df):floor(14/df),2:4),5,1)}, colors)

a = f_temp(ceil(1/df):floor(14/df));
sample_idx = 1:30:length(a);

data_pre = movmean(Z1_pre(ceil(1/df):floor(14/df),:),5,1);
data_post = movmean(Z1_post(ceil(1/df):floor(14/df),:),5,1);
[data_table, within_design] = gen_table_for_ranova(a(sample_idx), {data_pre(sample_idx,:), data_post(sample_idx,:)});
rm = fitrm(data_table,'measurements1-measurements16 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);

data_pre = movmean(Zprop_pre(ceil(1/df):floor(14/df),:),5,1);
data_post = movmean(Zprop_post(ceil(1/df):floor(14/df),:),5,1);
[data_table, within_design] = gen_table_for_ranova(a(sample_idx), {data_pre(sample_idx,:), data_post(sample_idx,:)});
rm = fitrm(data_table,'measurements1-measurements16 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);

save(fullfile(save_path,'monomer_puff_zapch1_pool.mat'),'Z1_pre','Z1_post','Zprop_pre', 'Zprop_post', 'f','df', 'M')

%% excitability_dend_inhibition
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

M = readtable('E:\data\polymer\pool_monomer.xlsx', 'Sheet','excitability_dend_inhibition');

base_path = 'E:\data';
save_path = 'E:\data\polymer\sub&supra_dend_inhibition';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end

f_con = [];
f_5nS = [];
f_10nS = [];
f_20nS = [];

I = [];
for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder(i));
    if ~isnan(M.control(i))
        [~, ~, f_pre_temp, I_temp, ~, ~] = sub_and_supra_classify(data_path, M.date{i}, M.cell(i), M.control(i),0,1,save_path);
        f_con = [f_con, f_pre_temp'];
    end
    if ~isnan(M.nS5(i))
        [~, ~, f_5nS_temp, I_temp, ~, ~] = sub_and_supra_classify(data_path, M.date{i}, M.cell(i), M.nS5(i),0,1,save_path);
        f_5nS = [f_5nS, f_5nS_temp'];
    end
    if ~isnan(M.nS10(i))
        [~, ~, f_10nS_temp, I_temp, ~, ~] = sub_and_supra_classify(data_path, M.date{i}, M.cell(i), M.nS10(i),0,1,save_path);
        f_10nS = [f_10nS, f_10nS_temp'];
    end
    if ~isnan(M.nS20(i))
        [~, ~, f_20nS_temp, I_temp, ~, ~] = sub_and_supra_classify(data_path, M.date{i}, M.cell(i), M.nS20(i),0,1,save_path);
        f_20nS = [f_20nS, f_20nS_temp'];
    end
    I = [I, I_temp'];
end
colors = [[30,30,30];[250,172,144];[243,119,129];[187,107,132]]/255;
lineplot_with_shaded_errorbar(I_temp, {f_con,f_5nS, f_10nS, f_20nS}, colors), xlim([0,450]);

[data_table, within_design] = gen_table_for_ranova(I_temp((I_temp>0)&(I_temp<=450)), {f_con((I_temp>0)&(I_temp<=450),:),f_5nS((I_temp>0)&(I_temp<=450),:),f_10nS((I_temp>0)&(I_temp<=450),:),[nan(9,1),f_20nS((I_temp>0)&(I_temp<=450),:)]});

rm = fitrm(data_table,'measurements1-measurements36 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);

save(fullfile(save_path,'excitability_dend_inhibition_pool.mat'),'I_temp','f_con', 'f_5nS', 'f_10nS','f_20nS', 'M')
%% Halo optogenetic

clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

M = readtable('E:\data\polymer\pool_monomer.xlsx', 'Sheet','optogenetic');
base_path = 'E:\data';

save_path = 'E:\data\polymer\NpHR_optogenetic';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end
power = [];
dist = [];
attenuation = [];
FR_LED = zeros(size(M,1), 3);
FR_soma = zeros(size(M,1), 3);
FR_dend = zeros(size(M,1), 3);
FR_soma_high = zeros(size(M,1), 3);
FR_dend_high = zeros(size(M,1), 3);
idx_rmv_LED = [];
idx_rmv_soma = [];
idx_rmv_dend = [];
idx_rmv_soma_high = [];
idx_rmv_dend_high = [];
for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder(i));
    if ~isnan(M.idx_opto(i))
        [pre_temp,stim_temp, post_temp, dur_stim_temp, no_stim_temp] = stim_AP(data_path, M.date{i}, M.cell(i), M.idx_opto(i),0,1,save_path);
        FR_LED(i,:) = [pre_temp,stim_temp, post_temp];
    else
        idx_rmv_LED = [idx_rmv_LED, i];
    end
    if ~isnan(M.idx_soma(i))
        [pre_temp,stim_temp, post_temp, dur_stim_temp, no_stim_temp] = stim_AP(data_path, M.date{i}, M.cell(i), M.idx_soma(i),0,1,save_path);
        FR_soma(i,:) = [pre_temp,stim_temp, post_temp];
    else
        idx_rmv_soma = [idx_rmv_soma, i];
    end
    if ~isnan(M.idx_dend(i))
        [pre_temp,stim_temp, post_temp, dur_stim_temp, no_stim_temp] = stim_AP(data_path, M.date{i}, M.cell(i), M.idx_dend(i),0,1,save_path);
        FR_dend(i,:) = [pre_temp,stim_temp, post_temp];
    else
        idx_rmv_dend = [idx_rmv_dend, i];
    end
        if ~isnan(M.idx_soma_high(i))
        [pre_temp,stim_temp, post_temp, dur_stim_temp, no_stim_temp] = stim_AP(data_path, M.date{i}, M.cell(i), M.idx_soma_high(i),0,1,save_path);
        FR_soma_high(i,:) = [pre_temp,stim_temp, post_temp];
    else
        idx_rmv_soma_high = [idx_rmv_soma_high, i];
    end
    if ~isnan(M.idx_dend_high(i))
        [pre_temp,stim_temp, post_temp, dur_stim_temp, no_stim_temp] = stim_AP(data_path, M.date{i}, M.cell(i), M.idx_dend_high(i),0,1,save_path);
        FR_dend_high(i,:) = [pre_temp,stim_temp, post_temp];
    else
        idx_rmv_dend_high = [idx_rmv_dend_high, i];
    end
end

FR_LED(idx_rmv_LED,:) = [];
FR_soma(idx_rmv_soma,:) = [];
FR_dend(idx_rmv_dend,:) = [];
FR_soma_high(idx_rmv_soma_high,:) = [];
FR_dend_high(idx_rmv_dend_high,:) = [];

colors = [[239,153,75];[220,165,195];[175,178,207];[137,62,129];[124,135,181]]/255;
errorbar_with_lines([1:3], {FR_LED'}, colors(1,:))
hold on
for i = 1:size(FR_LED,1)
    plot([1:3], FR_LED(i,:),'Color',[0.5,0.5,0.5])
end
xlim([0.5,3.5])
ylim([0,12])

errorbar_with_lines([1:3], {FR_soma'}, colors(2,:))
hold on
for i = 1:size(FR_soma,1)
    plot([1:3], FR_soma(i,:),'Color',[0.5,0.5,0.5])
end
xlim([0.5,3.5])
ylim([0,12])

errorbar_with_lines([1:3], {FR_dend'}, colors(3,:))
hold on
for i = 1:size(FR_dend,1)
    plot([1:3], FR_dend(i,:),'Color',[0.5,0.5,0.5])
end
xlim([0.5,3.5])
ylim([0,12])

errorbar_with_lines([1:3], {FR_soma_high'}, colors(4,:))
hold on
for i = 1:size(FR_soma_high,1)
    plot([1:3], FR_soma_high(i,:),'Color',[0.5,0.5,0.5])
end
xlim([0.5,3.5])
ylim([0,12])

errorbar_with_lines([1:3], {FR_dend_high'}, colors(5,:))
hold on
for i = 1:size(FR_dend_high,1)
    plot([1:3], FR_dend_high(i,:),'Color',[0.5,0.5,0.5])
end
xlim([0.5,3.5])
ylim([0,12])
save(fullfile(save_path,'NpHR_optogenetic_pool.mat'),'FR_LED','FR_soma', 'FR_dend', 'FR_soma_high','FR_dend_high', 'M')
%% uncaging
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

M = readtable('E:\data\polymer\pool_monomer.xlsx', 'Sheet','uncaging');
base_path = 'E:\data';

save_path = 'E:\data\polymer\uncaging';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end

amp_single = cell(1,2);
halfwidth_single = cell(1,2);
risetime_single = cell(1,2);

amp_sum_pool = cell(1,2);

for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder(i));
    [amp, halfwidth, risetime, maxdv] = single_spine_serise(data_path, M.date{i}, M.cell(i), M.branch(i), M.single_idx(i),0,1,save_path);
    amp_single{M.group(i)} = [amp_single{M.group(i)},amp];
    halfwidth_single{M.group(i)} = [halfwidth_single{M.group(i)},halfwidth];
    risetime_single{M.group(i)} = [risetime_single{M.group(i)},risetime];
    [amp_sum, halfwidth_sum, risetime_sum, maxdv_sum] = summation(data_path, M.date{i}, M.cell(i), M.branch(i), M.sum_idx(i),0,1,save_path);
    amp_sum_pool{M.group(i)} = [amp_sum_pool{M.group(i)},amp_sum'];
end

idx_rmv{1} = find((amp_single{1}>1));
idx_rmv{2} = find((amp_single{2}>1));

for i = 1:2
    halfwidth_single1{i} = halfwidth_single{i};
    halfwidth_single1{i}(idx_rmv{i}) = [];
    risetime_single1{i} = risetime_single{i};
    risetime_single1{i}(idx_rmv{i}) = [];
end
violinplot_with_datapoint(amp_single)
violinplot_with_datapoint(halfwidth_single1)

idx_rmv = cell(1,2);
for i = 1:2
    for j = 1:size(amp_sum_pool{i},2)
    amp_norm{i}(:,j) = amp_sum_pool{i}(:,j)/max(amp_sum_pool{i}(:,j));
    end
    idx_rmv{i} = find(amp_norm{i}(1,:)>=0.35);
    amp_norm1{i} = amp_norm{i};
    amp_norm1{i}(:,idx_rmv{i}) = [];
end


ft = fittype('d+b/(1+exp(-(x-a)/k))', 'independent', 'x');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [4, 1, 0,1];
opts.Lower = [0, -inf, 0, -inf];
opts.Upper = [8, inf, 1, inf];
curve_con = fit([1:8]',mean(amp_norm1{1},2), ft, opts);
curve_mono = fit([1:8]',mean(amp_norm1{2},2), ft, opts);

colors = [[0,0,0];[119,176,203]/255]; % color for Na rich polymer

errorbar_with_fitcurve([1:8], amp_norm1, {curve_con,curve_mono},colors),box on
[measurements, gp1, gp2] = gen_table_for_unbalanced_anova([1:8], amp_norm1);
[p,~,stats] = anovan(measurements,{gp1,gp2},'model','interaction');

save(fullfile(save_path, 'unaging_pool.mat'),'amp_norm', 'amp_norm1', 'halfwidth_single1', 'amp_single', 'amp_sum_pool')