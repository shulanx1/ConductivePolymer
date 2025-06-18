%% excitability
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

M = readtable('E:\data\polymer\pool_opto.xlsx', 'Sheet','excitability');
% M = readtable('E:\data\dendritic patch\pool\pass_filter.xlsx', 'Sheet','morph');
base_path = 'E:\data';
save_path = 'E:\data\polymer\sub&supra_opto';
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
    [params_pre_temp, ~, f_pre_temp, I_temp, AP_amp_pre_temp, AP_width_pre_temp] = sub_and_supra_classify(data_path, M.date{i}, M.cell(i), M.pre(i),0,1,save_path, 0);
    [params_post_temp, ~, f_post_temp, I_temp, AP_amp_post_temp, AP_width_post_temp] = sub_and_supra_classify(data_path, M.date{i}, M.cell(i), M.post(i),0,1,save_path, 1);
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



colors = [[0,0,0];[119,176,203]/255]; % color for Na rich polymer
% colors = [[0,0,0];[159,186,149]/255]; % color for K rich polymer
boxplot_pairwise(Rm, colors)
% boxplot_pairwise(Rm_fit.*Cm/1e3)
boxplot_pairwise(Cm, colors)
barplot_pairwise(AP_amp((sum(isnan(AP_amp),2)==0),:), colors), ylim([50,100])
barplot_pairwise(AP_width((sum(isnan(AP_width),2)==0),:), colors), ylim([1,3.5])
barplot_pairwise(-Ihold, colors)
lineplot_with_shaded_errorbar(I_temp, {f_pre,f_post}, colors), xlim([0,350])
[data_table, within_design] = gen_table_for_ranova(I_temp, {f_pre',f_post'});

rm = fitrm(data_table,'measurements1-measurements24 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);
p = zeros(1, length(I_temp));
for i = 1:length(I_temp)
    [p(i),~] = signrank(f_pre(i,:)',f_post(i,:)');
end
save(fullfile(save_path,'excitability_opto_pool.mat'),'Rm','Cm','AP_amp','AP_width', 'Ihold','I_temp','f_pre', 'f_post', 'anova_table', 'p', 'M')
%% na activation
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))
M = readtable('E:\data\polymer\pool_opto.xlsx', 'Sheet','na_activation');
base_path = 'E:\data';
save_path = 'E:\data\polymer\na_activation_opto';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end

Na_amp_pre = zeros(size(M,1),11);
Na_amp_post = zeros(size(M,1),11);
Vstep = zeros(size(M,1),11);
for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder(i));
    [Na_amp_pre(i,:), ~, ~] = na_activation(data_path, M.date{i}, M.cell(i), M.pre(i), nan, 0, 1, save_path, 0);
    [Na_amp_post(i,:), ~, Vstep(i,:)] = na_activation(data_path, M.date{i}, M.cell(i), M.post(i), nan, 0,1,save_path, 1);
end
G_pre = Na_amp_pre./(Vstep-53.9);
G_post = Na_amp_post./(Vstep-53.9);
ft = fittype('G*(1/(1+exp(-(x-Vh)/k)))', 'independent', 'x');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [60,-30,3];
opts.Lower = [0,-70,0];
opts.Upper = [Inf, 0, Inf];
opts.Exclude = [11];
curve_pre = fit(Vstep(1,:)',mean(G_pre,1)', ft, opts);
curve_post = fit(Vstep(1,:)',mean(G_post,1)', ft, opts);

colors = [[0,0,0];[119,176,203]/255]; % color for Na rich polymer
% colors = [[0,0,0];[159,186,149]/255]; % color for K rich polymer
% colors = [[0,0,0];[195,129,168]/255]; % color for pedotPSS
errorbar_with_fitcurve(Vstep(1,:), {G_pre', G_post'}, {curve_pre,curve_post},colors),xlim([-90,0]),box on
boxplot_pairwise(-[reshape(Na_amp_pre(:,6:11),[],1)/1e3,reshape(Na_amp_post(:,6:11)/1e3,[],1)], colors)
[data_table, within_design] = gen_table_for_ranova(Vstep(1,:), {G_pre, G_post});

rm = fitrm(data_table,'measurements1-measurements22 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);
postHoc = multcompare(rm,'treatment');
postHoc_g = postHoc.pValue(1);
p = zeros(1, length(Vstep));
for i = 1:length(Vstep)
    [p(i),~] = signrank(G_pre(:,i), G_post(:,i));
end

save(fullfile(save_path,'na_activation_opto_pool.mat'),'G_pre','G_post','Vstep','anova_table','p','M')
%% ca activation (wo ttx)
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))
M = readtable('E:\data\polymer\pool_opto.xlsx', 'Sheet','ca_activation');
base_path = 'E:\data';
save_path = 'E:\data\polymer\ca_activation_opto';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end

Na_amp_pre = zeros(size(M,1),15);
Na_amp_post = zeros(size(M,1),15);
Vstep = zeros(size(M,1),15);
for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder(i));
    [Na_amp_pre(i,:), ~, ~] = ca_activation(data_path, M.date{i}, M.cell(i), M.pre(i), nan, 0, 1, save_path, 0);
    [Na_amp_post(i,:), ~, Vstep(i,:)] = ca_activation(data_path, M.date{i}, M.cell(i), M.post(i), nan, 0,1,save_path, 1);
end
G_pre = Na_amp_pre./(Vstep-53.9);
G_post = Na_amp_post./(Vstep-53.9);
ft = fittype('G*(1/(1+exp(-(x-Vh)/k)))', 'independent', 'x');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [60,-30,3];
opts.Lower = [0,-70,0];
opts.Upper = [Inf, 0, Inf];
opts.Exclude = [11:15];
curve_pre = fit(Vstep(1,:)',mean(G_pre,1)', ft, opts);
curve_post = fit(Vstep(1,:)',mean(G_post,1)', ft, opts);

colors = [[0,0,0];[119,176,203]/255]; % color for Na rich polymer
% colors = [[0,0,0];[159,186,149]/255]; % color for K rich polymer
% colors = [[0,0,0];[195,129,168]/255]; % color for pedotPSS
errorbar_with_fitcurve(Vstep(1,:), {G_pre', G_post'}, {curve_pre,curve_post},colors),xlim([-90,0]),box on
boxplot_pairwise(-[reshape(Na_amp_pre(:,6:11),[],1)/1e3,reshape(Na_amp_post(:,6:11)/1e3,[],1)], colors)
[data_table, within_design] = gen_table_for_ranova(Vstep(1,1:11), {G_pre(:,1:11), G_post(:,1:11)});

rm = fitrm(data_table,'measurements1-measurements22 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);
postHoc = multcompare(rm,'treatment');
postHoc_g = postHoc.pValue(1);
p = zeros(1, length(Vstep));
for i = 1:length(Vstep)
    [p(i),~] = signrank(G_pre(:,i), G_post(:,i));
end

save(fullfile(save_path,'ca_activation_opto_pool.mat'),'G_pre','G_post','Vstep','anova_table','p','M')
%% na inactivation
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

M = readtable('E:\data\polymer\pool_opto.xlsx', 'Sheet','na_inactivation');
% M = readtable('E:\data\dendritic patch\pool\pass_filter.xlsx', 'Sheet','morph');
base_path = 'E:\data';
save_path = 'E:\data\polymer\na_inactivation_opto';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end
Na_act_pre = zeros(size(M,1),11);
Na_act_post = zeros(size(M,1),11);
Na_in_pre = zeros(size(M,1),11);
Na_in_post = zeros(size(M,1),11);
K_pre = zeros(size(M,1),11);
K_post = zeros(size(M,1),11);
Vstep = zeros(size(M,1),11);
for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder(i));
    [Na_act_pre(i,:),~,Na_in_pre(i,:), ~, K_pre(i,:), Vstep(i,:)] = na_inactivation(data_path, M.date{i}, M.cell(i), M.pre(i), nan, 0, 1, save_path, 0);
    [Na_act_post(i,:),~,Na_in_post(i,:), ~, K_post(i,:), Vstep(i,:)] = na_inactivation(data_path, M.date{i}, M.cell(i), M.post(i), nan, 0,1,save_path, 1);
end

GNaact_pre = Na_act_pre./(Vstep-53.9);
GNaact_post = Na_act_post./(Vstep-53.9);
ft = fittype('G*(1/(1+exp(-(x-Vh)/k)))', 'independent', 'x');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [60,-30,3];
opts.Lower = [0,-70,0];
opts.Upper = [Inf, 0, Inf];
opts.Exclude = [11];
curveNaact_pre = fit(Vstep(1,:)',mean(GNaact_pre,1)', ft, opts);
curveNaact_post = fit(Vstep(1,:)',mean(GNaact_post,1)', ft, opts);

colors = [[0,0,0];[119,176,203]/255]; % color for Na rich polymer
% colors = [[0,0,0];[159,186,149]/255]; % color for K rich polymer
% colors = [[0,0,0];[195,129,168]/255]; % color for pedotPSS
% errorbar_with_fitcurve(Vstep(1,:), {GNaact_pre', GNaact_post'}, {curveNaact_pre,curveNaact_post},colors)

GNain_pre = Na_in_pre./(0-53.9);
GNain_post = Na_in_post./(0-53.9);
ft = fittype('G*(1/(1+exp(-(-x+Vh)/k)))', 'independent', 'x');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [60,-30,3];
opts.Lower = [0,-70,0];
opts.Upper = [Inf, 0, Inf];
curveNain_pre = fit(Vstep(1,:)',mean(GNain_pre(:,:),1)', ft, opts);
curveNain_post = fit(Vstep(1,:)',mean(GNain_post(:,:),1)', ft, opts);

errorbar_with_fitcurve(Vstep(1,:), {GNain_pre(:,:)', GNain_post(:,:)'}, {curveNain_pre,curveNain_post},colors)
% boxplot_pairwise(-[reshape(Na_in_pre(:,1:6),[],1)/1e3,reshape(Na_in_post(:,1:6)/1e3,[],1)], colors)
[data_table, within_design] = gen_table_for_ranova(Vstep(1,:), {GNain_pre, GNain_post});

rm = fitrm(data_table,'measurements1-measurements22 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);
% postHoc = multcompare(rm,'treatment');
% postHoc_g = postHoc.pValue(1);
p = zeros(1, length(Vstep));
for i = 1:length(Vstep)
    [p(i),~] = signrank(GNain_pre(:,i), GNain_post(:,i));
end

GK_pre = K_pre./(Vstep+85);
GK_post = K_post./(Vstep+85);
for i = 1:size(M,1)
    GK_pre(i,:) = GK_pre(i,:)-min(GK_pre(i,2:end));
    GK_post(i,:) = GK_post(i,:)-min(GK_post(i,2:end));
end
% ft = fittype('G*(1/(1+exp(-(x-Vh)/k)))', 'independent', 'x');
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.StartPoint = [60,-30,3];
% opts.Lower = [0,-50,0];
% opts.Upper = [Inf, 20, Inf];
% curveK_pre = fit(Vstep(1,:)',mean(GK_pre,1)', ft, opts);
% curveK_post = fit(Vstep(1,:)',mean(GK_post,1)', ft, opts);
% 
% % colors = [[0,0,0];[119,176,203]/255]; % color for Na rich polymer
% colors = [[0,0,0];[159,186,149]/255]; % color for K rich polymer
% errorbar_with_fitcurve(Vstep(1,:), {GK_pre', GK_post'}, {curveK_pre,curveK_post},colors)
save(fullfile(save_path,'na_inactivation_pool_pedotPSS.mat'),'Vstep','GNain_pre','GNain_post','GNaact_pre','GNaact_post','GK_pre','GK_post','anova_table','p','M')
%% na recovery
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))


M = readtable('E:\data\polymer\pool_opto.xlsx', 'Sheet','na_recovery');
% M = readtable('E:\data\dendritic patch\pool\pass_filter.xlsx', 'Sheet','morph');
base_path = 'E:\data';
save_path = 'E:\data\polymer\na_recovery_opto';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end
tstep = [1, 5, 10, 20, 50, 100, 200, 500, 1000];
Na_rec_pre = zeros(size(M,1),length(tstep));
Na_rec_post = zeros(size(M,1),length(tstep));


for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder(i));
    [Na_rec_pre(i,:),~,Vstep(i,:)] = na_recovery(data_path, M.date{i}, M.cell(i), M.pre(i), nan, tstep, 0, 1,save_path);
    [Na_rec_post(i,:),~,Vstep(i,:)] = na_recovery(data_path, M.date{i}, M.cell(i), M.post(i), nan, tstep, 0,1, save_path);
end
for i = 1:size(M, 1)
    Na_rec_pre_norm(i,:) = Na_rec_pre(i,:)/min(Na_rec_pre(i,:));
    Na_rec_post_norm(i,:) = Na_rec_post(i,:)/min(Na_rec_post(i,:));
end

ft = fittype('b-a*exp(-x/tau)', 'independent', 'x');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [1, 1, 10];
opts.Lower = [0,0,0];
opts.Upper = [Inf, Inf, Inf];
curveNarec_pre = fit(tstep',mean(Na_rec_pre_norm,1)', ft, opts);
curveNarec_post = fit(tstep',mean(Na_rec_post_norm,1)', ft, opts);
[data_table, within_design] = gen_table_for_ranova(tstep, {Na_rec_pre_norm, Na_rec_post_norm});

rm = fitrm(data_table,'measurements1-measurements18 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);
% postHoc = multcompare(rm,'treatment');
% postHoc_g = postHoc.pValue(1);
p = zeros(1, length(tstep));
for i = 1:length(tstep)
    [p(i),~] = signrank(Na_rec_pre_norm(:,i),Na_rec_post_norm(:,i));
end

colors = [[0,0,0];[119,176,203]/255]; % color for Na rich polymer
% colors = [[0,0,0];[159,186,149]/255]; % color for K rich polymer
% colors = [[0,0,0];[195,129,168]/255]; % color for pedotPSS
errorbar_with_fitcurve(tstep, {Na_rec_pre_norm', Na_rec_post_norm'}, {curveNarec_pre,curveNarec_post},colors), set(gca, 'XScale','log'), box on
save(fullfile(save_path,'na_recovery_pool_pedotPSS.mat'),'tstep','Na_rec_post_norm','Na_rec_pre_norm','anova_table','p','M')