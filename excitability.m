%% excitability
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

M = readtable('E:\data\polymer\pool.xlsx', 'Sheet','excitability_Krich');
% M = readtable('E:\data\dendritic patch\pool\pass_filter.xlsx', 'Sheet','morph');
base_path = 'E:\data';
save_path = 'E:\data\polymer\sub&supra_Krich';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end
Rm = [];
RMP = [];
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
    [params_pre_temp, ~, f_pre_temp, I_temp, AP_amp_pre_temp, AP_width_pre_temp] = sub_and_supra_classify(data_path, M.date{i}, M.cell(i), M.pre(i),0,1,save_path);
    [params_post_temp, ~, f_post_temp, I_temp, AP_amp_post_temp, AP_width_post_temp] = sub_and_supra_classify(data_path, M.date{i}, M.cell(i), M.post(i),0,1,save_path);
    Rm = [Rm;[params_pre_temp(1), params_post_temp(1)]];
    Cm = [Cm;[params_pre_temp(5), params_post_temp(5)]];
    RMP = [RMP;[params_pre_temp(8), params_post_temp(8)]];
    AP_amp = [AP_amp; [median(AP_amp_pre_temp), median(AP_amp_post_temp)]];
    AP_width = [AP_width; [median(AP_width_pre_temp), median(AP_width_post_temp)]];
    Rm_fit = [Rm_fit;[params_pre_temp(6), params_post_temp(6)]];
    Ihold = [Ihold;[params_pre_temp(7), params_post_temp(7)]];
    f_pre = [f_pre,f_pre_temp'];
    f_post = [f_post,f_post_temp'];
    I = [I, I_temp'];
end

% idx_rmv = find(abs(diff(Ihold'))>200);
% Rm(idx_rmv,:) = [];
% Cm(idx_rmv,:) = [];
% Rm_fit(idx_rmv,:) = [];
% Ihold(idx_rmv,:) = [];
% f_pre(:,idx_rmv) = [];
% f_post(:,idx_rmv) = [];
% AP_amp(idx_rmv,:) = [];

% colors = [[0,0,0];[119,176,203]/255]; % color for Na rich polymer
colors = [[0,0,0];[159,186,149]/255]; % color for K rich polymer
boxplot_pairwise(Rm, colors)
% boxplot_pairwise(Rm_fit.*Cm/1e3)
boxplot_pairwise(Cm, colors)
barplot_pairwise(AP_amp((M.withttx==0)&(sum(isnan(AP_amp),2)==0),:), colors), ylim([50,100])
barplot_pairwise(AP_width((M.withttx==0)&(sum(isnan(AP_width),2)==0),:), colors), ylim([1,3.5])
barplot_pairwise(-Ihold, colors)
lineplot_with_shaded_errorbar(I_temp, {f_pre(:,find(M.withttx==0)),f_post(:,find(M.withttx==0))}, colors), xlim([0,350])
[data_table, within_design] = gen_table_for_ranova(I_temp((I_temp>0)&(I_temp<=350)), {f_pre((I_temp>0)&(I_temp<=350),find(M.withttx==0)),f_post((I_temp>0)&(I_temp<=350),find(M.withttx==0))});

rm = fitrm(data_table,'measurements1-measurements14 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);
% postHoc = multcompare(rm,'treatment');
% postHoc_g = postHoc.pValue(1);
p = zeros(1, length(I_temp));
for i = 1:length(I_temp)
    [p(i),~] = signrank(f_pre(i,find(M.withttx==0))',f_post(i,find(M.withttx==0))');
end
save(fullfile(save_path,'excitability_pool_Krich.mat'),'RMP','Rm','Cm','AP_amp','AP_width', 'Ihold','I_temp','f_pre', 'f_post', 'anova_table', 'p', 'M')
%% excitability_pdot_PSS
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

M = readtable('E:\data\polymer\pool.xlsx', 'Sheet','pedotPSS');
base_path = 'E:\data';
save_path = 'E:\data\polymer\sub&supra_PEDOT';
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
    [params_pre_temp, ~, f_pre_temp, I_temp, AP_amp_pre_temp,AP_width_pre_temp] = sub_and_supra_classify(data_path, M.date{i}, M.cell(i), M.pre(i), 0,1,save_path);
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


colors = [[0,0,0];[195,129,168]/255];
boxplot_pairwise(Rm, colors)
% boxplot_pairwise(Rm_fit.*Cm/1e3)
boxplot_pairwise(Cm, colors)
barplot_pairwise(AP_amp((M.withttx==0)&(sum(isnan(AP_amp),2)==0),:), colors), ylim([50,100])
barplot_pairwise(AP_width((M.withttx==0)&(sum(isnan(AP_amp),2)==0),:), colors), ylim([1,3.5])
barplot_pairwise(-Ihold, colors)
lineplot_with_shaded_errorbar(I_temp, {f_pre(:,find(M.withttx==0)),f_post(:,find(M.withttx==0))}, colors), xlim([0,350])
[data_table, within_design] = gen_table_for_ranova(I_temp((I_temp>0)&(I_temp<=350)), {f_pre((I_temp>0)&(I_temp<=350),find(M.withttx==0)),f_post((I_temp>0)&(I_temp<=350),find(M.withttx==0))});

rm = fitrm(data_table,'measurements1-measurements14 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);
% postHoc = multcompare(rm,'treatment');
% postHoc_g = postHoc.pValue(1);
p = zeros(1, length(I_temp));
for i = 1:length(I_temp)
    [p(i),~] = signrank(f_pre(i,find(M.withttx==0))',f_post(i,find(M.withttx==0))');
end
save(fullfile(save_path,'excitability_Pedot_pool.mat'),'Rm','Cm','AP_amp', 'AP_width','Ihold','I_temp','f_pre', 'f_post', 'anova_table', 'p', 'M')
%% na activation
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))
M = readtable('E:\data\polymer\pool.xlsx', 'Sheet','na_activation');
% M = readtable('E:\data\dendritic patch\pool\pass_filter.xlsx', 'Sheet','morph');
base_path = 'E:\data';
save_path = 'E:\data\polymer\na_activation';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end

Na_amp_pre = zeros(size(M,1),11);
Na_amp_post = zeros(size(M,1),11);
Vstep = zeros(size(M,1),11);
for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder(i));
    [Na_amp_pre(i,:), ~, ~] = na_activation(data_path, M.date{i}, M.cell(i), M.pre(i), M.ttx(i), 0, 0, save_path);
    [Na_amp_post(i,:), ~, Vstep(i,:)] = na_activation(data_path, M.date{i}, M.cell(i), M.post(i), M.ttx(i), 0,0,save_path);
%     [Na_amp_pre(i,:), ~, ~] = na_activation(data_path, M.date{i}, M.cell(i), M.pre(i), nan, 0, 1, save_path);
%     [Na_amp_post(i,:), ~, Vstep(i,:)] = na_activation(data_path, M.date{i}, M.cell(i), M.post(i), nan, 0,1,save_path);
end
%%
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

% colors = [[0,0,0];[119,176,203]/255]; % color for Na rich polymer
% colors = [[0,0,0];[159,186,149]/255]; % color for K rich polymer
colors = [[0,0,0];[195,129,168]/255]; % color for pedotPSS
errorbar_with_fitcurve(Vstep(1,:), {G_pre', G_post'}, {curve_pre,curve_post},colors),xlim([-90,0]),box on
boxplot_pairwise(-[reshape(Na_amp_pre(:,6:11),[],1)/1e3,reshape(Na_amp_post(:,6:11)/1e3,[],1)], colors)
[data_table, within_design] = gen_table_for_ranova(Vstep(1,:), {G_pre', G_post'});

rm = fitrm(data_table,'measurements1-measurements22 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);
% postHoc = multcompare(rm,'treatment');
% postHoc_g = postHoc.pValue(1);
p = zeros(1, length(Vstep));
for i = 1:length(Vstep)
    [p(i),~] = signrank(G_pre(:,i), G_post(:,i));
end

save(fullfile(save_path,'na_activation_pool_pedotPSS.mat'),'G_pre','G_post','Vstep','anova_table','p','M')
%% na inactivation
% clear all
% close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

M = readtable('E:\data\polymer\pool.xlsx', 'Sheet','na_inactivation');
% M = readtable('E:\data\dendritic patch\pool\pass_filter.xlsx', 'Sheet','morph');
base_path = 'E:\data';
save_path = 'E:\data\polymer\na_inactivation';
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
    [Na_act_pre(i,:),~,Na_in_pre(i,:), ~, K_pre(i,:), Vstep(i,:)] = na_inactivation(data_path, M.date{i}, M.cell(i), M.pre(i), M.ttx(i), 0, 1, save_path);
    [Na_act_post(i,:),~,Na_in_post(i,:), ~, K_post(i,:), Vstep(i,:)] = na_inactivation(data_path, M.date{i}, M.cell(i), M.post(i), M.ttx(i), 0,1,save_path);
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

errorbar_with_fitcurve(Vstep(1,:), {GNain_pre(:,:)', GNain_post(:,:)'}, {curveNain_pre,curveNain_post},colors),box on
% boxplot_pairwise(-[reshape(Na_in_pre(:,1:6),[],1)/1e3,reshape(Na_in_post(:,1:6)/1e3,[],1)], colors)
[data_table, within_design] = gen_table_for_ranova(Vstep(1,:), {GNain_pre', GNain_post'});

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
save(fullfile(save_path,'na_inactivation_pool.mat'),'Vstep','GNain_pre','GNain_post','GNaact_pre','GNaact_post','GK_pre','GK_post','anova_table','p','M')
%% na recovery
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))


M = readtable('E:\data\polymer\pool.xlsx', 'Sheet','na_recovery');
% M = readtable('E:\data\dendritic patch\pool\pass_filter.xlsx', 'Sheet','morph');
base_path = 'E:\data';
save_path = 'E:\data\polymer\na_recovery';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end
tstep = [1, 5, 10, 20, 50, 100, 200, 500, 1000];
Na_rec_pre = zeros(size(M,1),length(tstep));
Na_rec_post = zeros(size(M,1),length(tstep));


for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder(i));
    [Na_rec_pre(i,:),~,Vstep(i,:)] = na_recovery(data_path, M.date{i}, M.cell(i), M.pre(i), M.ttx(i), tstep, 0, 1,save_path);
    [Na_rec_post(i,:),~,Vstep(i,:)] = na_recovery(data_path, M.date{i}, M.cell(i), M.post(i), M.ttx(i), tstep, 0,1, save_path);
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
[data_table, within_design] = gen_table_for_ranova(tstep, {Na_rec_pre_norm', Na_rec_post_norm'});

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
save(fullfile(save_path,'na_recovery_pool.mat'),'tstep','Na_rec_post_norm','Na_rec_pre_norm','anova_table','p','M')
%% nap
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))



M = readtable('E:\data\polymer\pool.xlsx', 'Sheet','inap_Krich');
% M = readtable('E:\data\dendritic patch\pool\pass_filter.xlsx', 'Sheet','morph');
base_path = 'E:\data';
save_path = 'E:\data\polymer\inap_Krich';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end
Vstep = [-90:10:0];
Nap_pre = zeros(size(M,1),length(Vstep));
Nap_post = zeros(size(M,1),length(Vstep));


for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder(i));
    Nap_pre(i,:) = nap(data_path, M.date{i}, M.cell(i), M.pre(i), M.ttx(i), Vstep, 0,1,save_path);
    Nap_post(i,:) = nap(data_path, M.date{i}, M.cell(i), M.post(i), M.ttx(i),Vstep, 0,1,save_path);
end

GNap_pre = Nap_pre./(Vstep-53.9);
GNap_post = Nap_post./(Vstep-53.9);
ft = fittype('G*(1/(1+exp(-(x-Vh)/k)))', 'independent', 'x');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [60,-30,3];
opts.Lower = [0,-70,0];
opts.Upper = [Inf, 0, Inf];
curveNap_pre = fit(Vstep(1,:)',mean(GNap_pre,1)', ft, opts);
curveNap_post = fit(Vstep(1,:)',mean(GNap_post,1)', ft, opts);
[data_table, within_design] = gen_table_for_ranova(Vstep, {GNap_pre', GNap_post'});

rm = fitrm(data_table,'measurements1-measurements20 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);
% postHoc = multcompare(rm,'treatment');
% postHoc_g = postHoc.pValue(1);
p = zeros(1, length(Vstep));
for i = 1:length(Vstep)
    [p(i),~] = signrank(GNap_pre(:,i),GNap_post(:,i));
end

% colors = [[0,0,0];[119,176,203]/255]; % color for Na rich polymer
colors = [[0,0,0];[159,186,149]/255]; % color for K rich polymer
errorbar_with_fitcurve(Vstep(1,:), {GNap_pre', GNap_post'}, {curveNap_pre,curveNap_post},colors)
save(fullfile(save_path,'inap_Krich_pool.mat'),'Vstep','GNap_pre','GNap_post','anova_table','p','M')
%% ca activation
clear all
% close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

M = readtable('E:\data\polymer\pool.xlsx', 'Sheet','ca_activation_ttx_pedotPSS');
% M = readtable('E:\data\dendritic patch\pool\pass_filter.xlsx', 'Sheet','morph');
base_path = 'E:\data';
save_path = 'E:\data\polymer\ca_activation_pedotPSS';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end

K_pre = zeros(size(M,1),15);
K_post = zeros(size(M,1),15);
Ca_pre = zeros(size(M,1),15);
Ca_post = zeros(size(M,1),15);
Vstep = zeros(size(M,1),15);
for i = 1:size(M, 1)
%     if M.with_ttx(i) == 1
        data_path = fullfile(base_path,M.folder(i));
        [K_pre(i,:), Ca_pre(i,:), Vstep(i,:)] = ca_activation_ttx(data_path, M.date{i}, M.cell(i), M.pre(i), 0, 1,save_path);
        [K_post(i,:), Ca_post(i,:), Vstep(i,:)] = ca_activation_ttx(data_path, M.date{i}, M.cell(i), M.post(i), 0,1,save_path);
%     end
end

% idx_rmv = find([M.with_ttx]==0);
% K_pre(idx_rmv,:) = [];
% K_post(idx_rmv,:) = [];
% Ca_pre(idx_rmv,:) = [];
% Ca_post(idx_rmv,:) = [];
% Vstep(idx_rmv,:) = [];

GCa_pre = Ca_pre./(Vstep-120);
GCa_post = Ca_post./(Vstep-120);
ft = fittype('G*(1/(1+exp(-(x-Vh)/k)))', 'independent', 'x');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [10,0,3];
opts.Lower = [0,-70,0];
opts.Upper = [Inf, 50, Inf];
curveCa_pre = fit(Vstep(1,:)',mean(GCa_pre)', ft, opts);
curveCa_post = fit(Vstep(1,:)',mean(GCa_post)', ft, opts);

colors = [[0,0,0];[119,176,203]/255]; % color for Na rich polymer
% colors = [[0,0,0];[159,186,149]/255]; % color for K rich polymer
% colors = [[0,0,0];[195,129,168]/255]; % color for pedot PSS
errorbar_with_fitcurve(Vstep(1,:), {GCa_pre', GCa_post'}, {curveCa_pre,curveCa_post},colors), box on
[data_table, within_design] = gen_table_for_ranova(Vstep(1,:), {GCa_pre', GCa_post'});
rm = fitrm(data_table,'measurements1-measurements30 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table_Ca = anovaTable(AT, 'DV');
disp(anova_table_Ca);
% postHoc = multcompare(rm,'treatment');
% postHoc_g = postHoc.pValue(1);
pCa = zeros(1, length(Vstep(1,:)));
for i = 1:length(Vstep(1,:))
    [pCa(i),~] = signrank(GCa_pre(:,i),GCa_post(:,i));
end



GK_pre = K_pre./(Vstep+85);
GK_post = K_post./(Vstep+85);
ft = fittype('G*(1/(1+exp(-(x-Vh)/k)))', 'independent', 'x');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [60,-30,3];
opts.Lower = [0,-50,0];
opts.Upper = [Inf, 20, Inf];
opts.Exclude = [1];
curveK_pre = fit(Vstep(1,:)',mean(GK_pre)', ft, opts);
curveK_post = fit(Vstep(1,:)',mean(GK_post)', ft, opts);

% colors = [[0,0,0];[119,176,203]/255]; % color for Na rich polymer
% colors = [[0,0,0];[159,186,149]/255]; % color for K rich polymer
colors = [[0,0,0];[195,129,168]/255]; % color for pedot PSS
errorbar_with_fitcurve(Vstep(1,:), {GK_pre', GK_post'}, {curveK_pre,curveK_post},colors), box on
[data_table, within_design] = gen_table_for_ranova(Vstep(1,:), {GK_pre', GK_post'});
rm = fitrm(data_table,'measurements1-measurements30 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table_K = anovaTable(AT, 'DV');
disp(anova_table_K);
% postHoc = multcompare(rm,'treatment');
% postHoc_g = postHoc.pValue(1);
pK = zeros(1, length(Vstep(1,:)));
for i = 1:length(Vstep(1,:))
    [pK(i),~] = signrank(GK_pre(:,i),GK_post(:,i));
end
save(fullfile(save_path,'ca_activation_pool_pedotPSS.mat'),'Vstep','GCa_pre','GCa_post','anova_table_Ca','pCa','GK_pre','GK_post','anova_table_K','pK','M')

%% ca recovery
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

M = readtable('Z:\data\shulan\polymer\pool.xlsx', 'Sheet','ca_recovery_Krich');
% M = readtable('E:\data\dendritic patch\pool\pass_filter.xlsx', 'Sheet','morph');
base_path = 'Z:\data\shulan';
save_path = 'Z:\data\shulan\polymer\ca_recovery_Krich';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end

tstep = [1, 5, 10, 20, 50, 100, 200, 500, 1000];
Ca_rec_pre = zeros(size(M,1),length(tstep));
Ca_rec_post = zeros(size(M,1),length(tstep));


for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder(i));
    [Ca_rec_pre(i,:),~] = na_recovery_ttx(data_path, M.date{i}, M.cell(i), M.pre(i), tstep, 0, 1,save_path);
    [Ca_rec_post(i,:),~] = na_recovery_ttx(data_path, M.date{i}, M.cell(i), M.post(i), tstep, 0, 1, save_path);
end
for i = 1:size(M, 1)
    Ca_rec_pre_norm(i,:) = Ca_rec_pre(i,:)/min(Ca_rec_pre(i,:));
    Ca_rec_post_norm(i,:) = Ca_rec_post(i,:)/min(Ca_rec_post(i,:));
end

ft = fittype('b-a*exp(-x/tau)', 'independent', 'x');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [1, 1, 10];
opts.Lower = [0,0,0];
opts.Upper = [Inf, Inf, Inf];
curveCarec_pre = fit(tstep',mean(Ca_rec_pre_norm)', ft, opts);
curveCarec_post = fit(tstep',mean(Ca_rec_post_norm)', ft, opts);

% colors = [[0,0,0];[119,176,203]/255]; % color for Na rich polymer
colors = [[0,0,0];[159,186,149]/255]; % color for K rich polymer
errorbar_with_fitcurve(tstep, {Ca_rec_pre_norm', Ca_rec_post_norm'}, {curveCarec_pre,curveCarec_post},colors), set(gca, 'XScale','log'), box on
[data_table, within_design] = gen_table_for_ranova(tstep, {Ca_rec_pre_norm', Ca_rec_post_norm'});
rm = fitrm(data_table,'measurements1-measurements18 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);
% postHoc = multcompare(rm,'treatment');
% postHoc_g = postHoc.pValue(1);
p = zeros(1, length(tstep));
for i = 1:length(tstep)
    [p(i),~] = signrank(Ca_rec_pre_norm(:,i),Ca_rec_post_norm(:,i));
end

save(fullfile(save_path,'ca_recovery_Krich.mat'),'tstep','Ca_rec_post_norm','Ca_rec_pre_norm','anova_table','p')

%% capacitance_VC & capacitance_IC
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))


M = readtable('E:\data\polymer\pool.xlsx', 'Sheet','capacitance_VC');
base_path = 'E:\data';
Ra = zeros(size(M,1),2);
Cm = zeros(size(M,1),2);
Ileak = zeros(size(M,1),2);
for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder(i));
    [Ra(i,1), Cm(i,1), Ileak(i,1)] = capacitance_VC(data_path, M.date{i}, M.cell(i), M.pre(i));
    [Ra(i,2), Cm(i,2), Ileak(i,2)] = capacitance_VC(data_path, M.date{i}, M.cell(i), M.post(i));
end

% idx_rmv = find(abs(diff(Ileak'))>200);
% Ra(:,idx_rmv) = [];
% Cm(:,idx_rmv) = [];


colors = [[0,0,0];[119,176,203]/255];
boxplot_pairwise(Ra, colors)
boxplot_pairwise(Cm, colors)
boxplot_pairwise(Ileak,colors)


%% capacitance_IC
clear all
close all
M = readtable('E:\data\polymer\pool.xlsx', 'Sheet','capacitance_CC');
base_path = 'E:\data';
Rm = zeros(size(M,1),2);
Cm = zeros(size(M,1),2);
Ihold = zeros(size(M,1),2);
RMP = zeros(size(M,1),2);
for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder(i));
    [Rm(i,1), Cm(i,1), Ihold(i,1), RMP(i,1),~] = capacitance_CC(data_path, M.date{i}, M.cell(i), M.pre(i));
    [Rm(i,2), Cm(i,2), Ihold(i,2), RMP(i,1),~] = capacitance_CC(data_path, M.date{i}, M.cell(i), M.post(i));
end


colors = [[0,0,0];[119,176,203]/255];
boxplot_pairwise(Rm, colors)
boxplot_pairwise(Cm, colors)
boxplot_pairwise(Ihold,colors)

%% capacitance_VC & capacitance_IC, Pdot PSS
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))


M = readtable('Z:\data\shulan\polymer\pool.xlsx', 'Sheet','capacitance_VC_pdotPSS');
base_path = 'Z:\data\shulan';
Ra = zeros(size(M,1),2);
Cm = zeros(size(M,1),2);
Ileak = zeros(size(M,1),2);
for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder(i));
    [Ra(i,1), Cm(i,1), Ileak(i,1)] = capacitance_VC(data_path, M.date{i}, M.cell(i), M.pre(i));
    [Ra(i,2), Cm(i,2), Ileak(i,2)] = capacitance_VC(data_path, M.date{i}, M.cell(i), M.post(i));
end

% idx_rmv = find(abs(diff(Ileak'))>200);
% Ra(:,idx_rmv) = [];
% Cm(:,idx_rmv) = [];


colors = [[0,0,0];[195,129,168]/255];
boxplot_pairwise(Ra, colors)
boxplot_pairwise(Cm, colors)
boxplot_pairwise(Ileak,colors)


%% capacitance_IC, Pdot PSS
clear all
close all
M = readtable('Z:\data\shulan\polymer\pool.xlsx', 'Sheet','capacitance_CC_pdotPSS');
base_path = 'Z:\data\shulan';
Rm = zeros(size(M,1),2);
Cm = zeros(size(M,1),2);
Ihold = zeros(size(M,1),2);
RMP = zeros(size(M,1),2);
for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder(i));
    [Rm(i,1), Cm(i,1), Ihold(i,1), RMP(i,1),~] = capacitance_CC(data_path, M.date{i}, M.cell(i), M.pre(i));
    [Rm(i,2), Cm(i,2), Ihold(i,2), RMP(i,1),~] = capacitance_CC(data_path, M.date{i}, M.cell(i), M.post(i));
end


colors = [[0,0,0];[195,129,168]/255];
boxplot_pairwise(Rm, colors)
boxplot_pairwise(Cm, colors)
boxplot_pairwise(Ihold,colors)