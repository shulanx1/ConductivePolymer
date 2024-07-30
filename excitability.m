%% excitability
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

M = readtable('E:\data\polymer\pool.xlsx', 'Sheet','excitability');
base_path = 'E:\data';
Rm = [];
Ihold = [];
Rm_fit = [];
Cm = [];
f_pre = [];
f_post = [];
AP_amp = [];
I = [];
for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder(i));
    [params_pre_temp, ~, f_pre_temp, I_temp, AP_amp_pre_temp] = sub_and_supra_classify(data_path, M.date{i}, M.cell(i), M.pre(i));
    [params_post_temp, ~, f_post_temp, I_temp, AP_amp_post_temp] = sub_and_supra_classify(data_path, M.date{i}, M.cell(i), M.post(i));
    Rm = [Rm;[params_pre_temp(1), params_post_temp(1)]];
    Cm = [Cm;[params_pre_temp(5), params_post_temp(5)]];
    AP_amp = [AP_amp; [median(AP_amp_pre_temp), median(AP_amp_post_temp)]];
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
lineplot_with_shaded_errorbar(I_temp, {f_pre(find(M.withttx==0),:),f_post(find(M.withttx==0),:)}, colors), xlim([0,350])

%% excitability_pdot_PSS
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

M = readtable('E:\data\polymer\pool.xlsx', 'Sheet','pdotPSS');
base_path = 'E:\data';
Rm = [];
Ihold = [];
Rm_fit = [];
Cm = [];
f_pre = [];
f_post = [];
AP_amp = [];
I = [];
for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder(i));
    [params_pre_temp, ~, f_pre_temp, I_temp, AP_amp_pre_temp] = sub_and_supra_classify(data_path, M.date{i}, M.cell(i), M.pre(i));
    [params_post_temp, ~, f_post_temp, I_temp, AP_amp_post_temp] = sub_and_supra_classify(data_path, M.date{i}, M.cell(i), M.post(i));
    Rm = [Rm;[params_pre_temp(1), params_post_temp(1)]];
    Cm = [Cm;[params_pre_temp(5), params_post_temp(5)]];
    AP_amp = [AP_amp; [median(AP_amp_pre_temp), median(AP_amp_post_temp)]];
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
barplot_pairwise(AP_amp, colors), ylim([80,140])
barplot_pairwise(-Ihold, colors)
lineplot_with_shaded_errorbar(I_temp, {f_pre(:,find(M.withttx==0)),f_post(:,find(M.withttx==0))}, colors), xlim([0,350])
%% na activation
clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

M = readtable('E:\data\polymer\pool.xlsx', 'Sheet','na_activation');
base_path = 'E:\data';
Na_amp_pre = zeros(size(M,1),11);
Na_amp_post = zeros(size(M,1),11);
Vstep = zeros(size(M,1),11);
for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder(i));
    [Na_amp_pre(i,:), ~, ~] = na_activation(data_path, M.date{i}, M.cell(i), M.pre(i), M.ttx(i), 1, 0, 'gray');
    [Na_amp_post(i,:), ~, Vstep(i,:)] = na_activation(data_path, M.date{i}, M.cell(i), M.post(i), M.ttx(i), 1);
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
curve_pre = fit(Vstep(1,:)',mean(G_pre)', ft, opts);
curve_post = fit(Vstep(1,:)',mean(G_post)', ft, opts);

colors = [[0,0,0];[119,176,203]/255];
errorbar_with_fitcurve(Vstep(1,:), {G_pre', G_post'}, {curve_pre,curve_post},colors),xlim([-90,0]),box on
boxplot_pairwise(-[reshape(Na_amp_pre(:,6:11),[],1)/1e3,reshape(Na_amp_post(:,6:11)/1e3,[],1)], colors)


%% ca activation
clear all
% close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

M = readtable('E:\data\polymer\pool.xlsx', 'Sheet','ca_activation_ttx');
base_path = 'E:\data';
K_pre = zeros(size(M,1),15);
K_post = zeros(size(M,1),15);
Vstep = zeros(size(M,1),15);
for i = 1:size(M, 1)
    data_path = fullfile(base_path,M.folder(i));
    [K_pre(i,:), ~, Vstep(i,:)] = ca_activation_ttx(data_path, M.date{i}, M.cell(i), M.pre(i), 0, 0, 'gray');
    [K_post(i,:), ~, Vstep(i,:)] = ca_activation_ttx(data_path, M.date{i}, M.cell(i), M.post(i), 0);
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

colors = [[0,0,0];[119,176,203]/255];
errorbar_with_fitcurve(Vstep(1,:), {GK_pre', GK_post'}, {curveK_pre,curveK_post},colors), box on

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