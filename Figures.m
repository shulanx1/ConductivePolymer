close all
clear all
clc

datapath = 'C:\Users\xiao208\Box\data_backup\polymer\data_backup'; % path to the dataset downloaded from zenodo
addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))
%% Figure 3e, Extended figure 2a, Figure S22(a&b)
load(fullfile(datapath,'sub&supra','excitability_pool.mat'))
M = readtable(fullfile(datapath,'pool.xlsx'), 'Sheet','excitability');

colors = [[0,0,0];[119,176,203]/255]; % color for Na rich polymer
boxplot_pairwise(Rm(1:11,:), colors), ylabel('R_m [MOhm]')
boxplot_pairwise(Cm(1:11,:), colors), ylabel('C_m [pF]')
barplot_pairwise(AP_amp((M.withttx==0)&(sum(isnan(AP_amp),2)==0),:), colors), ylim([50,100]), ylabel('AP amp [mV]')
barplot_pairwise(AP_width((M.withttx==0)&(sum(isnan(AP_width),2)==0),:), colors), ylim([1,3.5]), ylabel('AP width [ms]')
barplot_pairwise(-Ihold, colors),ylabel('leak [pA]')
lineplot_with_shaded_errorbar(I_temp, {f_pre(:,find(M.withttx==0)),f_post(:,find(M.withttx==0))}, colors), xlim([0,350]), xlabel('I [pA]'); ylabel('f [Hz]')
[data_table, within_design] = gen_table_for_ranova(I_temp, {f_pre(:,find(M.withttx==0)),f_post(:,find(M.withttx==0))});

rm = fitrm(data_table,'measurements1-measurements24 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);
%% Extended figure 2c
load(fullfile(datapath,'na_activation','na_activation_pool_w&wo_ttx.mat'))
colors = [[0,0,0];[119,176,203]/255]; % color for Na rich polymer
boxplot_pairwise(-[G_pre(1:6,10)*(-53.9), G_post(1:6,10)*(-53.9)]/1000, colors), ylabel('INa [nA]')
load(fullfile(datapath,'na_activation','ca_activation_pool_w&wo_ttx.mat'))
boxplot_pairwise([GK_pre(:,13)*(30+85), GK_post(:,13)*(30+85)]/1000, colors), ylabel('IK [nA]')
%% Figure 3h, Extended figure 2d, Figure S22(c&d)
load(fullfile(datapath,'sub&supra_PEDOT','excitability_pool.mat'))
M = readtable(fullfile(datapath,'pool.xlsx'), 'Sheet','pedotPSS');

colors = [[0,0,0];[195,129,168]/255]; % color for pedotPSS
boxplot_pairwise(Rm(:,:), colors), ylabel('R_m [MOhm]')
boxplot_pairwise(Cm(:,:), colors), ylabel('C_m [pF]')
barplot_pairwise(AP_amp((M.withttx==0)&(sum(isnan(AP_amp),2)==0),:), colors), ylim([50,100]), ylabel('AP amp [mV]')
barplot_pairwise(AP_width((M.withttx==0)&(sum(isnan(AP_width),2)==0),:), colors), ylim([1,3.5]), ylabel('AP width [ms]')
barplot_pairwise(-Ihold, colors),ylabel('leak [pA]')
lineplot_with_shaded_errorbar(I_temp, {f_pre(:,find(M.withttx==0)),f_post(:,find(M.withttx==0))}, colors), xlim([0,350]), xlabel('I [pA]'); ylabel('f [Hz]')
[data_table, within_design] = gen_table_for_ranova(I_temp, {f_pre(:,find(M.withttx==0)),f_post(:,find(M.withttx==0))});

rm = fitrm(data_table,'measurements1-measurements24 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);
%% Extended figure 2f
load(fullfile(datapath,'na_activation_pedotPSS','na_activation_pedotPSS_pool.mat'))
colors = [[0,0,0];[195,129,168]/255]; % color for pedotPSS
boxplot_pairwise(-[G_pre(:,10)*(-53.9), G_post(:,10)*(-53.9)]/1000, colors), ylabel('INa [nA]')
load(fullfile(datapath,'ca_activation','ca_activation_pedotPSS_pool.mat'))
boxplot_pairwise([GK_pre(:,13)*(30+85), GK_post(:,13)*(30+85)]/1000, colors), ylabel('IK [nA]')
%% extended figure 3
colors = [[0,0,0];[119,176,203]/255]; % color for Na rich polymer

% Na activation
load(fullfile(datapath,'na_activation','na_activation_pool.mat'))
ft = fittype('G*(1/(1+exp(-(x-Vh)/k)))', 'independent', 'x');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [60,-30,3];
opts.Lower = [0,-70,0];
opts.Upper = [Inf, 0, Inf];
opts.Exclude = [11];
curve_pre = fit(Vstep(1,:)',mean(G_pre,1)', ft, opts);
curve_post = fit(Vstep(1,:)',mean(G_post,1)', ft, opts);
errorbar_with_fitcurve(Vstep(1,:), {G_pre', G_post'}, {curve_pre,curve_post},colors),xlim([-90,0]),box on
title('Na activation'), xlabel('Vcmd [mV]'), ylabel('GNa [nS]')

[data_table, within_design] = gen_table_for_ranova(Vstep(1,:), {G_pre', G_post'});
rm = fitrm(data_table,'measurements1-measurements22 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);

% Na inactivation
load(fullfile(datapath,'na_inactivation','na_inactivation_pool.mat'))
ft = fittype('G*(1/(1+exp(-(-x+Vh)/k)))', 'independent', 'x');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [60,-30,3];
opts.Lower = [0,-70,0];
opts.Upper = [Inf, 0, Inf];
curveNain_pre = fit(Vstep(1,:)',mean(GNain_pre(:,:),1)', ft, opts);
curveNain_post = fit(Vstep(1,:)',mean(GNain_post(:,:),1)', ft, opts);
errorbar_with_fitcurve(Vstep(1,:), {GNain_pre(:,:)', GNain_post(:,:)'}, {curveNain_pre,curveNain_post},colors),xlim([-90,0]),box on
title('Na inactivation'), xlabel('Vcmd [mV]'), ylabel('GNa [nS]')

[data_table, within_design] = gen_table_for_ranova(Vstep(1,:), {GNain_pre', GNain_post'});
rm = fitrm(data_table,'measurements1-measurements22 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);

% Na recovery
load(fullfile(datapath, 'na_recovery', 'na_recovery_pool.mat'))
ft = fittype('b-a*exp(-x/tau)', 'independent', 'x');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [1, 1, 10];
opts.Lower = [0,0,0];
opts.Upper = [Inf, Inf, Inf];
curveNarec_pre = fit(tstep',mean(Na_rec_pre_norm,1)', ft, opts);
curveNarec_post = fit(tstep',mean(Na_rec_post_norm,1)', ft, opts);
errorbar_with_fitcurve(tstep, {Na_rec_pre_norm', Na_rec_post_norm'}, {curveNarec_pre,curveNarec_post},colors), set(gca, 'XScale','log'), box on
title('Na inactivation'), xlabel('Vcmd [mV]'), ylabel('GNa [AU]')

[data_table, within_design] = gen_table_for_ranova(tstep, {Na_rec_pre_norm', Na_rec_post_norm'});
rm = fitrm(data_table,'measurements1-measurements18 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);

% K activation
load(fullfile(datapath, 'ca_activation', 'ca_activation_pool.mat'))
ft = fittype('G*(1/(1+exp(-(x-Vh)/k)))', 'independent', 'x');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [60,-30,3];
opts.Lower = [0,-50,0];
opts.Upper = [Inf, 20, Inf];
opts.Exclude = [1];
curveK_pre = fit(Vstep(1,:)',mean(GK_pre)', ft, opts);
curveK_post = fit(Vstep(1,:)',mean(GK_post)', ft, opts);
errorbar_with_fitcurve(Vstep(1,:), {GK_pre', GK_post'}, {curveK_pre,curveK_post},colors), box on
title('K activation'), xlabel('Vcmd [mV]'), ylabel('GK [nS]')

[data_table, within_design] = gen_table_for_ranova(Vstep(1,:), {GK_pre', GK_post'});
rm = fitrm(data_table,'measurements1-measurements30 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table_K = anovaTable(AT, 'DV');
disp(anova_table_K);
%% figure 4c
load(fullfile(datapath,'sub&supra_dend','excitability_pool.mat'))
M = readtable(fullfile(datapath,'pool_dstim.xlsx'), 'Sheet','excitability');

colors = [[0,0,0];[119,176,203]/255]; % color for Na rich polymer
lineplot_with_shaded_errorbar(I_temp, {f_pre(:,find(M.withttx==0)),f_post(:,find(M.withttx==0))}, colors), xlim([0,350]), xlabel('I [pA]'); ylabel('f [Hz]')
[data_table, within_design] = gen_table_for_ranova(I_temp, {f_pre(:,find(M.withttx==0)),f_post(:,find(M.withttx==0))});

rm = fitrm(data_table,'measurements1-measurements24 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);

%% figure 4f & 4h
load(fullfile(datapath, 'dstim','dstim.mat'))
colors = [[0,0,0];[119,176,203]/255];

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


a = padcat(attenuation_pool{:})';
b = padcat(attenuation_pool_control{:})';
errorbar_with_lines(power_pool([1:9,12]), {b([1:9,12],:), a([1:9,12],:)}, colors), box on
xlabel('Laser power [mW]'), ylabel('FR modulation [AU]')

%% figure S26c
load(fullfile(datapath, 'dstim','dstim.mat'))
colors = [[0,0,0];[119,176,203]/255];

idx2 = find(power_control==149);
idx2(2) = 23;

FR_pool_control_hightemp = FR_control(idx2,:);
errorbar_with_lines([1:3], {FR_pool_control_hightemp'}, colors(1,:))
hold on
for i = 1:size(FR_pool_control_hightemp,1)
    plot([1:3], FR_pool_control_hightemp(i,:),'Color',[0.5,0.5,0.5])
end
xlim([0.5,3.5])
%% figure S23c
load(fullfile(datapath, 'dstim','dstim_hyperpolarization.mat'))
colors = [[0,0,0];[119,176,203]/255];

errorbar_with_lines([1:3], { RMP_control'}, colors(1,:))
hold on
for i = 1:size( RMP_control,1)
    plot([1:3],  RMP_control(i,:),'Color',[0.5,0.5,0.5])
end
xlim([0.5,3.5])
%% figure S25c
load(fullfile(datapath, 'dstim','dstim_hyperpolarization.mat'))
colors = [[0,0,0];[119,176,203]/255];

dV = RMP(:,2)-RMP(:,1);
dV_control = RMP_control(:,2)-RMP_control(:,1);
barplot_with_datapoint({dV_control,dV}, colors), ylabel('DeltaV [mV]')

%% extended figure 4c
load(fullfile(datapath, 'temp_calibration','laser_temp_change.mat'))
colors = [[0,0,0];[119,176,203]/255];

lineplot_with_shaded_errorbar(power_mW, {padcat(dTemp_tissue{:})',padcat(dTemp_polymer{:})'}), xlim([0,25])
dT1 = padcat(dTemp_tissue{:});
dT2 = padcat(dTemp_polymer{:});
[data_table, within_design] = gen_table_for_ranova(power_mW(1:8), {dT2(:,1:8),dT1(:,1:8)});
rm = fitrm(data_table,'measurements1-measurements16 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);
%% figure S27c
load(fullfile(datapath, 'temp_calibration','temp_calibration_log.mat'))

% figure,plot(mdl);
figure,scatter(reshape(R_sample_acsf,[],1),reshape(temp_acsf_norm,[],1),'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'None')
hold on
x_range = xlim();
x_plot = x_range(1):0.01:x_range(2);
plot(x_plot, predict(mdl, x_plot'), 'r')
box on
%% figure 5i
load(fullfile(datapath, 'invivo','allMiceCombined.mat'))

[~,TFBaselineReject] = rmoutliers(ISIBaseline);
[~,TFPolymerReject] = rmoutliers(ISIPolymer);

ISIBaselineReject = ISIBaseline(~(TFPolymerReject | TFBaselineReject));
ISIPolymerReject = ISIPolymer(~(TFPolymerReject | TFBaselineReject));

temp = [ISIBaselineReject';ISIPolymerReject'];
templabels = cellstr([repmat('Baseline',size(ISIBaselineReject,2),1);repmat('Polymer ',size(ISIPolymerReject,2),1)]);
colors = [166/255 14/255 90/255;53/255 189/255 206/255];
figure,violinplot(temp,templabels,'ShowData',true,'ShowWhiskers',false,'ShowBox',false,'MarkerSize',5,'ViolinColor',colors);
ylabel('ISI (s)'); title('Polymer nPBDF');
set(gca,'TickDir','out','fontsize',14');box off;