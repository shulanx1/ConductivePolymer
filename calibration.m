clear all
close all

addpath(genpath(fullfile(pwd,'main')))
addpath(genpath(fullfile(pwd,'plotting')))

M = readtable('E:\data\polymer\pool_temp_calibration.xlsx', 'Sheet','calibration_acsf');
M = M(4:6,:);
base_path = 'E:\data';
save_path = 'E:\data\polymer\temp_calibration';
if ~exist(save_path, 'dir')
    mkdir(save_path)
end

temp_sample = 1:17;
R_sample_acsf = zeros(size(M,1), length(temp_sample));
baseline_acsf = zeros(size(M,1), 1);
temp_acsf_norm = zeros(size(M,1), length(temp_sample));
temp_acsf = zeros(size(M,1), length(temp_sample));

for i = 1:size(M,1)
    base_folder = fullfile(base_path, M.folder{i});
    [baseline, sample_temp, R_sample_temp] = resistance_temperature_calibration(base_folder, M.date{i}, M.idx_f(i),temp_sample, 0, 1, save_path);
    baseline_acsf(i) = baseline;
%     R_sample_acsf(i,:) = R_sample_temp/R_sample_temp(1);
%     temp_acsf_norm(i,:) = sample_temp/baseline;
    R_sample_acsf(i,:) = log(R_sample_temp(sample_temp==baseline)./R_sample_temp);
    temp_acsf_norm(i,:) =1/(baseline+273.15)- 1./(sample_temp+273.15);
    temp_acsf(i,:) = sample_temp-baseline;
end



mdl = fitlm(reshape(R_sample_acsf,[],1), reshape(temp_acsf_norm,[],1));
figure,plot(mdl);
figure,scatter(reshape(R_sample_acsf,[],1),reshape(temp_acsf_norm,[],1),'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'None')
hold on
x_range = xlim();
x_plot = x_range(1):0.01:x_range(2);
plot(x_plot, predict(mdl, x_plot'), 'r')
box on

% figure,plot(mdl);
% ft = fittype('a*x+b', 'independent', 'x');
% curve_acsf = fit(temp_acsf(1,:)', nanmean(R_sample_acsf, 1)', ft, 'Start', [0, 8]);
% errorbar_with_fitcurve(temp_acsf(1,:), {R_sample_acsf'}, {curve_acsf}); box on


% save(fullfile(save_path, 'temp_calibration.mat'),'mdl', 'M','baseline_acsf', 'R_sample_acsf', 'temp_acsf_norm', 'temp_acsf')





save(fullfile(save_path, 'temp_calibration_log.mat'),'mdl', 'M','baseline_acsf', 'R_sample_acsf', 'temp_acsf_norm', 'temp_acsf')
%% agar
% 
% 
% 
% M = readtable('Z:\data\Meisam\Heat shock\pool\pool_files.xlsx', 'Sheet','calibration_agar');
% base_path = 'Z:\data\Meisam';
% temp_sample = [23:38];
% R_sample = zeros(size(M,1), length(temp_sample));
% 
% for i = 1:size(M,1)
%     base_folder = fullfile(base_path, M.folder{i});
%     R_sample_temp = resistance_temperature_calibration(base_folder, M.date{i}, M.idx_f(i),temp_sample, 1, 0);
%     a = R_sample_temp(~isnan(R_sample_temp));
%     R_sample(i,:) = R_sample_temp/a(1);
% end
% ft = fittype('a*x+b', 'independent', 'x');
% curve = fit(temp_sample(~isnan(mean(R_sample, 1)))'-temp_baseline, mean(R_sample(~isnan(R_sample)), 1)', ft, 'Start', [0, 8]);
% errorbar_with_fitcurve(temp_sample-temp_baseline, {R_sample_acsf',R_sample'}, {curve_acsf,curve}); box on