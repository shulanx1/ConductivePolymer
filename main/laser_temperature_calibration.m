function [dT_laser, temp, t_down] = laser_temperature_calibration(base_path, date, idx_f, mdl, temp_baseline, if_plot, if_save, save_path)
if nargin < 4, load('E:\data\polymer\temp_calibration\temp_calibration.mat', 'mdl'); end
if nargin < 5, temp_baseline = 34; end
if nargin < 6, if_plot = 0; end
if nargin < 7, if_save = 0; end
if nargin < 8, save_path = base_path; end


preflix = strcat(strcat(date(5:6), date(1:4)),'_001');

filename = sprintf('%s.power_long_resistance.%d.wcp', preflix,idx_f);

try
    out=import_wcp(fullfile(base_path, filename),'debug');
catch
    date_plus_one = next_date_string(date);
    preflix = strcat(strcat(date_plus_one(5:6),date_plus_one(1:4)),'_001');
    filename = sprintf('%s.power_long_resistance.%d.wcp', preflix,idx_f);
    out=import_wcp(fullfile(base_path, filename),'debug');
end

dt = out.T(2)-out.T(1); % sample time in s
n_recording = length(out.rec_index);
   
Im = reshape(out.S{3}, [], 1);
Vm = reshape(out.S{4}, [], 1);
power = reshape(out.S{7}, [], 1);
t = [];
for n = 1:n_recording
    t = [t, out.T+out.T(end)*(n-1)];
end

[R_down,t_down, idx_down] = calc_R(Vm, Im,t,  dt);
R_down = R_down/10;
dt_down = (idx_down(2)-idx_down(1))*dt;
R_down_orig = R_down;
% R_down = R_down/mean(R_down(end-10:end));
% temp = predict(mdl,R_down')*temp_baseline;
R_baseline = mean(R_down(1:10));



% idx_rmv = find(abs(diff(temp))>= 50);
% for i = 1:length(idx_rmv)
%     temp(idx_rmv(i)+1) = mean(temp([idx_rmv(i)-4: idx_rmv(i)+6]));
% end


power_down = power(idx_down);

idx_laser_start = find(diff(power)>=10);
idx_laser_start(find(diff(idx_laser_start)<=1)+1) = [];
idx_laser_end = find(diff(power)<-10);
idx_laser_end(find(diff(idx_laser_end)<=1)+1) = [];

if isempty(idx_laser_start)
    flag = 1;
    idx_laser_start = [42, 82, 122, 162, 202]*1e3;
    idx_laser_end = [62, 102, 142, 182, 222]*1e3 + 1;
else
    flag = 0;
end
t_laser_start = t(idx_laser_start);
t_laser_end = t(idx_laser_end);

temp = real(1./(1/(temp_baseline+273.15)-predict(mdl,log(R_baseline./R_down'))))-273.15;
Fs = 1/dt_down;
Fc = 7;
Wn = Fc./(Fs/2);
[b1,a1] = butter(6, Wn, 'low');
temp = filtfilt(b1, a1, temp);

n_smooth = 1;
for n = 1:n_smooth
idx_rmv = find(abs(diff(temp))>=15)+1;
if ~isempty(idx_rmv)
    for i = 1:length(idx_laser_start)
        idx_rmv(find(abs(idx_rmv-idx_laser_start(i))<=2)) = [];
        idx_rmv(find(abs(idx_rmv-idx_laser_end(i))<=2)) = [];
    end
    for i = 1:length(idx_rmv)
        idx_smooth = max(1,idx_rmv(i)-5):min(idx_rmv(i)+5, length(temp));
        idx_smooth = setdiff(idx_smooth,idx_rmv);
        temp(idx_rmv(i)) = mean(temp(idx_smooth));
    end
end
end
% [yupper,ylower] = envelope(temp-temp_baseline, 1,'peak');
% temp = ylower+ temp_baseline;

a = [42, 82, 122, 162, 202]*1e3;
b = [62, 102, 142, 182, 222]*1e3 + 1;
idx_mid = floor((a+b)/2);

if flag==0
    idx_on = find(power(idx_mid)>10);
else
    idx_on = 1:5;
end

temp_laseroff = mean(temp(end-10:end));
temp_laseron = zeros(1, length(idx_laser_start));
for i = 1:length(idx_laser_start)
    idx_current = find(t_down>=t_laser_start(i)&t_down<=t_laser_end(i));
    if isempty(idx_current)
        idx_current = find(abs(t_down-t_laser_start(i))<=dt_down);
    end
    temp_laseron(i) = median(temp(idx_current));
    idx_laser_start_down(i) = idx_current(1);
    idx_laser_end_down(i) = idx_current(end);
end

dT_laser = temp_laseron-temp_laseroff;
temp_truncate = reshape(temp(65:564),100,[]);
temp_truncate = temp_truncate(:,idx_on);

if if_plot
    figure,subplot(211),
%     yyaxis left, 
    plot(t_down, temp),ylabel('T [degree C]');
%     yyaxis right, plot(t_down, R_down_orig),ylim([0,20]), ylabel('R [MOhm]'),
    x1 = xlim();
    subplot(212), plot(t, power), xlim(x1);
    
    lineplot_with_shaded_errorbar(t_down(1:100), temp_truncate)
end
if if_save
    save(fullfile(save_path, sprintf('laser_heating_%d.mat', idx_f)), 't_down', 'idx_down', 'R_down', 'temp', 'power_down', 'idx_laser_start', 'idx_laser_end','idx_laser_start_down','idx_laser_end_down','temp_truncate', 'dT_laser')
end
fprintf('finished file %s\n',filename)
end