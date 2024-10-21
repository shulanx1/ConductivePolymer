function [baseline, sample_ranged, R_sample_ranged] = resistance_temperature_calibration(base_path, date, idx_f, temp_sample_shift, if_plot, if_save, save_path)
if nargin < 4, temp_sample_shift = [1:15];end
if nargin < 5, if_plot = 0; end
if nargin < 6, if_save = 0; end

preflix = strcat(strcat(date(5:6), date(1:4)),'_001');
path = fullfile(base_path, date);
if nargin < 7, save_path = path; end

filename = sprintf('%s.resistance.%d.wcp', preflix,idx_f);

try
    out=import_wcp(fullfile(path, filename),'debug');
catch
    date_plus_one = next_date_string(date);
    preflix = strcat(strcat(date_plus_one(5:6),date_plus_one(1:4)),'_001');
    filename = sprintf('%s.resistance.%d.wcp', preflix,idx_f);
    out=import_wcp(fullfile(path, filename),'debug');
end

n_channel = out.channel_no;
n_recording = length(out.rec_index);

dt = out.T(2)-out.T(1); % sample time in s

   
Im = reshape(out.S{3}, [], 1);
Vm = reshape(out.S{4}, [], 1);
I_temp = reshape(out.S{5}, [], 1);
t = [];
for n = 1:n_recording
    t = [t, out.T+out.T(end)*(n-1)];
end


idx_rmv = find(diff(Vm)==0); % remove the data points after the recording was terminated
idx_rmv(find(diff(idx_rmv)==1)+1) = [];
Im(idx_rmv(end):end) = [];
Vm(idx_rmv(end):end) = [];
I_temp(idx_rmv(end):end) = [];
t(idx_rmv(end):end) = [];
[R_down,t_down, idx_down] = calc_R(Vm, Im,t,  dt);
dt_down = (idx_down(2)-idx_down(1))*dt;

Fs = 1/dt;
Fc = 60;
Wn = Fc./(Fs/2);
QF = Wn/35;
[d,c] = iirnotch(Wn,QF);

Fc = 120;
Wn = Fc./(Fs/2);
QF = Wn/15;
[d1,c1] = iirnotch(Wn,QF);


I_temp1 = filtfilt(d1, c1, I_temp);
I_temp2 = filtfilt(d, c, I_temp1);
I_temp_down = I_temp2(idx_down);
dt_down = t_down(2)-t_down(1);

temp = calc_temp(I_temp_down, dt_down);
temp_sample = [min(ceil(temp)):1:max(floor(temp))];

R_sample = nan(1, length(temp_sample));
[~,idx_start_sample] = max(temp);
for i = 1:length(temp_sample)
    idx_current_sample = find(abs(temp(idx_start_sample:end)-temp_sample(i))<=0.1);
    if ~isempty(idx_current_sample)
        R_sample(i) = mean(R_down(idx_current_sample));
    end
end
if if_plot
    figure,yyaxis left, plot(t_down, temp), yyaxis right, plot(t_down, R_down);
    figure, scatter(temp_sample, R_sample)
end

if length(temp_sample_shift)>length(temp_sample)
    sprintf('only %d temperature datapoints from %d to %d, choose a smaller range', length(temp_sample), temp_sample(1), temp_sample(2));
    return
else
    R2 = zeros(length(temp_sample)+1-length(temp_sample_shift), 1);
    for i = 1:length(temp_sample)+1-length(temp_sample_shift)
        X = temp_sample(i-1+temp_sample_shift)-temp_sample(i);
        y = R_sample(i-1+temp_sample_shift)/R_sample(i);
        mdl = fitlm(X,y);
        R2(i) = mdl.Rsquared.Adjusted;
    end
    [~,idx] = max(R2);
    baseline = temp_sample(i);
    sample_ranged = temp_sample(i-1+temp_sample_shift);
    R_sample_ranged = R_sample(i-1+temp_sample_shift);
end

if if_save
    save(fullfile(save_path, sprintf('R_temp_calibration_%d.mat', idx_f)), 't_down', 'R_down', 'temp', 'temp_sample', 'R_sample', 'sample_ranged', 'R_sample_ranged', 'baseline')
end
fprintf('finished date %s, file %s\n',date, filename)
end