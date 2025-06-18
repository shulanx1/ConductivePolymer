function [Na_amp_net,Na_amp,Vstep] = na_recovery(base_path, date, cell_idx, idx_f, idx_ttx, tstep, if_plot, if_save, save_path,color)
if nargin < 7, if_plot = 0; end
if nargin < 8, if_save = 0; end
if nargin < 10, color = 'darkblue'; end
preflix = strcat(strcat(date(5:6), date(1:4)),'_001');
path = fullfile(base_path, date, strcat('cell', string(cell_idx)));
if nargin < 9, save_path = path; end
for i = 1:length(tstep)
    filename{i} = sprintf('%s.na_recovery_%dms.%d.wcp', preflix,tstep(i), idx_f);
end

try
    for i = 1:length(tstep)
        out{i}=import_wcp(fullfile(path, filename{i}),'debug');
    end
catch
    date_plus_one = next_date_string(date);
    preflix = strcat(strcat(date_plus_one(5:6),date_plus_one(1:4)),'_001');
    for i = 1:length(tstep)
        filename{i} = sprintf('%s.na_recovery_%dms.%d.wcp', preflix,tstep(i),idx_f);
        out{i} = import_wcp(fullfile(path, filename{i}),'debug');
    end
end

if ~isnan(idx_ttx)
    for i = 1:length(tstep)
        filename_ttx{i} = sprintf('%s.na_recovery_%dms.%d.wcp', preflix,tstep(i),idx_ttx);
        out_ttx{i}=import_wcp(fullfile(path, filename_ttx{i}),'debug');
    end
end
n_channel = out{1}.channel_no;
n_recording = length(out);

dt = out{1}.T(2)-out{1}.T(1); % sample time in s

for i = 1:length(tstep)
    Im(:,i) = out{i}.S{3};
    if ~isnan(idx_ttx)
        Im_ttx(:,i) = out_ttx{i}.S{3};
    end
    Vm(:,i) = out{i}.S{4};
    
end
t = 0:dt:dt*(size(Vm,1)-1);

Na_amp_net = zeros(1, n_recording);
Na_amp = zeros(1, n_recording);

Vstep = [-90:10:10];




for i = 1:n_recording
    pulse_start = find(abs(diff(Vm(:,i)))>5);
    pulse_start(find(diff(pulse_start)<=2)+1) = [];
    pulse_end = pulse_start(2);
    pulse_start1 = pulse_start(1);
    pulse_start = pulse_start(3);
    I_baseline = mean(Im(1:pulse_start1,i));
    Im(:,i) = Im(:,i)-I_baseline;
    
%     I_baseline_ttx = mean(Im_ttx(1:pulse_start1,i));
%     Im_ttx(:,i) = Im_ttx(:,i)-I_baseline_ttx;
%     Im_net(:,i) = Im(:,i)-Im_ttx(:,i);
    if ~isnan(idx_ttx)
        I_baseline_ttx = mean(Im_ttx(1:pulse_start,i));
        Im_ttx(:,i) = Im_ttx(:,i)-I_baseline_ttx;
        Im_net(:,i) = Im(:,i)-Im_ttx(:,i);
    else
        Im_net(:,i) = Im(:,i)-Im(pulse_end,i);
    end
    
    [Na_amp_net(i),] = min(Im_net(pulse_start+10:pulse_start + 500, i));
    Na_amp(i) = min(Im(pulse_start+10:pulse_start + 500, i));

end

if if_plot
    addpath(genpath(fullfile(pwd,'GeneSetAnalysisMatlab')))
    if strcmp(color,'gray')
        cmap = gray(100);
        
    else
        cmap = custom_cmap(color);
        cmap = flip(cmap, 1);
    end
    colors = cmap(ceil(linspace(1, size(cmap,1),n_recording)),:);
    figure
    for i = 1:7
        plot(t, Im_net(:,i)-250*(i-1), 'Color', colors(i,:))
        hold on
    end
    xlim([0.95,1.3])
%     ylim([-3000,500])
end

if if_save
    if ~isnan(idx_ttx)
        save(fullfile(save_path,sprintf('%s_cell%s_na_recovery_%d.mat', date, string(cell_idx),idx_f)),'Im','Vm', 'Im_ttx', 'Vstep','Na_amp','Na_amp_net','Im_net', 'tstep');
    else
        save(fullfile(save_path,sprintf('%s_cell%s_na_recovery_%d.mat', date, string(cell_idx),idx_f)),'Im','Vm', 'Vstep','Na_amp','Na_amp_net','Im_net', 'tstep');
    end
end
fprintf('finished cell%s, from date %s, file %s\n',string(cell_idx), date, filename{1})
end