function [Ca_amp,tstep] = na_recovery_ttx(base_path, date, cell_idx, idx_f, tstep, if_plot, if_save, save_path,color)
if nargin < 6, if_plot = 0; end
if nargin < 7, if_save = 0; end
if nargin < 9, color = 'darkblue'; end
preflix = strcat(strcat(date(5:6), date(1:4)),'_001');
path = fullfile(base_path, date, strcat('cell', string(cell_idx)));
if nargin < 8, save_path = path; end
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

n_channel = out{1}.channel_no;
n_recording = length(out);

dt = out{1}.T(2)-out{1}.T(1); % sample time in s

for i = 1:length(tstep)
    Im(:,i) = out{i}.S{3};
    Vm(:,i) = out{i}.S{4};
    
end
t = 0:dt:dt*(size(Vm,1)-1);

Ca_amp = zeros(1, n_recording);

Im_orig = Im;



for i = 1:n_recording
    pulse_start = find(abs(diff(Vm(:,i)))>5);
    pulse_start(find(diff(pulse_start)<=2)+1) = [];
    pulse_end1 = pulse_start(2);
    pulse_end = pulse_start(4);
    pulse_start1 = pulse_start(1);
    pulse_start = pulse_start(3);
    I_baseline = mean(Im(1:pulse_start1,i));
    Im(:,i) = Im(:,i)-I_baseline;
    Im(:,i) = Im(:,i)-mean(Im(pulse_end1-100:pulse_end1,i));
    Ca_amp(i) = min(Im(pulse_start+200:pulse_end , i));

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
        plot(t, Im_orig(:,i)-250*(i-1), 'Color', colors(i,:))
        hold on
    end
    xlim([0.95,1.3])
%     ylim([-3000,500])
end

if if_save
    save(fullfile(save_path,sprintf('%s_cell%s_na_recovery_%d.mat', date, string(cell_idx),idx_f)),'Im','Vm','Ca_amp','Im_orig', 'tstep');
end
fprintf('finished cell%s, from date %s, file %s\n',string(cell_idx), date, filename{1})
end