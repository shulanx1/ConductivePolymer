function [Na_amp_net,Na_amp, Vstep] = na_activation(base_path, date, cell_idx, idx_f, idx_ttx, if_plot, if_save, color)
if nargin < 6, if_plot = 0; end
if nargin < 7, if_save = 0; end
if nargin < 8, color = 'darkblue'; end
preflix = strcat(strcat(date(5:6), date(1:4)),'_001');
path = fullfile(base_path, date, strcat('cell', string(cell_idx)));
filename = sprintf('%s.na_activation.%d.wcp', preflix,idx_f);
filename_ttx = sprintf('%s.na_activation.%d.wcp', preflix,idx_ttx);

try
    out=import_wcp(fullfile(path, filename),'debug');
    out_ttx=import_wcp(fullfile(path, filename_ttx),'debug');
catch
    date_plus_one = next_date_string(date);
    preflix = strcat(strcat(date_plus_one(5:6),date_plus_one(1:4)),'_001');
    filename = sprintf('%s.na_activation.%d.wcp', preflix,idx_f);
    filename_ttx = sprintf('%s.na_activation.%d.wcp', preflix,idx_ttx);
    out=import_wcp(fullfile(path, filename),'debug');
    out_ttx=import_wcp(fullfile(path, filename_ttx),'debug');
end

n_channel = out.channel_no;
n_recording = length(out.rec_index);

dt = out.T(2)-out.T(1); % sample time in s

   
Im = out.S{3};
Im_ttx = out_ttx.S{3};
Im_net = zeros(size(Im));
Vm = out.S{4};
t = 0:dt:dt*(size(Vm,1)-1);

Na_amp_net = zeros(1, n_recording);
Na_amp = zeros(1, n_recording);

Vstep = [-90:10:10];
pulse_start = find(abs(diff(Vm(:,1)))>5);
pulse_start = pulse_start(1);
for i = 1:n_recording
    I_baseline = mean(Im(1:pulse_start,i));
    Im(:,i) = Im(:,i)-I_baseline;
    I_baseline_ttx = mean(Im_ttx(1:pulse_start,i));
    Im_ttx(:,i) = Im_ttx(:,i)-I_baseline_ttx;
%     if i == 1
%         control_step = mean(Im(pulse_start+1500:pulse_start+2000, 1)); % adjust for change in Rin
%         ttx_step = mean(Im_ttx(pulse_start+1500:pulse_start+2000, 1));
%         step_diff = (control_step-ttx_step)/(-90+75);
%     end
    Im_net(:,i) = Im(:,i)-Im_ttx(:,i);%-step_diff*(Vm(:,i)+75);
    [Na_amp_net(i),] = min(Im_net(pulse_start:pulse_start + 500, i));
    if i <= 4
        Na_amp(i) = min(Im(pulse_start+50:pulse_start + 500, i)); % to get rid of the Cm transient
    else
        Na_amp(i) = min(Im(pulse_start:pulse_start + 500, i));
    end
end

if if_plot
    addpath(genpath(fullfile(pwd,'GeneSetAnalysisMatlab')))
    if strcmp(color,'gray')
        cmap = gray(100);
        cmap = flip(cmap, 1);
    else
        cmap = custom_cmap(color);
    end
    colors = cmap(ceil(linspace(1, size(cmap,1),n_recording)),:);
    figure
    for i = 1:n_recording
        plot(out.T, Im_net(:,i), 'Color', colors(i,:))
        hold on
    end
    xlim([0.498,0.508])
    ylim([-3500,1000])
end

if if_save
    save(fullfile(path,sprintf('na_activation_%d.mat', idx_f)),'Im','Vm', 'Im_ttx', 'Vstep','Na_amp','Na_amp_net','Im_net');
end
fprintf('finished cell%s, from date %s, file %s\n',string(cell_idx), date, filename)
end