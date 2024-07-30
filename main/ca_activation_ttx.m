function [K_amp_adj,Ca_amp, Vstep] = ca_activation_ttx(base_path, date, cell_idx, idx_f, if_plot, if_save, color)
if nargin < 5, if_plot = 0; end
if nargin < 6, if_save = 0; end
if nargin < 7, color = 'darkblue'; end
preflix = strcat(strcat(date(5:6), date(1:4)),'_001');
path = fullfile(base_path, date, strcat('cell', string(cell_idx)));
filename = sprintf('%s.ca_activation.%d.wcp', preflix,idx_f);

try
    out=import_wcp(fullfile(path, filename),'debug');
catch
    date_plus_one = next_date_string(date);
    preflix = strcat(strcat(date_plus_one(5:6),date_plus_one(1:4)),'_001');
    filename = sprintf('%s.ca_activation.%d.wcp', preflix,idx_f);
    out=import_wcp(fullfile(path, filename),'debug');
end

n_channel = out.channel_no;
n_recording = length(out.rec_index);

dt = out.T(2)-out.T(1); % sample time in s

   
Im = out.S{3};
Vm = out.S{4};
t = 0:dt:dt*(size(Vm,1)-1);

K_amp = zeros(1, n_recording);
Ca_amp = zeros(1, n_recording);

Vstep = [-90:10:50];
pulse_start = find(abs(diff(Vm(:,1)))>5);
pulse_start(find(diff(pulse_start)<=2)+1) = [];
pulse_end = pulse_start(2);
pulse_start = pulse_start(1);

for i = 1:n_recording
    I_baseline = mean(out.S{3}(1:pulse_start,i));
    Im(:,i) = out.S{3}(:,i)-I_baseline;
    Im(:,i) = Im(:,i)-mean(Im(pulse_end-100:pulse_end, i));

    K_amp(i) = mean(out.S{3}(pulse_end-100:pulse_end,i)-I_baseline);
    Ca_amp(i) = min(Im(pulse_start+200:pulse_start + 2000, i));
end
Im_orig = out.S{3};
P = polyfit(Vstep(1:4), K_amp(1:4), 1); % leak K channels
K_amp_adj = K_amp-P(1)*Vstep-P(2);

if if_plot
    addpath(genpath(fullfile(pwd,'GeneSetAnalysisMatlab')))
    if strcmp(color,'gray')
        cmap = gray(100);
        cmap = flip(cmap, 1);
    else
        cmap = custom_cmap(color);
    end
    colors = cmap(ceil(linspace(1, size(cmap,1),n_recording)),:);
%     figure
%     for i = 1:n_recording
%         plot(out.T, Im_net(:,i), 'Color', colors(i,:))
%         hold on
%     end
%     xlim([0.498,0.508])
%     ylim([-3500,1000])
    figure
    for i = 1:n_recording
        plot(out.T, Im(:,i), 'Color', colors(i,:))
        hold on
    end
    xlim([0.5, 0.6])
    ylim([-1000,1000])
end

if if_save
    save(fullfile(path,sprintf('na_activation_%d.mat', idx_f)),'Im','Vm','Vstep','Im_orig','K_amp_adj', 'Ca_amp');
end
fprintf('finished cell%s, from date %s, file %s\n',string(cell_idx), date, filename)
end