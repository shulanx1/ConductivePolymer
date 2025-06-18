function [Na_act_amp_net,Na_act_amp,Na_in_amp_net, Na_in_amp, K_amp_adj, Vstep] = na_inactivation(base_path, date, cell_idx, idx_f, idx_ttx, if_plot, if_save, save_path, if_opto, color)
if nargin < 6, if_plot = 0; end
if nargin < 7, if_save = 0; end
if nargin < 9, if_opto = 0; end
if nargin < 10, color = 'darkblue'; end
preflix = strcat(strcat(date(5:6), date(1:4)),'_001');
path = fullfile(base_path, date, strcat('cell', string(cell_idx)));
if nargin < 8, save_path = path; end
filename = sprintf('%s.na_inactivation.%d.wcp', preflix,idx_f);
if if_opto, filename = sprintf('%s.na_inactivation_opto.%d.wcp', preflix,idx_f);end

try
    out=import_wcp(fullfile(path, filename),'debug');
catch
    date_plus_one = next_date_string(date);
    preflix = strcat(strcat(date_plus_one(5:6),date_plus_one(1:4)),'_001');
    filename = sprintf('%s.na_inactivation.%d.wcp', preflix,idx_f);
    if if_opto, filename = sprintf('%s.na_inactivation_opto.%d.wcp', preflix,idx_f);end
    out=import_wcp(fullfile(path, filename),'debug');
end

if ~isnan(idx_ttx)
    filename_ttx = sprintf('%s.na_inactivation.%d.wcp', preflix,idx_ttx);
    if if_opto, filename_ttx = sprintf('%s.na_inactivation_opto.%d.wcp', preflix,idx_ttx);end
    out_ttx=import_wcp(fullfile(path, filename_ttx),'debug');
end

n_channel = out.channel_no;
n_recording = length(out.rec_index);

dt = out.T(2)-out.T(1); % sample time in s

   
Im = out.S{3};
Im_net = zeros(size(Im));
Vm = out.S{4};
t = 0:dt:dt*(size(Vm,1)-1);

Na_act_amp_net = zeros(1, n_recording);
Na_act_amp = zeros(1, n_recording);

Na_in_amp_net = zeros(1, n_recording);
Na_in_amp = zeros(1, n_recording);

Vstep = [-90:10:10];
pulse_start = find(abs(diff(Vm(:,1)))>5);
pulse_start(find(diff(pulse_start)<=2)+1) = [];
pulse_end = pulse_start(2);
pulse_end2 = pulse_start(3);
pulse_start = pulse_start(1);
% pulse_end = 20001;
% pulse_end2 = 21001;
% pulse_start = 10001;

for i = 1:n_recording
    I_baseline = mean(Im(1:pulse_start,i));
    Im(:,i) = Im(:,i)-I_baseline;
    Im(:,i) = Im(:,i)-mean(Im(pulse_end-100:pulse_end, i));
%     Im_ttx(:,i) = out_ttx.S{3}(:,i);
%     I_baseline_ttx = mean(Im_ttx(1:pulse_start,i));
%     Im_ttx(:,i) = Im_ttx(:,i)-I_baseline_ttx;
%     Im_ttx(:,i) = Im_ttx(:,i)-mean(Im_ttx(pulse_end-100:pulse_end, i));
% %     if i == 1
% %         control_step = mean(Im(pulse_start+1500:pulse_start+2000, 1)); % adjust for change in Rin
% %         ttx_step = mean(Im_ttx(pulse_start+1500:pulse_start+2000, 1));
% %         step_diff = (control_step-ttx_step)/(-90+75);
% %     end
%     Im_net(:,i) = Im(:,i)-Im_ttx(:,i);%-step_diff*(Vm(:,i)+75);
    if ~isnan(idx_ttx)
        Im_ttx(:,i) = out_ttx.S{3}(:,i);
        I_baseline_ttx = mean(Im_ttx(1:pulse_start,i));
        Im_ttx(:,i) = Im_ttx(:,i)-I_baseline_ttx;
        Im_ttx(:,i) = Im_ttx(:,i)-mean(Im_ttx(pulse_end-100:pulse_end, i));
        Im_net(:,i) = Im(:,i)-Im_ttx(:,i);%-step_diff*(Vm(:,i)+75);
    else
        Im_net(:,i) = Im(:,i)-Im(pulse_end,i);%-step_diff*(Vm(:,i)+75);
    end
    Na_act_amp_net(i) = min(Im_net(pulse_start:pulse_start + 500, i));
    if i <= 4
        Na_act_amp(i) = min(Im(pulse_start+50:pulse_start + 500, i)); % to get rid of the Cm transient
    else
        Na_act_amp(i) = min(Im(pulse_start:pulse_start + 500, i));
    end
    Na_in_amp_net(i) = min(Im_net(pulse_end:pulse_end + 500, i))-mean(Im_net(pulse_end-100:pulse_end, Vstep==0));
    Na_in_amp(i) = min(Im(pulse_end:pulse_end + 500, i))-mean(Im(pulse_end-100:pulse_end, Vstep==0)); 
    if ~isnan(idx_ttx)
        K_amp(i) = mean(out.S{3}(pulse_end-100:pulse_end-1,i)-I_baseline);
        K_amp_ttx(i) = mean(out_ttx.S{3}(pulse_end-100:pulse_end-1,i)-I_baseline_ttx);
    end
end
if ~isnan(idx_ttx)
    P = polyfit(Vstep(1:4), K_amp_ttx(1:4), 1); % leak K channels
    K_amp_adj = K_amp-P(1)*Vstep-P(2);
else
    K_amp_adj = nan(1, size(Im,2));
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
%     figure
%     for i = 1:n_recording
%         plot(out.T, Im_net(:,i), 'Color', colors(i,:))
%         hold on
%     end
%     xlim([0.498,0.508])
%     ylim([-3500,1000])
    figure
    for i = 1:n_recording
        plot(out.T, Im_net(:,i), 'Color', colors(i,:))
        hold on
    end
    xlim([0.998,1.008])
    ylim([-3500,1000])
end

if if_save
    if ~isnan(idx_ttx)
        save(fullfile(save_path,sprintf('%s_cell%s_na_activation_%d.mat', date, string(cell_idx),idx_f)),'Im','Vm', 'Im_ttx', 'Vstep','Na_act_amp','Na_act_amp_net','Na_in_amp','Na_in_amp_net','Im_net');
    else
        save(fullfile(save_path,sprintf('%s_cell%s_na_activation_%d.mat', date, string(cell_idx),idx_f)),'Im','Vm', 'Vstep','Na_act_amp','Na_act_amp_net','Na_in_amp','Na_in_amp_net','Im_net');
    end
end
fprintf('finished cell%s, from date %s, file %s\n',string(cell_idx), date, filename)
end