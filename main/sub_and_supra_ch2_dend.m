function [decay_R, Rin_pooled1,Rin_pooled2] = sub_and_supra_ch2_dend(base_path, date, cell, file_name, idx_f, if_plot, if_save, save_path)

if nargin < 6, if_plot = 0; end
if nargin < 7, if_save = 1;end


preflix = strcat(strcat(date(5:6), date(1:4)),'_001');
path = fullfile(base_path, date, strcat('cell', string(cell)));
if nargin < 8, save_path = path; end

filename = sprintf('%s.%s.%d.wcp', preflix,file_name, idx_f);
try
    out=import_wcp(fullfile(path, filename),'debug');
catch
    date_plus_one = next_date_string(date);
    preflix = strcat(strcat(date_plus_one(5:6),date_plus_one(1:4)),'_001');
    filename = sprintf('%s.%s.%d.wcp', preflix,file_name, idx_f);
    out=import_wcp(fullfile(path, filename),'debug');
end

n_channel = out.channel_no;
n_recording = length(out.rec_index);

dt = out.T(2)-out.T(1); % sample time in s
t = (-500*dt:dt:2500*dt)*1e3;

spike_time1 = {}; % spike time sorted as the rising phase when cross threshold
spike_time2 = {};
thr1 = -40;
thr2 = -10;
spike_waveform1 = [];
spike_waveform2 = [];
idx_plot1 = -30:30;
idx_plot2 = -40:40;
for i = 1:n_recording
    Vm1(:,i) = out.S{1}(:,i);
    Im1(:,i) = out.S{2}(:,i);
    Vm2(:,i) = out.S{3}(:,i);
    Im2(:,i) = out.S{4}(:,i);
    spike_time1{i} = [];
    spike_time2{i} = [];
    d_Vm1 = diff(Vm1(:,i));
    d_Vm2 = diff(Vm2(:,i));
    for n = floor(2/dt)+1:floor(2.5/dt)
        if Vm2(n, i)>thr2 && Vm2(n-1, i)<=thr2          
            [~, peak] = max(Vm2(n-10:n+10,i));
            spike_time2{i} = [spike_time2{i}, (n-9+peak)*dt]; % spike time in s
            spike_idx_temp = (n-9+peak);
            spike_waveform2 = [spike_waveform2, Vm2(spike_idx_temp+idx_plot2,i)];  
            bAP_idx_temp = find_bAP(spike_idx_temp, Vm1(:,i), dt, thr1);
            if isempty(bAP_idx_temp)
                thr1 = quantile(Vm1(:,i),0.95);
                bAP_idx_temp = find_bAP(spike_idx_temp, Vm1(:,i), dt, thr1);
            end
            if ~isempty(bAP_idx_temp)
                spike_waveform1 = [spike_waveform1, Vm1(bAP_idx_temp+idx_plot1,i)];  
                spike_time1{i} = [spike_time1{i}, bAP_idx_temp*dt];
            else
                spike_waveform1 = [spike_waveform1, zeros(length(idx_plot1),1)];
                spike_time1{i} = [spike_time1{i}, nan];
            end
        end
    end
    
    Ihold1 = mean(Im1(1:floor(0.5/dt),i)); % holding current, in pA
    V_sub1(i) = mean(Vm1(floor(0.7/dt)+1:floor(1/dt),i));
    V_supra1(i) = mean(Vm1(floor(2.2/dt)+1:floor(2.5/dt),i));
    Vrest1(i) = mean(Vm1(1:floor(0.5/dt) ,i));
    FR(i) = length(spike_time2{i})/0.5; % firing rate, Hz
    dv1(i) = d_Vm1(floor(2/dt)+1);
    RMP1(i) = quantile(Vm1(floor(2.2/dt)+1:floor(2.5/dt),i), 0.3);
    

    Ihold2 = mean(Im2(1:floor(0.5/dt),i)); % holding current, in pA
    V_sub2(i) = mean(Vm2(floor(0.7/dt)+1:floor(1/dt),i));
    V_supra2(i) = mean(Vm2(floor(2.2/dt)+1:floor(2.5/dt),i));
    Vrest2(i) = mean(Vm2(1:floor(0.5/dt) ,i));
    dv2(i) = d_Vm2(floor(2/dt)+1);
end
Istep_sub = -50:10:10*(n_recording-1)-50;
Istep_supra = -200:50:50*(n_recording-1)-200;
% P1 = polyfit(Istep_sub, dv1,1);
% P2 = polyfit(Istep_sub, dv2,1);
% if P1(1)>0.05
%     V_sub1 = V_sub1-Vm1(floor(0.5/dt)+5,:);
% end
% if P2(1)>0.05
%     V_sub2 = V_sub2-Vm2(floor(0.5/dt)+5,:);
% end
Ihold1 = sum(mean(Im1(1:floor(0.5/dt),:))/n_recording); % holding current, in pA
Vrest1 = mean(Vm1(1:floor(0.5/dt),:));
Ihold2 = sum(mean(Im2(1:floor(0.5/dt),:))/n_recording); % holding current, in pA                                                                                                                                                          
Vrest2 = mean(Vm2(1:floor(0.5/dt),:));

% P = polyfit(Istep_sub(1:8), V_sub1(1:8),1)*1e3;   % input resistance, MOhm
try
    P = polyfit(Istep_sub(1:8), V_sub1(1:8)-Vrest1(1:8),1)*1e3;   % input resistance, MOhm
    Rin_pooled1 = P(1);
    P = polyfit(Istep_sub(1:8), V_sub2(1:8)-Vrest2(1:8),1)*1e3;   % input resistance, MOhm
    Rin_pooled2 = P(1);
catch
    Rin_pooled1 = (V_sub1(1)-Vrest1(1))/Istep_sub(1)*1e3;
    Rin_pooled2 = (V_sub2(1)-Vrest2(1))/Istep_sub(1)*1e3;
end
decay_R = Rin_pooled1/Rin_pooled2;

AP_amp1 = [];
AP_width1 = [];
peak = find(idx_plot1==0);
rec_int = zeros(1, n_recording+1);

for i = 2:n_recording+1
    rec_int(i) =rec_int(i-1)+length(spike_time2{i-1});
end
    
if ~isempty(spike_waveform1)
for i = 1:size(spike_waveform1,2)
    rec_idx = 0;
    for j = 2:length(rec_int)
        if (i > rec_int(j-1))&&(i <= rec_int(j))
            rec_idx = j-1;
            break
        end
    end
    AP_amp1 = [AP_amp1,spike_waveform1(peak,i)-RMP1(rec_idx)];
    [~,locs,width,~]  = findpeaks(spike_waveform1(:,i), 'MinPeakProminence',1, 'MinPeakDistance',floor(0.005/dt));
    if isempty(locs)
        AP_width1 = [AP_width1, nan];
    else
        if length(locs)>1
            [~,idx_locs] = min(locs-peak);
            AP_width1 = [AP_width1, width(idx_locs)*dt];
        else
            AP_width1 = [AP_width1, width*dt];
        end
    end
%     spike_w = spike_waveform1(:,i)-min(spike_waveform1(:,i));
%     a = abs(spike_w-spike_w(peak)/2);
%     [~, idx1] = min(a(1:peak));
%     [~, idx2] = min(a(peak:length(spike_w)));
%     AP_width1 = [AP_width1, ((peak + idx2)-idx1)*dt];
end
end

AP_amp2 = [];
AP_width2 = [];
for i = 1:size(spike_waveform2,2)
    spike_w = spike_waveform2(:,i)-min(spike_waveform2(1:round(size(spike_waveform2,1)/2),i));
    AP_amp2 =[AP_amp2, max(spike_w)];
    [~,peak] = max(spike_w);
    a = abs(spike_w-max(spike_w)/2);
    [~, idx1] = min(a(1:peak));
    [~, idx2] = min(a(peak:length(spike_w)));
    AP_width2 = [AP_width2, ((peak + idx2)-idx1)*dt];
end
AP_ratio = median(AP_amp1)/median(AP_amp2);
steady_ratio = Rin_pooled1/Rin_pooled2;

if if_plot
figure,plot(out.T, out.S{1}(:,1)-Vrest1(1))
hold on, plot(out.T, out.S{3}(:,1)-Vrest2(1))
xlim([0.25,1.25])
end
if if_save
    save(fullfile(save_path,sprintf('%s_cell%s_%s_%d.mat', date, string(cell),file_name, idx_f)), 'decay_R','Rin_pooled1', 'Rin_pooled2', 'Vm1', 'Im1', 'Vm2', 'Im2', 'Istep_supra','spike_time1', 'spike_time2', 'FR', 't','Ihold1', 'Ihold2','Vrest1','Vrest2', 'AP_amp1', 'AP_width1', 'AP_amp2', 'AP_width2', 'AP_ratio', 'steady_ratio')
end
fprintf('finished cell%s, from date %s, file %s\n',string(cell), date, filename)
end
