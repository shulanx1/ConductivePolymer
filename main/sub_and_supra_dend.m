function [params,steady_ratio,AP_ratio, f, I,Rin_pooled1,Rin_pooled2] = sub_and_supra_dend(base_path, date, cell, file_name, idx_f, if_plot, if_save, save_path)

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
    V_sub1(i) = mean(Vm1(floor(0.8/dt)+1:floor(1/dt),i));
    V_supra1(i) = mean(Vm1(floor(2.2/dt)+1:floor(2.5/dt),i));
    Vrest1(i) = mean(Vm1(1:floor(0.5/dt) ,i));
    FR2(i) = length(spike_time2{i})/0.5; % firing rate, Hz
    dv1(i) = d_Vm1(floor(2/dt)+1);
    RMP1(i) = quantile(Vm1(floor(2.2/dt)+1:floor(2.5/dt),i), 0.3);
    

    Ihold2 = mean(Im2(1:floor(0.5/dt),i)); % holding current, in pA
    V_sub2(i) = mean(Vm2(floor(0.8/dt)+1:floor(1/dt),i));
    V_supra2(i) = mean(Vm2(floor(2.2/dt)+1:floor(2.5/dt),i));
    Vrest2(i) = mean(Vm2(1:floor(0.5/dt) ,i));
    dv2(i) = d_Vm2(floor(2/dt)+1);
end
Istep_sub = -50:10:10*(n_recording-1)-50;
Istep_supra = -200:50:50*(n_recording-1)-200;
try
P1 = flip(robustfit(Istep_sub, dv1));
P2 = flip(robustfit(Istep_sub, dv2));
catch
    P1 = polyfit(dv1, Istep_sub,1); 
    P2 = polyfit(dv2, Istep_sub,1); 
end
if P1(1)>0.05
    V_sub1 = V_sub1-Vm1(floor(0.5/dt)+5,:)+Vm1(floor(0.5/dt),:);
    V_supra1 = V_supra1-Vm1(floor(0.5/dt)+5,:)+Vm1(floor(0.5/dt),:);
end
if abs(P2(1))>0.05
    
    Fs = 1/dt;
    Fc = [0 1000];
    Wn = Fc./(Fs/2);
    %b = fir1(10,Wn,'bandpass');
    [b,a] = butter(6, Wn(2), 'low');
    dV2_t = filtfilt(b,a,diff(Vm2(floor(0.5/dt)+3:floor(0.5/dt)+500, 1)));
    for ii = 1:length(dV2_t)-1
        if (abs(dV2_t(ii+1))<std(dV2_t(300:end))/2)&&(abs(dV2_t(ii))>std(dV2_t(300:end))/2)
            break
        end
    end
    if (ii == length(dV2_t)-1)||(ii<=8)
        d_bridge = 5;
    else
        d_bridge = ii-3;
    end
    V_sub2 = V_sub2-Vm2(floor(0.5/dt)+d_bridge,:)+Vm2(floor(0.5/dt),:);
    V_supra2 = V_supra2-Vm2(floor(0.5/dt)+d_bridge,:)+Vm2(floor(0.5/dt),:);
    if if_plot
        Vm2_temp = Vm2;
        for i = 1:size(Vm2,2)
            Vm2_temp(floor(0.5/dt):floor(1/dt),i) = Vm2_temp(floor(0.5/dt):floor(1/dt),i)-Vm2_temp(floor(0.5/dt)+d_bridge,i)+Vm2_temp(floor(0.5/dt),i);
            Vm2_temp(floor(2/dt):floor(2.5/dt),i) = Vm2_temp(floor(2/dt):floor(2.5/dt),i)-Vm2_temp(floor(2/dt)+d_bridge,i)+Vm2_temp(floor(2/dt),i);
        end
        figure,plot(out.T, Vm2_temp(:,1)-Vrest2(1))
        hold on
        plot(out.T, Vm2(:,1)-Vrest2(1))
    end
end

Ihold1 = sum(mean(Im1(1:floor(0.5/dt),:))/n_recording); % holding current, in pA
Vrest1 = mean(Vm1(1:floor(0.5/dt),:));
Ihold2 = sum(mean(Im2(1:floor(0.5/dt),:))/n_recording); % holding current, in pA                                                                                                                                                          
Vrest2 = mean(Vm2(1:floor(0.5/dt),:));

% P = polyfit(Istep_sub(1:8), V_sub1(1:8),1)*1e3;   % input resistance, MOhm
try
P = flip(robustfit(Istep_supra(1:6), V_supra1(1:6)))*1e3;   % input resistance, MOhm
Rin_pooled1 = P(1);
P = flip(robustfit(Istep_supra(1:6), V_supra2(1:6)))*1e3;   % input resistance, MOhm
Rin_pooled2 = P(1);
catch
    Rin_pooled1 = (V_sub1(1)-Vrest1(1))/Istep_sub(1)*1e3;
    Rin_pooled2 = (V_sub2(1)-Vrest2(1))/Istep_sub(1)*1e3;
end
Vsag = Vrest2(1)-Vm2(floor(2/dt)+1:floor(2.5/dt), 1);
params(2) = max(Vsag(1: floor(0.2/dt)))-mean(Vsag(floor(0.4/dt):end)); % sag ratio
params(1) = Rin_pooled2;  % input resistance, MOhm
params(7) = mean(Ihold2);  % holding current, pA
params(8) = mean(Vrest2); % RMP, mV


decay_R = Rin_pooled1/Rin_pooled2;

if length(FR2)>=12
    f = FR2(1:12);
    I = Istep_supra(1:12);
else
    f = zeros(1, 12);
    I = zeros(1, 12);
    f(1:length(FR2)) = FR2;
    I(1:length(FR2)) = Istep_supra;
end
idx_s = find(FR2>2);
idx_s(idx_s<5) = [];
idx_s1 = find(FR2>4);
idx_s1(idx_s1<5) = [];
isi = [];
for i = 1:length(idx_s1)
    isi = [isi,diff(spike_time2{idx_s1(i)})/(0.5/length(spike_time2{idx_s1(i)})) ];
end
edges = 0:0.1:3;
isi_hist = histcounts(isi, edges);
params(3) = sum(isi_hist(edges<=0.4))/sum(isi_hist); % burstyness
try
    params(4) = min(diff(spike_time2{min(idx_s1)}));  % isi @ rhobase
catch
    params(4) = nan;
end

% fit Rin and Cm
% i_fit_range = [1:4,8:11]; % subthreshold current that will be used for curve fitting
i_fit_range = [11]; % subthreshold current that will be used for curve fitting
i_fit_range(i_fit_range>size(Vm2,2)) = [];
Rin_fit = zeros(length(i_fit_range),1); % fit for two components
Cm_fit = zeros(length(i_fit_range),1);
R_square = zeros(length(i_fit_range), 1);
for i = 1:length(i_fit_range)
    ft = fittype('a1*(1-exp(-x/tau1))+b', 'independent', 'x');
    V_temp = Vm2(floor(0.5/dt):floor(0.65/dt), i_fit_range(i))-mean(Vm2(1:floor(0.5/dt),i_fit_range(i)));
    t_temp = [0:dt:dt*(length(V_temp)-1)]*1e3; %ms
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [Istep_sub(i_fit_range(i))*params(1)/1e3,0,5];
    opts.Lower = [min([Istep_sub(i_fit_range(i))*200/1e3,Istep_sub(i_fit_range(i))/1e3]), -5,0];
    opts.Upper = [max([Istep_sub(i_fit_range(i))/1e3,Istep_sub(i_fit_range(i))*200/1e3]),0,500]; 
    [curve,gof,~] = fit(t_temp', V_temp, ft, opts);    
    R_square(i) = gof.adjrsquare;
    Rin_fit(i) = curve.a1/Istep_sub(i_fit_range(i))*1e3; % MOhm
    Cm_fit(i, 1) = curve.tau1/Rin_fit(i, 1)*1e3; %pF
end
% Cm_fit = [];
% Rin_fit = [];
% if ~isempty(Cm_fit)&&(~isnan(Cm_fit))
%     params(5) = Cm_fit(1);
%     params(6) = Rin_fit(1);
% else
%     Rin_fit = zeros(length(i_fit_range),1); % fit for two components
%     Cm_fit = zeros(length(i_fit_range),1);
%     R_square = zeros(length(i_fit_range), 1);
%     for i = 1:length(i_fit_range)
%         ft = fittype('a1*(1-exp(-x/tau1))+a2*(1-exp(-x/tau2))', 'independent', 'x');
%         V_temp = Vm2(floor(0.5/dt):floor(0.65/dt), i_fit_range(i))-mean(Vm2(1:floor(0.5/dt),i_fit_range(i)));
%         t_temp = [0:dt:dt*(length(V_temp)-1)]*1e3; %ms
%         opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%         opts.Display = 'Off';
%         opts.StartPoint = [Istep_sub(i_fit_range(i))*params(1)/1e3,Istep_sub(i_fit_range(i))*params(1)/1e3, 5,5];
%         opts.Lower = [min([Istep_sub(i_fit_range(i))*200/1e3,Istep_sub(i_fit_range(i))/1e3]), min([Istep_sub(i_fit_range(i))*200/1e3,Istep_sub(i_fit_range(i))/1e3]), 0, 0];
%         opts.Upper = [max([Istep_sub(i_fit_range(i))/1e3,Istep_sub(i_fit_range(i))*200/1e3]), max([Istep_sub(i_fit_range(i))/1e3,Istep_sub(i_fit_range(i))*200/1e3]), 500,500];
%         [curve,gof,~] = fit(t_temp', V_temp, ft, opts);
%         R_square(i) = gof.adjrsquare;
%         [~,idx_larger] = max([curve.tau1,curve.tau2]);
%         if idx_larger==1
%             Rin_fit(i, 1) = curve.a1/Istep_sub(i_fit_range(i))*1e3; % MOhm
%             Rin_fit(i, 2) = curve.a2/Istep_sub(i_fit_range(i))*1e3;
%             Cm_fit(i, 1) = curve.tau1/Rin_fit(i, 1)*1e3; %pF
%             Cm_fit(i, 2) = curve.tau2/Rin_fit(i, 2)*1e3;
%         else
%             Rin_fit(i, 2) = curve.a1/Istep_sub(i_fit_range(i))*1e3; % MOhm
%             Rin_fit(i, 1) = curve.a2/Istep_sub(i_fit_range(i))*1e3;
%             Cm_fit(i, 2) = curve.tau1/Rin_fit(i, 1)*1e3; %pF
%             Cm_fit(i, 1) = curve.tau2/Rin_fit(i, 2)*1e3;
%         end
%         Rin_fit(i) = curve.a1/Istep_sub(i_fit_range(i))*1e3; % MOhm
%         Cm_fit(i, 1) = curve.tau1/Rin_fit(i, 1)*1e3; %pF
%     end
    Rin_fit = median(Rin_fit(R_square>0.7,:));
    Cm_fit = median(Cm_fit(R_square>0.7,:));
    if ~isempty(Cm_fit)
        params(5) = Cm_fit(1);
        params(6) = Rin_fit(1);
    else
        params(5) = nan;
        params(6) = nan;
    end
% end
AP_amp1 = [];
AP_amp1_first = [];
AP_width1 = [];
peak = find(idx_plot1==0);
rec_int = zeros(1, n_recording+1);
spike_n = FR2/2;
spike_n_cum = cumsum(spike_n);
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
    if ~isempty(find(spike_n_cum+1 == i))
        AP_amp1_first = [AP_amp1_first,spike_waveform1(peak,i)-RMP1(rec_idx)];
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
AP_amp2_first = [];
AP_width2 = [];

for i = 1:size(spike_waveform2,2)
    spike_w = spike_waveform2(:,i)-min(spike_waveform2(1:round(size(spike_waveform2,1)/2),i));
    if ~isempty(find(spike_n_cum+1 == i))
        AP_amp2_first = [AP_amp2_first,max(spike_w)];
    end
    AP_amp2 =[AP_amp2, max(spike_w)];
    [~,peak] = max(spike_w);
    a = abs(spike_w-max(spike_w)/2);
    [~, idx1] = min(a(1:peak));
    [~, idx2] = min(a(peak:length(spike_w)));
    AP_width2 = [AP_width2, ((peak + idx2)-idx1)*dt];
end
AP_ratio = median(AP_amp1_first)/median(AP_amp2_first);

try
dV1 = V_supra1([1:5,7])-Vrest1([1:5,7]);
dV2 = V_supra2([1:5,7])-Vrest2([1:5,7]);

P = flip(robustfit(dV2,dV1));
steady_ratio = P(1);
catch
    dV1 = V_sub1-Vrest1;
    dV2 = V_sub2-Vrest2;

    P = polyfit(dV1, dV2, 1);
    steady_ratio = P(1);
end

if if_plot
    if exist('d_bridge', 'var')
        Vm2_temp = Vm2;
        for i = 1:size(Vm2,2)
            Vm2_temp(floor(0.5/dt):floor(1/dt),i) = Vm2_temp(floor(0.5/dt):floor(1/dt),i)-Vm2_temp(floor(0.5/dt)+d_bridge,i)+Vm2_temp(floor(0.5/dt),i);
            Vm2_temp(floor(2/dt):floor(2.5/dt),i) = Vm2_temp(floor(2/dt):floor(2.5/dt),i)-Vm2_temp(floor(2/dt)+d_bridge,i)+Vm2_temp(floor(2/dt),i);
        end
        figure,plot(out.T, Vm1(:,1)-Vrest1(1))
        hold on, plot(out.T, Vm2_temp(:,1)-Vrest2(1))
        xlim([0.25,1.25])
        fprintf('dbridge: %d \n', d_bridge)
    else
        figure,plot(out.T, Vm1(:,1)-Vrest1(1))
        hold on, plot(out.T, Vm2(:,1)-Vrest2(1))
        xlim([0.25,1.25])
    end
end
if if_save
    save(fullfile(save_path,sprintf('%s_cell%s_%s_%d.mat', date, string(cell),file_name, idx_f)),'params', 'steady_ratio','AP_ratio','Rin_pooled1', 'Rin_pooled2', 'Vm1', 'Im1', 'Vm2', 'Im2', 'Istep_supra','spike_time1', 'spike_time2', 'FR2', 't','Ihold1', 'Ihold2','Vrest1','Vrest2', 'AP_amp1', 'AP_width1', 'AP_amp2', 'AP_width2', 'AP_ratio', 'steady_ratio')
end
fprintf('finished cell%s, from date %s, file %s\n',string(cell), date, filename)
end
