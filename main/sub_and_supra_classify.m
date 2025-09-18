function [params,Vsupra, f, I, AP_amp, AP_width] = sub_and_supra_classify(base_path, date, cell, idx_f, if_plot, if_save, save_path, if_opto)
% Vsupra: representive spike waveform (normalized)
% params: [input resistance, sag ratio, burstyness, min. isi at rhobase, Cm(fitted), Rm(fitted)]
if nargin < 5, if_plot = 0; end
if nargin < 6, if_save = 1; end
if nargin < 8, if_opto = 0; end


preflix = strcat(strcat(date(5:6), date(1:4)),'_001');
path = fullfile(base_path, date, strcat('cell', string(cell)));
if nargin < 7, save_path = path; end 
filename = sprintf('%s.sub and supra_Ch2.%d.wcp', preflix,idx_f);
if if_opto, filename = sprintf('%s.sub and supra_Ch2_opto.%d.wcp', preflix,idx_f); end
try
    out=import_wcp(fullfile(path, filename),'debug');
catch
    try
        date_plus_one = next_date_string(date);
        preflix = strcat(strcat(date_plus_one(5:6),date_plus_one(1:4)),'_001');
        filename = sprintf('%s.sub and supra_Ch2.%d.wcp', preflix,idx_f);
        if if_opto, filename = sprintf('%s.sub and supra_Ch2_opto.%d.wcp', preflix,idx_f); end
        out=import_wcp(fullfile(path, filename),'debug');
    catch
        try
            preflix = strcat(strcat(date(5:6), date(1:4)),'_001');
            filename = sprintf('%s.sub and supra.%d.wcp', preflix,idx_f);
            if if_opto, filename = sprintf('%s.sub and supra_opto.%d.wcp', preflix,idx_f); end
            try
                out=import_wcp(fullfile(path, filename),'debug');
            catch
                date_plus_one = next_date_string(date);
                preflix = strcat(strcat(date_plus_one(5:6),date_plus_one(1:4)),'_001');
                filename = sprintf('%s.sub and supra.%d.wcp', preflix,idx_f);
                if if_opto, filename = sprintf('%s.sub and supra_opto.%d.wcp', preflix,idx_f); end
                out=import_wcp(fullfile(path, filename),'debug');
            end
        catch
            preflix = strcat(strcat(date(5:6), date(1:4)),'_001');
            filename = sprintf('%s.sub and supra_Ch1.%d.wcp', preflix,idx_f);
            if if_opto, filename = sprintf('%s.sub and supra_Ch1_opto.%d.wcp', preflix,idx_f); end
            try
                out=import_wcp(fullfile(path, filename),'debug');
            catch
                date_plus_one = next_date_string(date);
                preflix = strcat(strcat(date_plus_one(5:6),date_plus_one(1:4)),'_001');
                filename = sprintf('%s.sub and supra_Ch1.%d.wcp', preflix,idx_f);
                if if_opto, filename = sprintf('%s.sub and supra_Ch1_opto.%d.wcp', preflix,idx_f); end
                out=import_wcp(fullfile(path, filename),'debug');
            end
        end
    end
end

n_channel = out.channel_no;
n_recording = length(out.rec_index);

dt = out.T(2)-out.T(1); % sample time in s
t = (-500*dt:dt:2500*dt)*1e3;


spike_time2 = {};
thr2 = 0;%-10;
spike_waveform2 = {};
FR2 = zeros(1, n_recording);
for i = 1:n_recording
    if  contains(filename,'Ch1')
        Vm2(:,i) = out.S{1}(:,i);
        Im2(:,i) = out.S{2}(:,i);
    else
        Vm2(:,i) = out.S{3}(:,i);
        Im2(:,i) = out.S{4}(:,i);
    end
    spike_time2{i} = [];
    spike_waveform2{i} = [];
    d_Vm2 = diff(Vm2(:,i));
    dv2(i) = d_Vm2(floor(2/dt)+1);
    V_sub2(i) = mean(Vm2(floor(0.7/dt)+1:floor(1/dt),i));
    Vrest2(i) = mean(Vm2(1:floor(0.5/dt) ,i));
    Ihold(i) = mean(Im2(1:floor(0.5/dt),i)); % holding current, in pA
    for n = floor(2/dt)+1:floor(2.5/dt)
        if Vm2(n, i)>thr2 && Vm2(n-1, i)<=thr2          
            [~, peak] = max(Vm2(n-10:n+10,i));
            spike_time2{i} = [spike_time2{i}, (n-9+peak)*dt]; % spike time in s
            spike_idx_temp = (n-9+peak);
            spike_waveform2{i} = [spike_waveform2{i}, Vm2(spike_idx_temp-10:spike_idx_temp+100,i)];             
        end
    end
    FR2(i) = length(spike_time2{i})/0.5;
end

Istep_sub = -50:10:10*(n_recording-1)-50;
P2 = polyfit(Istep_sub, dv2,1);
if P2(1)>0.05
    V_sub2 = V_sub2-Vm2(floor(0.5/dt)+5,:);
end

P = polyfit(Istep_sub(1:6), V_sub2(1:6),1)*1e3;  
params(1) = P(1);  % input resistance, MOhm
params(7) = mean(Ihold);  % holding current, pA
params(8) = mean(Vrest2); % RMP, mV

Vsag = Vrest2(1)-Vm2(floor(2/dt)+1:floor(2.5/dt), 1);
params(2) = max(Vsag(1: floor(0.2/dt)))-mean(Vsag(floor(0.4/dt):end)); % sag ratio
Istep_supra = -200:50:50*(n_recording-1)-200;
len = 14;

if length(FR2)>=len
    f = FR2(1:len);
    I = Istep_supra(1:len);
else
    f = zeros(1, len);
    I = zeros(1, len);
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
i_fit_range = [1:4,8:11]; % subthreshold current that will be used for curve fitting
i_fit_range(i_fit_range>size(Vm2,2)) = [];
Rin_fit = zeros(length(i_fit_range),1); % fit for two components
Cm_fit = zeros(length(i_fit_range),1);
R_square = zeros(length(i_fit_range), 1);
for i = 1:length(i_fit_range)
%     ft = fittype('a1*(1-exp(-x/tau1))+a2*(1-exp(-x/tau2))', 'independent', 'x');
%     V_temp = Vm2(floor(0.5/dt):floor(0.65/dt), i_fit_range(i))-mean(Vm2(1:floor(0.5/dt),i_fit_range(i)));
%     t_temp = [0:dt:dt*(length(V_temp)-1)]*1e3; %ms
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     opts.StartPoint = [Istep_sub(i_fit_range(i))*params(1)/1e3,Istep_sub(i_fit_range(i))*params(1)/1e3, 5,5];
%     opts.Lower = [min([Istep_sub(i_fit_range(i))*200/1e3,Istep_sub(i_fit_range(i))/1e3]), min([Istep_sub(i_fit_range(i))*200/1e3,Istep_sub(i_fit_range(i))/1e3]), 0, 0];
%     opts.Upper = [max([Istep_sub(i_fit_range(i))/1e3,Istep_sub(i_fit_range(i))*200/1e3]), max([Istep_sub(i_fit_range(i))/1e3,Istep_sub(i_fit_range(i))*200/1e3]), 500,500];
%     [curve,gof,~] = fit(t_temp', V_temp, ft, opts);
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
%     [~,idx_larger] = max([curve.tau1,curve.tau2]);
%     if idx_larger==1
%         Rin_fit(i, 1) = curve.a1/Istep_sub(i_fit_range(i))*1e3; % MOhm
%         Rin_fit(i, 2) = curve.a2/Istep_sub(i_fit_range(i))*1e3;
%         Cm_fit(i, 1) = curve.tau1/Rin_fit(i, 1)*1e3; %pF
%         Cm_fit(i, 2) = curve.tau2/Rin_fit(i, 2)*1e3;
%     else
%         Rin_fit(i, 2) = curve.a1/Istep_sub(i_fit_range(i))*1e3; % MOhm
%         Rin_fit(i, 1) = curve.a2/Istep_sub(i_fit_range(i))*1e3;
%         Cm_fit(i, 2) = curve.tau1/Rin_fit(i, 1)*1e3; %pF
%         Cm_fit(i, 1) = curve.tau2/Rin_fit(i, 2)*1e3;
%     end
    Rin_fit(i) = curve.a1/Istep_sub(i_fit_range(i))*1e3; % MOhm
    Cm_fit(i, 1) = curve.tau1/Rin_fit(i, 1)*1e3; %pF
end

% Rin_fit = mean(Rin_fit(R_square>0.85,:));
% Cm_fit = mean(Cm_fit(R_square>0.85,:));
Cm_fit = [];
Rin_fit = [];
if ~isempty(Cm_fit)&&(~isnan(Cm_fit))
    params(5) = Cm_fit(1);
    params(6) = Rin_fit(1);
else
    Rin_fit = zeros(length(i_fit_range),1); % fit for two components
    Cm_fit = zeros(length(i_fit_range),1);
    R_square = zeros(length(i_fit_range), 1);
    for i = 1:length(i_fit_range)
        ft = fittype('a1*(1-exp(-x/tau1))+a2*(1-exp(-x/tau2))', 'independent', 'x');
        V_temp = Vm2(floor(0.5/dt):floor(0.65/dt), i_fit_range(i))-mean(Vm2(1:floor(0.5/dt),i_fit_range(i)));
        t_temp = [0:dt:dt*(length(V_temp)-1)]*1e3; %ms
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [Istep_sub(i_fit_range(i))*params(1)/1e3,Istep_sub(i_fit_range(i))*params(1)/1e3, 5,5];
        opts.Lower = [min([Istep_sub(i_fit_range(i))*200/1e3,Istep_sub(i_fit_range(i))/1e3]), min([Istep_sub(i_fit_range(i))*200/1e3,Istep_sub(i_fit_range(i))/1e3]), 0, 0];
        opts.Upper = [max([Istep_sub(i_fit_range(i))/1e3,Istep_sub(i_fit_range(i))*200/1e3]), max([Istep_sub(i_fit_range(i))/1e3,Istep_sub(i_fit_range(i))*200/1e3]), 500,500];
        [curve,gof,~] = fit(t_temp', V_temp, ft, opts);
        R_square(i) = gof.adjrsquare;
        %     ft = fittype('a1*(1-exp(-x/tau1))+b', 'independent', 'x');
        %     V_temp = Vm2(floor(0.5/dt):floor(0.65/dt), i_fit_range(i))-mean(Vm2(1:floor(0.5/dt),i_fit_range(i)));
        %     t_temp = [0:dt:dt*(length(V_temp)-1)]*1e3; %ms
        %     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        %     opts.Display = 'Off';
        %     opts.StartPoint = [Istep_sub(i_fit_range(i))*params(1)/1e3,0,5];
        %     opts.Lower = [min([Istep_sub(i_fit_range(i))*200/1e3,Istep_sub(i_fit_range(i))/1e3]), -5,0];
        %     opts.Upper = [max([Istep_sub(i_fit_range(i))/1e3,Istep_sub(i_fit_range(i))*200/1e3]),0,500];
        %     [curve,gof,~] = fit(t_temp', V_temp, ft, opts);
        %     R_square(i) = gof.adjrsquare;
        [~,idx_larger] = max([curve.tau1,curve.tau2]);
        if idx_larger==1
            Rin_fit(i, 1) = curve.a1/Istep_sub(i_fit_range(i))*1e3; % MOhm
            Rin_fit(i, 2) = curve.a2/Istep_sub(i_fit_range(i))*1e3;
            Cm_fit(i, 1) = curve.tau1/Rin_fit(i, 1)*1e3; %pF
            Cm_fit(i, 2) = curve.tau2/Rin_fit(i, 2)*1e3;
        else
            Rin_fit(i, 2) = curve.a1/Istep_sub(i_fit_range(i))*1e3; % MOhm
            Rin_fit(i, 1) = curve.a2/Istep_sub(i_fit_range(i))*1e3;
            Cm_fit(i, 2) = curve.tau1/Rin_fit(i, 1)*1e3; %pF
            Cm_fit(i, 1) = curve.tau2/Rin_fit(i, 2)*1e3;
        end
        Rin_fit(i) = curve.a1/Istep_sub(i_fit_range(i))*1e3; % MOhm
        Cm_fit(i, 1) = curve.tau1/Rin_fit(i, 1)*1e3; %pF
    end
    Rin_fit = mean(Rin_fit(R_square>0.75,:));
    Cm_fit = mean(Cm_fit(R_square>0.75,:));
    if ~isempty(Cm_fit)
        params(5) = Cm_fit(1);
        params(6) = Rin_fit(1);
    else
        params(5) = nan;
        params(6) = nan;
    end
end

Vsupra = [];
AP_amp = [];
AP_width = [];
if ~isempty(idx_s)
    if length(spike_time2)>=12
        idx_range = linspace(min(idx_s), 12, 4);
    else
        idx_range = linspace(min(idx_s), length(spike_time2), 4);
    end
    for i = round(idx_range(3))
        idx_waveform = find(diff(spike_time2{i})>100*dt);
        while (length(idx_waveform)==0)||(length(idx_waveform)<length(spike_time2{i})/2)
            if i==n_recording
                break
            end
            i = i + 1;
            idx_waveform = find(diff(spike_time2{i})>100*dt);
        end
        if isempty(idx_waveform)
            i = n_recording;
            idx_waveform = find(diff(spike_time2{i})>100*dt);
            while (length(idx_waveform)==0)||(length(idx_waveform)<length(spike_time2{i})/2)
                if i==1
                    break
                end
                i = i - 1;
                idx_waveform = find(diff(spike_time2{i})>100*dt);
            end
        end
        waveform_temp = (spike_waveform2{i}(:,idx_waveform)-mean(spike_waveform2{i}(:,idx_waveform)))./std(spike_waveform2{i}(:,idx_waveform));
        [U,S,V] = svd(waveform_temp);
        T = waveform_temp*V(:,1);
        l = waveform_temp'*T;
        [~, loc] = max(l);
        spike_temp = waveform_temp(:,loc)-mean(Vm2(floor(2/dt)+1:floor(2.5/dt), i));
        Vsupra = [Vsupra;(spike_temp-min(spike_temp))/(max(spike_temp)-min(spike_temp))];
        spike_waveform = spike_waveform2{i}(:,idx_waveform)-spike_waveform2{i}(1,idx_waveform);
        for ii = 1:size(spike_waveform,2)
            [amp_temp,max_amp] = max(spike_waveform(:,ii));
            [~,II] = sort(abs(spike_waveform(:,ii)-amp_temp/2));
            b = find(II<max_amp);
            a = find(II>max_amp);
            try
                halfwidth = (II(a(1))-II(b(1)))*dt*1000; %in ms
                AP_width = [AP_width,halfwidth];
            catch
            end
        end
%         AP_amp = [AP_amp, max(spike_waveform2{i}(:,idx_waveform))-Vrest2(round(idx_range(3)))];
        AP_amp = [AP_amp, max(spike_waveform2{i}(:,idx_waveform))-spike_waveform2{i}(1,idx_waveform)];
    end
end

if if_plot
    t = (-0.5:dt:1)*1e3;
    idx_plot = [1,12];

    figure
    hold on
    len = size(Im2,2);
    red = [0,0,0];
    pink = [0.8,0.8,0.8];%[255, 192, 203]/255;
    colors_p = [linspace(pink(1),red(1),len)', linspace(pink(2),red(2),len)', linspace(pink(3),red(3),len)'];

    for i =  1:length(idx_plot)
    plot(t, Vm2(floor(1.5/dt):floor(3/dt),idx_plot(i)),'Color', colors_p(idx_plot(i),:,:))
    end
    xlim([-500, 1000])
end
if if_save
    save(fullfile(save_path,sprintf('%s_cell%s_sub_and_supra_classify_%d.mat',date, string(cell), idx_f)),'params', 'Vsupra', 'Istep_supra','spike_time2', 'f','I', 't','Vrest2','isi_hist', 'edges','AP_amp','AP_width')
end

fprintf('finished cell%s, from date %s, file %s\n',string(cell), date, filename)
end
