function [f, Z1_abs, Z2_abs] = zap(base_path, date, cell, idx_f, if_plot, if_save, save_path)

if nargin < 5, if_plot = 0; end
if nargin < 6, if_save = 1;end

file_name = 'zap';
preflix = strcat(strcat(date(5:6), date(1:4)),'_001');
path = fullfile(base_path, date, strcat('cell', string(cell)));
if nargin < 7, save_path = path; end

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
thr2 = 0;
if n_recording == 1 
    Vm1 = movmean(out.S{1},100);
    Vm2 = movmean(out.S{3},100);
else
    spike_num = zeros(1, n_recording);
    for i = 1:n_recording
    spike_time2 = [];
    for n = 2:size(out.S{3},1)
        if out.S{3}(n, i)>thr2 && out.S{3}(n-1, i)<=thr2          
            spike_time2 = [spike_time2, n*dt]; % spike time in s
        end
    end
    spike_num(i) = length(spike_time2);
    end
    
    Vm1 = mean(out.S{1}(:,find(spike_num==0)),2);
    Vm2 = mean(out.S{3}(:,find(spike_num==0)),2);
end
Im = mean(out.S{4},2);
t = out.T;

amp = 10*round((max(Im)-min(Im))/20);
f0 = 0;
f1 = 15;
T = 15;
a = (2*pi*f1-2*pi*f0)/2/T;
b = 2*pi*f0;
Izap = amp*sin(a*t.^2+b*t);
Izap = Izap';

L = length(t);
n = 2^nextpow2(L);
fs = 1/dt;

f = fs*(0:(n/2))/n;
df = f(2)-f(1);
f_V2 = fft(Vm2-mean(Vm2),n);
f_I = fft(Izap, n);
Z2 = f_V2(1:n/2+1)./f_I(1:n/2+1);
Z2_abs = abs(Z2)*1e3; %Mohm
phy2 = angle(Z2);

f_V1 = fft(Vm1-mean(Vm1),n);
f_I = fft(Izap, n);
Z1 = f_V1(1:n/2+1)./f_I(1:n/2+1);
Z1_abs = abs(Z1)*1e3; %Mohm
phy1 = angle(Z1);

if if_plot
    figure
    subplot(211)
    plot(t, Vm2-mean(Vm2), 'k', 'Linewidth', 1)
    hold on
    plot(t, Vm1-mean(Vm1), 'r', 'Linewidth', 1)
    xlabel('t (s)')
    ylabel('V (mV)')
    subplot(212)
    plot(f(ceil(1/df):floor(14/df)), movmean(Z2_abs(ceil(1/df):floor(14/df)),5), 'k', 'Linewidth', 1)
    hold on
    plot(f(ceil(1/df):floor(14/df)), movmean(Z1_abs(ceil(1/df):floor(14/df)),5), 'r', 'Linewidth', 1)
    xlabel('f (Hz)')
    ylabel('|Z| Mohm')
%     subplot(313)
%     plot(f(ceil(1/df):floor(14/df)), movmean(phy2(ceil(1/df):floor(14/df)),5), 'k', 'Linewidth', 1)
%     hold on
%     plot(f(ceil(1/df):floor(14/df)), movmean(phy2(ceil(1/df):floor(14/df)),5), 'r', 'Linewidth', 1)
%     xlabel('f (Hz)')
%     ylabel('phase rad')
end

if if_save
    save(fullfile(save_path, sprintf('%s_cell%d_zap_%d.mat', date, cell, idx_f)), 'f', 'Z2', 'Z2_abs', 'phy2', 'Z1', 'Z1_abs', 'phy1','t', 'Vm2', 'Vm1','Izap')
end
fprintf('finished cell%s, from date %s, file %s\n',string(cell), date, filename)
end