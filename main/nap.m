function Na_amp_net = nap(base_path, date, cell_idx, idx_f, idx_ttx, Vstep,if_plot, if_save,save_path, color)
if nargin < 7, if_plot = 0; end
if nargin < 8, if_save = 0; end
if nargin < 10, color = [119,176,203]/255; end
preflix = strcat(strcat(date(5:6), date(1:4)),'_001');
path = fullfile(base_path, date, strcat('cell', string(cell_idx)));
if nargin < 9, save_path = path; end
filename = sprintf('%s.inap.%d.wcp', preflix,idx_f);
filename_ttx = sprintf('%s.inap.%d.wcp', preflix,idx_ttx);


try
    out=import_wcp(fullfile(path, filename),'debug');
    out_ttx=import_wcp(fullfile(path, filename_ttx),'debug');
catch
    date_plus_one = next_date_string(date);
    preflix = strcat(strcat(date_plus_one(5:6),date_plus_one(1:4)),'_001');
    filename = sprintf('%s.inap.%d.wcp', preflix,idx_f);
    out=import_wcp(fullfile(path, filename),'debug');
    filename_ttx = sprintf('%s.inap.%d.wcp', preflix,idx_ttx);
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
pulse_start = find(abs(diff(Vm(:,1)))>5);
pulse_start(find(diff(pulse_start)<=2)+1) = [];
pulse_start(2) = (pulse_start(2)-pulse_start(1))/2+pulse_start(1);
for i = 1:n_recording
    I_baseline = mean(Im(1:pulse_start,i));
    Im(:,i) = Im(:,i)-I_baseline;
    I_baseline_ttx = mean(Im_ttx(1:pulse_start,i));
    Im_ttx(:,i) = Im_ttx(:,i)-I_baseline_ttx;
    Im_net(:,i) = Im(:,i)-Im_ttx(:,i);%-step_diff*(Vm(:,i)+75);
end


% remove APs
spike_time = {}; 
spike_idx = {}; 
thr = -500;
for i = 1:n_recording
    spike_time{i} = [];
    spike_idx{i} = [];
    for n = 2:size(Im_net(:,i),1)
        if Im_net(n,i)<thr && Im_net(n-1,i)>=thr
            spike_time{i} = [spike_time{i}, n*dt]; % spike time in s
            spike_idx{i} = [spike_idx{i}, n]; % spike index
        end
    end
end

for i = 1:n_recording
    idx_rmv = [];
    for n = 1:length(spike_idx{i})
        idx_rmv = [idx_rmv,max(1,spike_idx{i}(n)-50):min(spike_idx{i}(n) +500, size(Im,1))];
    end
    idx_rmv = sort(unique(idx_rmv));
    if ~isempty(idx_rmv)
        Im(idx_rmv) = [];
        Im_net(idx_rmv) = [];
        Vm(idx_rmv) = [];
        Im_ttx(idx_rmv) = [];
        t(idx_rmv) = [];
    end
end

Na_amp_net = zeros(1, length(Vstep));


Ena = 53.9;
for i = 1:length(Vstep)
    V_idx = find(abs(Vm(pulse_start(1):pulse_start(2),1)-Vstep(i))<=0.25);
    Na_amp_net(i) = mean(Im_net(V_idx,1));
end

if if_plot
    figure
    for i = 1:n_recording
        plot(t, Im_net(:,i), 'Color', color)
        hold on
        xlim([t(pulse_start(1)), t(pulse_start(2))])
    end
end

if if_save
    save(fullfile(save_path,sprintf('%s_cell%s_nap_%d.mat', date,string(cell_idx),idx_f)),'Im','Vm', 'Im_ttx', 'Vstep','Na_amp_net','Im_net');
end
fprintf('finished cell%s, from date %s, file %s\n',string(cell_idx), date, filename)
end