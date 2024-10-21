function [pre, stim, post] = stim_sub(base_path, date, cell_idx,file_name, idx_f, if_plot, if_save, save_path)
if nargin < 6, if_plot = 0; end
if nargin < 7, if_save = 0; end

preflix = strcat(strcat(date(5:6), date(1:4)),'_001');
path = fullfile(base_path, date, strcat('cell', string(cell_idx)));
if nargin < 8, save_path = path; end

filename = sprintf('%s.%s.%d.wcp', preflix,file_name, idx_f);
try
    out=import_wcp(fullfile(path, filename),'debug');
catch
    date_plus_one = next_date_string(date);
    preflix = strcat(strcat(date_plus_one(5:6),date_plus_one(1:4)),'_001');
    filename = sprintf('%s.%s.%d.wcp', preflix,file_name,idx_f);
    out=import_wcp(fullfile(path, filename),'debug');
end

n_channel = out.channel_no;
n_recording = 1;%length(out.rec_index);

dt = out.T(2)-out.T(1); % sample time in s
t = (-500*dt:dt:2500*dt)*1e3;


spike_time2 = {};
thr2 = -10;%-10;
spike_waveform2 = {};
FR2 = zeros(1, n_recording);
stim_time_idx = {};
stim_end_idx = {};
stim_time = {};
stim_end = {};

for i = 1:n_recording
    Vm2(:,i) = out.S{3}(:,i);
    Im2(:,i) = out.S{4}(:,i);
    power(:,i) = out.S{7}(:,i);
    stim_time_idx{i} = find(diff(power(:,i))>15);
    stim_time_idx{i}(find(diff(stim_time_idx{i})<=1)+1) = [];
    stim_end_idx{i} = find(diff(power(:,i))<-15);
    stim_end_idx{i}(find(diff(stim_time_idx{i})<=1)+1) = [];
    
    stim_time{i} = out.T(stim_time_idx{i});
    stim_end{i} = out.T(stim_end_idx{i});
    spike_time2{i} = [];
%     spike_waveform2{i} = [];
    for n = 2:length(Vm2(:,i))
        if Vm2(n, i)>thr2 && Vm2(n-1, i)<=thr2   
            if (n-10>0)&&(n+10<=size(Vm2,1))
            [~, peak] = max(Vm2(n-10:n+10,i));
            spike_time2{i} = [spike_time2{i}, (n-9+peak)*dt]; % spike time in s
            spike_idx_temp = (n-9+peak);
%             spike_waveform2{i} = [spike_waveform2{i}, Vm2(spike_idx_temp-10:spike_idx_temp+100,i)];        
            end
        end
    end
    if ~isempty(spike_time2{i})% remove spikes
        for n = 1:length(spike_time2{i})
            idx_temp = [max(11,floor(spike_time2{i}(n)/dt)-100):1:min(size(Vm2,1)-10,floor(spike_time2{i}(n)/dt)+1200)];
            Vm2(idx_temp,i) = mean(Vm2([idx_temp(1)-10:idx_temp(1)-1,idx_temp(end)+1:idx_temp(end)+10],i));
        end
    end
end


pre = nan(n_recording, length(stim_time{1}));
stim = nan(n_recording, length(stim_time{1}));
post = nan(n_recording, length(stim_time{1}));

for i = 1:n_recording
    for n = 1:length(stim_time{i})
        pre(i, n) = mean(Vm2(stim_time_idx{i}(n)-10:stim_time_idx{i}(n),i));
        stim(i, n) = mean(Vm2(stim_end_idx{i}(n)-10:stim_end_idx{i}(n),i));
        post(i, n) = mean(Vm2(stim_end_idx{i}(n)+2000:stim_end_idx{i}(n)+2000+10,i));
    end

end


if if_plot
    figure
    for n = 1:length(stim_time{1})
        plot(out.T-stim_time{1}(n), Vm2(:,1))
        hold on
    end
    xlim([-0.2,stim_end{1}(1)-stim_time{1}(1)+0.2])
end

if if_save
    save(fullfile(save_path,sprintf('%s_cell%s_power_sub_%d.mat', date, string(cell_idx), idx_f)),'idx_f','Vm2','spike_time2', 'stim_time','stim_end','stim_time_idx','stim_end_idx','pre','stim','post')

end
fprintf('finished cell%s, from date %s, file %s\n',string(cell_idx), date, filename)
end