function [pre_stim_FR, stim_FR, post_stim_FR, dur_stim, no_stim] = stim_AP(base_path, date, cell_idx, idx_f, if_plot, if_save, save_path)
if nargin < 5, if_plot = 0; end
if nargin < 6, if_save = 0; end

preflix = strcat(strcat(date(5:6), date(1:4)),'_001');
path = fullfile(base_path, date, strcat('cell', string(cell_idx)));
if nargin < 7, save_path = path; end

filename = sprintf('%s.power_long.%d.wcp', preflix,idx_f);
try
    out=import_wcp(fullfile(path, filename),'debug');
catch
    date_plus_one = next_date_string(date);
    preflix = strcat(strcat(date_plus_one(5:6),date_plus_one(1:4)),'_001');
    filename = sprintf('%s.power_long.%d.wcp', preflix,idx_f);
    out=import_wcp(fullfile(path, filename),'debug');
end

n_channel = out.channel_no;
n_recording = 1;%length(out.rec_index);

dt = out.T(2)-out.T(1); % sample time in s
t = (-500*dt:dt:2500*dt)*1e3;


spike_time2 = {};
%thr2 = -10;%-10;
thr2 = 0; % idx = 122:127 nPBDF
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
    if isempty(stim_time_idx{i})
        stim_time_idx{i} = [42, 82, 122, 162, 202]'*1e3;
        stim_end_idx{i} = [62, 102, 142, 182, 222]'*1e3 + 1;
    end
    stim_end_idx{i} = [1;stim_end_idx{i}];
    
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
end

stim_spike_time = cell(n_recording, length(stim_time{1}));
stim_FR = nan(n_recording, length(stim_time{1}));
no_stim_spike_time = cell(n_recording, length(stim_time{1})+1);
no_stim_FR = nan(n_recording, length(stim_time{1})+1);

for i = 1:n_recording
    for n = 1:length(stim_time{i})
        stim_power(i,n) = mean(power(floor(stim_time{i}(n)/dt):floor(stim_end{i}(n+1)/dt),i));
        stim_spike_time{i,n} = spike_time2{i}(find((spike_time2{i}>=stim_time{i}(n))&(spike_time2{i}<stim_end{i}(n+1))));
        stim_FR(i,n) = length(stim_spike_time{i,n})/(stim_end{i}(n+1)-stim_time{i}(n));
        no_stim_spike_time{i,n} = spike_time2{i}(find((spike_time2{i}>=stim_end{i}(n))&(spike_time2{i}<stim_time{i}(n))));
        no_stim_FR(i,n) = length(no_stim_spike_time{i,n})/(stim_time{i}(n)-stim_end{i}(n));
    end
    no_stim_spike_time{i,length(stim_time{i})+1} = spike_time2{i}(find(spike_time2{i}>=stim_end{i}(length(stim_time{i})+1)));
    no_stim_FR(i,length(stim_time{i})+1) = length(no_stim_spike_time{i,length(stim_time{i})+1})/(out.T(end)-stim_end{i}(length(stim_time{i})+1));
end



if if_plot
    figure
    for i = 1:length(stim_time{1})
        plot(out.T-stim_time{1}(i),Vm2(:,1),'k')
        hold on
    end
    xlim([-1,2])
    y1 = ylim();
    yplot = y1(1):0.01:y1(2);

    plot(0*ones(size(yplot)),yplot,'Color',[0.5,0.5,0.5])
    plot((stim_end{1}(1+1)-stim_time{1}(1))*ones(size(yplot)),yplot,'Color',[0.5,0.5,0.5])

%     addpath(genpath(fullfile(pwd, 'plotSpikeRaster_v1.2')));
%     figure,[xPoints, yPoints] = plotSpikeRaster(spike_time2,'PlotType','vertline');
%     y1 = ylim();
%     yplot = y1(1):0.01:y1(2);
%     hold on
%     for i = 1:length(stim_time{1})
%         plot(stim_time{1}(i)*ones(size(yplot)),yplot,'Color',[0.5,0.5,0.5])
%         plot(stim_end{1}(i+1)*ones(size(yplot)),yplot,'Color',[0.5,0.5,0.5])
%     end
    xlim
end


dur_stim = stim_FR;
no_stim = no_stim_FR;
pre_stim_FR = no_stim_FR(1);
post_stim_FR = nanmean(no_stim_FR(2:end));
stim_FR = nanmean(stim_FR);

if (length(find(power<=5))<3)&(length(find(power<=3))>0)
    idx_rmv = find(power<=3);
    dur_stim(idx_rmv) = [];
    no_stim(idx_rmv) = [];
    temp1 = no_stim_FR(2:end);
    temp1(idx_rmv) = [];
    post_stim_FR = nanmean(temp1);
    stim_FR = nanmean(dur_stim);
end

if if_save
    save(fullfile(save_path,sprintf('%s_cell%s_power_long_%d.mat', date, string(cell_idx), idx_f)),'idx_f','Vm2','spike_time2', 'stim_time','stim_end','pre_stim_FR','stim_FR','post_stim_FR','dur_stim','no_stim','stim_spike_time','no_stim_spike_time')

end
fprintf('finished cell%s, from date %s, file %s\n',string(cell_idx), date, filename)
end