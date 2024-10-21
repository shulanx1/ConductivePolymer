function [pre, stim, post, Vstep] = stim_VC(base_path, date, cell_idx, idx_f, if_plot, if_save, save_path)
if nargin < 5, if_plot = 0; end
if nargin < 6, if_save = 0; end

preflix = strcat(strcat(date(5:6), date(1:4)),'_001');
path = fullfile(base_path, date, strcat('cell', string(cell_idx)));
if nargin < 7, save_path = path; end

filename = sprintf('%s.power_VC.%d.wcp', preflix,idx_f);
try
    out=import_wcp(fullfile(path, filename),'debug');
catch
    date_plus_one = next_date_string(date);
    preflix = strcat(strcat(date_plus_one(5:6),date_plus_one(1:4)),'_001');
    filename = sprintf('%s.power_VC.%d.wcp', preflix,idx_f);
    out=import_wcp(fullfile(path, filename),'debug');
end

n_channel = out.channel_no;
n_recording = length(out.rec_index);

dt = out.T(2)-out.T(1); % sample time in s
t = (-500*dt:dt:2500*dt)*1e3;


stim_time_idx = {};
stim_end_idx = {};
stim_time = {};
stim_end = {};
Vstep = zeros(1, length(n_recording));

for i = 1:n_recording
    Im2(:,i) = out.S{3}(:,i);
    Vm2(:,i) = out.S{4}(:,i);
    power(:,i) = out.S{7}(:,i);
    stim_time_idx{i} = find(diff(power(:,i))>15);
    stim_time_idx{i}(find(diff(stim_time_idx{i})<=1)+1) = [];
    stim_end_idx{i} = find(diff(power(:,i))<-15);
    stim_end_idx{i}(find(diff(stim_time_idx{i})<=1)+1) = [];
    Vstep(i) = mean(Vm2(:,i));
    stim_time{i} = out.T(stim_time_idx{i});
    stim_end{i} = out.T(stim_end_idx{i});
end
pre = nan(n_recording, length(stim_time{1}));
stim = nan(n_recording, length(stim_time{1}));
post = nan(n_recording, length(stim_time{1}));
for i = 1:n_recording
    for n = 1:length(stim_time{i})
        pre(i, n) = mean(Im2(stim_time_idx{i}(n)-10:stim_time_idx{i}(n),i));
        [~,I] = max(abs(Im2(stim_end_idx{i}(n)-100:stim_end_idx{i}(n)+200,i)-pre(i,n)));
        stim(i, n) = mean(Im2(stim_end_idx{i}(n)-100+I-11:stim_end_idx{i}(n)-100+I-1,i));
        post(i, n) = mean(Im2(stim_end_idx{i}(n)+2000:stim_end_idx{i}(n)+2000+10,i));
    end
end

if if_plot
    len = n_recording;
    red = [0,0,0];
    pink = [0.8,0.8,0.8];%[255, 192, 203]/255;
    colors_p = [linspace(pink(1),red(1),len)', linspace(pink(2),red(2),len)', linspace(pink(3),red(3),len)'];

    figure
    for i = 1:n_recording
     plot(out.T, Im2(:,i),'Color',colors_p(i,:))
     hold on
    end
    xlim([stim_time{1}(end)-0.2,stim_end{1}(end)+0.2])
    hold on
    y1 = ylim();
    yplot = y1(1):y1(2);
    plot(stim_end{1}(end)*ones(size(yplot)),yplot, 'r')
    plot(stim_time{1}(end)*ones(size(yplot)),yplot, 'r')
    
    figure, scatter(Vstep,stim(:,end)-pre(:,end))
end

if if_save
    save(fullfile(save_path,sprintf('%s_cell%s_power_VC_%d.mat', date, string(cell_idx), idx_f)),'idx_f','Vm2','Vstep', 'stim_time','stim_end','stim_time_idx','stim_end_idx','pre','stim','post')

end


fprintf('finished cell%s, from date %s, file %s\n',string(cell_idx), date, filename)
end