function [amp, halfwidth, risetime, maxdv]= single_spine_serise(base_path, date, cell_idx, branch_no, idx_single, if_plot, if_save, save_path)
if nargin < 6, if_plot = 0; end
if nargin < 7, if_save = 0; end



preflix = strcat(strcat(date(5:6), date(1:4)),'_001');
path = fullfile(base_path, date, strcat('cell', string(cell_idx)));
if nargin < 8, save_path = path; end

filename = sprintf('%s.single_spine_series.%d.wcp', preflix,idx_single);
try
    out=import_wcp(fullfile(path, filename),'debug');
catch
    date_plus_one = next_date_string(date);
    preflix = strcat(strcat(date_plus_one(5:6),date_plus_one(1:4)),'_001');
    filename = sprintf('%s.single_spine_series.%d.wcp', preflix,idx_single);
    out=import_wcp(fullfile(path, filename),'debug');
end

n_channel = out.channel_no;
n_recording = length(out.rec_index);

stim_start = cell(1,n_recording);
for i = 1
    Vm(:,i) = out.S{3}(:,i);
    PC(:,i) = out.S{7}(:,i);
    stim_start{i} = find(diff(PC(:,i))>90);
    stim_start{i} = stim_start{i}([1; find(diff(stim_start{i})>40)+1]);
    if size(stim_start{i}, 1) == 0
        continue
    end
    num_spines = length(stim_start{i});
end
T = out.T;
dt = T(2) - T(1);
t = (-499*dt:dt:4501*dt)*1e3;

EPSP = [];

% Fs = 1/dt;
% Fc = [0 200];
% Wn = Fc./(Fs/2);
% %b = fir1(10,Wn,'bandpass');
% [b1,a1] = butter(6, Wn(2), 'low');


if length(idx_single) == 1
    for i = 1:n_recording
        Vm(:,i) = out.S{3}(:,i);
        PC(:,i) = out.S{7}(:,i);
        shutter(:,i) = out.S{8}(:,i);
        stim_start{i} = find(diff(PC(:,i))>50);
        if size(stim_start{i}, 1) == 0
            continue
        end
        stim_start{i} = stim_start{i}([1; find(diff(stim_start{i})>50)+1]);
        for j = 1:length(stim_start{i})
            if shutter(stim_start{i}(j), i)< 4
                continue
            else
                stim_start_new = 500;
                EPSP(:,j) = Vm(stim_start{i}(j)-stim_start_new+1:stim_start{i}(j)+4501,i)-Vm(stim_start{i}(j), i);
                %EPSP = [EPSP, Vm(stim_start{i}(j)-500:stim_start{i}(j)+3500,i)];
                EPSP_temp = sgolayfilt(EPSP(:,j),2,61)';
                EPSP_temp = EPSP_temp - EPSP_temp(stim_start_new);
                [amp(j),~] = max(abs(EPSP(stim_start_new:stim_start_new+500,j)-EPSP(stim_start_new,j)));
                [amp_temp,max_amp(j)] = max(abs(EPSP_temp(stim_start_new:stim_start_new+500)));
                max_amp(j) = max_amp(j)+stim_start_new;
                [half,I] = sort(abs(EPSP_temp-amp_temp/2));
                b = find(I<max_amp(j) & I>stim_start_new);
                a = find(I>max_amp(j) & I<=stim_start_new+1000);
                halfwidth(j) = (I(a(1))-I(b(1)))*dt*1000; %in ms
                risetime(j) = (max_amp(j)-stim_start_new)*dt*1000; %in ms
                maxdv(j) = max(diff(movmean(EPSP(:,j), 10)))/dt/1000; %in V/s
                
            end
        end
    end
else
    for m = 1:length(idx_single)
        filename_single = sprintf('%s.single_spine_series.%d.wcp',preflix, idx_single(m));
        out=import_wcp(fullfile(path, filename_single),'debug');
        for i = 1:n_recording
            Vm(:,i) = out.S{3}(:,i);
            PC(:,i) = out.S{7}(:,i);
            shutter(:,i) = out.S{8}(:,i);
            stim_start{i} = find(diff(PC(:,i))>50);
            if size(stim_start{i}, 1) == 0
                continue
            end
            stim_start{i} = stim_start{i}([1; find(diff(stim_start{i})>50)+1]);
            for j = spine{m}%1:length(stim_start{i})
                if shutter(stim_start{i}(j), i)< 4
                    continue
                else
                    stim_start_new = 500;
                    EPSP(:,j) = Vm(stim_start{i}(j)-stim_start_new+1:stim_start{i}(j)+4001,i)-Vm(stim_start{i}(j), i);
                    %EPSP = [EPSP, Vm(stim_start{i}(j)-500:stim_start{i}(j)+3500,i)];
                    EPSP_temp = filtfilt(b1, a1, EPSP(:,j));
                    [amp(j),~] = max(abs(EPSP(stim_start_new:stim_start_new+500,j)-EPSP(stim_start_new,j)));
                    [amp_temp,max_amp(j)] = max(abs(EPSP_temp(stim_start_new:stim_start_new+500)-EPSP_temp(stim_start_new)));
                    max_amp(j) = max_amp(j)+stim_start_new;
                    [half,I] = sort(abs(EPSP(:,j)-amp_temp/2));
                    b = find(I<max_amp(j) & I>stim_start_new);
                    a = find(I>max_amp(j));
                    halfwidth(j) = (I(a(1))-I(b(1)))*dt*1000; %in ms
                    risetime(j) = (max_amp(j)-stim_start_new)*dt*1000; %in ms
                    maxdv(j) = max(diff(movmean(EPSP(:,j), 10)))/dt/1000; %in V/s
                end
            end
        end
        
    end
end
if if_plot
    figure
    for j = 1:size(EPSP,2)
        subplot(2,4,j)
        plot(t, EPSP(:,j), 'Color', [0.5,0.5,0.5])
        xlim([-20,150])
        %             ylim([-0.5,1.5])
        hold on
        scatter(t([I(a(1)),I(b(1)),max_amp(j)]), EPSP([I(a(1)),I(b(1)),max_amp(j)],j), 'O', 'MarkerEdgeColor', [0, 0.4470, 0.7410])
        xlabel('t (ms)')
        ylabel('Vm (mV)')
        %             xlim([min(t), max(t)])
        box off
    end
end

if if_save
    save(fullfile(save_path, sprintf('%s_cell%s_branch%s_single_spine_series_%d.mat', date, string(cell_idx),  string(branch_no), idx_single)),'EPSP', 'Vm', 'amp', 'risetime','halfwidth', 'maxdv','idx_single', 't', 'dt')
end

fprintf('finished branch %s, cell%s, from date %s, file %s\n',string(branch_no), string(cell_idx), date, filename)
end