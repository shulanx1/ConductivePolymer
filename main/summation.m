function [amp_sum, halfwidth_sum, risetime_sum, maxdv_sum]= summation(base_path, date, cell_idx, branch_no, idx_sum, if_plot, if_save, save_path)
if nargin < 6, if_plot = 0; end
if nargin < 7, if_save = 0; end



preflix = strcat(strcat(date(5:6), date(1:4)),'_001');
path = fullfile(base_path, date, strcat('cell', string(cell_idx)));
if nargin < 8, save_path = path; end

filename = sprintf('%s.summation.%d.wcp', preflix,idx_sum);
try
    out=import_wcp(fullfile(path, filename),'debug');
catch
    date_plus_one = next_date_string(date);
    preflix = strcat(strcat(date_plus_one(5:6),date_plus_one(1:4)),'_001');
    filename = sprintf('%s.summation.%d.wcp', preflix,idx_sum);
    out=import_wcp(fullfile(path, filename),'debug');
end

n_channel = out.channel_no;
n_recording = length(out.rec_index);

stim_start = cell(1,n_recording);
T = out.T;
dt = T(2) - T(1);
t = (-499*dt:dt:4501*dt)*1e3;

EPSP_sum = [];

for i = 1:n_recording
    Vm_sum(:,i) = out.S{3}(:,i);
    Vm1_sum(:,i) = out.S{1}(:,i);
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
            EPSP_sum(:,j) = Vm_sum(stim_start{i}(j)-stim_start_new+1:stim_start{i}(j)+4501,i)-Vm_sum(stim_start{i}(j), i);
%             EPSP_sum1(:,j) = Vm1_sum(stim_start{i}(j)-stim_start_new+1:stim_start{i}(j)+5501,i)-Vm1_sum(stim_start{i}(j), i);
            %EPSP = [EPSP, Vm(stim_start{i}(j)-500:stim_start{i}(j)+3500,i)];
            [amp_sum(j),max_amp_sum(j)] = max(abs(EPSP_sum(stim_start_new:stim_start_new+1200,j)-EPSP_sum(stim_start_new,j)));
            max_amp_sum(j) = max_amp_sum(j)+stim_start_new;
            [half,I] = sort(abs(EPSP_sum(:,j)-amp_sum(j)/2));
            b = find(I<max_amp_sum(j) & I>stim_start_new);
            a = find(I>max_amp_sum(j));
            halfwidth_sum(j) = (I(a(1))-I(b(1)))*dt*1000; %in ms
            risetime_sum(j) = (max_amp_sum(j)-stim_start_new)*dt*1000; %in ms
            maxdv_sum(j) = max(diff(movmean(EPSP_sum(:,j), 10)))/dt/1000; %in V/s
            [~,c] = max(diff(movmean(EPSP_sum(:,j), 10)));

        end
    end
end
if if_plot
    figure
    for j = 1:size(EPSP_sum,2)
        subplot(2,4,j)
        plot(t, EPSP_sum(:,j), 'Color', [0.5,0.5,0.5])
        hold on
        scatter(t([I(a(1)),I(b(1)),max_amp_sum(j)]), EPSP_sum([I(a(1)),I(b(1)),max_amp_sum(j)],j), 'O', 'MarkerEdgeColor', [0, 0.4470, 0.7410])
        scatter(t(c), EPSP_sum(c, j))
        xlabel('t (ms)')
        ylabel('Vm (mV)')
        xlim([min(t), max(t)])
        box off
    end
end

if if_save

    save(fullfile(save_path, sprintf('%s_cell%s_branch%s_sum_galvo_%d.mat',  date, string(cell_idx),  string(branch_no), idx_sum)),'EPSP_sum', 'Vm_sum', 'amp_sum', 'risetime_sum','halfwidth_sum', 'maxdv_sum', 'idx_sum', 't','dt')
end
fprintf('finished branch %s, cell%s, from date %s, file %s\n',string(branch_no), string(cell_idx), date, filename)
end