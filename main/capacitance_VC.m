function [Ra,Cm, Ileak] = capacitance_VC(base_path, date, cell_idx, idx_f, if_plot, if_save)
if nargin < 5, if_plot = 0; end
if nargin < 6, if_save = 0; end


preflix = strcat(strcat(date(5:6), date(1:4)),'_001');
path = fullfile(base_path, date, strcat('cell', string(cell_idx)));
filename = sprintf('%s.capacitance_VCstep.%d.wcp', preflix,idx_f);

try
    out=import_wcp(fullfile(path, filename),'debug');
catch
    try
        date_plus_one = next_date_string(date);
        preflix = strcat(strcat(date_plus_one(5:6),date_plus_one(1:4)),'_001');
        filename = sprintf('%s.capacitance_VCstep.%d.wcp', preflix,idx_f);
        out=import_wcp(fullfile(path, filename),'debug');
    catch
        try
            preflix = strcat(strcat(date(5:6), date(1:4)),'_001');
            filename = sprintf('%s.sealtest.%d.wcp', preflix,idx_f);
            out=import_wcp(fullfile(path, filename),'debug');
        catch
            date_plus_one = next_date_string(date);
            preflix = strcat(strcat(date_plus_one(5:6),date_plus_one(1:4)),'_001');
            filename = sprintf('%s.sealtest.%d.wcp', preflix,idx_f);
            out=import_wcp(fullfile(path, filename),'debug');
        end
    end
end


n_channel = out.channel_no;
if contains(filename, 'sealtest')
    n_recording = min([10, length(out.rec_index)-1]);
else
    n_recording = length(out.rec_index);
end

dt = out.T(2)-out.T(1); % sample time in s

   
Im = out.S{3};
Vm = out.S{4};
t = 0:dt:dt*(size(Vm,1)-1);

Cm = zeros(1, n_recording);
Ra = zeros(1, n_recording);
Ileak = zeros(1, n_recording);

for i = 1:n_recording
    pulse_start = find(diff(Vm(:,i))<-1);
    pulse_start = pulse_start(1);
    I_baseline = mean(Im(1:pulse_start,i));
    Ileak(i) = I_baseline;
    Im(:,i) = Im(:,i) - Ileak(i);
end

Im_mean = mean(Im,2);
I_step = mean(Im_mean(pulse_start+1000:pulse_start+6000));
I_trans = min(Im_mean(1:pulse_start+1000));
[~,start_point] = min(Im_mean(1:pulse_start+1000));
I_temp = (Im_mean(start_point:start_point + 500))-I_step;
t = 0:dt:dt*(size(I_temp,1)-1);
I_temp(I_temp>0) = 0;
Cm = -trapz(t, I_temp)*1e3/10; %pF
Ra = -10/(I_trans/1e3);  %in Mohm
Ileak = mean(Ileak);

if if_plot
    figure
    plot(0:dt:dt*(size(Im,1)-1), Im_mean)
end

if if_save
    save(fullfile(path,sprintf('capacitance_VCstep_%d.mat', idx_f)),'Cm','Ra', 'Im', 'Vm','Ileak','Im_mean');
end

fprintf('finished cell%s, from date %s, file %s\n',string(cell_idx), date, filename)
end