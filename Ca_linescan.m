if exist('r_fft', 'var')
    clearvars -except r_fft
else    clear all
end
close all

base_folder = 'E:\data\polymer\042125\cell2';

output_folder = fullfile(base_folder, 'Ca_output');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
addpath([pwd,'\main_Caimage']);
basefolder = pwd;

image_foldername = uigetdir(base_folder, 'Pick a image file');
% base_name = regexp(image_foldername,'\','split');
% base_name = base_name{end};
[filename, ~] = uigetfile({fullfile(base_folder, '*.wcp')}, 'Pick a ephys file');
base_name_temp = strsplit(filename,'.');
recording_idx = str2num(base_name_temp{3});
recording_type = base_name_temp{2};
base_name = sprintf('%s_%d', recording_type, recording_idx);


Ch2_dir = dir(fullfile(image_foldername, '*_Ch2_000001.ome.tif'));
Ch1_dir = dir(fullfile(image_foldername, '*_Ch1_000001.ome.tif'));
filename_meta = dir(fullfile(image_foldername, '*.xml'));

[img_Ch1,~] = imread(fullfile(image_foldername, Ch1_dir(1).name)); % load image
[img,~] = imread(fullfile(image_foldername, Ch2_dir(1).name)); % load image
if length(Ch2_dir) > 1 % stitch the files
    for i = 2:length(Ch2_dir)
        [img_Ch1_temp,~] = imread(fullfile(image_foldername, Ch1_dir(i).name)); % load image
        [img_temp,~] = imread(fullfile(image_foldername, Ch2_dir(i).name)); % load image
        img_Ch1 = [img_Ch1;img_Ch1_temp];
        img = [img; img_temp];
        clear img_Ch1_temp img_temp
    end
end
%%
meta_data = parseXML(fullfile(image_foldername, filename_meta(1).name));  % load metadata
img = double(img); % reformat data from uint to double for FFT

if ~exist('r_fft', 'var')
    [denoise_img,  r_fft] = custom_fft(img, 1);
    [denoise_Ch1, ~] = custom_fft(img_Ch1,  0, r_fft);
else
    [denoise_img, ~] = custom_fft(img, 0, r_fft);
    [denoise_Ch1, ~] = custom_fft(img_Ch1, 0, r_fft);
end
denoise_img = imrotate(denoise_img,90); % rotate image in 90 degree 
img = imrotate(img,90); % rotate image in 90 degree 
img_Ch1 = imrotate(denoise_Ch1,90); % rotate image in 90 degree 

%dwell_time = str2num(meta_data.Children(4).Children(10).Attributes(2).Value)*1e-6;  % read dwell time from metadata (s)
dwell_time = str2num(meta_data.Children(4).Children(12).Attributes(2).Value)*1e-6;  % read dwell time from metadata (s)
% dt_Ca = (dwell_time+2.81e-7)*size(img, 1);   % line rate
dt_Ca = (dwell_time+5e-7)*size(img, 1);   % line rate
t_Ca = dt_Ca*[0:(size(img, 2)-1)]; 

% extract dF/F
raw_f = extract_raw_trace(denoise_img, base_name, output_folder, img_Ch1);
raw_df = extract_df(raw_f, dt_Ca,  [0.05, 0.45]);
filt_df= sgolayfilt(raw_df',2,61)';


%% time align with ephys
h = figure(); 
ax1 = subplot(211);
plot(t_Ca, raw_df, 'Color', [0.5,0.5,0.5])
hold on, plot(t_Ca, filt_df, 'LineWidth', 2)
xlabel('T [s]')
ylabel('dF/F [A.U.]')


out=import_wcp(fullfile(base_folder, filename),'debug');

n_channel = out.channel_no;
n_recording = length(out.rec_index);

dt = out.T(2)-out.T(1); % sample time in s
t = out.T;

Im1 = out.S{2};
Vm1 = out.S{1};   
Im2 = out.S{4};
Vm2 = out.S{3};
thr2 = 0;
spike_time2 = {};
for n = 2:max(size(Vm2))
    for i = 1:min(size(Vm2))
        spike_time2{i} = [];
        if Vm2(n, i)>thr2 && Vm2(n-1, i)<=thr2
            [~, peak] = max(Vm2(n-10:n+10,i));
            spike_time2{i} = [spike_time2{i}, (n-9+peak)*dt]; % spike time in s, aligned to peak
        end
    end
end
ax2 = subplot(212);
plot(out.T, Vm2(:,1), 'k', 'LineWidth', 1)
linkprop([ax1,ax2], 'xlim')
savefig(h, fullfile(output_folder, sprintf('%s_traces.fig', base_name)));
%% save data
save_linescan;



