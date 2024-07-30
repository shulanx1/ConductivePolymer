function raster_plot(spike_time, color)
% spike_time: cells of time of spikes, each cell is a trial
addpath(genpath(fullfile(pwd, 'plotSpikeRaster_v1.2')))

if nargin < 2, color = [0,0,0]; end
LineFormat = struct();
LineFormat.Color = color;
for i = 1:length(spike_time)
    if isempty(spike_time{i})
        spike_time{i} = 1000; % add a fake spike so no error will be evoked in the plot function
    end
end
plotSpikeRaster(spike_time,'PlotType','vertline', 'LineFormat', LineFormat);
xlim([-0.01,0.1])
end