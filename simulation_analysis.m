
% load('E:\Code\simpl5pn_inhibition\results\shunt\shunt_051225_144138.mat')
% load('E:\Code\simpl5pn_inhibition\results\shunt_w_ionic\shunt_051225_171930.mat')
load('E:\Code\simpl5pn_inhibition\results\ionic\shunt_051225_221738.mat')
a = find(loc==0);
b = find(loc==2);
c = find(loc==3);

amp = S(a);

va = reshape(v(a,1,:),length(a),[]);
vb = reshape(v(b,1,:),length(b),[]);
vc = reshape(v(c,1,:),length(c),[]);

t = 0:dt:double(T);
stim_idx = t(find(diff(g)>1));
stimend_idx = t(find(diff(g)<-1));
if isempty(stim_idx)
    stim_idx = [199.900000000000,599.900000000000,999.900000000000,1399.90000000000];
    stimend_idx = [399.900000000000,799.900000000000,1199.90000000000];
end
t_spike = cell(3, length(amp));
t_spike(1,:) = spike_times(dt, va);
t_spike(2,:) = spike_times(dt, vb);
t_spike(3,:) = spike_times(dt, vc);

mod = zeros(3, length(amp));
for i = 1:length(amp)
    for j = 1:3
        spike_stim = find(((t_spike{j,i}>=stim_idx(1))&(t_spike{j,i}<stimend_idx(1)))|((t_spike{j,i}>=stim_idx(2))&(t_spike{j,i}<stimend_idx(2)))|((t_spike{j,i}>=stim_idx(3))&(t_spike{j,i}<stimend_idx(3)))|((t_spike{j,i}>=stim_idx(4))));
        spike_nostim = find(((t_spike{j,i}<stim_idx(1)))|((t_spike{j,i}>=stimend_idx(1))&(t_spike{j,i}<stim_idx(2)))|((t_spike{j,i}>=stimend_idx(2))&(t_spike{j,i}<stim_idx(3)))|((t_spike{j,i}>=stimend_idx(3))&(t_spike{j,i}<stim_idx(4))));
        mod(j,i) = (length(spike_stim)/0.3)/(length(spike_nostim)/0.4);
    end
end

function t_spike = spike_times(dt, v)
    % Get spike times from voltage trace.
    %
    % Parameters
    % ----------
    % dt : float
    %     simulation timestep
    % v : matrix
    %     compartment voltages v=v(compartment, time)
    %
    % Returns
    % -------
    % t_spike : array
    %     spike times
    t_spike = cell(1, size(v,1));
    for i = 1:size(v,1)
        thresh_cross = find(v(i, :) > 0); % Find indices where voltage crosses threshold

        if ~isempty(thresh_cross)
            spikes = find(diff(thresh_cross) > 1) + 1; % Detect gaps in threshold crossings
            spikes = [1, spikes]; % Include the first spike
            spikes = thresh_cross(spikes); % Map to original indices

            spikes_temp = spikes;
            for k = 1:length(spikes) % Detect the max amplitude point as spike time
                spike = spikes(k);
                spike_w_temp = v(1, spike:min(spike + floor(2/dt), size(v, 2)));
                [~, max_idx] = max(spike_w_temp);
                spikes_temp(k) = spike + max_idx - 1;
            end

            t_spike{i} = spikes_temp * dt; % Convert indices to time
        else
            t_spike{i} = []; % Return empty array if no spikes detected
        end
    end
end

