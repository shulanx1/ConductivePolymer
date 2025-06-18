function bAP_idx_temp = find_bAP(spike_idx_temp, Vm1, dt, thr)
    if nargin < 3, dt = 5e-5; end
    if nargin < 4, thr = -40;end
    [pks,locs,~,~]  = findpeaks(Vm1(spike_idx_temp:spike_idx_temp+floor(0.015/dt)), 'MinPeakProminence',1, 'MinPeakDistance',floor(0.005/dt));
    bAP_idx_temp = spike_idx_temp + min(locs(pks>=thr));
    if isempty(bAP_idx_temp)
        [pks,locs,~,~]  = findpeaks(Vm1(spike_idx_temp:spike_idx_temp+floor(0.015/dt)), 'MinPeakDistance',floor(0.005/dt));
        idx_rmv = locs<floor(0.002/dt);
        locs(idx_rmv) = [];
        pks(idx_rmv) = [];
        bAP_idx_temp = spike_idx_temp + min(locs(pks>=thr));
    end
end