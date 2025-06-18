
if (strcmp(recording_type,'bAP_incremental'))||(strcmp(recording_type,'bAP'))||(strcmp(recording_type,'noopto_somastim'))
    stim_idx = find(diff(Im2(:,1))>300);
    stim_idx(find(diff(stim_idx)<=15e-3/dt)+1) = [];
elseif strcmp(recording_type, 'paired_stim_Ca')
    stim_idx = find(diff(Im2)>300);
    stim_idx(find(diff(stim_idx)<=15e-3/dt)) = [];
    stim_idx = [stim_idx; 93501];  % add the apical stim
elseif strcmp(recording_type, 'noopto_dendstim')
    stim_idx = find(diff(Im1)>300);
    stim_idx(find(diff(stim_idx)<=15e-3/dt)) = [];
elseif strcmp(recording_type, 'uEPSP_ca')
    stim_idx = find(diff(out.S{7})>300);
    stim_idx(find(diff(stim_idx)<=15e-3/dt)) = [];
end


stim_time = t(stim_idx);
stim_idx_Ca = zeros(size(stim_time));
Ca_amp = zeros(size(filt_df, 1), length(stim_time));
for i = 1:length(stim_time)
    [~, stim_idx_Ca(i)] = min(abs(t_Ca-stim_time(i)));
    if max(stim_idx_Ca(i)- floor(0.1/dt_Ca) + floor(0.8/dt_Ca))>size(filt_df, 2)
        break
    end
    Ca_amp(:,i) = max(filt_df(:,stim_idx_Ca(i):stim_idx_Ca(i)- floor(0.1/dt_Ca) + floor(0.8/dt_Ca)), [], 2);
end

for i = 1:size(Ca_amp, 1)
    mean(Ca_amp(i, Ca_amp(i,:)~=0))
end
save(fullfile(output_folder, sprintf('data_%s.mat', base_name)), 'raw_f', 'raw_df', 'filt_df', 't_Ca', 'Vm1', 'Im1', 'Vm2', 'Im2', 't', 'Ca_amp', 'stim_idx', 'stim_idx_Ca');
