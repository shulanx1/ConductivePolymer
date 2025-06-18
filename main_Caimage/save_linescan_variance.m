dark_var = 0.0228^2;

Ca_amp = zeros(size(filt_df, 1), size(filt_df, 3));
GR_amp = zeros(size(filt_GR, 1), size(filt_GR, 3));
time_range_Ca_var = [0, 0.1];
idx_range = -time_range_idx_Ca(1) + floor(time_range_Ca_var(1)/dt_Ca)+1:-time_range_idx_Ca(1) + floor(time_range_Ca_var(2)/dt_Ca);
Ca_amp= max(filt_df(:,idx_range,:), [], 2)- min(filt_df(:,idx_range,:), [], 2);
GR_amp= max(filt_GR(:,idx_range,:), [], 2)-min(filt_GR(:,idx_range,:), [], 2);
Ca_amp = reshape(Ca_amp ,size(filt_df, 1), size(filt_df, 3));
GR_amp = reshape(Ca_amp ,size(filt_GR, 1), size(filt_GR, 3));
std_thr_df = std(filt_df_temp, [],2);
std_thr_GR = std(filt_GR_temp, [],2);
invade_df = zeros(size(raw_df, 1), size(filt_df, 3));
invade_GR = zeros(size(raw_df, 1), size(filt_GR, 3));
for i = 1:size(filt_df, 3)
    invade_df(:,i) = Ca_amp(:,i)>=1.5*std_thr_df;
    invade_GR(:,i) = GR_amp(:,i)>=1.5*std_thr_GR;
end

Ca_mean = zeros(size(raw_df, 1), size(filt_df, 2));
Ca_var = zeros(size(raw_df, 1), size(filt_df, 2));
time_range_Ca_var = [0.1, 0.2];
idx_range = -time_range_idx_Ca(1) + floor(time_range_Ca_var(1)/dt_Ca)+1:-time_range_idx_Ca(1) + floor(time_range_Ca_var(2)/dt_Ca);
CV = zeros(size(raw_df, 1)-1, length(idx_range));
for j = 1:size(raw_df, 1)
    Ca_mean(j,:) = mean(reshape(filt_f(j,:,find(invade_df(j,:))), [], length(find(invade_df(j,:))))-F0(j,find(invade_df(j,:))), 2);
    Ca_var(j,:) = var(filt_f(j,:,find(invade_df(j,:))), [], 3);
    a = Ca_var(j,idx_range)-dark_var;
    a(a<0)=0;
    CV(j,:) = sqrt(a)./Ca_mean(j,idx_range);
end
p_invade(1) = length(find(invade_df(1,:)))/size(filt_df, 3);
p_invade(2) = length(find(invade_df(2,:)))/size(filt_df, 3);



syms N
P_0 = 1-p_invade(1);
CV_s = nanmax(CV(1,:));
eqn1 = P_0 == (CV_s^2/(1/N+CV_s^2))^N;
N_soln = vpasolve(eqn1,N);

stim_idx_Ca = stim_time_idx_Ca;
save(fullfile(output_folder, sprintf('data_variance_%s.mat', base_name)), 'raw_f', 'raw_df', 'filt_df','raw_R', 'raw_GR', 'filt_GR', 't_Ca','t_Ca_plot','t_plot','Vm', 'Vm1', 'Im1', 'Vm2', 'Im2', 't', 'Ca_amp','GR_amp','std_thr_df','std_thr_GR','invade_GR','invade_df','p_invade', 'stim_idx', 'stim_idx_Ca', 'Ca_mean', 'Ca_var', 'CV');
