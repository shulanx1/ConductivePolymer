function [R_down,t_down,idx_down] = calc_R(Vm, Im, t, dt)
if nargin <3, dt = t(2)-t(1); end


idx_rising = find(diff(Vm)>1);
idx_falling = find(diff(Vm)<-1);
idx_rising(find(diff(idx_rising)<=1)+1) = [];
idx_falling(find(diff(idx_falling)<=1)+1) = [];

if isempty(idx_rising)
%     idx_rising = find(diff(Im)>15);
%     idx_falling = find(diff(Im)<-15);
%     idx_rising(find(diff(idx_rising)<=1)+1) = [];
%     idx_falling(find(diff(idx_falling)<=1)+1) = [];
    idx_rising = [204:400:length(Vm)];
    idx_falling = [4:400:length(Vm)];
end
if length(idx_falling)<length(idx_rising)
   idx_rising =  idx_rising(1:length(idx_falling));
elseif    length(idx_falling)>length(idx_rising)
   idx_falling =  idx_falling(1:length(idx_rising));
end

% isi = idx_falling(2:end)-idx_rising(1:end-1);
% idx_rmv = find(isi<quantile(isi,0.1)); % remove artifacts at the end of each trials
% idx_rising(idx_rmv) = [];
% idx_falling(idx_rmv+1) = [];


I_fall = zeros(1, length(idx_falling)-1);
I_rise = zeros(1, length(idx_falling)-1);
idx_down = idx_rising(1:end-1);
for i = 1:length(idx_falling)-1
    I_fall(i) = mean(Im(idx_falling(i)+20:idx_rising(i)-20));
    I_rise(i) = mean(Im(idx_rising(i)+20:idx_falling(i+1)-20));
end

if isempty(find(~isnan(I_fall)))||isempty(find(~isnan(I_rise)))
    I_fall = zeros(1, length(idx_falling)-1);
    I_rise = zeros(1, length(idx_falling)-1);
    idx_down = idx_falling(1:end-1);
    for i = 1:length(idx_falling)-1
        I_fall(i) = mean(Im(idx_falling(i)+20:idx_rising(i+1)-20));
        I_rise(i) = mean(Im(idx_rising(i)+20:idx_falling(i)-20));
    end
end

R_down = 5000./(I_rise-I_fall);
% remove scanning artifact
idx_rmv = find(abs(diff(R_down))>50)+1;
for i = 1:length(idx_rmv)
    idx_smooth = max(1,idx_rmv(i)-10):min(idx_rmv(i)+10, length(R_down));
    idx_smooth = setdiff(idx_smooth,idx_rmv);
    R_down(idx_rmv(i)) = mean(R_down(idx_smooth));
end

R_baseline = mean(R_down(1:10));
idx_rmv = find(R_down<R_baseline*0.3);
for i = 1:length(idx_rmv)
    idx_smooth = max(1,idx_rmv(i)-10):min(idx_rmv(i)+10, length(R_down));
    idx_smooth = setdiff(idx_smooth,idx_rmv);
    R_down(idx_rmv(i)) = mean(R_down(idx_smooth));
end
t_down = t(idx_down);
end