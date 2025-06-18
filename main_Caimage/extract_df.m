function [df,F0] = extract_df(f, dt, win)
% compute df/f from the raw fluorescence trace
% f: raw fluorescence
% dt: line rate (s)
% win: time window of baseline fluoresence (s), default to be the first 10 ms 
if nargin < 3, win = [0, 0.01]; end
if iscell(f)
    df = zeros(size(f{1}));
    for i = 1:size(f{1}, 1)
        baseline = mean(f{1}(i, ceil(win(1)/dt)+1:ceil(win(2)/dt)));
        baseline_r = mean(f{2}(i, ceil(win(1)/dt)+1:ceil(win(2)/dt)));
        df(i,:) = (f{1}(i,:)-baseline)/baseline_r;
    end
else
    df = zeros(size(f));
    for i = 1:size(f, 1)
        baseline = mean(f(i, ceil(win(1)/dt)+1:ceil(win(2)/dt)));
        df(i,:) = (f(i,:)-baseline)/baseline;
        F0(i) = baseline;
    end
end
end