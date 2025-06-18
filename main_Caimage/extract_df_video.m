function dF = extract_df_video(F, manual_flag, win)
% extract dF/F from the raw fluorescence from a video
% F: raw fluorescence trace in rows
% manual_flag: 1 for defining the baseline window manually, 0 for useing
% deconvolution method to extract the baseline and denoise
% win: datapoint index of the baseline. 
if nargin < 2, manual_flag = 0;end
if nargin < 3, win = [1, 10]; end
if manual_flag == 0
    [F_deconvolved, s, options] = deconvolveCa(double(F), 'foopsi', 'ar1', 'optimize_pars', true, 'optimize_b', true);
    b = options.b;
    dF = (F-b)/b;
else
    b = mean(F(:, win), 2);
    for i = 1:length(b), dF(i,:) = (F(i,:)-b(i))/b(i); end
end