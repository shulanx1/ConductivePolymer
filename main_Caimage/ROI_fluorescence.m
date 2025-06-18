function F = ROI_fluorescence(img, a)
% find out the average/median of fluorescence trace inside the ROI

% F: mean or median fluorescenc trace in the ROI, column vector

% data: MatFile with downsampled fluorescence dynamic
% a: spatial matrix, rearranged in column vector
% method: 'mean' or 'median', defalut is median

if nargin < 3, method = 'median'; end

data = reshape(img, [], size(img, 3));

F_raw = nonzeros(a*data)/sum(a);
F = sgolayfilt(F_raw,2,5);
end
