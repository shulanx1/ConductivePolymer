function [df, raw_f, A, Cn]= extract_raw_trace_video(img, thre, sen)
Cn = std(img, [], 3);
Cn1 = imgaussfilt(Cn,1.2);
FOV = [size(Cn, 1), size(Cn, 2)];
if nargin < 3, sen = 0.1; end
if nargin < 2, thre = 0.5; end
h = figure(1); imshow(stretch(Cn));
hold on
title(sprintf('ROIs'))
x = 1:FOV(1);
y = 1:FOV(2);
BW = imbinarize(stretch(Cn1),"adaptive","Sensitivity",sen);
a = bwboundaries(BW);
coor = cell(0);
count = 1;
for i = 1:length(a)
    if size(a{i}, 1)>20
        coor{count}(:,1) = a{i}(:,2);
        coor{count}(:,2) = a{i}(:,1);
        count = count + 1;
    end
end
raw_f = zeros(length(coor), size(img, 3));
A = sparse(length(coor), FOV(1)*FOV(2));
for i = 1:length(coor)
    A(i,:) = spatial_contour(coor{i}, FOV, Cn, 0.2);
    try
        raw_f(i,:) = ROI_fluorescence(img, A(i,:));
        baseline = mean(raw_f(i,1:10), 2);
        df(i,:) = (raw_f(i,:)-baseline)/baseline;
    catch
        continue
    end
end
% df = extract_df_video(raw_f);


idx = [];
for i = 1:size(df, 1)
    if max(df(i,:)) > thre
        idx = [idx, i];
        plot(coor{i}(:,1), coor{i}(:,2), 'm', 'Linewidth', 2)
    end
end
df = df(idx,:);
raw_f = raw_f(idx,:);
A = A(idx, :);

end