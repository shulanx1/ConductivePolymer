function raw_f = extract_raw_trace(img, base_name, output_folder, img_Ch1)

if nargin < 2, base_name = []; end
if nargin < 3 output_folder = []; end
if nargin < 4, img_Ch1 = []; end

h = figure(); imshow(stretch(img(:,1:min(size(img, 1)*20, size(img, 2)), 1)));
hold on
count = 0;
pos = {};
while (ishandle(h))
    p = drawrectangle('Color','m');
    try
        ydim = p.Position(4);
    catch 
        break
    end
    x1 = 1;
    x2 = size(img, 2);
    count = count + 1;
    pos{count} = [x1, p.Position(2), x2, p.Position(4)];
    rectangle('Position',[1, p.Position(2), min(size(img, 1)*20, size(img, 2)), p.Position(4)])
    if (~isempty(output_folder)) && (~isempty(base_name))
    savefig(h, fullfile(output_folder, [base_name, '_ROI.fig']))
    end
end
fprintf('%d ROIs are extracted\n', count);

raw_f = zeros(count, size(img, 2), size(img, 3));
for n = 1:size(img, 3)
    for j = 1:count
        x1 = pos{j}(1);
        y1 = pos{j}(2);
        x2 = pos{j}(3);
        ydim = pos{j}(4);
        raw_f_temp = zeros(ceil(ydim), size(img, 2));
        for i = 1:ceil(ydim)    
            x = [x1, x2]; 
            y = [y1+i-1, y1+i-1];                
            raw_f_temp(i, :) = improfile(img(:,:,n),x,y);    
        end
        raw_f(j, :, n) = mean(raw_f_temp, 1);
    end
end

if ~isempty(img_Ch1)
    raw_f = {raw_f};
    raw_f{2} = zeros(count, size(img_Ch1, 2), size(img_Ch1, 3));
    for n = 1:size(img, 3)
        for j = 1:count
            x1 = pos{j}(1);
            y1 = pos{j}(2);
            x2 = pos{j}(3);
            ydim = pos{j}(4);
            raw_f_temp = zeros(ceil(ydim), size(img, 2));
            for i = 1:ceil(ydim)
                x = [x1, x2];
                y = [y1+i-1, y1+i-1];
                raw_f_temp(i, :) = improfile(img_Ch1(:,:,n),x,y);
            end
            raw_f{2}(j, :, n) = mean(raw_f_temp, 1);
        end
    end
end

end