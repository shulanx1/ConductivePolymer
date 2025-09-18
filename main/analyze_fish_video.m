function data = analyze_fish_video(data, if_plot)
    if nargin < 2, if_plot = 0; end
    Fs = 750;
    coor = data.raw_coor;
    head_col = [1:2];
    spine_col = [3:4];
    tail_col = [5:6];
    cols = [head_col;spine_col;tail_col];
    
    ini_vect = [coor(1,cols(1,:))-coor(1,cols(2,:))]';
    theta = acos(dot(ini_vect, [0;1]) / (norm(ini_vect) * norm([0;1])));  
    if ini_vect(1)<0
        theta = 2*pi-theta;
    end
    data.M_rotate = [[cos(theta), -sin(theta)]; [sin(theta), cos(theta)]];
    
    data.rotated_coor = zeros(size(data.raw_coor));
    data.rotated_coor(:,cols(1,:)) = coor(:,cols(1,:))-repmat(coor(1,cols(1,:)),[size(coor,1),1]);
    for i = 2:size(cols,1)
        data.rotated_coor(:,cols(i,:)) = data.rotated_coor(:,cols(1,:)) - [data.M_rotate*[coor(:,cols(1,:))-coor(:,cols(i,:))]']';
    end
    
    data.directions = zeros(size(data.raw_coor,1),size(cols,1)-1);
    for j = 1:size(cols,1)-1
        for i = 1:size(data.raw_coor,1)
            data.directions(i,j) = acos(dot([coor(i,cols(j,:))-coor(i,cols(j+1,:))]', [coor(1,cols(j,:))-coor(i,cols(j+1,:))]') / (norm([coor(i,cols(j,:))-coor(i,cols(j+1,:))]') * norm([coor(1,cols(j,:))-coor(i,cols(j+1,:))]'))); 
        end
    end


    data.curvature = zeros(1,size(coor,1));
    for i = 1:size(data.raw_coor,1)
        data.curvature(i) = calc_curvature(coor(i,:));
    end
    
    data.max_curve  = max(data.curvature(i));
    idx1 = find(abs(diff(data.directions(:,1)))>=0.02);
    idx2 = find(data.directions(:,1)>=2.5);
    if isempty(idx2)
        [~,idx2] = max(data.directions(:,1));
    end
    data.turn_time = [idx2(1)-idx1(1)]/Fs;
    if if_plot
        colors = flip(jet(size(coor,1)),2);
        figure
        for i = 1:30
            hold on,plot(data.rotated_coor(i,[1,3]),data.rotated_coor(i,[2,4]), 'Color',colors(i,:), 'Linewidth', 1.5)
            plot(data.rotated_coor(i,[3,5]),data.rotated_coor(i,[4,6]), 'Color',colors(i,:), 'Linewidth', 0.5)
        end
    end
end

