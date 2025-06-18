function coor = find_contour(coor_orig, Cn, thre)
    if nargin < 3, thre = 0.2;end
    FOV = [size(Cn, 1), size(Cn, 2)];
    Cn1 = zeros(FOV + [2, 2]);  % add dark boarders
    Cn1(2:FOV(1) + 1, 2:FOV(2) + 1) = Cn;
    coor = coor_orig + [1, 1];
    value = Cn1(coor_orig);
    value_new = value;
    while isempty(find(value_new < value*thre))
        for i = 1:size(coor, 1)
            if (coor(1)<=1) || (coor(1)>=FOV(1) + 1) || (coor(2)<=1) || (coor(2)>=FOV(1) + 1) 
                continue;
            end
            coor_temp = find_2D_neighbors(coor(i,:));
            for j = 1:size(coor_temp, 1)
                if Cn1(coor_temp(j,:)) > value*thre
                    if isempty(ismember(coor_temp(j,:),coor,'rows'))
                        coor = [coor;coor_temp(j,:)];
                        value_new = [value_new, Cn1(coor_temp(j,:))];
                    end
                end
            end
        end
    end
    

end