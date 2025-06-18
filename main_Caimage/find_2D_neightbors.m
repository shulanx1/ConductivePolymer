function coor_list = find_2D_neightbors(coor)
coor_list = [coor+[-1, -1];
    coor+[-1,0];
    coor+[-1,1];
    coor + [0, -1];
    coor + [0,1];
    coor + [1, -1];
    coor + [1, 0];
    coor + [1, 1]];

end