function a = spatial_contour(coor,FOV, Cn, thre)
% generate a spatial map with manually selected x and y coordinates

% a: weight of the spetial components

% coor: coordinates, x and y in column vector
% Cn: reference image
% thre: threshold (normalized brightness in the reference image)

if nargin <4, thre = 1;end

coor = [coor; coor(1,:)]';
x = coor(1,:);
y = coor(2,:);


left = min(x); right = max(x);  top = min(y); bottom = max(y);
[gridX, gridY] = meshgrid(left:1:right, top:1:bottom); % generate the coordinates for the search grid
gridX = reshape(gridX,1,[]);
gridY = reshape(gridY,1,[]);

in = inpolygon(gridX,gridY,x,y);
Cn1 = reshape(Cn(uint16(top:1:bottom), uint16(left:1:right)),1,[]);
idx_corr = Cn1 > thre;

idx = in & idx_corr;

a = sparse(FOV(1), FOV(2));
b = uint16(find(idx));
for i = 1:length(b)
    a(uint16(gridY(b(i))),uint16(gridX(b(i)))) = 1;
end

a = reshape(a, [] ,1);
end