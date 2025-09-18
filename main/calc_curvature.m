function K = calc_curvature(p)
x1 = p(1); y1 = p(2); % Point 1
x2 = p(3); y2 = p(4); % Point 2 (middle)
x3 = p(5); y3 = p(6);% Point 3

% Calculate the lengths of the triangle sides
a = sqrt((x2-x1)^2 + (y2-y1)^2); % Side opposite point 3
b = sqrt((x3-x2)^2 + (y3-y2)^2); % Side opposite point 1
c = sqrt((x3-x1)^2 + (y3-y1)^2); % Side opposite point 2

% Calculate the area of the triangle using Heron's formula
s = (a+b+c)/2;
Area = sqrt(s*(s-a)*(s-b)*(s-c));

% Calculate the curvature (Menger curvature)
if abs(Area) < eps || abs(a*b*c) < eps
    K = 0; % Collinear points
else
    K = 4 * Area / (a * b * c);
end