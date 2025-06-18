function Y = stretch(X, T1, T2)
% contrast enhancement by stretching the LUT to fit the maximum display range

% X: imput image, can be in uint8, uint16 or double format
% T1: minimum display range, default to be the minimum pixel value
% T2: maximum display range, default to be 90% of the maximum pixel value

range = [min(reshape(X, [], 1)), quantile(reshape(X, [], 1),0.99)];

if nargin <2 , T1 = range(1) + 0*(range(2)-range(1)); end
if nargin <3 , T2 = range(1) + 0.75*(range(2)-range(1)); end

if isa(X,'uint8')
    L = 2^8;%256;
    f_diff = L/abs(T2 - T1);
    % compute histogram of X
    h = hist(X(:), 0:L-1);
    % initialize and assign Y
    nrow = size(X, 1);
    ncol = size(X, 2);
    Y = zeros(nrow, ncol, 'uint8');
    for i = 1:nrow
        for j = 1:ncol
            if X(i, j)>max([T1, T2])
                Y(i, j) = L-1;
            elseif X(i, j)<min([T1, T2])
                Y(i, j) = 0;
            else
                Y(i, j) = (X(i, j) - T1)*f_diff;
            end           
        end
    end   
elseif isa(X, 'uint16')
    L = 2^16;
    f_diff = L/abs(T2 - T1);
    % compute histogram of X
    h = hist(X(:), 0:L-1);
    % initialize and assign Y
    nrow = size(X, 1);
    ncol = size(X, 2);
    Y = zeros(nrow, ncol, 'uint16');
    for i = 1:nrow
        for j = 1:ncol
            if X(i, j)>max([T1, T2])
                Y(i, j) = L-1;
            elseif X(i, j)<min([T1, T2])
                Y(i, j) = 0;
            else
                Y(i, j) = (X(i, j) - T1)*f_diff;
            end           
        end
    end 
elseif isa(X, 'double')
    L = 1;
    f_diff = L/abs(T2 - T1);
    % compute histogram of X
    h = hist(X(:), 0:L-1);
    % initialize and assign Y
    nrow = size(X, 1);
    ncol = size(X, 2);
    Y = zeros(nrow, ncol, 'double');
    for i = 1:nrow
        for j = 1:ncol
            if X(i, j)>max([T1, T2])
                Y(i, j) = L;
            elseif X(i, j)<min([T1, T2])
                Y(i, j) = 0;
            else
                Y(i, j) = (X(i, j) - T1)*f_diff;
            end           
        end
    end 
    elseif isa(X, 'single')
    L = 1;
    f_diff = L/abs(T2 - T1);
    % compute histogram of X
    h = hist(X(:), 0:L-1);
    % initialize and assign Y
    nrow = size(X, 1);
    ncol = size(X, 2);
    Y = zeros(nrow, ncol, 'single');
    for i = 1:nrow
        for j = 1:ncol
            if X(i, j)>max([T1, T2])
                Y(i, j) = L;
            elseif X(i, j)<min([T1, T2])
                Y(i, j) = 0;
            else
                Y(i, j) = (X(i, j) - T1)*f_diff;
            end           
        end
    end 
else
    Y=X;
end