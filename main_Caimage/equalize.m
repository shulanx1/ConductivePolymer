function Y = equalize(X)
% contrast enhancement by equalizing the LUT (fit the LUT with the same
% cdf of pixel values), the background noise is more obvious compared to
% the stretch method

% X: input image

X = double(X)/max(reshape(X, [], 1));
L = 1;
res = 256;  % assume uint8 resolution
% compute estimated cdf
F_est = zeros(1, res);
h = hist(X(:), 0:1/res:1 - 1/res);
for i = 1:res
    F_est(i) = sum(h(1:i))/sum(h);
end
% assign Y_s based on F_est
nrow = size(X, 1);
ncol = size(X, 2);
Y_s = zeros(nrow, ncol);
for i = 1: nrow
    for j = 1:ncol
        Y_s(i, j) = F_est(max([ceil(X(i, j)*res),1]));
    end
end
% normalize Y_s to get Y
Y = (Y_s - min(Y_s(:)))/(max(Y_s(:)) - min(Y_s(:)));
%plot F_est
% figure
% plot(0:L-1, F_est)
% xlim([0, L-1]);
% xlabel('pixel intensity')
% ylabel('estimated cumulative propability')
end