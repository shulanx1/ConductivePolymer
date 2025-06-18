function [Y, r] = custom_fft(X, manual_flag, r, plot_flag)
% filter the linescan with FFT to remove the PMT linenoise
% X: imput image
% manual_flag: if the range to remove will be selected manualy. default to
% be 0
% plot_flag: plot the raw and filtered image, default to be 0

if nargin < 2, manual_flag = 0; end
if nargin < 3, r = [0.02, 0.45]; end % width/2 & hight of the blocked region in Fourier plane
if nargin < 4, plot_flag = 0; end

dim = size(X);
X1 = fft2(X, 2^nextpow2(max(dim)), 2^nextpow2(max(dim)));
X1_shift = fftshift(X1);

if manual_flag
    h = figure();
    imagesc(abs(X1_shift));
    colormap('gray');
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    imcontrast;
    p = drawrectangle('Linewidth', 1,'Color','Cyan');
    y1 = p.Position(2);
    y2 = p.Position(2)+p.Position(4);
    x1 = p.Position(1);
    x2 = p.Position(1)+p.Position(3);
    X1_shift(round(y1):round(y2), round(x1):round(x2)) = 0;
    X1_shift(size(X1_shift, 2)-[round(y2):-1:round(y1)], round(x1):round(x2)) = 0;
    imagesc(abs(X1_shift));
    close(h)
else
    x1 = round(size(X1_shift, 1)/2 - r(1)*size(X1_shift, 1));
    x2 = round(size(X1_shift, 1)/2 + r(1)*size(X1_shift, 1));
    y1 = 1;
    y2 = round(r(2)*size(X1_shift, 2));
    X1_shift(round(y1):round(y2), round(x1):round(x2)) = 0;
    X1_shift(size(X1_shift, 2)-[round(y2):-1:round(y1)], round(x1):round(x2)) = 0;
end

Y = ifft2(fftshift(X1_shift));
Y = abs(Y(1:dim(1), 1:dim(2)));

if plot_flag
    figure
    subplot(211), imshow(stretch(imrotate(X,90)));
    subplot(212), imshow(stretch(imrotate(Y,90)));
end
end