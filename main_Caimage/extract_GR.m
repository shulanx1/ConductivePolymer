function df = extract_GR(f, R, dt, win)
if nargin < 3, win = [0, 0.01]; end
df = zeros(size(f));
for i = 1:size(f, 1)
    baseline = mean(f(i, ceil(win(1)/dt)+1:ceil(win(2)/dt)));
    df(i,:) = (f(i,:)-baseline)./R(i,:);
end
end