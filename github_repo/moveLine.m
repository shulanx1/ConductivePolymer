
clc;
clear;
clf;
figure(1);
grid on;

l = plot([0 4],[3 4],'-b');

for i=0:10
    set(l,'XData',[i i+4]);
    axis([0 20 0 8])
    grid on;
    pause(0.2);
end
hold on;
p = plot(5,5,'ob');
hold off;
