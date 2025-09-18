

% CP = [-3.5 -3 -2 -1 0 1 3 3.5;
%     0    2   1 -1 1 2 2.5 1];

% cp1 = [1 3 4; 4 5 1];
% for i=1:20
%    cp1 = [cp1 [4;1]]; 
% end
% cp1 = [cp1 [8 3 1; 4 2 9]];

CP = [1 2 3 4; 2 4 2 4];

main = DeCasteObj(CP);

% TODO
% grafico della curvatura dinamico
% implementare suddivisione 
% 

%----------------------- 
% %curvature 
% h = 1/400;
% xt = main.X(1,:)
% yt = main.X(2,:);
% n = size(xt,2);
% %calcolo le derivate prima
% xt1 = (xt(2:n) - xt(1:n-1))/h;
% yt1 = (yt(2:n) - yt(1:n-1))/h;
% 
% %calcolo derivata seconda
% xt2 = ( xt(3:n) - 2*xt(2:n-1) + xt(1:n-2) ) / (h^2);
% yt2 = ((yt(3:n) - 2*yt(2:n-1) + yt(1:n-2))) / (h^2);
% 
% xt1 = xt1(2:size(xt1,2)); %non esiste per i = 0 la d''
% yt1 = yt1(2:size(yt1,2)); %non esiste per i = 0 la d''
% 
% size(xt1,2)
% size(xt2,2)
% 
% kt = ( abs((xt1.*yt2) - (yt1.*xt2) )) ./ ( (xt1.^2 + yt1.^2).^(3/2) )
% curv = figure(2);
% plot(linspace(0,1,398),kt,'-b')
% % axis equal
% axis([0 1 min(kt)-0.1 max(kt)+0.1])
% grid on
% set(curv,'units','normalized')
% set(curv,'color','w');
% set(curv, 'Position', [0 0 0.2 0.2]);








