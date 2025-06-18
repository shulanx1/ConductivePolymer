close all
% Load TIFF image
tif_dir = 'E:\data\polymer\IHC\61589 CMP S1 NeuN.tif';
[img, map] = imread(tif_dir);  % Replace with your file path[1][5]

% load txt file
pix_per_micron = 0.8045;

txt_dir = 'E:\data\polymer\IHC\NeuN\61589_nPBDF.txt';
s = readtable(txt_dir,'Delimiter','\t','HeaderLines',1);
x_coords = s.Var8*pix_per_micron;
y_coords = s.Var9*pix_per_micron;

txt_dir = 'E:\data\polymer\IHC\GFAP\61589_nPBDF.txt';
s = readtable(txt_dir,'Delimiter','\t','HeaderLines',1);
x_coords_G = s.Var8*pix_per_micron;
y_coords_G = s.Var9*pix_per_micron;

selected_points_N = cell(1,5);
selected_points_G = cell(1,5);
total_area = zeros(1,5);
polygons = cell(1,5);
for i = 1:5
    figure(i)
    imshow(img);
    hold on;
    % Draw polygon interactively
    h = drawpolygon('Color','r','LineWidth',1.5);  % Creates ROI object[1]
    title('Draw polygon: Left-click vertices, Double-click to finish');

    % Get polygon coordinates
    polygonVertices = h.Position;  % Nx2 matrix [x,y][1]
    polygons{i} = polygonVertices;
    polyObj = polyshape(polygonVertices);
    total_area(i) = area(polyObj)/(pix_per_micron^2);

    % Detect points inside polygon
    in = inpolygon(x_coords, y_coords, polygonVertices(:,1), polygonVertices(:,2));
    selected_points_N{i} = [x_coords(in), y_coords(in)];
    in_G = inpolygon(x_coords_G, y_coords_G, polygonVertices(:,1), polygonVertices(:,2));
    selected_points_G{i} = [x_coords_G(in_G), y_coords_G(in_G)];

    % Visualize results
    plot(polygonVertices(:,1), polygonVertices(:,2), 'r-', 'LineWidth',1.5);  % Polygon outline
    plot(selected_points_N{i}(:,1), selected_points_N{i}(:,2), 'bo', 'MarkerSize',4);  % Detected points
    hold off;
end
[filedir,filename,~] = fileparts(txt_dir);
save(fullfile(filedir, sprintf('%s.mat', filename)), 'selected_points_G', 'selected_points_N', 'total_area', 'polygons')
%%
M = readtable('E:\data\polymer\pool_IHC.xlsx', 'Sheet','Sheet1');
nPBDF_NeuN = zeros(3, 5);
control_NeuN = zeros(3, 5);
nPBDF_GFAP = zeros(3, 5);
control_GFAP = zeros(3, 5);
c_control = 1;
c_nPBDF = 1;
for i = 1:6
    if M.type(i)==1
        load(fullfile(M.folder{i}, M.filename{i}))
        for ii = 1:length(selected_points_N)
            nPBDF_NeuN(c_nPBDF, ii) = length(selected_points_N{ii})/(total_area(ii)/1e6);
            nPBDF_GFAP(c_nPBDF, ii) = length(selected_points_G{ii})/(total_area(ii)/1e6);
        end
        c_nPBDF = c_nPBDF+1;
    elseif M.type(i)==2
        load(fullfile(M.folder{i}, M.filename{i}))
        for ii = 1:length(selected_points_N)
            control_NeuN(c_control, ii) = length(selected_points_N{ii})/(total_area(ii)/1e6);
            control_GFAP(c_control, ii) = length(selected_points_G{ii})/(total_area(ii)/1e6);
        end
        c_control = c_control+1;
    end
end
colors = [[119,176,203]/255;[0,0,0];[0.5,0.5,0.5]]; % color for Na rich polymer
NeuN = zeros(size(nPBDF_NeuN,1), size(nPBDF_NeuN,2)+size(control_NeuN,2));
NeuN(:,2*[0:size(nPBDF_NeuN,2)-1]+1) = nPBDF_NeuN;
NeuN(:,2*[0:size(nPBDF_NeuN,2)-1]+2) = control_NeuN;
barplot_with_datapoint(NeuN, colors(1:2,:));

GFAP = zeros(size(nPBDF_GFAP,1), size(nPBDF_GFAP,2)+size(control_GFAP,2));
GFAP(:,2*[0:size(nPBDF_GFAP,2)-1]+1) = nPBDF_GFAP;
GFAP(:,2*[0:size(nPBDF_GFAP,2)-1]+2) = control_GFAP;
barplot_with_datapoint(GFAP, colors(1:2,:));

[data_table, within_design] = gen_table_for_ranova([1:5], {nPBDF_NeuN', control_NeuN'});
rm = fitrm(data_table,'measurements1-measurements10 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);

[data_table, within_design] = gen_table_for_ranova([1:5], {nPBDF_GFAP', control_GFAP'});
rm = fitrm(data_table,'measurements1-measurements10 ~ 1', 'WithinDesign', within_design);
AT = ranova(rm, 'WithinModel','treatment*voltage');
anova_table = anovaTable(AT, 'DV');
disp(anova_table);