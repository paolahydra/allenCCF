%% use centroids from csv file to make an ROI image
% d_red are the images that define who gets processeed here

mouseNum = 993031;
excelFileNameG = '/Users/galileo/dati/Jan2021/anatomy/mouse_993031/ShrinkToObjectCenters_green.csv';
excelFileNameR = '/Users/galileo/dati/Jan2021/anatomy/mouse_993031/ShrinkToObjectCenters_red.csv';
cellProfilerOutputFolder = '/Users/galileo/dati/Jan2021/anatomy/mouse_993031';
target_img = '03.015';

%% extract centroids for the specific image
R = readtable(excelFileNameR);
filenames = unique(R.FileName_red, 'stable');
imageNumbers = unique(R.ImageNumber, 'stable');
imgN = imageNumbers(contains(filenames, target_img)); %this is valid across colors. I will use it in the green channel.
assert(contains(unique(R.FileName_red(R.ImageNumber == imgN)), target_img))
xR = R.Location_Center_X(R.ImageNumber==imgN);
yR = R.Location_Center_Y(R.ImageNumber==imgN);
% get freen centroids
G = readtable(excelFileNameG);
x = G.Location_Center_X(G.ImageNumber==imgN);
y = G.Location_Center_Y(G.ImageNumber==imgN);
%%
    im = imread(fullfile(cellProfilerOutputFolder, filenames{imgN}));
    
    %% get centroids for red image
    ax = axes; hold on 
    axis image
    axis on
    ax.YDir = 'reverse';
    ax.XLim = [0.5, size(im,2)+0.5];
    ax.YLim = [0.5, size(im,1)+0.5];
    %with axis lines:
%     ax.Box = 'on';
%     ax.XAxis.Visible = 'on';
%     ax.YAxis.Visible = 'on';
%     ax.XTick = [];
%     ax.YTick = [];
    %without axis lines:
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    scatter(x, y, 4, [0 1 0], 'filled');
   
    %% make ROI image instead
    roi_array = zeros(size(im),'uint8');
    for r = 1:length(x)
        roi_array(y(r)+1, x(r)+1) = 250;
    end
    figure; imshow(roi_array, [])
    hold on
    scatter(x+0.5, y+0.5, 8, [0 1 0], 'filled');
    
    
    
    %% scale the coordinates first, than make an ROI image already in the DS dimensions
    microns_per_pixel = 3.8852;
    microns_per_pixel_after_downsampling = 10;
    original_image_size = size(im);
    disp('downsampling coordinates...')
    xDS = round((x+1) * microns_per_pixel/microns_per_pixel_after_downsampling);
    yDS = round((y+1) * microns_per_pixel/microns_per_pixel_after_downsampling);
    image = imread('/Users/galileo/dati/Jan2021/anatomy/mouse_993031/mouse_993031_03.015.tif (green).png');
    imageDS = imresize(image, [round(original_image_size(1)*microns_per_pixel/microns_per_pixel_after_downsampling)  NaN]);
%     %check
%     figure; hold on
%     imshow(imageDS, []);
%     scatter(xDS, yDS, 6, [0 1 0], 'filled');
    % now make and save the ROI image
    roi_array = zeros(size(imageDS),'uint8');
    for r = 1:length(x)
        roi_array(yDS(r), xDS(r)) = 255;
    end
    imwrite(roi_array,fullfile(cellProfilerOutputFolder, 'ROI_G_DS.png'))
    %this would work if I did not rotate the images in HistologyBrowser...
    
    