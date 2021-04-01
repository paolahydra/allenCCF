% #1: 35235 x 13365 %pixel dimensions of one image, zen series 1, acquired with 5x axioscan
% #2: 11745 x 4455 %corresponding pixel dimensions of the same image in series 2
% #3: 3915 x 1485
series(1).pxsize = 1.29507612; %um, with the 5x axioscan
series(2).pxsize =35235*series(1).pxsize/11745; %scaling to get pixel size in um for series 2
series(3).pxsize =35235*series(1).pxsize/3915;
%%  SET FILE AND PARAMETERS -- run only this now

% * remember to run one cell at a time, instead of the whole script at once *

% directory of histology images
image_folder = 'D:\axioscan_processing\992234';
image_tag = 'mouse_992234_';
% directory to save the processed images -- can be the same as the above image_folder
% results will be put inside a new folder called 'processed' inside of this image_folder
save_folder = image_folder;

% name of images, in order anterior to posterior or vice versa
% once these are downsampled they will be named ['original name' '_processed.tif']
image_file_names = dir([image_folder filesep image_tag '*.tif']); % get the contents of the image_folder
image_file_names = natsortfiles({image_file_names.name});

% image_file_names = {'slide no 2_RGB.tif','slide no 3_RGB.tif','slide no 4_RGB.tif'}; % alternatively, list each image in order
folder_processed_images = fullfile(save_folder, 'processed');
if ~exist(folder_processed_images)
    mkdir(folder_processed_images)
end

% atlas_reference_size_um = [8000 11400];  %used for 993031 -> scale to 34
% mm wide in illustrator
atlas_reference_size_um = [7400 10729]; %this should scale to 32 mm wide in illustrator
% atlas_reference_size_um = [6500 8000]; % smaller slices at lambda level
% atlas_reference_size_um = [6500 10000]; % 7N level
save_file_name = image_tag;

%fix the fact that multichannel is not kept except for the first file in the list
% processImageNum = 6;
warning('off')

for f = 1:length(image_file_names)
    ordinal = regexp(image_file_names{f}, '\_', 'split');
    ordinal = str2double(ordinal{3});

    histology_figure = figure('Name','Histology Viewer');
    HistologyCropper_PP_2021Jan3(histology_figure, save_folder, image_file_names(f), atlas_reference_size_um, save_file_name, ordinal)
end

