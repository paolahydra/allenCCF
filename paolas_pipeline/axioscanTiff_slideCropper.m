function axioscanTiff_slideCropper(image_folder, image_tag, save_folder, imResolution_umperpixel, wait2confirmROI)

save_file_name = image_tag;     
if ~exist(save_folder)
    mkdir(save_folder)
end

%% name of images, in order anterior to posterior or vice versa
% once these are downsampled they will be named ['original name' '_processed.tif']
image_file_names = dir([image_folder filesep image_tag '*.tif']); % get the contents of the image_folder
image_file_names = natsortfiles({image_file_names.name});
% image_file_names = {'slide no 2_RGB.tif','slide no 3_RGB.tif','slide no 4_RGB.tif'}; % alternatively, list each image in order


%% get the image resolution from the filename, or else manually input it.
% imResolution_umperpixel = regexp(image_file_names{1}, '\_', 'split');
% imResolution_umperpixel = regexp(imResolution_umperpixel{end}, 'umppx', 'split');
% imResolution_umperpixel = str2double(imResolution_umperpixel{1})*1e6; %in microns
% imResolution_umperpixel = 2.60; % in microns

%% do the cropping
atlas_reference_size_um = [8000 11400];
reference_size = round(atlas_reference_size_um/imResolution_umperpixel);

warning('off')
for f = 1:length(image_file_names)
    fprintf('Opening %s...\n',image_file_names{f})
    disp('Click on the center of each slice.')
    disp('When done, close the figure to move on to the next one.')
    histology_figure = figure('Name','Histology Viewer');
    HistologyCropper_PP(histology_figure, image_folder, save_folder, image_file_names(f), reference_size, save_file_name, f, wait2confirmROI)
end
disp('Done.')

