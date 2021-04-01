function axioscanTiff_slideCropper(image_folder, image_tag)

% directory to save the processed images -- can be the same as the above image_folder
% results will be put inside a new folder called 'processed' inside of this image_folder
save_folder = fullfile(image_folder, 'startingSingleSlices');    %change it if you want to
save_file_name = image_tag;     %change it if you want to

if ~exist(save_folder)
    mkdir(save_folder)
end

%% name of images, in order anterior to posterior or vice versa
% once these are downsampled they will be named ['original name' '_processed.tif']
image_file_names = dir([image_folder filesep image_tag '*.tif']); % get the contents of the image_folder
image_file_names = natsortfiles({image_file_names.name});
% image_file_names = {'slide no 2_RGB.tif','slide no 3_RGB.tif','slide no 4_RGB.tif'}; % alternatively, list each image in order


%% get the image resolution from the filename, or else manually input it.
imResolution_umperpixel = regexp(image_file_names{1}, '\_', 'split');
imResolution_umperpixel = regexp(imResolution_umperpixel{end}, 'umppx.tif', 'split');
imResolution_umperpixel = str2double(imResolution_umperpixel{1});
% imResolution_umperpixel = 2.40;

%% do the cropping
atlas_reference_size_um = [8000 11400];
reference_size = round(atlas_reference_size_um/imResolution_umperpixel);

warning('off')
for f = 1:length(image_file_names)
    fprintf('Opening %s...\n',image_file_names(f))
    histology_figure = figure('Name','Histology Viewer');
    HistologyCropper_PP(histology_figure, save_folder, image_file_names(f), reference_size, save_file_name, f)
end
disp('Done.')

