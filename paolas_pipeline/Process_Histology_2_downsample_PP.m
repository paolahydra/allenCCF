
%% LOAD AND PROCESS SLICE PLATE IMAGES
% PAOLA: here I just need to scroll through to save the dowsampled images.
% Made it automatic.

% finds or creates a folder location for processed images -- 
% a folder within save_folder called processed
folder_processed_images = fullfile(save_folder, 'processed'); %this just means 'downsampled', but I kept the old naming for consistency
if ~exist(folder_processed_images)
    mkdir(folder_processed_images)
end

% close all figures
close all

file_name_suffix = '_processed';
if ~use_already_downsampled_image
    % loop through images and downsample
    for f = 1: length(image_file_names)
        fname = fullfile(folder_preprocessed_images, image_file_names{f});
        % load histology image
        INFO = imfinfo(fname);
        nChannels = length(INFO);
        clear image
        for ch = 1:nChannels
            image(:,:,ch) = imread(fname, 'tif', ch); %this is the original image. Quality will be preserved.
        end
        
        % resize (downsample) image to reference atlas size
        original_image_size = size(image);
        image = imresize(image, [round(original_image_size(1)*microns_per_pixel/microns_per_pixel_after_downsampling)  NaN]);
        image = image*gain;
        if size(image,3) == 2
            image = cat(3, image, zeros(size(image(:,:,1))));
        end
        
        % to visualiye RGB
        image = image(:,:,1:3);
        imwrite(image, fullfile(folder_processed_images, [image_file_names{f}(1:end-4) file_name_suffix '.tif']))
        

%         %% to save all the channels
%         fname  = fullfile(folder_processed_images, [image_file_names{f}(1:end-4) file_name_suffix '.tif']);
%         imwrite(image(:,:,1), fname);
%         for i = 2:size(image, 3)
%             imwrite(image(:,:,i), fname, 'WriteMode', 'append');
%         end
    end
end
clear image original_image_size f   

% % if the images need to be downsampled to 10um pixels (use_already_downsampled_image = false), 
% % this will downsample and allow you to adjust contrast of each channel of each image from image_file_names
% %
% % if the images are already downsampled (use_already_downsampled_image = true), this will allow
% % you to adjust the contrast of each channel
% %
% % Open Histology Viewer figure
% try; figure(histology_figure);
% catch; histology_figure = figure('Name','Histology Viewer'); end
% warning('off', 'images:initSize:adjustingMag'); warning('off', 'MATLAB:colon:nonIntegerIndex');
% 
% % Function to downsample and adjust histology image
% HistologyBrowser(histology_figure, save_folder, image_folder, image_file_names, folder_processed_images, image_files_are_individual_slices, ...
%             use_already_downsampled_image, microns_per_pixel, microns_per_pixel_after_downsampling, gain)




