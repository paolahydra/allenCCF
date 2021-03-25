% ------------------------------------------------------------------------
%          Crop, Rotate, Adjust Contrast, and Downsample Histology
% ------------------------------------------------------------------------


%%  SET FILE AND PARAMETERS - PAOLA: this can be made slimmer
% moved to generalPipeline.m
% run this script from there now.



%% PAOLA: %% FIRST GO THROUGH TO FLIP HORIZONTAL SLICE ORIENTATION, ROTATE, CROP, SHARPEN (disabled), and CHANGE ORDER (untested)
% this is now done on the full resolution image, that will be used for
% further analysis. DO NOT save changes in the contrast or gain of these images. 
% close all figures
% make SliceFlipper such that it reads and saves original (not _processed) images
close all
folder_preprocessed_images = fullfile(save_folder, 'preprocessed');
if ~exist(folder_preprocessed_images)
    mkdir(folder_preprocessed_images)
    for f = 1: length(image_file_names)
        fname = fullfile(image_folder, image_file_names{f});
        copyfile(fname, folder_preprocessed_images)
    end
end            
% this takes images from folder_processed_images ([save_folder/processed]),
% and allows you to rotate, flip, sharpen, crop, and switch their order, so they
% are in anterior->posterior or posterior->anterior order, and aesthetically pleasing
% 
% it also pads images smaller than the reference_size and requests that you
% crop images larger than this size
%
% note -- presssing left or right arrow saves the modified image, so be
% sure to do this even after modifying the last slice in the folder
slice_figure = figure('Name','Slice Viewer');
SliceFlipper_PP(slice_figure, folder_preprocessed_images, atlas_reference_size, gain)  %I have edited this file to save at least some of the transformations (rotation angle and whether flipped) -PP
%PP: gain here is now used for visualization only -- not for saving!

