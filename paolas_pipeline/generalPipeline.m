addpath(genpath('/Users/galileo/GitHub/allenCCF'))
rmpath(genpath('/Users/galileo/GitHub/WangLab_Allen'))
addpath(genpath('/Users/galileo/GitHub/matlabUtilities/'))

%% PP's preprocessing of axioscan images
% 1. batch convert all the axioscans in series 2 tiff in axioscan. ADD NAME
% -- perhaps avoid saturating right tail
% 2. split in single figures using my matlab code: MODIFY!
edit axioscanTiff_slideCropper_PP % TO DO: fix input and output folder definition. Check frame size for cropping (make it standard)



%startingSingleSlices should be 2935x2060 or else registration errors will occour



%% put all the filesystem and parameter definition here up front (moved from Process_Histology)

% * remember to run one cell at a time, instead of the whole script at once *

% directory of histology images
image_folder = '/Users/galileo/dati/registered_brains_completed/992234'; %this has been fixed in the next version...
save_file_name = 'mouse992234_';  %check again this one

microns_per_pixel = 2.60;
wait2confirmROI = 0;    % if true, you will need to double-click to confirm each ROI. If false, a cropped image is automatically saved.
                        % wait2confirmROI = 0; is much faster -- IF you don't make mistakes!

% directory to save the processed images -- can be the same as the above image_folder
% results will be put inside a new folder called 'processed' inside of this image_folder
save_folder = fullfile(image_folder, 'startingSingleSlices');
cd(save_folder)

% if the images are cropped (image_file_are_individual_slices = false),
% name to save cropped slices as; e.g. the third cropped slice from the 2nd
% image containing many slices will be saved as: save_folder/processed/save_file_name02_003.tif


% if the images are individual slices (as opposed to images of multiple
% slices, which must be cropped using the cell CROP AND SAVE SLICES)
image_files_are_individual_slices = true;

% use images that are already at reference atlas resolution (here, 10um/pixel)
use_already_downsampled_image = false; 

% pixel size parameters: microns_per_pixel of large images in the image
% folder (if use_already_downsampled_images is set to false);
% microns_per_pixel_after_downsampling should typically be set to 10 to match the atlas
microns_per_pixel = 3.8852;
microns_per_pixel_after_downsampling = 10;

% ----------------------
% additional parameters
% ----------------------

% increase gain if for some reason the images are not bright enough
gain = 3.3;   % PP- for visualization only during cropping, and for atlas alignment

% size in pixels of reference atlas brain coronal slice, typically 800 x 1140
atlas_reference_size = [800 1140]; 
reference_size = [1320         800        1140]; %this is reloaded later as size(tv)

% -----------------------
% auto: naming definition
% -----------------------

% name of images, in order anterior to posterior or vice versa
% once these are downsampled they will be named ['original name' '_processed.tif']
image_file_names = dir(fullfile(save_folder, '*.tif')); % get the contents of the image_folder
image_file_names = natsortfiles({image_file_names.name});
% image_file_names = {'slide no 2_RGB.tif','slide no 3_RGB.tif','slide no 4_RGB.tif'}; % alternatively, list each image in order


%%
Process_Histology_1_PP; %this will interactively allow you to crop, flip, rotate (and permute - untested) slices
% this can be quite fast if you don't dwell too much on rotations. Count 30
% sec per slice, or less

% Images are saved in the 'preprocessed' folder.
% no furter manipulation should be done to the images after this stage.


%% you will need to do cell detection on the *preprocessed* images.
% Further transformations will be saved and it should be possible to
% realign new images or points later on. This needs a final test.

% step 1:
% run  batch_split_invertColor_savePNG.ijm script in the 'preprocessed'
% folder

% step 2:
% run cellprofiler pipeline

%%
% % consider closing the previous figure when you are done preprocessing:
% close all
Process_Histology_2_downsample_PP; %this will automatically downsample your *preprocessed* images and save them in the 'processed' folder for registration.
disp('Downsampled and boosted images were saved in the processed folder')
% It also increases the gain for better visualization during registration.

%%
Navigate_Atlas_and_Register_Slices_PP
% Manually align every slice starting from several interspersed slices to get a gist of the correct angles.
% After you have registered at least two slices, help yourself by calling
% the table T belo: in the workspace.
% Call it again to update it after saving more transformations, before
% looking at it in the workspace.
%%
T = saveTransformTable(fullfile(folder_processed_images, 'transformations'), image_file_names, reference_size);

% fig_table = tabulateImageInfo(image_file_names, save_folder, T) %check
% original size of images %not a general purpose function
%% load a reference table to check previous registration parameters
sf = '/Users/galileo/dati/registered_brains_completed/993031/startingSingleSlices/preprocessed';
image_file_names_31 = dir(fullfile(sf, '*.tif')); % get the contents of the image_folder
image_file_names_31 = natsortfiles({image_file_names_31.name});
% image_file_names_31'
T_31 = saveTransformTable(fullfile('/Users/galileo/dati/registered_brains_completed/993031/startingSingleSlices/processed', 'transformations'), image_file_names_31, reference_size);

%% When done, tabulate the ROI data (from cellprofiler analysis and registration)
%% red (TH)

% load the reference brain annotations
if ~exist('av','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
end
if ~exist('st','var')
    disp('loading structure tree...')
    st = loadStructureTree(structure_tree_location);
end
if ~exist('tv','var')
    tv = readNPY(template_volume_location);
%     tv_plot = tv; %for coronal slices
end


object_tag = 'red';
braincolor = 'r';
T_roi = Register_and_Tabulate_Rois(object_tag, save_folder, save_file_name, av, st, tv, microns_per_pixel, microns_per_pixel_after_downsampling, reference_size);

                                  

black_brain = false;
fwireframe = [];
fwireframe = plotWireFrame(T_roi, braincolor, black_brain, fwireframe, microns_per_pixel, microns_per_pixel_after_downsampling );


%% green (rabies)
object_tag = 'green';
T_roi = Register_and_Tabulate_Rois(object_tag, save_folder, save_file_name, av, st, tv, microns_per_pixel, microns_per_pixel_after_downsampling, reference_size);
braincolor = 'g';


% black_brain = false;
% fwireframe = plotWireFrame(T_roi, braincolor, black_brain, fwireframe );


%% script for further analysis of ROIs and plotting
edit analyzeDistributionOfCells

