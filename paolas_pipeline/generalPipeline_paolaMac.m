% This software is based on:
% https://github.com/cortex-lab/allenCCF/

% Arber lab mainteined repository (forked from cortex-lab):
% https://github.com/paolahydra/allenCCF/tree/sliceRegistration
%
% % The following are needed for full functionality:
% Images of mouse brain slices (individually cropped or with multiple slices per image; coronal, sagittal, or transverse)
% Know the resolution (in microns per pixel) of these images
% A computer mouse with a scroll wheel
% MATLAB (R2017 or above used for testing)
% This repository. Add all folders and subfolders to your MATLAB path. All user-oriented scripts are in the 'SHARP-Track' folder.
% The npy-matlab repository: http://github.com/kwikteam/npy-matlab
% The Allen Mouse Brain Atlas volume and annotations (download all 4 files
% from this link: http://data.cortexlab.net/allenCCF/ )



% * remember to run one section at a time, instead of the whole script at once *

%%  always run: general settings (set once)
addpath(genpath('/Users/galileo/GitHub/allenCCF'))
rmpath(genpath('/Users/galileo/GitHub/WangLab_Allen'))
addpath(genpath('/Users/galileo/GitHub/matlabUtilities/'))

% directory of reference atlas files
pathToAtlas = '/Users/galileo/Documents/MATLAB/codeArberLab/anatomyRegistration/cortexLabCode/allen brain template files';
annotation_volume_location = fullfile(pathToAtlas, 'annotation_volume_10um_by_index.npy');
structure_tree_location = fullfile(pathToAtlas, 'structure_tree_safe_2017.csv');
template_volume_location = fullfile(pathToAtlas, 'template_volume_10um.npy');

% other stable settings:
% plane to view ('coronal', 'sagittal', 'transverse')
plane = 'coronal';
% transformation to use for registration:
transformationType = 'pwl';     %use 'projective', or 'pwl' (piece-wise linear: more advanced).


%%  always run: specific settings for the brain to register:
% move your images to a local disk (SSD possibly) for much faster processing!
image_folder = '/Users/galileo/dati/registered_brains_completed/992232'; % change this
image_tag = 'mouse_992232_';                                             % change this
microns_per_pixel = 3.8852; %take this value from your tiff filename


cd(image_folder)
save_folder = fullfile(image_folder, 'startingSingleSlices');


%% 1. do once, then skip: PP's preprocessing of axioscan images in ImageJ
% 1. batch convert all the axioscans ito tiff in ImageJ, using the macro: 
% batch_convert2tiff_highestResSeries_general.ijm.  Depending on how
% your images were acquired, you may want to choose the highest resolution
% series, or the second-highest one (there is a script for this too). 
% For cell detection, I have had good results for cell detection starting 
% from an image with 3.6 um per pixel.
% -- avoid saturating the right tail of the histogram if you want to
% detect stuff.


%% 2. do once, then skip: split axioscans in single figures (one per slice)
wait2confirmROI = 0;    % if true, you will need to double-click to confirm each ROI. If false, a cropped image is automatically saved.
                        % wait2confirmROI = 0; is much faster -- IF you don't make mistakes!
axioscanTiff_slideCropper(image_folder, image_tag, save_folder, microns_per_pixel, wait2confirmROI);

 
%% always run: filesystem and parameter definition - don't need to change
% directory of single histology images
image_folder = save_folder;

% if the images are individual slices (as opposed to images of multiple
% slices, which must be cropped using the cell CROP AND SAVE SLICES)
image_files_are_individual_slices = true;

% use images that are already at reference atlas resolution (here, 10um/pixel)
use_already_downsampled_image = false; 

% pixel size parameters: microns_per_pixel of large images in the image
% folder (if use_already_downsampled_images is set to false);
% microns_per_pixel_after_downsampling should typically be set to 10 to match the atlas
microns_per_pixel_after_downsampling = 10;

% ----------------------
% additional parameters
% ----------------------

% increase gain if for some reason the images are not bright enough
gain = 4;   % PP- for visualization only during cropping, and for atlas alignment

% size in pixels of reference atlas brain coronal slice, typically 800 x 1140
atlas_reference_size = [800 1140]; 
reference_size = [1320         800        1140]; %this is reloaded later as size(tv)

% -----------------------
% auto: naming definition
% -----------------------

% name of images, in order anterior to posterior or vice versa
% once these are downsampled they will be named ['original name' '_processed.tif']
image_file_names = dir([image_folder filesep '*.tif']); % get the contents of the image_folder
image_file_names = natsortfiles({image_file_names.name});
% image_file_names = {'slide no 2_RGB.tif','slide no 3_RGB.tif','slide no 4_RGB.tif'}; % alternatively, list each image in order


%% 3. do once, then skip: check all images for some to flip or adjust
Process_Histology_1_PP; 
%this will interactively allow you to crop, flip, rotate (and permute - untested) slices

% this can be quite fast if you don't dwell too much on rotations. 
% NOTE May 2021: No need to rotate, nor crop, unless you want to.
% Just check every slice and flip if necessary.

% NOTE:
% Images are saved in the 'preprocessed' folder only if you scroll through
% each of them in this process. When you get to the last image, go back
% once to save the last image too.

% IMPORTANT:
% no furter manipulation should be done to the images after this stage.

%% you will need to do cell detection on the *preprocessed* images.
% step 1:
% run  batch_split_invertColor_savePNG.ijm script in the 'preprocessed'
% folder

% step 2:
% run cellprofiler pipeline

%% 4. do once, then skip: downsample images for atlas registration (to the folder 'processed') - automatic and fast...
% % consider closing the previous figure when you are done preprocessing:
% close all
folder_preprocessed_images = fullfile(save_folder, 'preprocessed');     
Process_Histology_2_downsample_PP; %this will automatically downsample your *preprocessed* images and save them in the 'processed' folder for registration.
disp('Downsampled and boosted images were saved in the processed folder')
% This also increases the gain for better visualization during registration.


%% 5. Register each slice to the reference atlas
set(0, 'DefaultFigureWindowStyle', 'docked')
Navigate_Atlas_and_Register_Slices_PP;

%% as you are registering new slices, run this to keep your table of transformations T up to date.
T = saveTransformTable(fullfile(folder_processed_images, 'transformations'), image_file_names, reference_size);


%%



