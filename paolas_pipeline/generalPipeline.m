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

% addpath(genpath('C:\GitHub\allenCCF')) %clone the repository from : https://github.com/paolahydra/allenCCF/tree/sliceRegistration_confocal     
addpath(genpath('/Users/galileo/GitHub/allenCCF'))

addpath(genpath('\\tungsten-nas.fmi.ch\tungsten\scratch\garber\BrainRegistration\code and atlas'))
% pathToAtlas = '\\tungsten-nas.fmi.ch\tungsten\scratch\garber\BrainRegistration\code and atlas\allen brain template files\';
pathToAtlas = '/Users/galileo/Documents/MATLAB/codeArberLab/anatomyRegistration/cortexLabCode/allen brain template files';




annotation_volume_location = fullfile(pathToAtlas, 'annotation_volume_10um_by_index.npy');
structure_tree_location = fullfile(pathToAtlas, 'structure_tree_safe_2017.csv');
template_volume_location = fullfile(pathToAtlas, 'template_volume_10um.npy');


% other stable settings:
% plane to view ('coronal', 'sagittal', 'transverse')
plane = 'coronal';
% transformation to use for registration:
transformationType = 'pwl';     %use 'projective', or 'pwl' (piece-wise linear: more advanced).


%%  set once, then always run: specify paths and settings for the specific brain to register
% move your images to a local disk (SSD possibly) for much faster processing!
input_folder = '/Users/galileo/dati/registered_brains_completed/Chiara';   %change this
image_tag = 'MAX_Rabies_Cerv_uni_1_';                                               %change this - use an unequivocal tag for your experiment
microns_per_pixel = 1.2980; %take this value from your tile metadata

% increase gain if for some reason the images are not bright enough
gain = 1;   % for visualization only: during cropping or atlas alignment

if ~strcmp( image_tag(end), '_')
    image_tag = cat(2, image_tag, '_');
end

cd(input_folder)
image_folder = fullfile(input_folder, 'startingSingleSlices');
if ~exist(image_folder, 'dir')
    mkdir(image_folder);
end
save_folder = fullfile(image_folder, 'processed');
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end
% folder_processed_images = fullfile(image_folder, 'processed');
% if ~exist(folder_processed_images, 'dir')
%     mkdir(folder_processed_images);
% end

%% do once, then skip: move your MAX_ full resolution images in the startingSingleSlices folder 
% (which was just created inside your main folder)

% first you need to run the imageJ script:
% GitHub/imageJ_batchProcessing/tiffs/batch_maxProjection_saturation_8bit.ijm

%% do once, then skip: save the script and then create a new version with specific parameters - continue with the new script.

% DO SAVE the script (generalPipeline.m)!!. Do not do ctrl-S: 
% *** you must click on the save button, otherwise this will not work. ***

originalscript = which('generalPipeline');
[a, b] = fileparts(originalscript);
scriptname = fullfile(a, sprintf('generalPipeline_%s.m',image_tag(1:end-1)));
copyfile(originalscript, scriptname)
edit(scriptname)


%% always run: filesystem and parameter definition - don't need to change

% if the images are individual slices (as opposed to images of multiple
% slices, which must be cropped using the cell CROP AND SAVE SLICES)
image_files_are_individual_slices = true;

% use images that are already at reference atlas resolution (here, 10um/pixel)
use_already_downsampled_image = false; 

% pixel size parameters: microns_per_pixel of large images in the image
% folder (if use_already_downsampled_images is set to false);
% microns_per_pixel_after_downsampling should typically be set to 10 to match the atlas
microns_per_pixel_after_downsampling = 10;


% additional parameters
% size in pixels of reference atlas brain coronal slice, typically 800 x 1140
atlas_reference_size = [800 1140]; 
reference_size = [1320 800 1140];


% auto: naming definition
% name of images, in order anterior to posterior or vice versa
% once these are downsampled they will be named ['original name' '_processed.tif']
image_file_names = dir([image_folder filesep '*.tif']); % get the contents of the image_folder
image_file_names = natsortfiles({image_file_names.name});
% image_file_names = {'slide no 2_RGB.tif','slide no 3_RGB.tif','slide no 4_RGB.tif'}; % alternatively, list each image in order


%% set and do once, then skip: tranform your full resolution image into a registration-ready image, keeping track of the transformation
%set only the rotation to a standard coronal orientation:
Transf.rotation = 90;  % specify angle (deg) 
% positive angle is counterclockwise rotation. 
% To rotate the image clockwise, specify a negative value for angle.


% nothing to set below:
atlas_reference_size_um = microns_per_pixel_after_downsampling * atlas_reference_size;
Transf.atlas_reference_size = atlas_reference_size; 
Transf.reference_originalImage_RowCol_size = round(atlas_reference_size_um/microns_per_pixel);
Transf.downsamplingFactor = microns_per_pixel/microns_per_pixel_after_downsampling;

confocalTiff_slideDilater(image_folder, image_tag, save_folder, Transf); %save_folder is 'processed'

%% now check for slice flipping and update transformation file
slice_figure = figure('Name','Slice Viewer');
SliceFlipper_PP_confocal(slice_figure, save_folder, atlas_reference_size, gain)  


% IMPORTANT:
% no furter manipulation should be done to the images after this stage.



%% Register each slice to the reference atlas
set(0, 'DefaultFigureWindowStyle', 'docked')
folder_processed_images = save_folder; %legacy
Navigate_Atlas_and_Register_Slices_PP;


%% as you are registering new slices, run this to keep your table of transformations T up to date.
T = saveTransformTable(fullfile(folder_processed_images, 'transformations'), image_file_names, reference_size);




%% 6. do once: when finished with the registr to atlas, do this to register and tabulate the detected cells too.
object_tag = 'green'; 
tabulateData_confocal;


%% plot?
braincolor = 'g';
fwireframe = [];
black_brain = false;
fwireframe = plotWireFrame(T_roi, braincolor, black_brain, fwireframe, microns_per_pixel, microns_per_pixel_after_downsampling );



%% post-registration (still evaluate all mandatory blocks above before starting)
edit analyzeDistributionOfCells




