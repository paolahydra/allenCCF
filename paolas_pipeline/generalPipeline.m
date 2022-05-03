% This software is based on:
% https://github.com/cortex-lab/allenCCF/

% Arber lab mainteined repository (forked from cortex-lab):
% https://github.com/paolahydra/allenCCF/tree/sliceRegistration_confocal
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


%%  set once, then always run: specify paths and settings for the specific brain to register
% move your images to a local disk (SSD possibly) for much faster processing!
input_folder = '/Users/galileo/dati/registered_brains_completed/';   %change this
image_tag = 'SC_';                                               %change this - use an unequivocal tag for your experiment
microns_per_pixel = 16; %take this value from your tile metadata

% increase gain if for some reason the images are not bright enough
gain = 1;   % for visualization only: during cropping or atlas alignment


% stable settings (no need to change):
repositoryTag = 'SRC'; %this is the SliceRegistrationConfocal pipeline - these scripts are added to gitignore, useful to know where they came from

generalPipelinesFolder = 'generalPipelineUsedScripts'; %leave empty if you want to save the 
                        % mouse-specific script in the same folder as the
                        % generalPipeline.m script.  Or else, specify a new
                        % folder here.
                        
plane = 'coronal';  % plane to view ('coronal', 'sagittal', 'transverse')
% transformation to use for registration:
transformationType = 'pwl';     %use 'projective', or 'pwl' (piece-wise linear: more advanced).                       
                          
%%  always run: (set once to make it more easily portable across the computers you use)
[~,localhost] = system('hostname');
% archstr = computer('arch');     % it could be of use too

if strcmp(localhost(1:end-1), 'Paolas-MacBook-Pro.local')
    rootFolderScripts = '/Users/galileo/GitHub/';
    addpath(genpath('/Users/galileo/GitHub/matlabUtilities/npy-matlab-master/npy-matlab'))
    pathToAtlas = '/Users/galileo/Documents/MATLAB/codeArberLab/anatomyRegistration/cortexLabCode/allen brain template files';
elseif strcmp(localhost(1:end-1), 'f452d-96617e') % put yours here
    rootFolderScripts = 'C:\GitHub\';    % % -------------------   put yours here
    addpath(genpath('\\tungsten-nas.fmi.ch\tungsten\scratch\garber\BrainRegistration\code and atlas')) %includes npy
    pathToAtlas = '\\tungsten-nas.fmi.ch\tungsten\scratch\garber\BrainRegistration\code and atlas\allen brain template files\';
else 
    rootFolderScripts = 'C:\Users\patepaol\Documents\GitHub';  % % -------------------   put yours here
    addpath(genpath('\\tungsten-nas.fmi.ch\tungsten\scratch\garber\BrainRegistration\code and atlas')) %includes npy
    pathToAtlas = '\\tungsten-nas.fmi.ch\tungsten\scratch\garber\BrainRegistration\code and atlas\allen brain template files\';
end
clear localhost
addpath(genpath(fullfile(rootFolderScripts, 'allenCCF')))


% no need to change below:

annotation_volume_location = fullfile(pathToAtlas, 'annotation_volume_10um_by_index.npy');
structure_tree_location = fullfile(pathToAtlas, 'structure_tree_safe_2017.csv');
template_volume_location = fullfile(pathToAtlas, 'template_volume_10um.npy');

A = matlab.desktop.editor.getActive;                                                
if ~isempty(generalPipelinesFolder)                       
    genPipScriptLibrary = fullfile(rootFolderScripts, generalPipelinesFolder);
    if ~exist(genPipScriptLibrary, 'dir')
        mkdir(genPipScriptLibrary)
    end
else
    [genPipScriptLibrary, ~] = fileparts(A.Filename);
end


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


%
A.save;
scriptname = fullfile(genPipScriptLibrary, sprintf('generalPipeline_%s_%s.m',repositoryTag, image_tag(1:end-1)));

if strcmp(A.Filename, scriptname)
    fprintf('Script %s already existing and in use.\n',sprintf('generalPipeline_%s_%s.m',repositoryTag, image_tag(1:end-1)))
else
    if exist(scriptname, 'file') %do not overwrite
        scriptname = fullfile(genPipScriptLibrary, sprintf('generalPipeline_%s_%s_s.m',repositoryTag, image_tag(1:end-1)), string(datetime('now', 'Format','yyyyMMdd_hhmmss')));
    end
    copyfile(A.Filename, scriptname)
    A.close
    clear A
    matlab.desktop.editor.openAndGoToLine(scriptname, getcurrentline+4);
end

%%


%% do once, then skip: move your MAX_ full resolution images in the startingSingleSlices folder 
% (which was just created inside your main folder)

% first you need to run the imageJ script:
% GitHub/imageJ_batchProcessing/tiffs/batch_maxProjection_saturation_8bit.ijm

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




