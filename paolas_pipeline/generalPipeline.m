
% https://github.com/paolahydra/allenCCF/tree/sliceRegistration

% This software is based on:
% https://github.com/cortex-lab/allenCCF/

%
% % The following are needed for full functionality:
% Images of mouse brain slices (individually cropped or with multiple slices per image; coronal, sagittal, or transverse)
% Know the approximate resolution (in microns per pixel) of these images
% A computer mouse with a scroll wheel
% MATLAB (R2017 or above used for testing)
% This repository. Add all folders and subfolders to your MATLAB path. All user-oriented scripts are in the 'SHARP-Track' folder.
% The npy-matlab repository: http://github.com/kwikteam/npy-matlab
% The Allen Mouse Brain Atlas volume and annotations (download all 4 files
% from this link: http://data.cortexlab.net/allenCCF/ )



%% HOW TO RUN THIS PIPELINE
% * remember to run one section at a time, instead of the whole script at once *

% follow the instructions in the section's header:

%   - "set once", means that the section contains experiment specific
%   settings that you need to set once when you start. Once set, you don't
%   need to change those anymore (they will be saved with your
%   expt-specific script).

%   - "always run", means that you need to evaluate the section every time
%   you re-start the program. More clearly, after closing matlab, or simply 
%   clearing the workspace, you can pick up from where you had left, since  
%   the intermediate output is always saved. However, when you restart you 
%   need to run all the "always-run" sections to reload the needed
%   variables). 

%   - "do once, then skip": this is a section that only saves intermediate
%   output to the disk, so you only need to evaluate once. If you later
%   restart the program and you have already done that step, you can skip it.



%%  set once, then always run: specify paths and settings for the specific brain to register
% move your images to a local disk (SSD possibly) for much faster processing!
input_folder = '/Users/galileo/dati/registered_brains_completed/mouse_1031707';   %where your images are
image_tag = 'mouse_1031707_';      %fixed root of your image filenames: use an unequivocal tag for your experiment
image_files_are_individual_slices = 1;     % Do your images contain MULTIPLE brain slices?
microns_per_pixel = 5.1803;        %best guess (no need for submicron precision. If you are way off, you will crop too much or too little)

% increase gain if for some reason the images are not bright enough
gain = 1;   % for visualization only: during cropping or atlas alignment







% stable settings (no need to change):
repositoryTag = 'SR'; %this is the SliceRegistrationConfocal pipeline - these scripts are added to gitignore, useful to know where they came from

generalPipelinesFolder = 'generalPipelineUsedScripts'; %leave empty if you want to save the 
                        % mouse-specific script in the same folder as the
                        % generalPipeline.m script.  Or else, specify a new
                        % folder here.

plane = 'coronal';  % plane to view ('coronal', 'sagittal', 'transverse')
% transformation to use for registration:
transformationType = 'pwl';     %use 'projective', or 'pwl' (piece-wise linear: more advanced).                       

if ~strcmp( image_tag(end), '_')
    image_tag = cat(2, image_tag, '_');
end

A = matlab.desktop.editor.getActive; 
[~, fn] = fileparts(A.Filename);
if ~strcmp(fn, 'generalPipeline')
    if ~contains(fn, image_tag(1:end-1))
        warning('this script will be saved with the new settings when you evaluate the next session.')
    end
end
clear fn


%%  always run: (set once to make it more easily portable across all the computers you use)

% single machine setting (change with your paths):
rootFolderScripts = '/Users/galileo/GitHub/'; %path to this repo
pathToAtlas = '/Users/galileo/Documents/MATLAB/codeArberLab/anatomyRegistration/cortexLabCode/allen brain template files'; % path to atlas data
addpath(genpath('/Users/galileo/Documents/MATLAB/codeArberLab/anatomyRegistration/cortexLabCode/npy-matlab-master/npy-matlab')) %path to npy-matlab
addpath(genpath(fullfile(rootFolderScripts, 'allenCCF')))

% % if you use the code across different machines (with different filesystem
% % and/or paths), you may find convenient to specify all of them at once
% % here:
% [~,localhost] = system('hostname');
% % archstr = computer('arch');     % it could be of use too
% 
% if strcmp(localhost(1:end-1), 'Paolas-MacBook-Pro.local')
%     rootFolderScripts = '/Users/galileo/GitHub/';
%     addpath(genpath('/Users/galileo/GitHub/matlabUtilities/npy-matlab-master/npy-matlab'))
%     pathToAtlas = '/Users/galileo/Documents/MATLAB/codeArberLab/anatomyRegistration/cortexLabCode/allen brain template files';
% elseif strcmp(localhost(1:end-1), 'f452d-96617e') % put yours here
%     rootFolderScripts = 'C:\GitHub\';    % % -------------------   put yours here
%     addpath(genpath('\\tungsten-nas.fmi.ch\tungsten\scratch\garber\BrainRegistration\code and atlas')) %includes npy
%     pathToAtlas = '\\tungsten-nas.fmi.ch\tungsten\scratch\garber\BrainRegistration\code and atlas\allen brain template files\';
% else 
%     rootFolderScripts = 'C:\Users\patepaol\Documents\GitHub';  % % -------------------   put yours here
%     addpath(genpath('\\tungsten-nas.fmi.ch\tungsten\scratch\garber\BrainRegistration\code and atlas')) %includes npy
%     pathToAtlas = '\\tungsten-nas.fmi.ch\tungsten\scratch\garber\BrainRegistration\code and atlas\allen brain template files\';
% end
% clear localhost
% addpath(genpath(fullfile(rootFolderScripts, 'allenCCF')))




% no need to change below:

annotation_volume_location = fullfile(pathToAtlas, 'annotation_volume_10um_by_index.npy');
structure_tree_location = fullfile(pathToAtlas, 'structure_tree_safe_2017.csv');
template_volume_location = fullfile(pathToAtlas, 'template_volume_10um.npy');

                                             
if ~isempty(generalPipelinesFolder)                       
    genPipScriptLibrary = fullfile(rootFolderScripts, generalPipelinesFolder);
    if ~exist(genPipScriptLibrary, 'dir')
        mkdir(genPipScriptLibrary)
    end
else
    [genPipScriptLibrary, ~] = fileparts(A.Filename);
end



cd(input_folder)
image_folder = fullfile(input_folder, 'startingSingleSlices');
if ~exist(image_folder, 'dir')
    mkdir(image_folder);
end
folder_preprocessed_images = fullfile(image_folder, 'preprocessed');  
if ~exist(folder_preprocessed_images, 'dir')
    mkdir(folder_preprocessed_images);
end
save_folder = fullfile(image_folder, 'processed');
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

% save and start using the experiment-specific script
A.save; %save the current script
scriptname = fullfile(genPipScriptLibrary, sprintf('generalPipeline_%s_%s.m',repositoryTag, image_tag(1:end-1)));

if strcmp(A.Filename, scriptname)
    fprintf('Script %s already existing and in use.\n',sprintf('generalPipeline_%s_%s.m',repositoryTag, image_tag(1:end-1)))
else
    if exist(scriptname, 'file') %do not overwrite
        scriptname = fullfile(genPipScriptLibrary, sprintf('generalPipeline_%s_%s_%s.m',repositoryTag, image_tag(1:end-1), string(datetime('now', 'Format','yyyyMMdd_hhmmss'))));
    end
    copyfile(A.Filename, scriptname)
    if ~strcmp(which('generalPipeline'), A.Filename)
        A.close
    end
    clear A
    matlab.desktop.editor.openAndGoToLine(scriptname, getcurrentline+4);
end


%% 1. do once, then skip: assess, preprocess and locate your starting images
if ~image_files_are_individual_slices

    %% 1.0 (optional) - PP's preprocessing of ZEISS axioscan images in ImageJ
    
    % This step only applies if you start with Zeiss .czi format (axioscan).
    % If your images are already .tiff, you can skip this (or may need
    % different preprocessing).
    
    % goal is to batch convert all the axioscans .czi files to tiff.
    % Download macros from:
    % https://github.com/paolahydra/imageJ_batchProcessing/tree/main/czi_ZeissAxioscans/aaa_bestToUse
    
    
    % Instructions.
    % Drag the following macro into FIJI, then Run it (select a parent folder,
    % and it will batch-convert all .czi files in all subfolders):
    % batch_convert2tiff_highestResSeries_general.ijm
    
    % Depending on how your images were acquired, you may want to choose the
    % highest resolution series, or the second-highest one (there is a script
    % for this too).
    % For cell detection, I have had good results starting from an image with
    % 3.6 um per pixel.
    % NOTE: avoid saturating the right tail of the histogram if you want to
    % detect objects.
    
    
    %% 1.1  do once, then skip: split multislice images in input_folder into single-slice images
    % PRO TIP: if it is the first time you run this, set wait2confirmROI = 1:
    % you will get a much better sense of what you are doing.
    % Then set it to 0 and you will run this step way faster!
    
    wait2confirmROI = 0;    % if true, you will need to double-click to confirm each ROI. If false, a cropped image is automatically saved.
    % wait2confirmROI = 0; is much faster -- IF you don't make mistakes!
    axioscanTiff_slideCropper(input_folder, image_tag(1:end-1), image_folder, microns_per_pixel, wait2confirmROI);
    
else
    
    image_file_names = dir([input_folder filesep image_tag(1:end-1) '*.tif']); % get the contents of the image_folder
    image_file_names = natsortfiles({image_file_names.name});
    assert(~isempty(image_file_names), 'check that images are located in the input_folder or your image name definition')
    for f = 1:length(image_file_names)
        movefile(image_file_names{f}, image_folder)
    end
end

%% 2. always run: filesystem and parameter definition - don't need to change

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


%% 3. do once, then skip: check all images for some to flip or adjust
Process_Histology_1_PP; 
%this will interactively allow you to crop, flip, rotate (and permute - untested) slices

% NOTE May 2021: No need to rotate, nor crop, unless you want to.
% Just check every slice and flip if necessary.
% this step can be quite fast if you don't dwell too much on rotations/cropping. 

% IMPORTANT:
% no furter manipulation should be done to the images after this stage.

%% you will need to do cell detection on the *preprocessed* images.
% e.g.:
% step 1:
% run  batch_split_invertColor_savePNG.ijm script in the 'preprocessed'
% folder

% step 2:
% run cellprofiler pipeline or trackmate or other detection software. Get
% csv with coordinates

%% 4. do once, then skip: downsample images for atlas registration (to the folder 'processed') - automatic and fast...
% % consider closing the previous figure when you are done preprocessing:
% close all

% gain = 1;  % reset the gain if you found a better value in the previous step

Process_Histology_2_downsample_PP; %this will automatically downsample your *preprocessed* images and save them in the 'processed' folder for registration.
disp('Downsampled and boosted images were saved in the processed folder')
% This also increases the gain for better visualization during
% registration. For some very dark images you may need to set a higher gain and
% re-run this block.


%% 5. Register each slice to the reference atlas
set(0, 'DefaultFigureWindowStyle', 'docked')
Navigate_Atlas_and_Register_Slices_PP;

%% as you are registering new slices, run this to keep your table of transformations T up to date.
% PRO TIP: Copy and paste the follOwing line in your command window and
% hit enter to run it.
% This way, it will be easily accessible in the command window as a
% recently run command, simply through the up arrow.
% Run the line to update the content of the T variable as often as you
% register new slices.
% Inspect the T table to have an overview of all your offsets and
% registration parameters.

% It is useful to first coarsely register a few (2-4) slices with good
% landmarks spanning a good extent of your volume, so to establish the best
% offsets for your volume. 
% NOTE however that the DV angle may change at the level of the 
% midbrain/hindbrain transition, so it is best to establish offsets
% separately for forebrain and hindbrain if needed.

% Then, you can more precisely register your specific slices of interest
% following the suggested values in the T table.

T = saveTransformTable(fullfile(save_folder, 'transformations'), image_file_names, reference_size);



%% 6. do once: when finished with the registr to atlas, do this to register and tabulate the detected cells too.
object_tag = 'green'; 
tabulateData;


%% demo: plot all detected cells in the 3D model
braincolor = 'g';
fwireframe = [];
black_brain = false;
fwireframe = plotWireFrame(T_roi, braincolor, black_brain, fwireframe, microns_per_pixel, microns_per_pixel_after_downsampling );



%% post-registration (still evaluate all mandatory blocks above before starting)
edit analyzeDistributionOfCells




