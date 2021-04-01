addpath(genpath('/Users/galileo/GitHub/allenCCF'))

%% PP's preprocessing of axioscan images
% 1. batch convert all the axioscans ito tiff in ImageJ, using the macro: 
% batch_convert2tiff_highestResSeries_general.ijm.  Depending on how
% your images were acquired, you may want to choose the highest resolution
% series, or the second-highest one. For cell detection, I have had good 
% results for cell detection starting from an image with 3.9 um per pixel.
% -- avoid saturating the right tail of the histogram.


% 2. split in single figures using my matlab code: MODIFY!
image_folder = 'D:\axioscan_processing\992234';
image_tag = 'mouse_992234_';
axioscanTiff_slideCropper(image_file_names(f)); % TO DO: fix input and output folder definition. Check frame size for cropping (make it standard)



%% put all the filesystem and parameter definition here up front (moved from Process_Histology)

% * remember to run one cell at a time, instead of the whole script at once *

% directory of histology images
image_folder = '/Users/galileo/dati/Mar2021/reTestCortexLabAlignmentToAllenAtlas/mouse993031';

% directory to save the processed images -- can be the same as the above image_folder
% results will be put inside a new folder called 'processed' inside of this image_folder
save_folder = image_folder;

% if the images are cropped (image_file_are_individual_slices = false),
% name to save cropped slices as; e.g. the third cropped slice from the 2nd
% image containing many slices will be saved as: save_folder/processed/save_file_name02_003.tif
save_file_name = 'mouse993031_';  %check again this one

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


% -----------------------
% auto: naming definition
% -----------------------

% name of images, in order anterior to posterior or vice versa
% once these are downsampled they will be named ['original name' '_processed.tif']
image_file_names = dir([image_folder filesep '*.tif']); % get the contents of the image_folder
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

edit Navigate_Atlas_and_Register_Slices_PP.m; % set paths if standalone running, and run it block by block (not as crucial in fact). This can be changed to be called from here with some more checks for standalone running too.
% otherwise the path is taken from here

% Manually align every slice (help yourself with an excel table, otherwise it is impossible)
% Save the transformations.

% Done. Extras below.




%% how to retrieve and re-apply transformations
% run cell detection and cell density algorithms with cellProfiler

%% apply the same transformation to other images

%change this one:
trasformanda_tag = '.tif (green).png'; %adapt here


transf_atlasreg_folder = fullfile(image_folder, 'processed/transformations');
trasformanda_folder = fullfile(image_folder, 'preprocessed');
d = dir(fullfile(transf_atlasreg_folder, '*processed_transform_data.mat'));
for t = 1:length(d)
    % extract image root with ordinal tag
    fileroot = regexp(d(t).name, '_processed_transform_data.mat', 'split');
    fileroot = fileroot{1};

    % load the transformation data
    transfMAT = fullfile(transf_atlasreg_folder, d(t).name);
    st = load(transfMAT); %save_transform
    transform_data = st.save_transform;
    clear st
    
    % load the image or data to be transformed
    trasformanda_PNG = fullfile(trasformanda_folder, [fileroot, trasformanda_tag]); % to be transformed by applying transformations
    im = imread(trasformanda_PNG);
    original_image_size = size(im);
    im = imresize(im, [round(original_image_size(1)*microns_per_pixel/microns_per_pixel_after_downsampling)  NaN]);
    
    % apply the transformation
%     if ~isempty(transform_data.transform_points{1}) && ~isempty(transform_data.transform_points{2})
        ud.current_pointList_for_transform = transform_data.transform_points{1};
        ud_slice.pointList = transform_data.transform_points{2};
%     end           %if it is empty, I prefer to receive an error

    % load allen ref location
    slice_num = transform_data.allen_location{1};
    slice_angle = transform_data.allen_location{2};
    % create transformed histology image
    ud.ref = uint8(squeeze(tv_plot(slice_num,:,:)));
    R = imref2d(size(ud.ref));
    ud.curr_slice_trans = imwarp(im, transform_data.transform, 'OutputView',R);
    
    % do something with it...
    % eg plot it with atlas boundaries?? <-- too complicated for now.
    % Skipping
    % add
    % save it
   
    
%     if ud.currentAngle(1) == 0 && ud.currentAngle(2) == 0
%         curr_annotation = squeeze(av_plot(ud.currentSlice,:,:));
%     else
%         curr_annotation = ud.im_annotation; %check line 783 of AtlasTransformBrowser
%     end
%     % too complicated for now.
%     
%     atlas_vert_1 = double(curr_annotation(1:end-2,:));
%     atlas_vert_2 = double(curr_annotation(3:end,:));
%     atlas_vert_offset = abs( atlas_vert_1 - atlas_vert_2 ) > 0;
%     shifted_atlas_vert1 = zeros(size(curr_annotation(:,:)));
%     shifted_atlas_vert1(3:end,:) = atlas_vert_offset;
%     shifted_atlas_vert2 = zeros(size(curr_annotation(:,:)));
%     shifted_atlas_vert2(1:end-2,:) = atlas_vert_offset;
% 
%     atlas_horz_1 = double(curr_annotation(:,1:end-2));
%     atlas_horz_2 = double(curr_annotation(:,3:end));
%     atlas_horz_offset = abs( atlas_horz_1 - atlas_horz_2 )>0;
%     shifted_atlas_horz1 = zeros(size(curr_annotation(:,:)));
%     shifted_atlas_horz1(:,3:end) = atlas_horz_offset;
%     shifted_atlas_horz2 = zeros(size(curr_annotation(:,:)));
%     shifted_atlas_horz2(:,1:end-2) = atlas_horz_offset;
% 
%     shifted_atlas = shifted_atlas_horz1 + shifted_atlas_horz2 + shifted_atlas_vert1 + shifted_atlas_vert2;
% 
%     atlas_boundaries = (shifted_atlas>0); ud.atlas_boundaries = atlas_boundaries;
% 
%     if ud.showAtlas
%         image_blend =  uint8( imfuse(ud.curr_im, atlas_boundaries/3.5*(1+.35*isa(ud.curr_im,'uint16')),'blend','Scaling','none') )* 2;
%         set(ud.im, 'CData', image_blend); 
%     end
end

%% apply the same transormation to ROIs 

%change this one:
% trasformanda_tag = '.tif (green).png'; %adapt here
object_tag = 'green';
trasformanda_tag = '.tif (green)_dilate.png'; %adapt here - FIX THIS THROUGHOUT

roiLocation_cat_matfile = matfile(fullfile(folder_processed_images, sprintf('%s_roiLocations_All.mat',save_file_name)), 'Writable', true); %beware of duplicates
roiTable_cat_matfile = matfile(fullfile(folder_processed_images, sprintf('%s_roiTable_All.mat',save_file_name)), 'Writable', true); %beware of duplicates

transf_atlasreg_folder = fullfile(image_folder, 'processed/transformations');
trasformanda_folder = fullfile(image_folder, 'preprocessed');
d = dir(fullfile(transf_atlasreg_folder, '*processed_transform_data.mat'));

% load the reference brain annotations
if ~exist('av','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
end
if ~exist('st','var')
    disp('loading structure tree...')
    st = loadStructureTree(structure_tree_location);
end

for t = 1:length(d)
    % extract image root with ordinal tag
    fileroot = regexp(d(t).name, '_processed_transform_data.mat', 'split');
    fileroot = fileroot{1};

    % load the transformation data
    transfMAT = fullfile(transf_atlasreg_folder, d(t).name);
    trmat = load(transfMAT); %save_transform
    transform_data = trmat.save_transform;
    clear trmat
    ud.current_pointList_for_transform = transform_data.transform_points{1};
    ud_slice.pointList = transform_data.transform_points{2};
    % load allen ref location
    slice_num = transform_data.allen_location{1};
    slice_angle = transform_data.allen_location{2};
    
        
    % load the image or data to be transformed
    trasformanda_PNG = fullfile(trasformanda_folder, sprintf('%s%s', fileroot, trasformanda_tag)); % to be transformed by applying transformations
    im = imread(trasformanda_PNG);
    im2 = imread(fullfile(trasformanda_folder, sprintf('%s%s', fileroot, '.tif (green).png')));
    original_image_size = size(im);
    im = imresize(im, [round(original_image_size(1)*microns_per_pixel/microns_per_pixel_after_downsampling)  NaN]);
    im2 = imresize(im2, [round(original_image_size(1)*microns_per_pixel/microns_per_pixel_after_downsampling)  NaN]);
    
    % create transformed histology image
    ud.ref = uint8(squeeze(tv_plot(slice_num,:,:)));
    R = imref2d(size(ud.ref));
    ud.curr_slice_trans = imwarp(im, transform_data.transform, 'OutputView',R);
    transformed_slice_image = imwarp(im2, transform_data.transform, 'OutputView',R);
    
    rois = uint8(imregionalmax(ud.curr_slice_trans));
    
    % do something with it...
    
    
    figure; imshow(imfuse(rois, transformed_slice_image));
    title('transformed slice image, fused with ROIs')
    
    % make sure the rois are in a properly size image
    assert(size(rois,1)==800&size(rois,2)==1140&size(rois,3)==1,'roi image is not the right size');

    
    % initialize array of locations (AP, DV, ML relative to bregma) in reference space
    % and the correponding region annotations
    roi_location = zeros(sum(rois(:)>0),3);
    roi_annotation = cell(sum(rois(:)>0),3);
    
    % get location and annotation for every roi pixel
    [pixels_row, pixels_column] = find(rois>0);
    
    % generate other necessary values
    bregma = allenCCFbregma(); % bregma position in reference data space
    atlas_resolution = 0.010; % mm
    offset_map = get_offset_map(slice_angle);
    
    
    % loop through every pixel to get ROI locations and region annotations
    for pixel = 1:length(pixels_row)
        
        % get the offset from the AP value at the centre of the slice, due to
        % off-from-coronal angling
        offset = offset_map(pixels_row(pixel),pixels_column(pixel));
        
        % use this and the slice number to get the AP, DV, and ML coordinates
        ap = -(slice_num-bregma(1)+offset)*atlas_resolution;
        dv = (pixels_row(pixel)-bregma(2))*atlas_resolution;
        ml = (pixels_column(pixel)-bregma(3))*atlas_resolution;
        
        roi_location(pixel,:) = [ap dv ml];                                             % this is the one!!
        
        % finally, find the annotation, name, and acronym of the current ROI pixel
        ann = av(slice_num+offset,pixels_row(pixel),pixels_column(pixel));
        name = st.safe_name{ann};
        acr = st.acronym{ann};
        
        roi_annotation{pixel,1} = ann;
        roi_annotation{pixel,2} = name;
        roi_annotation{pixel,3} = acr;
        
    end
    
    roi_table = table(roi_annotation(:,2),roi_annotation(:,3), ...
        roi_location(:,1),roi_location(:,2),roi_location(:,3), roi_annotation(:,1), ...
        'VariableNames', {'name', 'acronym', 'AP_location', 'DV_location', 'ML_location', 'avIndex'});
    
    disp(roi_table)
    
    save(fullfile(folder_processed_images, sprintf('%s_%s_rois.mat',fileroot, object_tag)), 'roi_table', 'roi_location')
    if t == 1
        roiLocation_cat_matfile.roi_location = roi_location;
    else
        roiLocation_cat_matfile.roi_location = cat(1, roiLocation_cat_matfile.roi_location, roi_location);
    end
    % eg plot it with atlas boundaries?? 
    % add
    % save it
   
end
%% make the ROI image from centroids of detected cells -- does not work
% 
%change this one:
trasformanda_tag = 'green'; %adapt here


transf_atlasreg_folder = fullfile(image_folder, 'processed/transformations');
trasformanda_folder = fullfile(image_folder, 'preprocessed');
d = dir(fullfile(transf_atlasreg_folder, '*processed_transform_data.mat'));

excelFileName = fullfile(trasformanda_folder, sprintf('ShrinkToObjectCenters_%s.csv', trasformanda_tag));
T = readtable(excelFileName);
fnameField = sprintf('FileName_%s', trasformanda_tag);
filenames = unique(T.(fnameField), 'stable');
imageNumbers = unique(T.ImageNumber, 'stable');

for t = 1:length(d)
    % extract image root with ordinal tag
    fileroot = regexp(d(t).name, '_processed_transform_data.mat', 'split');
    fileroot = fileroot{1};

    % load the transformation data
    transfMAT = fullfile(transf_atlasreg_folder, d(t).name);
    st = load(transfMAT); %save_transform
    transform_data = st.save_transform;
    clear st
    ud.current_pointList_for_transform = transform_data.transform_points{1};
    ud_slice.pointList = transform_data.transform_points{2};
    % load allen ref location
    ud.currentSlice = transform_data.allen_location{1}; 
    ud.currentAngle = transform_data.allen_location{2};
    
        
    % load the data to be transformed
    trasformanda_imgN = imageNumbers(contains(filenames, fileroot)); %this is valid across colors. I will use it in the green channel. 
    x = T.Location_Center_X(T.ImageNumber==trasformanda_imgN);
    y = T.Location_Center_Y(T.ImageNumber==trasformanda_imgN);
    
    %downsample them
    xDS = round((x+1) * microns_per_pixel/microns_per_pixel_after_downsampling);
    yDS = round((y+1) * microns_per_pixel/microns_per_pixel_after_downsampling);
    A = [xDS, yDS, zeros(size(x))];
    
    % apply the transformation
    At = transform_data.transform.T*A';
    
    % load images
    trasformanda_PNG = fullfile(trasformanda_folder, sprintf('%s.tif (%s).png', fileroot, trasformanda_tag));
    im = imread(trasformanda_PNG);
    original_image_size = size(im);
    imDS = imresize(im, [round(original_image_size(1)*microns_per_pixel/microns_per_pixel_after_downsampling)  NaN]);
    % create transformed histology image
    ud.ref = uint8(squeeze(tv_plot(ud.currentSlice,:,:)));
    R = imref2d(size(ud.ref));
    ud.curr_slice_trans = imwarp(imDS, transform_data.transform, 'OutputView',R);
    
    % plot centroids in original format
    figure;
    imshow(im, []);  hold on
    scatter(x+1, y+1, 4, [0 1 0], 'filled');
    
    % plot centroids in downsampled format
    figure; hold on
    imshow(imDS, []);
    scatter(xDS, yDS, 4, [0 1 0], 'filled'); %not very precise already, but OK...
    
    % plot centroids in registered format?
    figure; hold on
    imshow(ud.curr_slice_trans, []);
    scatter(At(1,:), At(2,:), 4, [0 1 0], 'filled'); %nope
    
end

%% draft plot points on fwireframe


ProbeColors = .75*[1.3 1.3 1.3; 1 .75 0;  .3 1 1; .4 .6 .2; 1 .35 .65; .7 .7 .9; .65 .4 .25; .7 .95 .3; .7 0 0; .6 0 .7; 1 .6 0]; 


% show a table of regions that the probe goes through, in the console
show_region_table = true;      
% black brain?
black_brain = true;
fwireframe = [];
% create a new figure with wireframe
fwireframe = plotBrainGrid([], [], fwireframe, black_brain);
hold on; 
fwireframe.InvertHardcopy = 'off';

% for selected_probe = probes
%     
% % get the probe points for the currently analyzed probe 
% if strcmp(plane,'coronal')
%     curr_probePoints = probePoints.pointList.pointList{selected_probe,1}(:, [3 2 1]);
% elseif strcmp(plane,'sagittal')
%     curr_probePoints = probePoints.pointList.pointList{selected_probe,1}(:, [1 2 3]);
% elseif strcmp(plane,'transverse')
%     curr_probePoints = probePoints.pointList.pointList{selected_probe,1}(:, [1 3 2]);
% end
curr_probePoints = roiLocation_cat_matfile.roi_location;
selected_probe = 1; % this should be a single brain...

% focus on wireframe plot
figure(fwireframe);

% plot probe points
hp = plot3(curr_probePoints(:,1), curr_probePoints(:,3), curr_probePoints(:,2), '.','linewidth',2, 'color',[ProbeColors(selected_probe,:) .2],'markers',10);

% end