% ------------------------------------------------------------------------
%          Run Allen Atlas Browser
% ------------------------------------------------------------------------


%% ENTER FILE LOCATION AND PROBE-SAVE-NAME
% directory of histology
if ~exist('folder_processed_images','var')
    folder_processed_images = fullfile(image_folder, 'processed');
end

% name the saved probe points, to avoid overwriting another set of probes going in the same folder
probe_save_name_suffix = ''; 






%% GET PROBE TRAJECTORY POINTS

% load the reference brain and region annotations
if ~exist('av','var') || ~exist('st','var') || ~exist('tv','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
    tv = readNPY(template_volume_location);
end

% select the plane for the viewer
if strcmp(plane,'coronal')
    av_plot = av;
    tv_plot = tv;
elseif strcmp(plane,'sagittal')
    av_plot = permute(av,[3 2 1]);
    tv_plot = permute(tv,[3 2 1]);
elseif strcmp(plane,'transverse')
    av_plot = permute(av,[2 3 1]);
    tv_plot = permute(tv,[2 3 1]);
end

% create Atlas viewer figure
f = figure('Name','Atlas Viewer'); 

% show histology in Slice Viewer
try; figure(slice_figure_browser); title('');
catch; slice_figure_browser = figure('Name','Slice Viewer'); end
reference_size = size(tv);
sliceBrowser(slice_figure_browser, folder_processed_images, f, reference_size);


% use application in Atlas Transform Viewer
% use this function if you have a processed_images_folder with appropriately processed .tif histology images
f = AtlasTransformBrowser(f, tv_plot, av_plot, st, slice_figure_browser, folder_processed_images, probe_save_name_suffix, plane, transformationType); 

% use the simpler version, which does not interface with processed slice images
% just run these two lines instead of the previous 5 lines of code
% 
% save_location = processed_images_folder;
% f = allenAtlasBrowser(tv_plot, av_plot, st, save_location, probe_save_name_suffix);

