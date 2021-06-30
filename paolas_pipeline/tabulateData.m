%% When done, tabulate the ROI data (from cellprofiler analysis and registration)
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

T_roi = Register_and_Tabulate_Rois(object_tag, save_folder, image_tag, av, st, tv, microns_per_pixel, microns_per_pixel_after_downsampling, reference_size);
