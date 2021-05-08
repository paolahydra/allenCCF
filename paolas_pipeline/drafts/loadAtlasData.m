%% load basic stuff
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