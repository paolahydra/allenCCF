function T_roi = Register_and_Tabulate_Rois(object_tag, image_folder, save_file_name, av, st, tv_plot, microns_per_pixel, microns_per_pixel_after_downsampling, reference_size )
% reapply the same transormations to ROIs to register detected cells
% it assumes it gets path definitions from the generalPipeline script
%
% %change this one as needed:
% object_tag = 'green';
% trasformanda_tag = '.tif (green)_dilate.png';

%% paths
transformanda_tag = sprintf('(%s)_dilate.png', object_tag); % if you get an error here check your filenames

transformanda_folder = fullfile(image_folder, 'preprocessed');
folder_processed_images = fullfile(image_folder, 'processed');
transf_atlasreg_folder = fullfile(image_folder, 'processed/transformations');
d = dir(fullfile(transf_atlasreg_folder, '*processed_transform_data.mat'));

roiTable_name = fullfile(folder_processed_images, sprintf('%s%s_roiTable_All.csv',save_file_name, object_tag));
if exist(roiTable_name, 'file')
    warning('roi_table file already exists and it will be overwritten')
end


% % load the reference brain annotations
% if ~exist('av','var')
%     disp('loading reference atlas...')
%     av = readNPY(annotation_volume_location);
% end
% if ~exist('st','var')
%     disp('loading structure tree...')
%     st = loadStructureTree(structure_tree_location);
% end

%%
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
    transformanda_PNG = dir(fullfile(transformanda_folder, sprintf('%s*%s', fileroot, transformanda_tag)));
    if ~isempty(transformanda_PNG)
        transformanda_PNG = fullfile(transformanda_folder, transformanda_PNG.name); % to be transformed by applying transformations
    else
        warning(sprintf('cell detection probably not applied to %s. Skipping...', fileroot))
        continue
    end
    
    im = imread(transformanda_PNG);
%     im2 = imread(fullfile(trasformanda_folder, sprintf('%s.tif (%s).png', fileroot, object_tag))); % just for checking overlay
    original_image_size = size(im);
    im = imresize(im, [round(original_image_size(1)*microns_per_pixel/microns_per_pixel_after_downsampling)  NaN]);
%     im2 = imresize(im2, [round(original_image_size(1)*microns_per_pixel/microns_per_pixel_after_downsampling)  NaN]);
    
    % create transformed histology image
    ud.ref = uint8(squeeze(tv_plot(slice_num,:,:)));
    R = imref2d(size(ud.ref));
    ud.curr_slice_trans = imwarp(im, transform_data.transform, 'OutputView',R);
%     transformed_slice_image = imwarp(im2, transform_data.transform, 'OutputView',R);
    
    if length(unique(ud.curr_slice_trans(:))) > 1 % at least one object was detected
        rois = uint8(imregionalmax(ud.curr_slice_trans));
    else
        rois = uint8(zeros(size((ud.curr_slice_trans))));
    end
    
    % do something with it...
%     figure; imshow(imfuse(rois, transformed_slice_image));
%     title('transformed slice image, fused with ROIs')
    % save images?
    
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
    offset_map = get_offset_map(slice_angle, reference_size);
    
    
    % loop through every pixel to get ROI locations and region annotations
    for pixel = 1:length(pixels_row) 
        % get the offset from the AP value at the centre of the slice, due to
        % off-from-coronal angling
        offset = offset_map(pixels_row(pixel),pixels_column(pixel));
        
        % use this and the slice number to get the AP, DV, and ML coordinates
        ap = -(slice_num-bregma(1)+offset)*atlas_resolution;
        dv = (pixels_row(pixel)-bregma(2))*atlas_resolution;
        ml = (pixels_column(pixel)-bregma(3))*atlas_resolution;
        
        roi_location(pixel,:) = [ap dv ml];   
        
        % finally, find the annotation, name, and acronym of the current ROI pixel
        ann = av(slice_num+offset,pixels_row(pixel),pixels_column(pixel));
        name = st.safe_name{ann};
        acr = st.acronym{ann};
        
        roi_annotation{pixel,1} = ann;
        roi_annotation{pixel,2} = name;
        roi_annotation{pixel,3} = acr;
    end
    
    filename_origin = repmat({sprintf('%s%s', fileroot, transformanda_tag)}, length(pixels_row),1);
    roi_table = table(roi_annotation(:,2),roi_annotation(:,3), ...
        roi_location(:,1),roi_location(:,2),roi_location(:,3), cat(1, roi_annotation{:,1}), filename_origin, ...
        'VariableNames', {'name', 'acronym', 'AP_location', 'DV_location', 'ML_location', 'avIndex', 'roiFIle'});
    
    if t == 1
        %         writetable(roi_table, roiTable_name, 'WriteVariableNames', true, 'WriteMode', 'overwrite') %only available matlab 2021
        writetable(roi_table, roiTable_name, 'WriteVariableNames', true)
    else
        T_roi = readtable(roiTable_name);
        T_roi = cat(1, T_roi, roi_table);
        writetable(T_roi, roiTable_name, 'WriteVariableNames', true)
        %         writetable(roi_table, roiTable_name, 'WriteVariableNames', false, 'WriteMode', 'append') %only available matlab 2021
    end
    
end


end
