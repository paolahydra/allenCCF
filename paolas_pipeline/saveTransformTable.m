% display which slices have already been registered, and their coordinates
% suggest coordinates for next registration
function T = saveTransformTable(folder_transformations, image_file_names, reference_size)
% folder_transformations = fullfile(folder_processed_images, 'transformations');
% reference_size =  [ 1320         800        1140 ];   %check where this comes from
reference_size = reference_size(2:3);   %[800        1140 ];
bregma = allenCCFbregma();

d = dir (fullfile(folder_transformations, '*_transform_data.mat'));
transfs = {d.name};
T = table;
%build it first with all the data available.
warning('off')
for t = 1:length(image_file_names)
    T.sliceNum(t) = t;
    T.sliceName(t) = {image_file_names{t}(1:end-4)};
    if sum(contains(transfs, T.sliceName(t)))
        TR = load(fullfile(folder_transformations, transfs{contains(transfs, T.sliceName(t))}));
%         T.allenSlice(t) = TR.save_transform.allen_location{1}; %bregma(1) - TR.save_transform.allen_location{1};
%         T.allenAngle_DV(t) = TR.save_transform.allen_location{2}(1); %round(atand(TR.save_transform.allen_location{2}(1)/(reference_size(1)/2)));
%         T.allenAngle_ML(t) = TR.save_transform.allen_location{2}(2); %round(atand(TR.save_transform.allen_location{2}(2)/(reference_size(2)/2)));
        T.allenSlice(t) = bregma(1) - TR.save_transform.allen_location{1};
        T.suggestedAllenSlice(t) = NaN;
        T.d_slice(t) = NaN;
        T.DV_deg(t) = round(atand(TR.save_transform.allen_location{2}(1)/(reference_size(1)/2))*10)/10;
        T.ML_deg(t) = round(atand(TR.save_transform.allen_location{2}(2)/(reference_size(2)/2))*10)/10;
        T.class(t) = {class(TR.save_transform.transform)};
    else
        T.allenSlice(t) = NaN;
        T.suggestedAllenSlice(t) = NaN;
        T.d_slice(t) = NaN;
        T.DV_deg(t) = NaN;
        T.ML_deg(t) = NaN;
        T.class(t) = {''};
    end
end

%fill in the blanks with estimates
asl = T.allenSlice;
if sum(~isnan(asl))>=2
    T.suggestedAllenSlice = round(interp1(find(~isnan(asl)), asl(~isnan(asl)), 1:length(asl), 'linear', 'extrap'))';
    T.d_slice = cat(1, NaN(1,1), diff(T.suggestedAllenSlice));
    asl = T.DV_deg;
    T.DV_deg = (interp1(find(~isnan(asl)), asl(~isnan(asl)), 1:length(asl), 'nearest', 'extrap'))';
    asl = T.ML_deg;
    T.ML_deg = (interp1(find(~isnan(asl)), asl(~isnan(asl)), 1:length(asl), 'nearest', 'extrap'))';
end

writetable(T, fullfile(folder_transformations, 'dataTable_transformations.csv'))
end


%% notes on how to convert allen_location in actual coordinates and degrees
% %atlasRes = 0.01; %ndp

% if strcmp(plane,'coronal')
%     ap = -(currentSlice-bregma(1))*atlasRes;
%     dv = (pixel(1)-bregma(2))*atlasRes;
%     ml = (pixel(2)-bregma(3))*atlasRes;
%     set(angleText, 'String', ['Slice ' num2str(bregma(1) - slice_num) ', DV angle ' num2str(round(atand(ap_angle/(ref_size(1)/2)),1)) '^{\circ}, ML angle ' num2str(round(atand(ml_angle/(ref_size(2)/2)),1)) '^{\circ}']);



% %%
% function updateStereotaxCoords(currentSlice, pixel, bregma, bregmaText, angleText, slice_num, ap_angle, ml_angle, ref_size, plane)
% atlasRes = 0.010; % mm
% if strcmp(plane,'coronal')
%     ap = -(currentSlice-bregma(1))*atlasRes;
%     dv = (pixel(1)-bregma(2))*atlasRes;
%     ml = (pixel(2)-bregma(3))*atlasRes;
%     set(angleText, 'String', ['Slice ' num2str(bregma(1) - slice_num) ', DV angle ' num2str(round(atand(ap_angle/(ref_size(1)/2)),1)) '^{\circ}, ML angle ' num2str(round(atand(ml_angle/(ref_size(2)/2)),1)) '^{\circ}']);
% elseif strcmp(plane,'sagittal')
%     ap = -(pixel(2)-bregma(1))*atlasRes;
%     dv = (pixel(1)-bregma(2))*atlasRes;
%     ml = -(currentSlice-bregma(3))*atlasRes;
%     set(angleText, 'String', ['Slice ' num2str(bregma(1) - slice_num) ', DV angle ' num2str(round(atand(ap_angle/(ref_size(1)/2)),1)) '^{\circ}, AP angle ' num2str(round(atand(ml_angle/(ref_size(2)/2)),1)) '^{\circ}']);
% elseif strcmp(plane,'transverse')
%     ap = -(pixel(2)-bregma(1))*atlasRes;
%     dv = (currentSlice-bregma(2))*atlasRes;
%     ml = -(pixel(1)-bregma(3))*atlasRes;
%     set(angleText, 'String', ['Slice ' num2str(bregma(1) - slice_num) ', ML angle ' num2str(round(atand(ap_angle/(ref_size(1)/2)),1)) '^{\circ}, AP angle ' num2str(round(atand(ml_angle/(ref_size(2)/2)),1)) '^{\circ}']);
% end
% set(bregmaText, 'String', sprintf('%.2f AP, %.2f DV, %.2f ML', ap, dv, ml));