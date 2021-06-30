%% DEMO: show a single coronal slice with cells (vectorial)



%% define a specific region of interest to highlight at the AP level considered (if any)
% ROI = find(strcmp(st.name, 'Facial motor nucleus'));
% ROI = find(strcmp(st.name, 'Parvicellular reticular nucleus'));


% % or leave it empty:
 ROI = [];


%% define the AP level of interest

% % retrieve slice info from a saved transformation, with offset

% trdata = load('/Users/galileo/dati/registered_brains_completed/992234/startingSingleSlices/processed/transformations/mouse_992234_05.018_processed_transform_data.mat');
% slice_num = trdata.save_transform.allen_location{1};
% slice_angle = trdata.save_transform.allen_location{2};
% offset_map = get_offset_map(slice_angle, reference_size);
% avSliceNum = trdata.save_transform.allen_location{1};
% avSlice = squeeze(offset_map(avSliceNum, :, :));
% figure; imshow(avSlice, [])


% % retrieve AP level info from the AtlasViewer or from saved coordinates
% (in mm) -  no ML, no DV offsets here - registered straightned data 

AP_mm = -5.34;  %-1.9; % mm: you get this number by hovering over your region of interest with the mouse
showVolume = 0.250; %thickness of pseudo-slice to show in mm

%% convert the AP coordinate in atlas slice number (do not change)
bregma = allenCCFbregma();
atlas_resolution = 0.010; % mm
avSliceNum = bregma(1) - AP_mm/atlas_resolution; 

avSlice = squeeze(av(avSliceNum, :, :));
% figure; imshow(avSlice, [])
% imwrite(avSlice, 'atlasBoundaries_7Nlevel.png')


%%  make a figure
figure;
axATL = axes;
[coords, coordsReg, h] = sliceOutlineWithRegionVec(avSlice, ROI, [1,0,0], axATL);
% export_fig('atlasBoundaries_7Nlevel.pdf')
hold on


%% define and add the cells!
for i = 1:length(S)
    if i == 1 && exist('p', 'var')
        delete(p)
    end
    
    S(i).pltIdx = S(i).T_roi.avIndex ~= 1 ...
        & S(i).T_roi.AP_location >= AP_mm - showVolume/2 ...
        & S(i).T_roi.AP_location < AP_mm + showVolume/2;
    
    % transform coordinates
    ap_pixel = bregma(1) - S(i).T_roi.AP_location(S(i).pltIdx)./atlas_resolution; %OK
    ml_pixel = bregma(3) + S(i).T_roi.ML_location(S(i).pltIdx)./atlas_resolution; %OK
    dv_pixel = bregma(2) + S(i).T_roi.DV_location(S(i).pltIdx)./atlas_resolution; %OK
    
    p(i) = plot( ml_pixel, dv_pixel, '.','linewidth',2, 'color', S(i).braincolor, 'markers',10);
end



%% save the figure in vectorial format (.pdf or .eps)
%export_fig(sprintf('atlasBoundaries_bregma%1.2f.pdf', AP_mm))
fig = gcf;
filename=sprintf('2_atlasBoundaries_bregma%1.2f.pdf', AP_mm);
saveas(fig,filename) 

