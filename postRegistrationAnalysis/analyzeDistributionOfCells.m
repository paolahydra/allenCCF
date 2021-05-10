st = loadStructureTree(structure_tree_location);
% avIndex IS THE ROW INDEX IN THIS TABLE! (NOT the field "index").

% 622	'Bed nuclei of the stria terminalis'
% 1260	'stria terminalis'
% 572	'Striatum'

% 598	'Central amygdalar nucleus lateral part'
% 599	'Central amygdalar nucleus medial part'
% 597	'Central amygdalar nucleus capsular part'


% 801   'Subthalamic nucleus'
% 798   'Parasubthalamic nucleus'
% 1222  'nigrostriatal tract'


% 795	'Lateral hypothalamic area'
% 716	'Hypothalamus'
% 793	'Posterior hypothalamic nucleus'
% 721	'Paraventricular hypothalamic nucleus'

% 803	'Zona incerta'
% 805	'Fields of Forel'
% 646	'Ventral medial nucleus of the thalamus'


% 823	'Substantia nigra reticular part'
% 867	'Substantia nigra compact part'

% 826	'Midbrain reticular nucleus retrorubral area'
% 827	'Midbrain reticular nucleus'
% 807	'Midbrain'

% 853	'Cuneiform nucleus'
% 868	'Pedunculopontine nucleus'


%% load data to analyze
folders2Include = uipickfiles('FilterSpec', '/Users/galileo/dati/registered_brains_completed');
plotFWire = 0;
object_tag = 'green'; %'green' for rabies cells

S = loadTabDataFromMultipleBrains(folders2Include, plotFWire, object_tag); %updated plotting to exclude point that are in 'root' (avIndex = 1)

TH = loadTabDataFromMultipleBrains(folders2Include, plotFWire, 'red'); %updated plotting to exclude point that are in 'root' (avIndex = 1)

folder2save = '/Users/galileo/dati/registered_brains_completed/figures_993030_993031_2nd';
mkdir(folder2save)
cd(folder2save)

%% find unique targets and sort them
for i = 1:length(S)
    
    all_Names = S(i).T_roi.name;
    all_Files = S(i).T_roi.roiFIle;
    all_acIndex = S(i).T_roi.avIndex;
    all_isPositiveML_location = S(i).T_roi.ML_location;
    all_isPositiveML_location = all_isPositiveML_location>0; %this should be on the right
    
    TargetAreas_R = unique(all_acIndex(all_isPositiveML_location));
    for n = 1:length(TargetAreas_R)
        TargetNames_R{n,1} = all_Names{find(all_acIndex==TargetAreas_R(n), 1)};
        % count how many registered slices contributed to any given area count
        TargetSliceCount_R{n,1} = length(unique(all_Files(all_acIndex==TargetAreas_R(n) & all_isPositiveML_location)));
    end
    
    TargetAreas_L = unique(all_acIndex(~all_isPositiveML_location));
    for n = 1:length(TargetAreas_L)
        TargetNames_L{n,1} = all_Names{find(all_acIndex==TargetAreas_L(n), 1)};
        % count how many registered slices contributed to any given area count
        TargetSliceCount_L{n,1} = length(unique(all_Files(all_acIndex==TargetAreas_L(n) & ~all_isPositiveML_location)));
    end
    
    % 
    tbl = tabulate(all_acIndex(all_isPositiveML_location));
    tbl(tbl(:,2)==0,:) = [];
    [~, i_sorting] = sort(tbl(:,2), 'descend');
    tbl = tbl(i_sorting,:);
    TargetNames_R = TargetNames_R(i_sorting, :);
    TargetSliceCount_R = TargetSliceCount_R(i_sorting, :);
    S(i).cellCounts_right = table(tbl(:,1), TargetNames_R, tbl(:,2), TargetSliceCount_R, 'VariableNames', {'avIndex','name','count','nSlicesDenominator'});
    
    tbl = tabulate(all_acIndex(~all_isPositiveML_location));
    tbl(tbl(:,2)==0,:) = [];
    [~, i_sorting] = sort(tbl(:,2), 'descend');
    tbl = tbl(i_sorting,:);
    TargetNames_L = TargetNames_L(i_sorting, :);
    TargetSliceCount_L = TargetSliceCount_L(i_sorting, :);
    S(i).cellCounts_left = table(tbl(:,1), TargetNames_L, tbl(:,2), TargetSliceCount_L, 'VariableNames', {'avIndex','name','count','nSlicesDenominator'});

end

% for a given region of interest that I want to compare, I need to reload
% each registered slice that contained that area, mask the area, and sum
% the extent of pixels across all slices (as a proxy of the volume
% considered).
%
% for now, the number of slices will do to have a very coarse sense (but
% it's actuall not even good because I would need the number of slices that
% intersect that region of interest, not the number of slices that have at
% least one cell

%% isolate a specific region of interest
% useful bits in:
% edit script_sliceMovie
atlas_resolution = 0.010; % mm

SNR = 823;  % SNR index
SNC = 867;  % SNC
RR = 826;   % 'Midbrain reticular nucleus retrorubral area'
MB = find(strcmp(st.name,'Midbrain'));
MBRN = find(strcmp(st.name,'Midbrain reticular nucleus'));
VTA = find(strcmp(st.name,'Ventral tegmental area'));

CUN = find(strcmp(st.name,'Cuneiform nucleus'));
PPN = find(strcmp(st.name,'Pedunculopontine nucleus'));

%% Create the plots
szAP = size(tv,1);
szDV = size(tv,2);
szLR = size(tv,3);

fh = figure; set(fh, 'Position', [-1896         -72        1479        1024]);
figcol = 0.7;
fh.Color = [figcol figcol figcol];
% subtightplot(2,2,1); 
bregma = allenCCFbregma();
isROI = av==SNR | av==RR; % >0 for original av, >1 for by_index
gridIn3D(double(isROI), 0.25, 10, bregma, 'w'); %ndp: function contourHands = gridIn3D(volData, contourHeight, sliceSpacing, origin, contourColor)
axis vis3d
set(gca, 'ZDir', 'reverse')
axis equal
axis off
% view([-30    25]);
view([-90, 8]);
hold on

ax = gca;
isthmus_grid_handles = ax.Children; 

totHandles = length(ax.Children);
isROI = av==1; % >0 for original av, >1 for by_index
gridIn3D(double(isROI), 0.25, 50, bregma, [0 0 0]); %ndp: function contourHands = gridIn3D(volData, contourHeight, sliceSpacing, origin, contourColor)
handles.root = ax.Children(1:end-totHandles);
%%
makeVideo = 0;
if makeVideo
    WriterObj = VideoWriter('allenCCF_isthmus_movie_Root.mp4', 'MPEG-4');
    WriterObj.FrameRate=30;
    open(WriterObj);
    for f = -90:270
        view([f, 8])
        drawnow;
        if makeVideo
            frame = getframe(fh);
            writeVideo(WriterObj,frame);
        end
    end
    close(WriterObj);
end
export_fig('SNR_RRF_outline_withRoot.png', '-nocrop', '-png', '-m5')
delete(handles.root)
export_fig('SNR_RRF_outline.png', '-nocrop', '-png', '-m5')
%% add extra regions
totHandles = length(ax.Children);
isROI = av==MBRN; % >0 for original av, >1 for by_index
gridIn3D(double(isROI), 0.25, 10, bregma, [0.85 0.4 0.85]); %ndp: function contourHands = gridIn3D(volData, contourHeight, sliceSpacing, origin, contourColor)
handles.MBRN = ax.Children(1:end-totHandles);

for h = 1:length(handles.MBRN)
    handles.MBRN(h).Color = [0.85 0.4 0.85];
end
for h = 1:length(handles.MBRN)
    handles.MBRN(h).Visible = 'off';
end
%%
totHandles = length(ax.Children);
isROI = av==CUN; % >0 for original av, >1 for by_index
gridIn3D(double(isROI), 0.25, 10, bregma, [0.4 0.8 0.8]); %ndp: function contourHands = gridIn3D(volData, contourHeight, sliceSpacing, origin, contourColor)
handles.CUN = ax.Children(1:end-totHandles);

for h = 1:length(handles.CUN)
    handles.CUN(h).Color = [0.4 0.8 0.8];
end
for h = 1:length(handles.CUN)
    handles.CUN(h).Visible = 'off';
end
%%
totHandles = length(ax.Children);
isROI = av==PPN; % >0 for original av, >1 for by_index
gridIn3D(double(isROI), 0.25, 10, bregma, [0.5 1 1]); %ndp: function contourHands = gridIn3D(volData, contourHeight, sliceSpacing, origin, contourColor)
handles.PPN = ax.Children(1:end-totHandles);

for h = 1:length(handles.PPN)
    handles.PPN(h).Visible = 'off';
end
%%
totHandles = length(ax.Children);
isROI = av==MB; % >0 for original av, >1 for by_index
gridIn3D(double(isROI), 0.25, 10, bregma, [0 0 0]); %ndp: function contourHands = gridIn3D(volData, contourHeight, sliceSpacing, origin, contourColor)
handles.MB = ax.Children(1:end-totHandles);

for h = 1:length(handles.MB)
    handles.MB(h).Visible = 'off';
end


%% make a complex video - does not work. Making more than one
% start with all handles off, except for SNR-RR
makeVideo = 0;
if makeVideo   
    F = -30;
    view([F, 8])
    drawnow;
    export_fig('SNR_RRF_outline.png', '-nocrop', '-png', '-m5')
    
    %% PPN
    for h = 1:length(handles.PPN)
        handles.PPN(h).Visible = 'on';
    end
    WriterObj = VideoWriter('allenCCF_allROIs_movie.mp4', 'MPEG-4');
    WriterObj.FrameRate=30;
    open(WriterObj);
    for f = -30:330
        view([f, 8])
        drawnow;
        frame = getframe(fh);
        writeVideo(WriterObj,frame);
    end
    close(WriterObj);
    
    %% CUN
    for h = 1:length(handles.CUN)
        handles.CUN(h).Visible = 'on';
    end
    WriterObj = VideoWriter('allenCCF_allROIs_CUN_movie.mp4', 'MPEG-4');
    WriterObj.FrameRate=30;
    open(WriterObj);
    for f = -30:330
        view([f, 8])
        drawnow;
        frame = getframe(fh);
        writeVideo(WriterObj,frame);
    end
    close(WriterObj);
    
    %% MBRN
    for h = 1:length(handles.MBRN)
        handles.MBRN(h).Visible = 'on';
    end
    WriterObj = VideoWriter('allenCCF_allROIs_MBRN_movie.mp4', 'MPEG-4');
    WriterObj.FrameRate=30;
    open(WriterObj);
    for f = -30:330
        view([f, 8])
        drawnow;
        frame = getframe(fh);
        writeVideo(WriterObj,frame);
    end
    close(WriterObj);
    
    %% MB
    for h = 1:length(handles.MB)
        handles.MB(h).Visible = 'on';
    end
    WriterObj = VideoWriter('allenCCF_allROIs_MB_movie.mp4', 'MPEG-4');
    WriterObj.FrameRate=30;
    open(WriterObj);
    for f = -30:-1:-90
        view([f, 8])
        drawnow;
        frame = getframe(fh);
        writeVideo(WriterObj,frame);
    end
    for f = 8:36
        view([-90, f])
        drawnow;
        frame = getframe(fh);
        writeVideo(WriterObj,frame);
    end
    close(WriterObj);
    
end



%% add TH+ cells
view([-90, 8]);
CT = cbrewer2('Reds', 7, 'cubic');
CT = flipud(CT);

for i = 1:length(S)
    % first add only SNR and SNC cells, consider both sides (there are
    % sometimes very few cells in the contra side too)
    pltIdx = TH(i).T_roi.avIndex == SNR ...
        | TH(i).T_roi.avIndex == SNC ...
        | TH(i).T_roi.avIndex == RR ...
        | TH(i).T_roi.avIndex == MB ...
        | TH(i).T_roi.avIndex == MBRN ...
        | TH(i).T_roi.avIndex == VTA ...
        ; % & TH(i).T_roi.ML_location <= 0;         %use brackets if you want to apply this to all regions and not just the last one

    
    % transform coordinates
    ap_pixel = bregma(1) - TH(i).T_roi.AP_location(pltIdx)./atlas_resolution; %OK
    ml_pixel = bregma(3) + TH(i).T_roi.ML_location(pltIdx)./atlas_resolution; %OK
    dv_pixel = bregma(2) + TH(i).T_roi.DV_location(pltIdx)./atlas_resolution; %OK
    
    p_th(i) = plot3(ap_pixel, ml_pixel, dv_pixel, '.','linewidth',2, 'color', CT(i,:), 'markers',3);
    
end

for h = 1:length(p_th)
    p_th(h).Visible = 'on';
end

makeVideo = 0;
if makeVideo
    WriterObj = VideoWriter('allenCCF_isthmus_movie_TH.mp4', 'MPEG-4');
    WriterObj.FrameRate=30;
    open(WriterObj);
    for f = -90:270
        view([f, 8])
        drawnow;
        if makeVideo
            frame = getframe(fh);
            writeVideo(WriterObj,frame);
        end
    end
    close(WriterObj);
end

%% define and add the cells!
for i = 1:length(S)
    % first add only SNR and SNC cells, consider both sides (there are
    % sometimes very few cells in the contra side too)
    delete(p)
    S(i).pltIdx = S(i).T_roi.avIndex == SNR ...
        | S(i).T_roi.avIndex == SNC ...
        | S(i).T_roi.avIndex == RR ...
        | S(i).T_roi.avIndex == MB ...
        | S(i).T_roi.avIndex == MBRN; ...
%         | S(i).T_roi.avIndex == PPN  ...
%         | S(i).T_roi.avIndex == CUN  ...
        ; %& S(i).T_roi.ML_location <= 0 ...
%         ; % & S(i).T_roi.ML_location <= 0;
    
    % transform coordinates
    ap_pixel = bregma(1) - S(i).T_roi.AP_location(S(i).pltIdx)./atlas_resolution; %OK
    ml_pixel = bregma(3) + S(i).T_roi.ML_location(S(i).pltIdx)./atlas_resolution; %OK
    dv_pixel = bregma(2) + S(i).T_roi.DV_location(S(i).pltIdx)./atlas_resolution; %OK
    
    p(i) = plot3(ap_pixel, ml_pixel, dv_pixel, '.','linewidth',2, 'color', S(i).braincolor, 'markers',12);
    
end
export_fig('SNR_RRF_outline_withCells.png', '-nocrop', '-png', '-m5')
makeVideo = 0;
if makeVideo
    view([-90, 8]);
    WriterObj = VideoWriter('allenCCF_isthmus_movie_withAllCells.mp4', 'MPEG-4');
    WriterObj.FrameRate=30;
    open(WriterObj);
    for f = 8:20
        view([-90 f])
        drawnow;
        if makeVideo
            frame = getframe(fh);
            writeVideo(WriterObj,frame);
        end
    end
    for f = -90:8
        view([f 20])
        drawnow;
        if makeVideo
            frame = getframe(fh);
            writeVideo(WriterObj,frame);
        end
    end
    close(WriterObj);
end
%% set some manual boundaries to the region of interest
folder_pointlist = '/Users/galileo/dati/registered_brains_completed/993030/startingSingleSlices/processed/';
load(fullfile(folder_pointlist, 'probe_points.mat'), 'pointList');
object_num = 3; % second row is first attempt with only 4 points for fitting a plane. It was too ventral

% pointList.pointList{2,1}; %the points
% pointList.pointList{2,3}; % the registerd slices they come from

bregma = allenCCFbregma(); % bregma position in reference data space
atlas_resolution = 0.010; % mm
 
curr_objectPoints = pointList.pointList{object_num,1}(:, [3 2 1]);
% I DO NOT NEED THE COORDINATES TO PLOT IT ON WIREFRAME!!!
% THIS IS ALL I NEED:
% hp = plot3(curr_objectPoints(:,1), curr_objectPoints(:,3), curr_objectPoints(:,2), '.','linewidth',2, 'color','g','markers',30);   

% for point = 1:size(curr_objectPoints,1)  % OK - this was just to
% double-check. I don't need this info here
%     % find the annotation, name, and acronym of the current ROI pixel
%     ann = av(curr_objectPoints(point,1),curr_objectPoints(point,2),curr_objectPoints(point,3));
%     name = st.safe_name{ann};
%     disp(name)
% %     acr = st.acronym{ann};
% end

% fit a plane to these points, to be used as separator. If I add more
% points in more than two slices, I can use a polynomial fit. 
% plot3(ap_pixel, ml_pixel, dv_pixel, '.','linewidth',2, 'color',braincolor,'markers',10);
hp = plot3(curr_objectPoints(:,1), curr_objectPoints(:,3), curr_objectPoints(:,2), '.','linewidth',2, 'color','g','markers',20);
sf = fit([curr_objectPoints(:,1), curr_objectPoints(:,3)], curr_objectPoints(:,2), 'poly32');
h_sf = plot(sf);
% h_sf = plot(sf, [curr_objectPoints(:,1), curr_objectPoints(:,3)], curr_objectPoints(:,2));    %does not do what I wanted
% delete(hp)
% delete(h_sf)
makeVideo = 0;
if makeVideo
    view([-90, 8]);
    WriterObj = VideoWriter('allenCCF_isthmus_movie_withAllCells_surface.mp4', 'MPEG-4');
    WriterObj.FrameRate=30;
    open(WriterObj);
    for f = 8:20
        view([-90 f])
        drawnow;
        if makeVideo
            frame = getframe(fh);
            writeVideo(WriterObj,frame);
        end
    end
    for f = -90:8
        view([f 20])
        drawnow;
        if makeVideo
            frame = getframe(fh);
            writeVideo(WriterObj,frame);
        end
    end
    close(WriterObj);
end
delete(hp)
delete(h_sf)
delete(p)
%% now use the surface to separate out the unwanted cells
for i = 1:length(S)
    % transform coordinates
    ap_pixel = bregma(1) - S(i).T_roi.AP_location(S(i).pltIdx)./atlas_resolution; %OK
    ml_pixel = bregma(3) + S(i).T_roi.ML_location(S(i).pltIdx)./atlas_resolution; %OK
    dv_pixel = bregma(2) + S(i).T_roi.DV_location(S(i).pltIdx)./atlas_resolution; %OK
    
    leftPX = S(i).T_roi.ML_location(S(i).pltIdx) <= 0;
    leftPX_ind = find(leftPX);
    rightPX_ind = find(~leftPX);
    % apply the surface only on the left side, where it was defined, then
    % mirror it to the right side and apply it there too.
    % z axis is inverted!
    cells_inc_L = leftPX_ind(dv_pixel(leftPX) > sf(ap_pixel(leftPX), ml_pixel(leftPX)));     % in pltIdx reference
    cells_inc_R = rightPX_ind(dv_pixel(~leftPX) > sf(ap_pixel(~leftPX), 2*bregma(3)-ml_pixel(~leftPX))); % in pltIdx reference
    
    % put left and right back together
    S(i).cells_inc = false(sum(S(i).pltIdx),1);
    S(i).cells_inc(cat(1, cells_inc_L(:), cells_inc_R(:))) = true;
    
    p_filt(i) = plot3(ap_pixel(S(i).cells_inc), ml_pixel(S(i).cells_inc), dv_pixel(S(i).cells_inc), '.','linewidth',2, 'color', S(i).braincolor, 'markers',12);

end

makeVideo = 0;
if makeVideo
    view([-90, 8]);
    WriterObj = VideoWriter('allenCCF_isthmus_movie_ventralCellsONLY_surface.mp4', 'MPEG-4');
    WriterObj.FrameRate=30;
    open(WriterObj);
%     for f = 8:20
%         view([-90 f])
%         drawnow;
%         if makeVideo
%             frame = getframe(fh);
%             writeVideo(WriterObj,frame);
%         end
%     end
    for f = -90:8
        view([f, 8])
        drawnow;
        if makeVideo
            frame = getframe(fh);
            writeVideo(WriterObj,frame);
        end
    end
    close(WriterObj);
end