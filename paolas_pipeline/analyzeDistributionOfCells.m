%% this is a demo. You will need to adapt it to your own purposes.
%this demo will show you how to:
% 1: plot all your detected cells in a 3D model of the brain
% 2. tabulate all the regions occupied by the detected cells
% 3. identify specific regions of interest, isolate cells within them and plot it.
% 4. create some manual boundaries and use them to split ROIs.


if ~exist('st','var')
    disp('loading structure tree...')
    st = loadStructureTree(structure_tree_location);
end
atlas_resolution = 0.010; % mm
bregma = allenCCFbregma(); % bregma position in reference data space

%% load data to analyze

folder2save = 'D:\afonanar\postRegistrationAnalyzedData';
mkdir(folder2save)
cd(folder2save)

folders2Include = uipickfiles();
plotFWire = 0;
object_tag = 'green';

S = loadTabDataFromMultipleBrains(folders2Include, plotFWire, object_tag, microns_per_pixel); %updated plotting to exclude point that are in 'root' (avIndex = 1)

% TH = loadTabDataFromMultipleBrains(folders2Include, plotFWire, 'red', microns_per_pixel); %updated plotting to exclude point that are in 'root' (avIndex = 1)


%% plot brain wire with all cells (much faster with brainwire precalculated!!)
black_brain = false;
fwireframe = [];
for i = 1:length(S)
    fwireframe = plotWireFrame(S(i).T_roi, S(i).braincolor, black_brain, fwireframe, microns_per_pixel, microns_per_pixel_after_downsampling );
end
makeVideo = 0;
if makeVideo  
    WriterObj = VideoWriter('allenCCF_allCells_Root_movie_full.mp4', 'MPEG-4');
    WriterObj.FrameRate=30;
    open(WriterObj);
    for f = -37:322
        view([f, 22])
        drawnow;
        frame = getframe(fwireframe);
        writeVideo(WriterObj,frame);
    end
    close(WriterObj);
end


%% find unique targets and sort them (make distribution tables)
for i = 1:length(S)
    
    all_Names = S(i).T_roi.name;
    all_Files = S(i).T_roi.roiFIle;
    all_acIndex = S(i).T_roi.avIndex;
    all_isPositiveML_location = S(i).T_roi.ML_location;       % left hemisphere only
    all_isPositiveML_location = all_isPositiveML_location>0;  % right hemisphere only
    
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
% it's actually not even good because I would need the number of slices that
% intersect that region of interest, not the number of slices that have at
% least one cell

%% define specific regions of interest
% examples:
MB = find(strcmp(st.name,'Midbrain'));
MBRN = find(strcmp(st.name,'Midbrain reticular nucleus'));
VTA = find(strcmp(st.name,'Ventral tegmental area'));

CUN = find(strcmp(st.name,'Cuneiform nucleus'));
PPN = find(strcmp(st.name,'Pedunculopontine nucleus'));

VIIN = find(strcmp(st.name, 'Facial motor nucleus'));
PCRT = find(strcmp(st.name, 'Parvicellular reticular nucleus'));
IRT = find(strcmp(st.name, 'Intermediate reticular nucleus'));

%% Create the wire plots
szAP = size(tv,1);
szDV = size(tv,2);
szLR = size(tv,3);

fh = figure; 
figcol = 0.7;
fh.Color = [figcol figcol figcol];

bregma = allenCCFbregma();
isROI = (av==CUN | av==PPN); 
gridIn3D(double(isROI), 0.25, 10, bregma, 'w');
axis vis3d
set(gca, 'ZDir', 'reverse')
axis equal
axis off

view([-90, 8]);
hold on

ax = gca;
mlr_grid_handles = ax.Children; 

totHandles = length(ax.Children);
isROI = av==1; % >0 for original av, >1 for by_index
gridIn3D(double(isROI), 0.25, 50, bregma, [0 0 0]); %ndp: function contourHands = gridIn3D(volData, contourHeight, sliceSpacing, origin, contourColor)
handles.root = ax.Children(1:end-totHandles);

% % to toggle line visibility:
% for rl = 1: length(handles.root)
%     handles.root(rl).Visible = 'off';
% end


%% define and add the cells!
for i = 1:length(S)
    % first add only SNR and SNC cells, consider both sides (there are
    % sometimes very few cells in the contra side too)
    if i == 1 && exist('p', 'var')
        delete(p)
    end
    S(i).pltIdx = S(i).T_roi.avIndex == CUN ...
        | S(i).T_roi.avIndex == PPN; 
    
%     % to only show cells in the ROIs on the left side:
%     S(i).pltIdx = ( S(i).T_roi.avIndex == CUN ...
%                     | S(i).T_roi.avIndex == PPN ) ...
%                   & S(i).T_roi.ML_location <= 0;       %  Note the ()s.

    
    % transform coordinates
    ap_pixel = bregma(1) - S(i).T_roi.AP_location(S(i).pltIdx)./atlas_resolution;
    ml_pixel = bregma(3) + S(i).T_roi.ML_location(S(i).pltIdx)./atlas_resolution;
    dv_pixel = bregma(2) + S(i).T_roi.DV_location(S(i).pltIdx)./atlas_resolution; 
    
    p(i) = plot3(ap_pixel, ml_pixel, dv_pixel, '.','linewidth',2, 'color', S(i).braincolor, 'markers',12); 
end


%% save figure, export a movie
export_fig('PPN_CUN_outline_withCells.png', '-nocrop', '-png', '-m5')
makeVideo = 0;
if makeVideo
    view([-90, 8]);
    WriterObj = VideoWriter('allenCCF_mlr_movie_withAllCells.mp4', 'MPEG-4');
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

% retrieve manually set boundary points:
folder_pointlist = '/Users/galileo/dati/registered_brains_completed/993030/startingSingleSlices/processed/';
load(fullfile(folder_pointlist, 'probe_points.mat'), 'pointList');
object_num = 3; % second row is first attempt with only 4 points for fitting a plane. It was too ventral
% pointList.pointList{2,1}; %the points
% pointList.pointList{2,3}; % the registerd slices they come from
curr_objectPoints = pointList.pointList{object_num,1}(:, [3 2 1]);


% fit a plane to these points, to be used as separator. If I add more
% points in more than two slices, I can use a polynomial fit. 

% plot3(ap_pixel, ml_pixel, dv_pixel, '.','linewidth',2, 'color',braincolor,'markers',10);
hp = plot3(curr_objectPoints(:,1), curr_objectPoints(:,3), curr_objectPoints(:,2), '.','linewidth',2, 'color','g','markers',20);
sf = fit([curr_objectPoints(:,1), curr_objectPoints(:,3)], curr_objectPoints(:,2), 'poly32');
h_sf = plot(sf);

makeVideo = 0;
if makeVideo
    view([-90, 8]);
    WriterObj = VideoWriter('allenCCF_mlr_movie_withAllCells_surface.mp4', 'MPEG-4');
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
    WriterObj = VideoWriter('allenCCF_mlr_movie_ventralCellsONLY_surface.mp4', 'MPEG-4');
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


