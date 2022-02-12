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