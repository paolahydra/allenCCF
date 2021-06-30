% demo to plot distribution of cells on different 2D planes


make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]); %now use subplot as you normally would! :)
if ~make_it_tight,  clear subplot;  end

bregma = allenCCFbregma();
isROI = av==SNR; % >0 for original av, >1 for by_index


newmap = colorcet('R1', 'N', 3+1);  % actually the old map - adding 1 to get rid of yellow for now
% newmap = brewermap(7, 'Dark2'); %too dark
% newmap = brewermap(7, 'Pastel1'); %too light
% newmap = brewermap(7, 'Accent'); % still too light
% newmap = colorcet('CBD1', 'N', 3); %nope
% newmap = brewermap(7, 'Set1'); % could do, but not good for colorblindness
for i = 1:length(S)
    S(i).braincolor = newmap(i,:);
end

%% isolate cells in ipsilateral SNR only, from three+ different brains
for i = 1:length(S)
    S(i).pltIdx_leftSNR = S(i).T_roi.avIndex == SNR & S(i).T_roi.ML_location <= 0;     
     % transform coordinates
    P(i).x = bregma(1) - S(i).T_roi.AP_location(S(i).pltIdx_leftSNR)./atlas_resolution; %OK
    P(i).y = bregma(3) + S(i).T_roi.ML_location(S(i).pltIdx_leftSNR)./atlas_resolution; %OK
    P(i).z = bregma(2) + S(i).T_roi.DV_location(S(i).pltIdx_leftSNR)./atlas_resolution; %OK
end
X = cat(1, P.x);
Y = cat(1, P.y);
Z = cat(1, P.z);
brain = []; %categorical label
for i = 1:length(S)
    brain = cat(1, brain, i*ones(size(P(i).x)));
end



%% plot brain wire with SNR grid and cells 
fwireframe = figure('Color','w');

isROIgrid = av==1; %root
gridIn3D(double(isROIgrid), 0.25, 15, bregma, [0.5 0.5 0.5]);
hold on

isROIgrid = av==SNR;
gridIn3D(double(isROIgrid), 0.25, 15, bregma, [0.7 0.7 0.7]); %ndp: function contourHands = gridIn3D(volData, contourHeight, sliceSpacing, origin, contourColor)
axis vis3d
set(gca, 'ZDir', 'reverse')
axis equal
axis off
% view([-30    25]);
view([-60, 12]);
hold on

view([90 -90])
l_H = line(c_H(2,:),c_H(1,:), 'Color', 'k', 'LineWidth',1);
% export_fig('horizontal_withBrainANDCells.pdf')
export_fig('horizontal_withBrainANDCells.png', '-m2')
delete(l_H)

view([-180 0])
l_S = line(c_S(1,:), zeros(size(c_S(2,:))), c_S(2,:), 'Color', 'k', 'LineWidth',2);
export_fig('sagittal_withBrainANDCells.png', '-m2')
% export_fig('sagittal_withBrainANDCells.pdf')
delete(l_S)

view([90 0])
l_F = line(zeros(size(c_F(2,:))), c_F(1,:), c_F(2,:), 'Color', 'k', 'LineWidth',2);
export_fig('coronal_withBrainANDCells.png', '-m2')
delete(l_F)

%% set up one figure at the time and save it, FGS

fH(1) = figure('Color','w');

isROIgrid = av==SNR | av==RR;
gridIn3D(double(isROIgrid), 0.25, 15, bregma, [0.7 0.7 0.7]); %ndp: function contourHands = gridIn3D(volData, contourHeight, sliceSpacing, origin, contourColor)
axis vis3d
set(gca, 'ZDir', 'reverse')
axis equal
axis off
% view([-30    25]);
view([-60, 12]);
hold on

for i = 1:length(S)
    if i == 1 && exist('p_snr', 'var')
        delete(p_snr)
    end
    p_snr(i) = plot3(P(i).x, P(i).y, P(i).z, '.','linewidth',2, 'color', S(i).braincolor, 'markers',10);
end
legend(p_snr, 'brain 1', 'brain 2', 'brain 3')
legend('boxoff')

export_fig('SNR_3Dmodel_snap.pdf')

makeVideo = 0;
if makeVideo
    WriterObj = VideoWriter('allenCCF_leftSNRcells_movie_Root.mp4', 'MPEG-4');
    WriterObj.FrameRate=30;
    open(WriterObj);
    for f = -60:299
        view([f, 12])
        drawnow;
        frame = getframe(fH(1));
        writeVideo(WriterObj,frame);
    end
    close(WriterObj);
end

delete(handles.root)

% % was this faster?
% fwireframe = plotWireFrame(T_roi, braincolor, black_brain, fwireframe, microns_per_pixel, microns_per_pixel_after_downsampling );

%% 2: HORIZONTAL VIEW (left)
% I can add one representative contour of root at the ~center of ROI
% (overlay in white!)

figN = 2;

fH(figN) = figure('Color','w');
isROI_H =squeeze(nansum(isROI, 2) > 0);
scatter_handles(figN).h = scatterhist(Y,X,'Group',brain,'Kernel','on', 'Location', 'NorthEast', 'Direction', 'out', 'Color', cat(1,S.braincolor), 'Parent', fH(figN));
scatter_handles(figN).h(1).NextPlot = 'add';


c = contourc(double(isROI_H), 0.25*[1 1]);    
%     function cout = parseContours(c)
startInd = 1;
cInd = 1;
cout = {};
while startInd<size(c,2)
    nC = c(2,startInd);
    cout{cInd} = c(:,startInd+1:startInd+nC);
    cInd = cInd+1;
    startInd = startInd+nC+1;
end
c_H = cout{1}; %left
l_H = line(c_H(1,:), c_H(2,:), 'Color', 'k', 'LineWidth',2);


% scatter_handles(figN).h(1).YLimMode = 'auto';
axis image
XLIM = scatter_handles(figN).h(1).XLim;
YLIM = scatter_handles(figN).h(1).YLim;
offset = 10;
scatter_handles(figN).h(1).XLim = [XLIM(1)-offset, XLIM(2)+offset];
scatter_handles(figN).h(1).YLim = [YLIM(1)-offset, YLIM(2)+offset];

YL = scatter_handles(figN).h(1).XTick; 
YL = (YL-bregma(3)).* atlas_resolution;
scatter_handles(figN).h(1).XTickLabel = YL(:);
scatter_handles(figN).h(1).XLabel.String = 'ML';


XL = scatter_handles(figN).h(1).YTick; 
XL = (bregma(1)-XL).* atlas_resolution;
scatter_handles(figN).h(1).YTickLabel = XL(:);
scatter_handles(figN).h(1).YLabel.String = 'AP';

scatter_handles(figN).h(1).Box = 'off';

hLegend = findobj(gcf,'Type','legend');
hLegend.String = {'brain 1', 'brain 2', 'brain 3'};
hLegend.Box = 'off';

figure(fH(figN))
export_fig('SNR_horizontalDistribution.pdf')

% 
%     P(i).x = bregma(1) - S(i).T_roi.AP_location(S(i).pltIdx_leftSNR)./atlas_resolution; %OK
%     P(i).y = bregma(3) + S(i).T_roi.ML_location(S(i).pltIdx_leftSNR)./atlas_resolution; %OK
%     P(i).z = bregma(2) + S(i).T_roi.DV_location(S(i).pltIdx_leftSNR)./atlas_resolution; %OK


%% 3: paraSAGITAL VIEW
% isROI_S =permute(nansum(isROI, 3) > 0, [2,1]);
% xlabel('x')
% ylabel('z')

figN = 3;

fH(figN) = figure('Color','w');
isROI_S =permute(nansum(isROI, 3) > 0, [2,1]);
scatter_handles(figN).h = scatterhist(X,Z,'Group',brain,'Kernel','on', 'Location', 'NorthEast', 'Direction', 'out', 'Color', cat(1,S.braincolor), 'Parent', fH(figN));
scatter_handles(figN).h(1).NextPlot = 'add';


c = contourc(double(isROI_S), 0.25*[1 1]);    
%     function cout = parseContours(c)
startInd = 1;
cInd = 1;
cout = {};
while startInd<size(c,2)
    nC = c(2,startInd);
    cout{cInd} = c(:,startInd+1:startInd+nC);
    cInd = cInd+1;
    startInd = startInd+nC+1;
end
c_S = cout{1}; %only one here actually
l_S = line(c_S(1,:), c_S(2,:), 'Color', 'k', 'LineWidth',2);


% scatter_handles(figN).h(1).YLimMode = 'auto';
axis image
XLIM = scatter_handles(figN).h(1).XLim;
YLIM = scatter_handles(figN).h(1).YLim;
offset = 10;
scatter_handles(figN).h(1).XLim = [XLIM(1)-offset, XLIM(2)+offset];
scatter_handles(figN).h(1).YLim = [YLIM(1)-offset, YLIM(2)+offset];

ZL = scatter_handles(figN).h(1).YTick; 
ZL = (ZL-bregma(2)).* atlas_resolution;
scatter_handles(figN).h(1).YTickLabel = ZL(:);
scatter_handles(figN).h(1).YLabel.String = 'DV';


XL = scatter_handles(figN).h(1).XTick; 
XL = (bregma(1)-XL).* atlas_resolution;
scatter_handles(figN).h(1).XTickLabel = XL(:);
scatter_handles(figN).h(1).XLabel.String = 'AP';

scatter_handles(figN).h(1).Box = 'off';

%need to revert the DV axis in this case to make it natural
scatter_handles(figN).h(1).YDir = 'reverse'; % scatter
scatter_handles(figN).h(3).XDir = 'reverse'; % histogram


hLegend = findobj(gcf,'Type','legend');
hLegend.String = {'brain 1', 'brain 2', 'brain 3'};
hLegend.Box = 'off';

figure(fH(figN))
export_fig('SNR_parasagittalDistribution.pdf')


%% 4: FRONTAL VIEW
% isROI_F =squeeze(nansum(isROI, 1) > 0);
% xlabel('y')
% ylabel('z')

figN = 4;

fH(figN) = figure('Color','w');
isROI_F =squeeze(nansum(isROI, 1) > 0);
scatter_handles(figN).h = scatterhist(Y,Z,'Group',brain,'Kernel','on', 'Location', 'NorthEast', 'Direction', 'out', 'Color', cat(1,S.braincolor), 'Parent', fH(figN));
scatter_handles(figN).h(1).NextPlot = 'add';


c = contourc(double(isROI_F), 0.25*[1 1]);    
%     function cout = parseContours(c)
startInd = 1;
cInd = 1;
cout = {};
while startInd<size(c,2)
    nC = c(2,startInd);
    cout{cInd} = c(:,startInd+1:startInd+nC);
    cInd = cInd+1;
    startInd = startInd+nC+1;
end
c_F = cout{1}; %LEFT?
l_F = line(c_F(1,:), c_F(2,:), 'Color', 'k', 'LineWidth',2);


% scatter_handles(figN).h(1).YLimMode = 'auto';
axis image
XLIM = scatter_handles(figN).h(1).XLim;
YLIM = scatter_handles(figN).h(1).YLim;
offset = 10;
scatter_handles(figN).h(1).XLim = [XLIM(1)-offset, XLIM(2)+offset];
scatter_handles(figN).h(1).YLim = [YLIM(1)-offset, YLIM(2)+offset];

ZL = scatter_handles(figN).h(1).YTick; 
ZL = (ZL-bregma(2)).* atlas_resolution;
scatter_handles(figN).h(1).YTickLabel = ZL(:);
scatter_handles(figN).h(1).YLabel.String = 'DV';


YL = scatter_handles(figN).h(1).XTick; 
YL = (YL-bregma(3)).* atlas_resolution;
scatter_handles(figN).h(1).XTickLabel = YL(:);
scatter_handles(figN).h(1).XLabel.String = 'ML';

scatter_handles(figN).h(1).Box = 'off';

%need to revert the DV axis in this case to make it natural
scatter_handles(figN).h(1).YDir = 'reverse'; % scatter
scatter_handles(figN).h(3).XDir = 'reverse'; % histogram


hLegend = findobj(gcf,'Type','legend');
hLegend.String = {'brain 1', 'brain 2', 'brain 3'};
hLegend.Box = 'off';

figure(fH(figN))
export_fig('SNR_frontalDistribution.pdf')


