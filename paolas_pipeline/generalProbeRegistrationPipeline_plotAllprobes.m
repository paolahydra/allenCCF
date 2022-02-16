% This software is based on:
% https://github.com/cortex-lab/allenCCF/

% Arber lab mainteined repository (forked from cortex-lab):
% https://github.com/paolahydra/allenCCF/tree/sliceRegistration
%
% % The following are needed for full functionality:
% Images of mouse brain slices (individually cropped or with multiple slices per image; coronal, sagittal, or transverse)
% Know the resolution (in microns per pixel) of these images
% A computer mouse with a scroll wheel
% MATLAB (R2017 or above used for testing)
% This repository. Add all folders and subfolders to your MATLAB path. All user-oriented scripts are in the 'SHARP-Track' folder.
% The npy-matlab repository: http://github.com/kwikteam/npy-matlab
% The Allen Mouse Brain Atlas volume and annotations (download all 4 files
% from this link: http://data.cortexlab.net/allenCCF/ )



% * remember to run one section at a time, instead of the whole script at once *
doinitialize = 0;
if doinitialize
    fwireframe = [];
    fwireframe = plotBrainGrid([], [], fwireframe, black_brain);  %-PP uncomment
    hold on; 
    fwireframe.InvertHardcopy = 'off';
end
cab(5)
makeVideo = 0;
%%  set once, then always run: specify paths and settings for the specific brain to register
% move your images to a local disk (SSD possibly) for much faster processing!
image_folder = '\\tungsten-nas.fmi.ch\tungsten\scratch\garber\patepaol\intan_recs\anatomy_probeMapping\1046411';   %change this
image_tag = 'mouse_1046411';                                               %change this - use an unequivocal tag for your experiment
probe_lengths =         4.250;  % how far into the brain did you go from the surface, either for each probe or just one number for all -- in mm
active_probe_length =   0.200;  % from the bottom tip, how much of the probe contained recording sites -- in mm
nSessions_down =        3;      % number of sessions you turned down by an active length

brainNum = 6;


microns_per_pixel = 5.1803; %3.4536; %3.8852; %take this value from your tiff filename
cab(5)
% increase gain if for some reason the images are not bright enough
gain = 3;   % for visualization only: during cropping or atlas alignment


% stable settings (no need to change):
repositoryTag = 'SR_GPR'; %this is the SliceRegistrationConfocal pipeline - these scripts are added to gitignore, useful to know where they came from

generalPipelinesFolder = 'generalPipelineUsedScripts'; %leave empty if you want to save the 
                        % mouse-specific script in the same folder as the
                        % generalPipeline.m script.  Or else, specify a new
                        % folder here.

plane = 'coronal';  % plane to view ('coronal', 'sagittal', 'transverse')
% transformation to use for registration:
transformationType = 'pwl';     %use 'projective', or 'pwl' (piece-wise linear: more advanced).                       
  

%%  always run: (set once to make it more easily portable across the computers you use)
[~,localhost] = system('hostname');
% archstr = computer('arch');     % it could be of use too

if strcmp(localhost(1:end-1), 'Paolas-MacBook-Pro.local')
    rootFolderScripts = '/Users/galileo/GitHub/';
    addpath(genpath('/Users/galileo/GitHub/matlabUtilities/npy-matlab-master/npy-matlab'))
    pathToAtlas = '/Users/galileo/Documents/MATLAB/codeArberLab/anatomyRegistration/cortexLabCode/allen brain template files';
elseif strcmp(localhost(1:end-1), 'f452d-96617e') % put yours here
    rootFolderScripts = 'C:\GitHub\';    % % -------------------   put yours here
    addpath(genpath('\\tungsten-nas.fmi.ch\tungsten\scratch\garber\BrainRegistration\code and atlas')) %includes npy
    pathToAtlas = '\\tungsten-nas.fmi.ch\tungsten\scratch\garber\BrainRegistration\code and atlas\allen brain template files\';
else 
    rootFolderScripts = 'C:\Users\patepaol\Documents\GitHub';  % % -------------------   put yours here
    addpath(genpath('\\tungsten-nas.fmi.ch\tungsten\scratch\garber\BrainRegistration\code and atlas')) %includes npy
    pathToAtlas = '\\tungsten-nas.fmi.ch\tungsten\scratch\garber\BrainRegistration\code and atlas\allen brain template files\';
end
clear localhost
addpath(genpath(fullfile(rootFolderScripts, 'allenCCF')))


% no need to change below:

annotation_volume_location = fullfile(pathToAtlas, 'annotation_volume_10um_by_index.npy');
structure_tree_location = fullfile(pathToAtlas, 'structure_tree_safe_2017.csv');
template_volume_location = fullfile(pathToAtlas, 'template_volume_10um.npy');

A = matlab.desktop.editor.getActive;                                                
if ~isempty(generalPipelinesFolder)                       
    genPipScriptLibrary = fullfile(rootFolderScripts, generalPipelinesFolder);
    if ~exist(genPipScriptLibrary, 'dir')
        mkdir(genPipScriptLibrary)
    end
else
    [genPipScriptLibrary, ~] = fileparts(A.Filename);
end


if ~strcmp( image_tag(end), '_')
    image_tag = cat(2, image_tag, '_');
end

cd(image_folder)
save_folder = fullfile(image_folder, 'startingSingleSlices');
folder_processed_images = fullfile(save_folder, 'processed');


%
A.save;
scriptname = fullfile(genPipScriptLibrary, sprintf('generalPipeline_%s_%s.m',repositoryTag, image_tag(1:end-1)));

if strcmp(A.Filename, scriptname)
    fprintf('Script %s already existing and in use.\n',sprintf('generalPipeline_%s_%s.m',repositoryTag, image_tag(1:end-1)))
else
    if exist(scriptname, 'file') %do not overwrite
        scriptname = fullfile(genPipScriptLibrary, sprintf('generalPipeline_%s_%s_%s.m',repositoryTag, image_tag(1:end-1), string(datetime('now', 'Format','yyyyMMdd_hhmmss'))));
    end
    copyfile(A.Filename, scriptname)
%     A.close
    clear A
    matlab.desktop.editor.openAndGoToLine(scriptname, getcurrentline+4);
end


%% always run: filesystem and parameter definition - don't need to change
% directory of single histology images
image_folder = save_folder;

% if the images are individual slices (as opposed to images of multiple
% slices, which must be cropped using the cell CROP AND SAVE SLICES)
image_files_are_individual_slices = true;

% use images that are already at reference atlas resolution (here, 10um/pixel)
use_already_downsampled_image = false; 

% pixel size parameters: microns_per_pixel of large images in the image
% folder (if use_already_downsampled_images is set to false);
% microns_per_pixel_after_downsampling should typically be set to 10 to match the atlas
microns_per_pixel_after_downsampling = 10;


% additional parameters
% size in pixels of reference atlas brain coronal slice, typically 800 x 1140
atlas_reference_size = [800 1140]; 
reference_size = [1320 800 1140];


% auto: naming definition


% % ENTER PARAMETERS AND FILE LOCATION

% % file location of probe points
% folder_processed_images = 'C:\Drive\Histology\brainX\processed';
% 
% % directory of reference atlas files
% annotation_volume_location = 'C:\Drive\Histology\for tutorial\annotation_volume_10um_by_index.npy';
% structure_tree_location = 'C:\Drive\Histology\for tutorial\structure_tree_safe_2017.csv';

% name of the saved probe points
% probe_save_name_suffix = 'electrode_track_1';
probe_save_name_suffix = '';

% either set to 'all' or a list of indices from the clicked probes in this file, e.g. [2,3]
probes_to_analyze = 'all';  % [1 2]

% --------------
% key parameters
% --------------

% distance queried for confidence metric -- in um
probe_radius = 70; 

% overlay the distance between parent regions in gray (this takes a while)
show_parent_category = false; 

% plot this far or to the bottom of the brain, whichever is shorter -- in mm
distance_past_tip_to_plot = 0.5;

% set scaling e.g. based on lining up the ephys with the atlas
% set to *false* to get scaling automatically from the clicked points
scaling_factor = false;    % 1.05-1.12 reasonable range


% ---------------------
% additional parameters
% ---------------------
% plane used to view when points were clicked ('coronal' -- most common, 'sagittal', 'transverse')
plane = 'coronal';

% probe insertion direction 'down' (i.e. from the dorsal surface, downward -- most common!) 
% or 'up' (from a ventral surface, upward)
probe_insertion_direction = 'down';

% show a table of regions that the probe goes through, in the console
show_region_table = true;
      
% black brain?
black_brain = true;


% close all




% % GET AND PLOT PROBE VECTOR IN ATLAS SPACE

% load the reference brain annotations
if ~exist('av','var') || ~exist('st','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
end

% select the plane for the viewer
if strcmp(plane,'coronal')
    av_plot = av;
elseif strcmp(plane,'sagittal')
    av_plot = permute(av,[3 2 1]);
elseif strcmp(plane,'transverse')
    av_plot = permute(av,[2 3 1]);
end

% load probe points
probePoints = load(fullfile(folder_processed_images, ['probe_points' probe_save_name_suffix]));
ProbeColors = .75*[1.3 1.3 1.3; 1 .75 0;  .3 1 1; .4 .6 .2; 1 .35 .65; .7 .7 .9; .65 .4 .25; .7 .95 .3; .7 0 0; .6 0 .7; 1 .6 0]; 
% order of colors: {'white','gold','turquoise','fern','bubble gum','overcast sky','rawhide', 'green apple','purple','orange','red'};
% fwireframe = []; %-PP uncomment

% scale active_probe_length appropriately
active_probe_length = active_probe_length*100;

% determine which probes to analyze
if strcmp(probes_to_analyze,'all')
    probes = 1:size(probePoints.pointList.pointList,1);
else
    probes = probes_to_analyze;
end 





% % PLOT EACH PROBE -- FIRST FIND ITS TRAJECTORY IN REFERENCE SPACE

% create a new figure with wireframe
% fwireframe = plotBrainGrid([], [], fwireframe, black_brain);  %-PP uncomment
% hold on; 
% fwireframe.InvertHardcopy = 'off';
clear probeRecap

for selected_probe = probes
    
% get the probe points for the currently analyzed probe 
if strcmp(plane,'coronal')
    curr_probePoints = probePoints.pointList.pointList{selected_probe,1}(:, [3 2 1]);
%     bregma = allenCCFbregma();
%     curr_probePoints(:,3) = 2*bregma(3) - curr_probePoints(:,3); %reflect about the midline because probes were assigned to the wrong hemisphere here
elseif strcmp(plane,'sagittal')
    curr_probePoints = probePoints.pointList.pointList{selected_probe,1}(:, [1 2 3]);
elseif strcmp(plane,'transverse')
    curr_probePoints = probePoints.pointList.pointList{selected_probe,1}(:, [1 3 2]);
end



% get user-defined probe length from experiment
if length(probe_lengths) > 1
    probe_length = probe_lengths(selected_probe);
else
    probe_length = probe_lengths;
end

% get the scaling-factor method to use
if scaling_factor
    use_tip_to_get_reference_probe_length = false;
    reference_probe_length = probe_length * scaling_factor;
    disp(['probe scaling of ' num2str(scaling_factor) ' determined by user input']);    
else
    use_tip_to_get_reference_probe_length = true;
    disp(['getting probe scaling from histology data...']);
end

% get line of best fit through points
% m is the mean value of each dimension; p is the eigenvector for largest eigenvalue
[m,p,s] = best_fit_line(curr_probePoints(:,1), curr_probePoints(:,2), curr_probePoints(:,3));
if isnan(m(1))
    disp(['no points found for probe ' num2str(selected_probe)])
    continue
end

% ensure proper orientation: want 0 at the top of the brain and positive distance goes down into the brain
if p(2)<0
    p = -p;
end

% determine "origin" at top of brain -- step upwards along tract direction until tip of brain / past cortex
ann = 10;
out_of_brain = false;
while ~(ann==1 && out_of_brain) % && distance_stepped > .5*active_probe_length)
    m = m-p; % step 10um, backwards up the track
    ann = av_plot(round(m(1)),round(m(2)),round(m(3))); %until hitting the top
    if strcmp(st.safe_name(ann), 'root')
        % make sure this isn't just a 'root' area within the brain
        m_further_up = m - p*25; % is there more brain 250 microns up along the track?
        ann_further_up = av_plot(round(max(1,m_further_up(1))),round(max(1,m_further_up(2))),round(max(1,m_further_up(3))));
        if strcmp(st.safe_name(ann_further_up), 'root')
            m_further_up = m - p*50; % is there more brain 500 total microns up along the track?
            ann_further_up = av_plot(round(max(1,m_further_up(1))),round(max(1,m_further_up(2))),round(max(1,m_further_up(3))));
            if strcmp(st.safe_name(ann_further_up), 'root')
                out_of_brain = true;
            end
        end
    end
end

% focus on wireframe plot
figure(fwireframe);

% plot probe points
% hp = plot3(curr_probePoints(:,1), curr_probePoints(:,3), curr_probePoints(:,2), '.','linewidth',2, 'color',[ProbeColors(selected_probe,:) .2],'markers',10);  %-PP uncomment

% plot brain entry point
% plot3(m(1), m(3), m(2), 'r*','linewidth',1)   %-PP uncomment

% use the deepest clicked point as the tip of the probe, if no scaling provided (scaling_factor = false)
if use_tip_to_get_reference_probe_length
    % find length of probe in reference atlas space
    if strcmp(probe_insertion_direction, 'down')
        [depth, tip_index] = max(curr_probePoints(:,2));
    elseif strcmp(probe_insertion_direction, 'up')
        [depth, tip_index] = min(curr_probePoints(:,2));    
    end
    reference_probe_length_tip = sqrt(sum((curr_probePoints(tip_index,:) - m).^2)); 
    
    % and the corresponding scaling factor
    shrinkage_factor = (reference_probe_length_tip / 100) / probe_length;
    
    % display the scaling
    disp(['probe length of ' num2str(reference_probe_length_tip/100) ' mm in reference atlas space compared to a reported ' num2str(probe_length) ' mm']);
    disp(['probe scaling of ' num2str(shrinkage_factor)]); disp(' ');
    
    % plot line the length of the probe in reference space
    probe_length_histo = round(reference_probe_length_tip);
    
% if scaling_factor is user-defined as some number, use it to plot the length of the probe
else 
    probe_length_histo = round(reference_probe_length * 100); 
end

% find the percent of the probe occupied by electrodes
percent_of_tract_with_active_sites = min([active_probe_length / (probe_length*100), 1.0]);
active_site_start = probe_length_histo*(1-percent_of_tract_with_active_sites);
active_probe_position = round([active_site_start  probe_length_histo]);

% % plot line the length of the active probe sites in reference space
% plot3(m(1)+p(1)*[active_probe_position(1) active_probe_position(2)], m(3)+p(3)*[active_probe_position(1) active_probe_position(2)], m(2)+p(2)*[active_probe_position(1) active_probe_position(2)], ...
%     'Color', ProbeColors(selected_probe,:), 'LineWidth', 1);   %-PP uncomment
% plot line the length of the entire probe in reference space
plot3(m(1)+p(1)*[1 probe_length_histo], m(3)+p(3)*[1 probe_length_histo], m(2)+p(2)*[1 probe_length_histo], ...
    'Color', ProbeColors(brainNum,:), 'LineWidth', 1);


%% ----------------------------------------------------------------
% Get and plot brain region labels along the extent of each probe
% ----------------------------------------------------------------

% convert error radius into mm
error_length = round(probe_radius / 10);

% find and regions the probe goes through, confidence in those regions, and plot them
probeRecap(selected_probe).borders_table = plotDistToNearestToTip(m, p, av_plot, st, probe_length_histo, error_length, active_site_start, distance_past_tip_to_plot, show_parent_category, show_region_table, plane, nSessions_down); % plots confidence score based on distance to nearest region along probe
title(['Probe ' num2str(selected_probe)],'color',ProbeColors(selected_probe,:))
savefig(fullfile(folder_processed_images, sprintf('bordersTable_probe%d.fig',selected_probe)))
export_fig(fullfile(folder_processed_images, sprintf('bordersTable_%s_probe%d.png',image_tag, selected_probe)), '-m5');

% add some more digested info to proberecap and save it
probeRecap(selected_probe).shrinkage_factor = shrinkage_factor;
probeRecap(selected_probe).active_probe_position = active_probe_position;
probeRecap(selected_probe).m = m;
probeRecap(selected_probe).p = p;
% % there would be some rounding issue creating a mismatch between the
% % figure and the info here, if using this method:
% probeRecap(selected_probe).DVstart = 10* (probeRecap(selected_probe).active_probe_position(1) - diff(probeRecap(selected_probe).active_probe_position)* (0:nSessions_down-1) );
% probeRecap(selected_probe).DVtip =   10* (probeRecap(selected_probe).active_probe_position(2) - diff(probeRecap(selected_probe).active_probe_position)* (0:nSessions_down-1) );
% probeRecap(selected_probe).adjustedActivelength = 10*diff(probeRecap(selected_probe).active_probe_position);

% this is the same method used in the figure, with no rounding to the next 10.
actL = 10*(probe_length_histo-active_site_start);
probeRecap(selected_probe).adjustedActivelength = actL;
probeRecap(selected_probe).DVstart = round(10*active_site_start - actL*(0:nSessions_down-1));
probeRecap(selected_probe).DVtip =   round(10*probe_length_histo - actL*(0:nSessions_down-1));

Probe = probeRecap(selected_probe);
save(fullfile(folder_processed_images, sprintf('Probe%d.mat',selected_probe)), '-struct', 'Probe');

end



%% memo:
% % plot line the length of the entire probe in reference space
% probe_length_histo = active_probe_position(2);
% plot3(m(1)+p(1)*[1 probe_length_histo], m(3)+p(3)*[1 probe_length_histo], m(2)+p(2)*[1 probe_length_histo], ...
%     'Color', ProbeColors(selected_probe,:), 'LineWidth', 1, 'LineStyle',':');
% % plot points along the shank in reference space
% plot3(m(1)+p(1)*[active_probe_position(1) active_probe_position(2)], m(3)+p(3)*[active_probe_position(1) active_probe_position(2)], m(2)+p(2)*[active_probe_position(1) active_probe_position(2)], ...
%     'Color', ProbeColors(selected_probe,:), 'LineWidth', 1);


%%

%% load cells (retro_fromMed in github)
folders2Include={'\\tungsten-nas.fmi.ch\tungsten\scratch\garber\patepaol\intan_recs\anatomy_probeMapping\reference\992234', ...
                '\\tungsten-nas.fmi.ch\tungsten\scratch\garber\patepaol\intan_recs\anatomy_probeMapping\reference\993030',  ...
                '\\tungsten-nas.fmi.ch\tungsten\scratch\garber\patepaol\intan_recs\anatomy_probeMapping\reference\993031'};
plotFWire = 0;
object_tag = 'green'; %'green' for rabies cells
S = loadTabDataFromMultipleBrains(folders2Include, plotFWire, object_tag); %updated plotting to exclude point that are in 'root' (avIndex = 1)

%% define and add the cells!
clear p
SNR = 823;
SNC = 867;
bregma = allenCCFbregma();
atlas_resolution = 0.010; % mm
for i = 1:length(S)
    % first add only SNR and SNC cells, consider both sides (there are
    % sometimes very few cells in the contra side too)
    if i == 1 & exist('p', 'var')
        delete(p)
    end
    S(i).pltIdx = S(i).T_roi.avIndex == SNR ...
        | S(i).T_roi.avIndex == SNC; 

    
    % transform coordinates
    ap_pixel = bregma(1) - S(i).T_roi.AP_location(S(i).pltIdx)./atlas_resolution; %OK
    ml_pixel = bregma(3) - S(i).T_roi.ML_location(S(i).pltIdx)./atlas_resolution; %OK
    dv_pixel = bregma(2) + S(i).T_roi.DV_location(S(i).pltIdx)./atlas_resolution; %OK
    %reflect about midline
    ml_pixel = 2*bregma(3) - ml_pixel;
    
    p(i) = plot3(ap_pixel, ml_pixel, dv_pixel, '.','linewidth',2, 'color', S(3).braincolor, 'markers',7);
    
end
%%

if makeVideo   
    F = -29.2;
    view([F, 17])
    drawnow;
    export_fig('allProbes_SNR_SNC_cells_probes.png', '-nocrop', '-png', '-m5')
    
    %%

    WriterObj = VideoWriter('allProbes_SNR_SNC_NOcells_probes.avi', 'Motion JPEG AVI');
    WriterObj.FrameRate=30;
    open(WriterObj);
    for f = -30:329
        view([f, 17])
        drawnow;
        frame = getframe(fwireframe);
        writeVideo(WriterObj,frame);
    end
    close(WriterObj);
end






