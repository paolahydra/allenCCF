function S = loadTabDataFromMultipleBrains(folders2Include, plotFWire)
% T_roi has already beed produced (Register_and_Tabulate_Rois). 
% Here T_rois will be loaded for multiple brains to be compared
% 
% specify inputs and settings:
% folders2Include = uipickfiles('FilterSpec', '/Users/galileo/dati/registered_brains_completed');
% plotFWire = 1;

%%
if plotFWire
    object_tag = 'green'; %'green' for rabies cells
    black_brain = false;
    microns_per_pixel = 3.8852;
    microns_per_pixel_after_downsampling = 10;
    fwireframe = []; %reset the figure
end
map = colorcet('R1', 'N', length(folders2Include));
% N = 5; map = colorcet('R1', 'N', N); image([1:N]), colormap(map) %helper to select a nice colormap

for i =1:length(folders2Include)
    folder2load = fullfile(folders2Include{i}, 'startingSingleSlices', 'processed');
    temp = dir(fullfile(folder2load, sprintf('*_%s_roiTable_All.csv', object_tag)));
    S(i).roiTable_name = fullfile(folder2load, temp.name);
    temp = regexp(temp.name, sprintf('%s_roiTable_All.csv', object_tag), 'split');
    S(i).save_file_name = temp{1};
    clear temp
    S(i).T_roi = readtable(S(i).roiTable_name);
    S(i).braincolor = map(i, :);
    if plotFWire
        fwireframe = plotWireFrame(S(i).T_roi, S(i).braincolor, black_brain, fwireframe, microns_per_pixel, microns_per_pixel_after_downsampling );
    end
end

end