function fig_table = tabulateImageInfo(image_file_names, save_folder, T)
%not a general purpose function

fig_table = table({''},{''}, {[0 0 ]}, 0, ...
    'VariableNames', {'name', 'dateModified', 'size', 'allenSlice'});

for fn = 1:length(image_file_names)
    figpath = fullfile(save_folder, image_file_names{fn});
    in = imfinfo(figpath);
    fig_table.name{fn} = image_file_names{fn};
    fig_table.dateModified{fn} = in.FileModDate;
    fig_table.size{fn} = [in.Width, in.Height];
    fig_table.allenSlice(fn) = T.allenSlice(fn);
end
end