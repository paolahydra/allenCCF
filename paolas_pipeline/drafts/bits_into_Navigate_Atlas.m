% bits
transformationType = 'pwl'; %'projective'; %can change to 'affine' or 'pwl'
f = AtlasTransformBrowser(f, tv_plot, av_plot, st, slice_figure_browser, folder_processed_images, probe_save_name_suffix, plane, transformationType); 

T = saveTransformTable(fullfile(folder_processed_images, 'transformations'), image_file_names, reference_size);

