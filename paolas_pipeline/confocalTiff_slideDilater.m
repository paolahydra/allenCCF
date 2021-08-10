function confocalTiff_slideDilater(image_folder, image_tag, save_folder, Transf)
    
if ~exist(save_folder, 'dir')
    mkdir(save_folder)
end
file_name_suffix = '_processed';

%% transf settings
% T.rotation = 90;
% T.reference_size_image = reference_size_originalImage;


%% name of images, in order anterior to posterior or vice versa
% once these are downsampled they will be named ['original name' '_processed.tif']
image_file_names = dir([image_folder filesep image_tag '*.tif']); % get the contents of the image_folder
image_file_names = natsortfiles({image_file_names.name});
% image_file_names = {'slide no 2_RGB.tif','slide no 3_RGB.tif','slide no 4_RGB.tif'}; % alternatively, list each image in order

%%
warning('off')
for f = 1:length(image_file_names)
    
    %% load the image - some channel transformation may occour, but no coordinate change here
    fprintf('Opening %s...\n',image_file_names{f})
    INFO = imfinfo(fullfile(image_folder, image_file_names{f}));
    nChannels = length(INFO);
    
    %PAOLA: for now only use the first three channels. Later add the
    %possibility to choose
    clear hist_image
    for ch = 1:min([3, nChannels])  % only use the first three channels (as the last one is often empty/noisy)
        hist_image(:,:,ch) = imread(fullfile(image_folder, image_file_names{f}), 'tif', ch); %this is the original image. Quality will be preserved.
    end
    if nChannels==2
        % add an extra channel to make an RGB image
        hist_image(:,:,3) = zeros(size(hist_image(:,:,1)));
    end

    % save some useful info for coordinate retrieval
    Transf.originalImage_RowCol_size = [INFO(1).Height, INFO(1).Width];
    
    %% rotate to standard coronal orientation (coordinate change!)
    J = imrotate(hist_image, Transf.rotation);

    %% do the dilating - we asume that the confocal image will always be smaller than the reference size, which will always be true unless wrong parameters
    assert(Transf.reference_originalImage_RowCol_size(1) > size(J, 1) && Transf.reference_originalImage_RowCol_size(2) > size(J, 2), 'check you pixel size!')
    % atlas_reference_size_um = microns_per_pixel_after_downsampling * atlas_reference_size;
    % T.reference_size_image = round(atlas_reference_size_um/microns_per_pixel);
    
    padSize_rows =    ceil( (Transf.reference_originalImage_RowCol_size(1) - size(J,1)) /2 );
    padSize_columns = ceil( (Transf.reference_originalImage_RowCol_size(2) - size(J,2)) /2 );

    J = padarray(J, [padSize_rows padSize_columns], 0);
        
    % now crop the excess out
    J = J(1:Transf.reference_originalImage_RowCol_size(1), 1:Transf.reference_originalImage_RowCol_size(2), :);
    
%     figure; 
%     imshow(J)

    %% downsample 
    
     J = imresize(J, [round(Transf.reference_originalImage_RowCol_size(1) * Transf.downsamplingFactor)  NaN]);
     % now crop the excess out
     J = J(1:Transf.atlas_reference_size(1), 1:Transf.atlas_reference_size(2), :);


    %% do the flipping separately and interactively 
    Transf.flipped = 0;
    Transf.furtherRotation = 0;
    
    %% save in the processed folder
    imwrite(J, fullfile(save_folder, [image_file_names{f}(1:end-4) file_name_suffix '.tif']))
    save(fullfile(save_folder,[image_file_names{f}(1:end-4) file_name_suffix '_transf.mat']), '-struct', 'Transf')
end
disp('Done.')

