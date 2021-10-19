function T_roi = Register_and_Tabulate_Rois_confocal(object_tag, image_tag, input_folder, av, st, tv_plot, microns_per_pixel, reference_size )
warning('check the output to see if the coordinates need flipping about a specific axis')

folder_processed_images = fullfile(input_folder, 'startingSingleSlices/processed');
transf_atlasreg_folder  = fullfile(input_folder, 'startingSingleSlices/processed/transformations');

%% set up INPUTS

dilateF = 50; %if too small it sometimes disappear during transformations. If too big I am afraid it might become less precise, but I am not sure about this.
        

%1. list all the registered images
d = dir (fullfile(transf_atlasreg_folder, '*_transform_data.mat'));
transfs = {d.name};
transfs = natsortfiles(transfs);

%2. list also all the preprocessing transformation files
d = dir (fullfile(folder_processed_images, '*_transf.mat'));
transfs_prepr = {d.name};
transfs_prepr = natsortfiles(transfs_prepr);

%3. list all csvs calculated on 'startingSingleSlices' images 
celldetection_csvs = dir([input_folder filesep '*csv']);
celldetection_csvs = natsortfiles({celldetection_csvs.name});
celldetection_csvs = natsortfiles(celldetection_csvs);


if (length(transfs) == length(transfs_prepr)) && (length(transfs) == length(celldetection_csvs))
    % all files have been processed. Matching them should be
    % straightforward...
    easymatching = 1;
else
    easymatching = 0;
end

k = zeros(3,1);
a = strfind(transfs{1}, image_tag);
if ~isempty(a), k(1) = 1; end
a = strfind(transfs_prepr{1}, image_tag);
if ~isempty(a), k(2) = 1; end
a = strfind(celldetection_csvs{1}, image_tag);
if ~isempty(a), k(3) = 1; end
if sum(k)<3
    warning('Your file naming system is inconsistent.')
    fprintf('WARNING. Your file naming system is inconsistent.\n  You should keep the same image tag in all files of the same brain.\n  You should use zeros before single-digit numbers so that your numbering is made of the same number of digits across files.\n\n  I will try to work with this anyways, but do check your output.\n')
    fprintf('\nLooking for image tag:\n  %s\n', image_tag);
    fprintf('\nHere are example filenames that I am trying to pair (It looks like you are not making my job easy):\n  %s\n  %s\n  %s\n',transfs{1},transfs_prepr{1},celldetection_csvs{1})
    fprintf('----------------------------------------------------\n\n')
    useImageTag = 0;
else
    useImageTag = 1;
end
%% either you really find the sorting variable part in the filename and use that one...
% find which number is the sorting one
clear NUMs colDouble
SORTERcol = 0;
for i = 1:length(transfs)
    NUMs(i,:) = regexpi(transfs{i}, '\d*', 'match');
end
for c = 1:size(NUMs,2)
    col = (NUMs(:,c));
    clear colDouble
    for ij = 1:length(col)
        colDouble(ij,1) = str2double(col(ij));
    end
    if length(unique((colDouble)))<length(colDouble)
        continue
    else
        [~,bi] = sort(colDouble);
        if isequal(bi(:)', 1:length(colDouble))
            %this is not the sorting element in natsortfiles (which we
            %assume we trust.)
            SORTERcol = c;
        else
            continue
        end
    end   
end
    



    
%% set up OUTPUTS
roiTable_name = fullfile(folder_processed_images, sprintf('%s%s_roiTable_All.csv',image_tag, object_tag));
if exist(roiTable_name, 'file')
    warning('roi_table file already exists and it will be overwritten')
end


%% loop through all the registered images and find the corresponding celldetection file, if any, and register ROIs
diary(fullfile(input_folder, sprintf('LOG_Register_and_Tabulate_Rois_%s', datestr(now,'YYMMDD_hhmmss'))))
fprintf('DilateF used:   %d \n\n', dilateF)
countIM = 0;
for i = 1:length(transfs)
    if useImageTag
        [~,filename,~] = fileparts(transfs{i});
        starts = regexpi(filename, image_tag, 'start'); %take also the image_tag together
        ends = regexpi(filename, '\d*', 'end');
        image_var_pos_i = starts : ends(SORTERcol)+1; %plus 1 deals with missing zeros in two-digit numbers.
        if starts > ends %what the heck of a filename did you make?
            ends = regexpi(filename, image_tag, 'end'); %take also the image_tag together
            starts = regexpi(filename, '\d*', 'start');
            image_var_pos_i = starts(SORTERcol) : ends;
        end
        vartag = transfs{i}(image_var_pos_i);

        %now check if you actually have matching files
        if sum(contains(transfs, vartag))==1 && sum(contains(celldetection_csvs, vartag))==1
            countIM = countIM+1;
            transf_file_prepr = fullfile(folder_processed_images, transfs_prepr{contains(transfs_prepr, vartag)});
            transf_file_atlas = fullfile(folder_processed_images, 'transformations', transfs{contains(transfs, vartag)});
            csv_file = fullfile(input_folder, celldetection_csvs{contains(celldetection_csvs, vartag)});
            
            fileroot = regexp(transfs{contains(transfs, vartag)}, '_preprocessed_transform_data.mat', 'split');
            fileroot = fileroot{1};
        else
            fprintf('I could not find matching transf and csv files for: \n   %s\n',transfs{i})
            disp('Skipping it')
            continue
        end
    elseif easymatching %just trust the sorting
        countIM = countIM+1;
        transf_file_prepr = fullfile(folder_processed_images, transfs_prepr{i});
        transf_file_atlas = fullfile(folder_processed_images, 'transformations', transfs{i});
        csv_file = fullfile(input_folder, celldetection_csvs{i});
        
        fileroot = regexp(transfs{i}, '_preprocessed_transform_data.mat', 'split');
        fileroot = fileroot{1};
    else %try without image_tag
        [~,filename,~] = fileparts(transfs{i});
        starts = regexpi(filename,  '\d*', 'start');
        ends = regexpi(filename, '\d*', 'end');
        image_var_pos_i = max(starts(SORTERcol)-3,1) : ends(SORTERcol)+1; %plus 1 deals with missing zeros in two-digit numbers.
        vartag = transfs{i}(image_var_pos_i);

        %now check if you actually have matching files
        if sum(contains(transfs, vartag))==1 && sum(contains(celldetection_csvs, vartag))==1
            countIM = countIM+1;
            transf_file_prepr = fullfile(folder_processed_images, transfs_prepr{contains(transfs_prepr, vartag)});
            transf_file_atlas = fullfile(folder_processed_images, 'transformations', transfs{contains(transfs, vartag)});
            csv_file = fullfile(input_folder, celldetection_csvs{contains(celldetection_csvs, vartag)});
            
            fileroot = regexp(transfs{contains(transfs, vartag)}, '_preprocessed_transform_data.mat', 'split');
            fileroot = fileroot{1};
        else
            fprintf('I could not find matching transf and csv files for: \n   %s\n',transfs{i})
            disp('Skipping it')
            continue
        end
    end
    % state the pairing
    [~, FP_trprep] = fileparts(transf_file_prepr);
    [~, FP_tratl] = fileparts(transf_file_atlas);
    [~, FP_csv] = fileparts(csv_file);
    fprintf('----------------------------------------------------\n')
    fprintf('registering files:\n    %s\n    %s\n    %s\n', FP_trprep, FP_tratl, FP_csv )
        
        %% read the preprocessing transformations in
        T1 = load(transf_file_prepr);
        %% read the csv file in
        TC = readtable(csv_file); %coordinates are in um
        x = TC.POSITION_X/microns_per_pixel; %in image coordinates
        y = TC.POSITION_Y/microns_per_pixel; %in image coordinates

        %% Kevin: you need to flip the x axis of your coordinates - this needs to be made more robust
%         warning('check the output to see if the coordinates need flipping about a specific axis')
%         x = T1.originalImage_RowCol_size(2) - x;
        
        
        %% read and parse the atlas transformation
        trmat = load(transf_file_atlas); %save_transform
        transform_data = trmat.save_transform;
        clear trmat
        
%         current_pointList_for_transform = transform_data.transform_points{1};
%         ud_slice.pointList = transform_data.transform_points{2};
        % load allen ref location
        slice_num = transform_data.allen_location{1};
        slice_angle = transform_data.allen_location{2};
        
        % set up transformed roi image
        ref = uint8(squeeze(tv_plot(slice_num,:,:)));
        R = imref2d(size(ref));
        
        % generate other necessary values
        bregma = allenCCFbregma(); % bregma position in reference data space
        atlas_resolution = 0.010; % mm
        offset_map = get_offset_map(slice_angle, reference_size);
         
        %% for each detected cell, make an image and reapply all the transformations, store new xy coords
        
        im0 = false(T1.originalImage_RowCol_size);
        roi_location = zeros(length(x),3);
        roi_annotation = cell(length(x),3);

        pixels_row = nan(length(x),1);
        pixels_column = nan(length(x),1);
        
        tic
        disp('    reapplying transformations to the detected cells...')
        parfor p_i = 1:length(x)
%              disp(p_i)
            X = round(x(p_i));
            Y = round(y(p_i));
            IM = im0;
            IM(Y-dilateF:Y+dilateF, X-dilateF:X+dilateF) = true;
            
            %% reapply all the preprocessing transformations
            % rotate to standard coronal orientation (coordinate change!)
            J = imrotate(IM, T1.rotation);
            
            % do the cropping/dilating
            if ( T1.reference_originalImage_RowCol_size(1) < size(J, 1) || T1.reference_originalImage_RowCol_size(2) < size(J, 2) )
                % first crop the image to the reference dimension
                
                % the following assumes that the amount to be cropped is small and that the
                % ROI is fairly centered within the image.
                
                cropSize_rows =    max([0, ceil( -1*(T1.reference_originalImage_RowCol_size(1) - size(J,1)) /2 )]);
                cropSize_columns = max([0, ceil( -1*(T1.reference_originalImage_RowCol_size(2) - size(J,2)) /2 )]);
                J(1:cropSize_rows, :, :) = [];
                J(end-cropSize_rows+1:end, :, :) = [];
                J(:,1:cropSize_columns, :) = [];
                J(:,end-cropSize_columns+1:end, :) = [];
            end
            % now dilate it if needed
            % you need to dilate the image to the reference size
            % atlas_reference_size_um = microns_per_pixel_after_downsampling * atlas_reference_size;
            % T.reference_size_image = round(atlas_reference_size_um/microns_per_pixel);
            
            padSize_rows =    ceil( (T1.reference_originalImage_RowCol_size(1) - size(J,1)) /2 );
            padSize_columns = ceil( (T1.reference_originalImage_RowCol_size(2) - size(J,2)) /2 ); 
            J = padarray(J, [padSize_rows padSize_columns], 0);
            
            % finally crop any excess out
            J = J(1:T1.reference_originalImage_RowCol_size(1), 1:T1.reference_originalImage_RowCol_size(2), :);

            
            % downsample
            J = imresize(J, [round(T1.reference_originalImage_RowCol_size(1) * T1.downsamplingFactor)  NaN]);
            % now crop the excess out
            J = J(1:T1.atlas_reference_size(1), 1:T1.atlas_reference_size(2), :);
            
            % flipping (and further rotation) done interactively
            if T1.flipped
                J = flip(J,2);
            end
            if T1.furtherRotation %this is currently disabled anyhow, bc I did not want to deal with cropping, especially from several now unsaved iterations of rotations
                J = imrotate(J,T1.furtherRotation,'nearest','crop'); %check that it works reliably
            end
            

            %% reapply the atlas reg transformations
            
            curr_slice_trans = imwarp(J, transform_data.transform, 'OutputView',R);
            if sum(curr_slice_trans(:))==0
                fprintf('cell #%d is being wrongly assigned to the center of the image. Increase dilateF.\n', p_i)
                error('Some ROIs are not correctly transformed. Increase the value assigned to dilateF.' )
            end
            rois = uint8(imregionalmax(curr_slice_trans));
            
            [p_row, p_column] = find(rois>0);
            %            pixels_row(p_i) = round(mean(p_row));
            %            pixels_column(p_i) = round(mean(p_column));
            pixels_row(p_i) = round(mean(p_row));
            pixels_column(p_i) = round(mean(p_column));
        end
        toc
        
        tic
        disp('    annotating transformed cells...') 
        % loop again through cells (this is faster outside of parloop)
        for p_i = 1:length(x)
            p_row = pixels_row(p_i);
            p_column = pixels_column(p_i);
            % get the offset from the AP value at the centre of the slice, due to
            % off-from-coronal angling
            offset = offset_map(p_row,p_column);
            
            % use this and the slice number to get the AP, DV, and ML coordinates
            ap = -(slice_num-bregma(1)+offset)*atlas_resolution;
            dv = (p_row-bregma(2))*atlas_resolution;
            ml = (p_column-bregma(3))*atlas_resolution;
            
            roi_location(p_i,:) = [ap dv ml];
            
            % finally, find the annotation, name, and acronym of the current ROI pixel
            ann = av(slice_num+offset,p_row,p_column);
            name = st.safe_name{ann};
            acr = st.acronym{ann};
            
            roi_annotation{p_i,1} = ann;
            roi_annotation{p_i,2} = name;
            roi_annotation{p_i,3} = acr;
            
        end
        toc
        
        disp('    saving results in ...roiTable_All.csv...')
        filename_origin = repmat(fileroot, length(x), 1);
        roi_table = table(roi_annotation(:,2),roi_annotation(:,3), ...
            roi_location(:,1),roi_location(:,2),roi_location(:,3), cat(1, roi_annotation{:,1}), cellstr(filename_origin), ...
            'VariableNames', {'name', 'acronym', 'AP_location', 'DV_location', 'ML_location', 'avIndex', 'roiFIle'});
        
        if countIM == 1
            %         writetable(roi_table, roiTable_name, 'WriteVariableNames', true, 'WriteMode', 'overwrite') %only available matlab 2021
            writetable(roi_table, roiTable_name, 'WriteVariableNames', true)
            T_roi = roi_table;
        else
            T_roi = readtable(roiTable_name);
            T_roi = cat(1, T_roi, roi_table);
            writetable(T_roi, roiTable_name, 'WriteVariableNames', true)
            %         writetable(roi_table, roiTable_name, 'WriteVariableNames', false, 'WriteMode', 'append') %only available matlab 2021
        end
        % check output 
        transfImage = imread(fullfile(transf_atlasreg_folder, sprintf('%s_preprocessed_transformed.tif',fileroot)));
        figure; 
        imshow(transfImage, []); title(fileroot, 'Interpreter', 'none');
        hold on
        scatter(pixels_column, pixels_row, 20, 'g','filled')
        savefig(fullfile(folder_processed_images, sprintf('%s_registeredCells',fileroot)))

end
diary off
end