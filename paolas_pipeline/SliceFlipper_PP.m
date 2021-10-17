function SliceFlipper_PP(slice_figure, folder_preprocessed_images, reference_size, gain)
% crop, sharpen, and flip slice images

processed_images = dir([folder_preprocessed_images filesep '*tif']);
ud.processed_image_names = natsortfiles({processed_images.name});
ud.total_num_files = size(processed_images,1); disp(['found ' num2str(ud.total_num_files) ' slice images']);

ud.slice_num = 1;
ud.rotate_angle = 0;
ud.flipped = 0;
ud.reference_size = reference_size;
ud.gain = gain; 

ud.grid = 0;
ud = load_next_slice(ud,folder_preprocessed_images);

ud.grid = zeros(size(ud.current_slice_image),class(ud.original_slice_image)); 
gridsize = round(size(ud.current_slice_image,1)/20);
ud.grid(1:gridsize:end,:,:) = 150 + 20000*(isa(ud.original_slice_image,'uint16')); 
ud.grid(:,1:gridsize:end,:) = 150 + 20000*(isa(ud.original_slice_image,'uint16')); 
if size(ud.current_slice_image,1)>3000 
    for xi = 1:5
        x1 = (1:gridsize:size(ud.current_slice_image,1)-xi)+xi;
        y1 = (1:gridsize:size(ud.current_slice_image,2)-xi)+xi;
        ud.grid(x1,:,:) = 150 + 20000*(isa(ud.original_slice_image,'uint16'));
        ud.grid(:,y1,:) = 150 + 20000*(isa(ud.original_slice_image,'uint16'));
    end
elseif size(ud.current_slice_image,1)>1500 
    for xi = 1:3
        x1 = (1:gridsize:size(ud.current_slice_image,1)-xi)+xi;
        y1 = (1:gridsize:size(ud.current_slice_image,2)-xi)+xi;
        ud.grid(x1,:,:) = 150 + 20000*(isa(ud.original_slice_image,'uint16'));
        ud.grid(:,y1,:) = 150 + 20000*(isa(ud.original_slice_image,'uint16'));
    end
end

if size(ud.current_slice_image*ud.gain,3) > 3
    figure
    for ch_id = 1:4
        subplot(2,2,ch_id), hold on
        title(sprintf('channel %d', ch_id))
        imshow(ud.current_slice_image(:,:,ch_id)*ud.gain);
    end
    answer = inputdlg('Please indicate max 3 channels to keep (e.g.: 1:3, or 1 2 4):');
    if contains(answer, ':')
       ud.channels2Keep = eval(answer{1});
    else
        ud.channels2Keep = str2num(answer{1});
    end
    image2show = ud.current_slice_image(:,:,ud.channels2Keep);
    ud.grid = ud.grid(:,:,1:length(ud.channels2Keep));
    close
else
    image2show = ud.current_slice_image;
end
imshow(image2show*ud.gain + ud.grid)
title(['Slice ' num2str(ud.slice_num) ' / ' num2str(ud.total_num_files)])
set(slice_figure, 'UserData', ud);

     

% key function for slice
set(slice_figure, 'KeyPressFcn', @(slice_figure,keydata)SliceCropHotkeyFcn(keydata, slice_figure, folder_preprocessed_images));
% scroll function for slice
set(slice_figure, 'WindowScrollWheelFcn', @(src,evt)SliceScrollFcn(slice_figure, evt))


fprintf(1, '\n Controls: \n \n');
fprintf(1, 'right: save and see next image \n');
fprintf(1, 'left: save and see previous image \n');
fprintf(1, 'scroll: rotate slice \n');
fprintf(1, 's: sharpen \n');
fprintf(1, 'g: toggle grid \n');
fprintf(1, 'c: crop slice further \n');
fprintf(1, 'f: flip horizontally \n');
fprintf(1, 'w: switch order (move image forward) \n');
fprintf(1, 'r: reset to original \n');
fprintf(1, 'n: go to slice num... \n');
fprintf(1, 'i: increase gain... \n');
fprintf(1, 'd: decrease gain... \n');
fprintf(1, 'v: change channel mode view \n'); % under construction
fprintf(1, 'delete: delete current image \n');

% --------------------
% respond to keypress
% --------------------
function SliceCropHotkeyFcn(keydata, slice_figure, folder_preprocessed_images)

ud = get(slice_figure, 'UserData');



switch lower(keydata.Key)   
    case 'leftarrow' % save and previous slice
%         imwrite(ud.current_slice_image, fullfile(folder_preprocessed_images, ud.processed_image_name))   
        fname = fullfile(folder_preprocessed_images, ud.processed_image_name);
        imwrite(ud.current_slice_image(:,:,1), fname);
        for i = 2:size(ud.current_slice_image, 3)
            imwrite(ud.current_slice_image(:,:,i), fname, 'WriteMode', 'append');
        end
        
        rotate_angle = ud.rotate_angle;
        flipped = ud.flipped;
        save(fullfile(folder_preprocessed_images, sprintf('%s_transf.mat',ud.processed_image_name)), 'rotate_angle', 'flipped');
        ud.slice_num = ud.slice_num - 1*(ud.slice_num>1);
        ud = load_next_slice(ud,folder_preprocessed_images);

    case 'rightarrow' % save and next slice      
%         imwrite(ud.current_slice_image, fullfile(folder_preprocessed_images, ud.processed_image_name))
        fname = fullfile(folder_preprocessed_images, ud.processed_image_name);
        imwrite(ud.current_slice_image(:,:,1), fname);
        for i = 2:size(ud.current_slice_image, 3)
            imwrite(ud.current_slice_image(:,:,i), fname, 'WriteMode', 'append');
        end
        
        rotate_angle = ud.rotate_angle;
        flipped = ud.flipped;
        save(fullfile(folder_preprocessed_images, sprintf('%s_transf.mat',ud.processed_image_name)), 'rotate_angle', 'flipped');
        ud.slice_num = ud.slice_num + 1*(ud.slice_num < length(ud.processed_image_names));
        ud = load_next_slice(ud,folder_preprocessed_images);
        
    case 'delete' % delete and previous slice
        delete(fullfile(folder_preprocessed_images, ud.processed_image_name));
        ud.slice_num = ud.slice_num - 1*(ud.slice_num>1);  
        
        processed_images = dir([folder_preprocessed_images filesep '*tif']);
        ud.processed_image_names = natsortfiles({processed_images.name});
        ud.total_num_files = size(processed_images,1); disp(['found ' num2str(ud.total_num_files) ' processed slice images']);
        ud = load_next_slice(ud,folder_preprocessed_images);
    
    case 'n' %save and go to slice num...
        fname = fullfile(folder_preprocessed_images, ud.processed_image_name);
        imwrite(ud.current_slice_image(:,:,1), fname);
        for i = 2:size(ud.current_slice_image, 3)
            imwrite(ud.current_slice_image(:,:,i), fname, 'WriteMode', 'append');
        end
        
        rotate_angle = ud.rotate_angle;
        flipped = ud.flipped;
        save(fullfile(folder_preprocessed_images, sprintf('%s_transf.mat',ud.processed_image_name)), 'rotate_angle', 'flipped');
        
        ud.slice_num = input('Go to slice num: ');
        ud = load_next_slice(ud,folder_preprocessed_images);
        
    case 'i' %increase gain
        ud.gain = ud.gain + 1;
        disp(ud.gain)
        ud = load_next_slice(ud,folder_preprocessed_images);
        
    case 'd' %decrease gain
        ud.gain = ud.gain - 1;
        disp(ud.gain)
        ud = load_next_slice(ud,folder_preprocessed_images);
        
    case 'g' % grid
        if sum(ud.grid(:)) == 0
            ud.grid = zeros(size(ud.current_slice_image),class(ud.original_slice_image));
            gridsize = round(size(ud.current_slice_image,1)/20);
            ud.grid(1:gridsize:end,:,:) = 150 + 20000*(isa(ud.original_slice_image,'uint16'));
            ud.grid(:,1:gridsize:end,:) = 150 + 20000*(isa(ud.original_slice_image,'uint16'));
            if size(ud.current_slice_image,1)>3000
                for xi = 1:5
                    x1 = (1:gridsize:size(ud.current_slice_image,1)-xi)+xi;
                    y1 = (1:gridsize:size(ud.current_slice_image,2)-xi)+xi;
                    ud.grid(x1,:,:) = 150 + 20000*(isa(ud.original_slice_image,'uint16'));
                    ud.grid(:,y1,:) = 150 + 20000*(isa(ud.original_slice_image,'uint16'));
                end
            elseif size(ud.current_slice_image,1)>1500
                for xi = 1:3
                    x1 = (1:gridsize:size(ud.current_slice_image,1)-xi)+xi;
                    y1 = (1:gridsize:size(ud.current_slice_image,2)-xi)+xi;
                    ud.grid(x1,:,:) = 150 + 20000*(isa(ud.original_slice_image,'uint16'));
                    ud.grid(:,y1,:) = 150 + 20000*(isa(ud.original_slice_image,'uint16'));
                end
            end
        else
            ud.grid = zeros(size(ud.current_slice_image),class(ud.original_slice_image)); 
        end
    case 'c' % 
        cropped_slice_rect = imrect;
        slice_position = cropped_slice_rect.getPosition;
        try
            ud.current_slice_image = ud.current_slice_image(slice_position(2):slice_position(2)+slice_position(4),slice_position(1):slice_position(1)+slice_position(3),:);          
        catch; disp('crop out of bounds'); 
        end        
        
        try; ud.current_slice_image = padarray(ud.current_slice_image, [...
                                floor((ud.size(1) - size(ud.current_slice_image,1)) / 2) + ...
                                mod(size(ud.current_slice_image,1),2),...
                                floor((ud.size(2) - size(ud.current_slice_image,2)) / 2) + ...
                                mod(size(ud.current_slice_image,2),2)], 0);
            ud.original_ish_slice_image = ud.current_slice_image;                            
        catch; disp('cropping failed');
        end              
        
        ud.size = size(ud.current_slice_image); 
        ud.grid = imresize(ud.grid, ud.size(1:2));         
        
	case 's' % sharpen
%         ud.current_slice_image = localcontrast(ud.current_slice_image);
%         ud.original_ish_slice_image = localcontrast(ud.original_ish_slice_image);
%         
    case 'w' % switch order
        if ud.slice_num < length(ud.processed_image_names)
            disp('switching order -- moving this image forward')
            next_processed_image_name = ud.processed_image_names{ud.slice_num+1};
%             next_slice_image = imread(fullfile(folder_preprocessed_images, next_processed_image_name));
            fname = fullfile(folder_preprocessed_images, next_processed_image_name);
            INFO = imfinfo(fname);
            nChannels = length(INFO);
            clear A
            for ch = 1:nChannels  % only use the first three channels (as the last one is often empty/noisy)
                A(:,:,ch) = imread(fname, 'tif', ch); %this is the original image. Quality will be preserved.
            end
            next_slice_image = A;
    
%             imwrite(next_slice_image, fullfile(folder_preprocessed_images, ud.processed_image_name))  
            fname = fullfile(folder_preprocessed_images, ud.processed_image_name);
            imwrite(ud.current_slice_image(:,:,1), fname);
            for i = 2:size(ud.current_slice_image, 3)
                imwrite(next_slice_image(:,:,i), fname, 'WriteMode', 'append');
            end
            
%             imwrite(ud.current_slice_image, fullfile(folder_preprocessed_images, next_processed_image_name))
            fname = fullfile(folder_preprocessed_images, next_processed_image_name);
            imwrite(ud.current_slice_image(:,:,1), fname);
            for i = 2:size(ud.current_slice_image, 3)
                imwrite(ud.current_slice_image(:,:,i), fname, 'WriteMode', 'append');
            end
            
            ud.current_slice_image = next_slice_image; 
            ud.size = size(ud.current_slice_image); 
            ud.grid = imresize(ud.grid, ud.size(1:2));             
        end
        
    case 'f' % flip horizontally
        ud.current_slice_image = flip(ud.current_slice_image,2);
        ud.original_ish_slice_image = flip(ud.original_ish_slice_image,2);
        ud.flipped = ~ud.flipped;
    case 'r' % return to original image
        ud.current_slice_image = ud.original_slice_image;
        ud.original_ish_slice_image = ud.original_slice_image;
        ud.size = size(ud.current_slice_image); 
        ud.grid = imresize(ud.grid, ud.size(1:2)); 
        ud.rotate_angle = 0;
        ud.flipped = 0;

end



% in all cases, update image and title
if size(ud.current_slice_image*ud.gain,3) > 3
    image2show = ud.current_slice_image(:,:,ud.channels2Keep);
    ud.grid = ud.grid(:,:,1:length(ud.channels2Keep));
else
    image2show = ud.current_slice_image;
end
imshow(image2show*ud.gain + ud.grid)
title(['Slice ' num2str(ud.slice_num) ' / ' num2str(ud.total_num_files)])


set(slice_figure, 'UserData', ud);

function ud = load_next_slice(ud,folder_processed_images)
    ud.processed_image_name = ud.processed_image_names{ud.slice_num};
%     ud.current_slice_image = imread(fullfile(folder_processed_images, ud.processed_image_name));
    fname = fullfile(folder_processed_images, ud.processed_image_name);
    INFO = imfinfo(fname);
    nChannels = length(INFO);
    
    
%     clear A
%     for ch = 1:nChannels  
%         A(:,:,ch) = imread(fname, 'tif', ch); %this is the original image. Quality will be preserved.
% %         A(ch).image = imread(fname, 'tif', ch); %this is the original image. Quality will be preserved.
%     end
% %     figure
% %     for ch_id = 1:nChannels
% %         subplot(2,2,ch_id), hold on
% %         title(sprintf('channel %d', ch_id))
% %         imshow(A(ch_id).image*ud.gain);
% %     end
%     ud.current_slice_image = A;
    
    
    
   
    
    %PAOLA: for now only use the first three channels. Later add the
    %possibility to choose
    clear hist_image
    if nChannels == 1
        hist_image = imread(fname, 'tif');
    else
        for ch = 1:min([3, nChannels])  % only use the first three channels (as the last one is often empty/noisy)
            hist_image(:,:,ch) = imread(fname, 'tif', ch); %this is the original image. Quality will be preserved.
        end
        if nChannels==2
            % add an extra channel to make an RGB image
            hist_image(:,:,3) = zeros(size(hist_image(:,:,1)));
        end
    end
    ud.current_slice_image = hist_image;
    
    
    
    
    
    disp(['loaded ' ud.processed_image_name])
    
    
    % pad if possible (if small enough)
    try; ud.current_slice_image = padarray(ud.current_slice_image, [floor((ud.reference_size(1) - size(ud.current_slice_image,1)) / 2) + mod(size(ud.current_slice_image,1),2) ...
                                                  floor((ud.reference_size(2) - size(ud.current_slice_image,2)) / 2) + mod(size(ud.current_slice_image,2),2)],0);
    end 
    
    % crop out any odd row or column - to prevent mismatches in case of
    % further cropping
    if mod(size(ud.current_slice_image,1),2)
        ud.current_slice_image(end,:,:) = [];
    end
    if mod(size(ud.current_slice_image,2),2)
        ud.current_slice_image(:,end,:) = [];
    end

    ud.original_slice_image = ud.current_slice_image;             
    ud.original_ish_slice_image = ud.current_slice_image;        

    ud.size = size(ud.current_slice_image); 
    if ud.size(1) > ud.reference_size(1)+1 || ud.size(2) > ud.reference_size(2)+2
        disp(['Slice ' num2str(ud.slice_num) ' / ' num2str(ud.total_num_files) ' is ' num2str(ud.size(1)) 'x' num2str(ud.size(2)) ' pixels:']);
%         disp(['I suggest you crop this image down to under ' num2str(ud.reference_size(1)) ' x ' num2str(ud.reference_size(2)) ' pxl'])
    end        
    ud.grid = imresize(ud.grid, ud.size(1:2)); 
    ud.rotate_angle = 0;
    
%     imwrite(ud.current_slice_image, fullfile(folder_processed_images, ud.processed_image_name)) 
    fname = fullfile(folder_processed_images, ud.processed_image_name);
    imwrite(ud.current_slice_image(:,:,1), fname);
    for i = 2:size(ud.current_slice_image, 3)
        imwrite(ud.current_slice_image(:,:,i), fname, 'WriteMode', 'append');
    end


% function to rotate slice by scrolling
function SliceScrollFcn(fig, evt)

ud = get(fig, 'UserData');

%modify based on scrolling
ud.rotate_angle = ud.rotate_angle + evt.VerticalScrollCount*.75;

ud.current_slice_image = imrotate(ud.original_ish_slice_image,ud.rotate_angle,'nearest','crop');
if size(ud.current_slice_image*ud.gain,3) > 3
    image2show = ud.current_slice_image(:,:,ud.channels2Keep);
    ud.grid = ud.grid(:,:,1:length(ud.channels2Keep));
else
    image2show = ud.current_slice_image;
end
imshow(image2show*ud.gain + ud.grid)
title(['Slice ' num2str(ud.slice_num) ' / ' num2str(ud.total_num_files)])

set(fig, 'UserData', ud);



