function HistologyCropper_PP_2021Jan3(histology_figure, save_folder, image_file_names, reference_size_um, save_file_name, ordinal)
global reference_size
% set up histology figure
ud_histology.file_num = 1;
ud_histology.num_files = length(image_file_names);
ud_histology.slice_num = ones(length(image_file_names),1);
ud_histology.save_file_name =save_file_name;
ud_histology.ax = axes;
ud_histology.ordinal = ordinal;

% ud_histology.hist_image = imread(fullfile(save_folder, [image_file_names{ud_histology.file_num}(1:end-4) '.tif']));
for ch = 1:3 %this doesn't work well for HB images yet. How to make it general?
    img(:,:,ch) = imread(fullfile(save_folder, [image_file_names{ud_histology.file_num}(1:end-4) '.tif']), 'tif', ch);
    ud_histology.hist_image(:,:,ch) = img(:,:,ch);
    img(:,:,ch) = imadjust(img(:,:,ch), stretchlim(img(:,:,ch))); % [0, 0.950], [0,1]);
    img(:,:,ch) = imadjust(img(:,:,ch), [0.005, 0.975]); %also saturate it a little bit
%      figure; imshow( img(:,:,ch), [])
end
imgAll = mean(img,3);
imgAll = cast(imgAll, 'uint8');
ud_histology.imgAll = imadjust(imgAll, stretchlim(imgAll));

ud_histology.seriesNumber = str2double(image_file_names{ud_histology.file_num}(end-4));
switch ud_histology.seriesNumber
    case 1
        reference_size = round(reference_size_um/1.2951);
    case 2
        reference_size = round(reference_size_um/3.8852);
    case 3
        reference_size = round(reference_size_um/11.6557);
end

figure(histology_figure); 
ud_histology.im = imshow(ud_histology.imgAll, []);
warning('off', 'MATLAB:colon:nonIntegerIndex');

ud_histology.cropped_slice_rect = {};

set(histology_figure, 'UserData', ud_histology);

% crop and switch image function for histology
set(histology_figure, 'KeyPressFcn', @(histology_figure,keydata)HistologyCropHotkeyFcn(histology_figure, keydata, save_folder, image_file_names, reference_size_um));

fprintf(1, '\n Controls: \n \n');
fprintf(1, 'click and drag: crop image \n');
fprintf(1, 'space: next image \n \n');

disp('select ROIs')
ud_histology = crop_and_save_image(ud_histology, histology_figure, save_folder, reference_size);

end

% --------------------------------------------------------------------------
% use imrect to crop, and then save the cropped image in 'processed' folder
% --------------------------------------------------------------------------
function ud_histology = crop_and_save_image(ud_histology, histology_figure, save_folder, reference_size)
global reference_size
try % get first slice ROI
    ud_histology = get(histology_figure, 'UserData');
    N_slicesInThisImage = ud_histology.slice_num(ud_histology.file_num);
%     ud_histology.file_num
    
    ud_histology.h(N_slicesInThisImage) = drawpoint(ud_histology.ax);
    pos = ud_histology.h(N_slicesInThisImage).Position;
    hroi = drawrectangle( 'Position', [pos(1) - floor(reference_size(2)/2) , pos(2) - floor(reference_size(1)/2), reference_size(2),  reference_size(1)], ...
        'DrawingArea', 'unlimited', 'InteractionsAllowed', 'translate', 'FaceSelectable', 1);
    pos = customWait(hroi);
    
    %% pad and crop - new 2021
    im =  ud_histology.hist_image;
    if pos(1) < 0
        % pad to the left and shift the position to 0 or 1
        im = padarray(im, [0 ceil(abs(pos(1))) ], 0, 'pre');
        pos(1) = 1; %not crazy precise -- it does not matter
    end
    if pos(2) < 0
        % pad above and shift up the position
        im = padarray(im, [ceil(abs(pos(2))) 0], 0, 'pre');
        pos(2) = 1; %not crazy precise -- it does not matter
    end
    % pad extra for convenience and that's it
    im = padarray(im, [ceil(pos(4)) ceil(pos(3))], 0, 'post'); 
    
    % now crop
    pos = round(pos);
    ud_histology.slice_image = im(pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3), :);
    
    
    imwrite(ud_histology.slice_image, fullfile(save_folder, 'processed', ...
        [ud_histology.save_file_name num2str(ud_histology.ordinal,'%.2d') '.' num2str(ud_histology.slice_num(ud_histology.file_num),'%.3d') '.tif']))
    disp([ud_histology.save_file_name num2str(ud_histology.ordinal,'%.2d') '.' num2str(ud_histology.slice_num(ud_histology.file_num),'%.3d') ' saved!'])
    
    ud_histology.slice_num(ud_histology.file_num) = ud_histology.slice_num(ud_histology.file_num) + 1;
    hroi.Visible = 'off';
    delete(hroi)
    
catch; ud_histology.slice_image = zeros(100,100,'uint8');
end

set(histology_figure, 'UserData', ud_histology);

try
    ud_histology = crop_and_save_image(ud_histology, histology_figure, save_folder, reference_size);
catch; disp('')
    end
end


function pos = customWait(hROI)

% Listen for mouse clicks on the ROI
l = addlistener(hROI,'ROIClicked',@clickCallback);
% Block program execution
uiwait;
% Remove listener
delete(l);
% Return the current position
pos = hROI.Position;
end

function clickCallback(~,evt)
if strcmp(evt.SelectionType,'double')
    uiresume;
end
end
% -------------------------------
% respond to keypress (space bar)
% -------------------------------
function HistologyCropHotkeyFcn(histology_figure, keydata, save_folder, image_file_names, reference_size_um)
global reference_size
ud_histology = get(histology_figure, 'UserData');

    switch lower(keydata.Key)    
        
    case 'space' % move onto next image
        
        if ud_histology.file_num + 1 <= ud_histology.num_files
            ud_histology.file_num = ud_histology.file_num + 1;
            set(histology_figure, 'UserData', ud_histology);

            ud_histology.hist_image = [];
            for ch = 1:3 
                img(:,:,ch) = imread(fullfile(save_folder, [image_file_names{ud_histology.file_num}(1:end-4) '.tif']), 'tif', ch);
                ud_histology.hist_image(:,:,ch) = img(:,:,ch);
                img(:,:,ch) = imadjust(img(:,:,ch), stretchlim(img(:,:,ch))); % [0, 0.950], [0,1]);
                img(:,:,ch) = imadjust(img(:,:,ch), [0.005, 0.975]); %also saturate it a little bit
                %      figure; imshow( img(:,:,ch), [])
            end
            imgAll = mean(img,3);
            imgAll = cast(imgAll, 'uint8');
            ud_histology.imgAll = imadjust(imgAll, stretchlim(imgAll));
            
            ud_histology.seriesNumber = str2double(image_file_names{ud_histology.file_num}(end-4));
            switch ud_histology.seriesNumber
                case 1
                    reference_size = round(reference_size_um/1.2951);
                case 2
                    reference_size = round(reference_size_um/3.8852);
                case 3
                    reference_size = round(reference_size_um/11.6557);
            end
            
            
            figure(histology_figure); 
            delete(ud_histology.ax);

            for i = 1:length(ud_histology.h)
                delete(ud_histology.h(i))
            end
            ud_histology.h = {};
            
            
            ud_histology.ax = axes;
            ud_histology.im = imshow(ud_histology.imgAll, []);
            
            set(histology_figure, 'UserData', ud_histology);
            
            ud_histology = crop_and_save_image(ud_histology, histology_figure, save_folder, reference_size);

        else
           disp('That was the last file -- close and move on to the next cell') 
        end
    end

try   
set(histology_figure, 'UserData', ud_histology); 
catch; disp('')
end


end
