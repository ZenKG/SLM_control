function [image, canvx, canvy, valm] = func_readImage

% possible image types to list
imagetypes = {'*.tif';'*.tiff'};
%imagetypes = [imagetypes ; upper(imagetypes)]
imtypstr = [imagetypes' ; repmat({';'},1,2)];
imtypstr = [imtypstr{:}];

% popup, ask for a file
[filename, filepath] = uigetfile({imtypstr,'TIF Image Files';'*.*','All Files' },'Select an TIF image');



% if cancel
if ~filename
    errordlg('Please select an image')   
else
    % read selected file
    image = double(imread(strcat(filepath, filename)));
    valm = max(max(image));
    [canvx,canvy] = size(image);    
end
