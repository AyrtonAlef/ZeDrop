function captImage(vidAdap,vidID,vidFmt,imFmt)
%CAPTIMAGE Capture an image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Capture an image of what is being captured by the camera. It offers the
%   user the option of choosing the directory and the name of the captured 
%   image.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT
% vidAdap - character vector that specifies the name of the adaptor used to communicate with the device
% vidID - numeric scalar value that identifies a particle device available through the specified adaptor
% vidFmt - character vector that specifies a particular video format supported by the device
% imFmt - format of the captured image

% Select the path where the user wants to save the captured image
fprintf('------------------------- CAPTURE IMAGE ---------------------------\n');
fprintf('Select the path to save the captured image... \n');
filePath = uigetdir('C:\'); %Selection of the path to save the captured image
name = input('Enter a name for the captured image: ','s');
fileName = strcat(name,imFmt);
fullfileName = fullfile(filePath,fileName);
%vidDiskFmt = '.avi'; %Video format ('.avi','.mp4','.mj2')
%imFmt = '.tif'; %Image format ('.jpeg','.png','.tif')
%vidComp = 'Motion JPEG AVI'; %Video compression ('Archival','Motion JPEG AVI','Motion JPEG 2000','MPEG-4','Uncompressed AVI','Indexed AVI' e 'Grayscale AVI')
%vidQual = 90; %Quaity of the video (related with the compression ratio)

%- Capture image 
fprintf('Please wait, capturing image...\n');
%-- Clean/Delete video object
close all %Close all open windows
%delete(vidobj) %Delete any video objects that may have been created
clear vidobj vidsrc %Clean objects related to the video
%-- Create video object
[vidobj, vidsrc] = createVidObj(vidAdap,vidID,vidFmt);
%-- Capture image
close all %Close all open windows
frame = getsnapshot(vidobj); %Capture an image
fprintf('Capture image completed!\n');

% Showing image
imshow(frame); %Show image
movegui(gcf,'center') %Center image window
set(gcf,'Name',fileName)

% Saving image
fprintf('Please wait, exporting image..\n');
imwrite(frame,fullfileName)
fprintf('Export image completed!\n');

pause(3);
end

