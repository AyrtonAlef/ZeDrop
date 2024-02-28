function vid2frame(filePath,vidName,imName,imFmt,tInterval)
%VID2FRAME Converts the video into n frames according to the provided
%interval (tInterval).
%   INPUT
% filePath - Path in which is the video file
% vidName - Video name
% imName - Name of the sequence of frames that will be recorded
% imFMt - Image format
% tInterval - Time interval in which frames should be taken (s)

fullfilevidPath = fullfile(filePath,vidName); %Full path of the video
vidobj = VideoReader(fullfilevidPath); %Read Video File
nFrames = vidobj.NumFrames; %Get total video number of frames
%duration = vidobj.Duration;
%nFrames = duration*frameRate;
frameRate = vidobj.FrameRate; %Get video framerate
frameInterval = tInterval * frameRate; %Frame interval between image acquisition
%ROI = [600 400 1200 300];
i = 1;
for k = 1: round(frameInterval): nFrames %Getting a image from every video second
    frame = read(vidobj,k);
    framegray = rgb2gray(frame);
    %framecrop = imcrop(framegray,ROI);
    fullimName = strcat(imName,'_',num2str(i),imFmt);
    fullfileimPath = fullfile(filePath,fullimName);
    imwrite(framegray,fullfileimPath);
    i = i + 1; %Interval of image name
end
end

vid2frame(filePath,fullvidName,fileName,imFmt,tCaptInterval);
