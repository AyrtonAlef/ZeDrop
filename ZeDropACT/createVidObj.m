function [vidobj,vidsrc] = createVidObj(vidAdap,vidID,vidFmt)
%CREATEVIDOBJ Create a video object
%   INPUT
% vidAdap - Video adaptor (name of the adaptor used to commnunicate with the device)
% vidID - Device ID (numeric scalar value that identifies a particular device available through the specified adaptor)
% vidFmt - Video format (video format or name of device configuration file)
% OUTPUT
% vidobj - Video object

vidobj = videoinput(vidAdap,vidID,vidFmt); %Initializing communication with the camera
vidsrc = getselectedsource(vidobj); %Summary of video properties
vidsrc.Brightness = 15; %Increase brightness
vidobj.Timeout = Inf;
vidobj.ReturnedColorSpace = 'grayscale';
vidobj.Timeout = Inf;
vidobj.LoggingMode = 'disk'; %Image recording destination ('memory','disk','disk&memory')

end

