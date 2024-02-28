function [incMesaTotal] = routineContIncDrop(vidAdap,vidID,vidFmt,imFmt,smInc,passoMotor,nEntradasSemFim,nDentesCoroa,wMin,wMax,framerate,limTotalFrames,incMesaTotal)
%ROUTINECONTINCDROP Execution of continuous inclined drop experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Execution of continuous inclined drop experiment. Performs a continuous
%   tilt of the stage until reaching a pre-established maximum tilt. During
%   tilting, drop images are captured allowing further analysis. The 
%   routine requires the user to define the initial stage tilt (incIni), 
%   the initial stage tilt rate (taxaIncIni), the maximum stage tilt 
%   (incMax), the stage tilt rate during continuous tilting (taxaIncEns), 
%   the waiting time for the drop equilibrium (tEsp) and the number of 
%   images to be captured (nImages). The routine allows choosing the 
%   directory where the images will be saved, as well as defining the name 
%   of the images. At the end, a text file is exported with all relevant 
%   parameters of the experiment. *Requires the connection of the tilting 
%   stage system.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT
% vidAdap - Character vector that specifies the name of the adaptor used to communicate with the device
% vidID - Numeric scalar value that identifies a particle device available through the specified adaptor 
% vidFmt - Character vector that specifies a particular video format supported by the device
% imFmt - Format of the captured image 
% smInc - Syringe motor control
% passoMotor - Number of steps required for one revolution in the motor [steps/rev]
% nEntradasSemFim - Worm number of starts
% nDentesCoroa - Number of worm gear teeth
% wMin - Minimum stable rotation speed (RPM) in which the motor can reach. Corresponds to a sm.RPM = 2 (Microstep)
% wMax - Maximum stable rotation speed (RPM) in which the motor can reach
% framerate - Actual camera frame rate
% incMesaTotal - Current stage tilt [°]
%   OUTPUT
% incMesaTotal - Current stage tilt [°]

fprintf('------------------- INCLINED DROP - CONTINUOUS EXPERIMENT --------------------\n');
% Experiment parameters (entered by the user)
fprintf('Experiment parameters: \n');
taxaIncMax = (wMax*360*nEntradasSemFim)/(nDentesCoroa*60); %Maximum stage tilt rate in °/s
taxaIncMin = (wMin*360*nEntradasSemFim)/(nDentesCoroa*60); %Minimum stage tilt rate in °/s
incStart = incMesaTotal; %Current stage tilt [°] 
incIni = input("- Initial tilt (°): "); %Initial stage tilt before starting the measurement series
prompttaxaIncIni = strcat("- Initial stage tilt rate (°/s) (max. ",num2str(taxaIncMax,'%.2f')," °/s): ");
while 1
    taxaIncIni = input(prompttaxaIncIni); %Initial stage tilt rate in °/s
    if taxaIncIni > taxaIncMax %Check that the stage tilt does not exceed the maximum motor rotation speed
        fprintf("Unable to proceed with the routine. Please decrease the initital tilt rate. \n")
    else
        break;
    end
end
incMax = input("- Maximum tilt (°): "); %Maximum tilt reached by the stage
prompttaxaInc = strcat("- Stage tilt rate (°/s) (max. ",num2str(taxaIncMax,'%.2f')," °/s): ");
while 1
    %taxaInc = input("- Stage tilt rate (°/s) (max. %.2f °/s): "); %Stage tilt rate in °/s
    taxaIncEns = input(prompttaxaInc); %Stage tilt rate in °/s
    if taxaIncEns > taxaIncMax %Check that the stage tilt rate does not exceed the maximum motor rotation speed
        fprintf("Unable to continue the routine. Please lower the stage tilt rate. \n")
    else
        break;
    end
end
sentidoInc = input("- Tilt direction (clockwise (CW)/ counterclockwise (CCW)) [CW]: ", "s");
if isempty(sentidoInc) || (sentidoInc ~= "CW") || (sentidoInc ~= "CCW")
    sentidoInc = "CW";
end
if sentidoInc == "CW"
    sentido = -1;
    sentidoRet = "CCW";
elseif sentidoInc == "CCW"
    sentido = 1;
    sentidoRet = "CW";
end
tEspera = input("- Waiting time for drop equilibrium (s): "); %Waiting time for drop equilibrium
nImages = input("- Number of captured images: "); %Number of captured images during tilting

% Experiment parameters (calculated)
incVarIni = abs(incIni - incMesaTotal); %Initial tilt variation
incVarEns = abs(incMax - incIni); %Tilt variation between each measurement
incVarRet = abs(incMax - incStart); %Tilt variation to return the stage to its initial tilt position
tEnsaio = incVarEns/taxaIncEns; %Total experiment time (increase in stage tilting) in s
%tCaptInterval = tEnsaio/nImages; %Time interval between captured images in s
%volCaptInterval = taxaDepGota*tCaptInterval; %Volume interval between captured images in uL
totalFrame = ceil(framerate*tEnsaio); %Total number of frames acquired during the experiment
frameInterval = round(totalFrame/nImages); %Interval between frames to acquire the desired number of images
tCaptInterval = frameInterval/framerate; %Time interval between captured images in s
incCaptInterval = taxaIncEns*tCaptInterval; %Inclination interval between captured images in uL
erroIncTotal = 0;

%- Drive motor settings
nPassosIncIni = inc2passo(incVarIni,passoMotor,nEntradasSemFim,nDentesCoroa); %No of steps corresponding to the initial stage tilt 
nPassosIncEns = inc2passo(incVarEns,passoMotor,nEntradasSemFim,nDentesCoroa); %No of steps corresponding to the stage tilt 

%- Image capture settings
fprintf('Path selection to export video and images... \n');
filePath = uigetdir('C:\'); %Path selection to export video and images
name = input('Enter a name for the video and image series: ','s');
vidDiskFmt = '.avi'; %Video format ('.avi','.mp4','.mj2')
%imFmt = '.tif'; %Images format ('.jpeg','.png','.tif')
vidComp = 'Motion JPEG AVI'; %Video compression ('Archival','Motion JPEG AVI','Motion JPEG 2000','MPEG-4','Uncompressed AVI','Indexed AVI' e 'Grayscale AVI')
vidQual = 90; %Video quality (related with the compression ratio)

%- Begin of experiment
fprintf('--------------------- EXECUTION OF TILTING ROUTINE ----------------------\n');
fprintf('Press any key to proceed with the tilting routine... \n');
pause;
tic
close all %Close all open windows
%delete(vidobj) %Delete any video objects that may have been created
clear vidobj vidsrc %Clear video related objects
clockEnsaio1 = toc;

%- Initial stage tilt
tic
fprintf('Please wait, initial stage tilt... \n');
if taxaIncIni <= taxaIncMin %Tilt rate less than the minimum motor rotation speed RPM. Motor will work in pulses
    smInc.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a tilt rate = 0,0818 uL/s
    wMotorEns = incs2RPM(taxaIncIni,nEntradasSemFim,nDentesCoroa); %Motor rotation speed (RPM) in order to achieved the desired tilt rate (°/s)
    tadpasso = tadicionalpasso(passoMotor,wMotorEns,wMin); %Time added to each step (in seconds) to reach the desired motor rotational speed
    for i = 1:floor(nPassosIncIni) %Execution of an additional tilt on the stage
        move(smInc,sentido*1); %Executes 1 rotation step on the motor
        pause(tadpasso);
    end
elseif (taxaIncMin < taxaIncIni) && (taxaIncIni <= taxaIncMax)
    wRPM = fincs2RPM(taxaIncIni,nEntradasSemFim,nDentesCoroa); %Calculation of the rotational speed [RPM] necessary to be sent to the motor to reach the desired tilt rate [°/s]
    smInc.RPM = floor(wRPM);
    move(smInc,sentido*floor(nPassosIncIni)); %Motor movement 
end
release(smInc); %Motor release
fprintf('Initial stage tilt completed! \n');
erroIncIni = passo2inc(nPassosIncIni - floor(nPassosIncIni),passoMotor,nEntradasSemFim,nDentesCoroa); %Error in stage tilt
erroIncTotal = erroIncTotal + erroIncIni;
incMesaTotal = incMesaTotal - sentido*incVarIni; %Update of total stage tilt
clockIncIni = toc;

%- Waiting time for drop equilibrium
tic
fprintf('Please wait, waiting time for drop reach equilibrium... \n');
pause(tEspera);
fprintf('Waiting time completed! \n');
clockEspGota = toc;

%- Create video object
tic
[vidobj, vidsrc] = createVidObj(vidAdap,vidID,vidFmt);

%- Update of image names
%fileName = strcat(name,'_',num2str(k),'_V(uL)_',num2str((volGotaIni+k*volGotaVar),'%.2f'));
fileName = name;
fullvidName = strcat(fileName,vidDiskFmt);
fullfilevidPath = fullfile(filePath,fullvidName);

%- Preparation of video object
diskLogger = VideoWriter(fullfilevidPath,vidComp);
diskLogger.FrameRate = framerate; %Matching the video frame rate to the actual camera frame rate
diskLogger.Quality = vidQual; %Definition of the video quality
vidobj.DiskLogger = diskLogger;
if totalFrame <= limTotalFrames %Continuous frame capture during video recording
    vidobj.FramesPerTrigger = totalFrame; %Total frames to acquire
    vidobj.FrameGrabInterval = 1; %Frame interval
else %Capturing spaced frames during video recording
    vidobj.FramesPerTrigger = nImages; %Total frames to acquire
    vidobj.FrameGrabInterval = frameInterval; %Frame interval
end
clockPrepVid = toc; %Time for creating and preparing the video object

%- Start video recording
tic
fprintf('Start video recording. \n');
start(vidobj) %Start video recording
    
%- Continuous increase of stage tilt
fprintf('Please wait, continuous increase of stage tilt... \n');
if taxaIncEns <= taxaIncMin %Stage tilt rate less than minimum motor rotation speed. Motor will work in pulses
    smInc.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a tilt rate = 0,0818 uL/s
    wMotorEns = incs2RPM(taxaIncEns,nEntradasSemFim,nDentesCoroa); %Motor rotation speed (RPM) in order to achieved the desired tilt rate (°/s)
    tadpasso = tadicionalpasso(passoMotor,wMotorEns,wMin); %Time added to each step (in seconds) to reach the desired motor rotational speed
    for i = 1:floor(nPassosIncEns) %Execution of an additional tilt on the stage
        move(smInc,sentido*1); %Executes 1 rotation step on the motor
        pause(tadpasso);
    end
elseif (taxaIncMin < taxaIncEns) && (taxaIncEns <= taxaIncMax)
    wRPM = fincs2RPM(taxaInc,nEntradasSemFim,nDentesCoroa); %Calculation of the rotational speed [RPM] necessary to be sent to the motor to reach the desired tilt rate [°/s]
    taxaIncEns = fRPM2incs(ceil(wRPM)); %Update of tilt rate
    % Update variables
    tEnsaio = incVarEns/taxaIncEns;
    %tCaptInterval = tEnsaio/nImages;
    %volCaptInterval = taxaDepGota*tCaptInterval;
    totalFrame = ceil(framerate*tEnsaio);
    frameInterval = round(totalFrame/nImages);
    tCaptInterval = frameInterval/framerate;
    incCaptInterval = taxaIncEns*tCaptInterval;
    %vidobj.FramesPerTrigger = totalFrame;
    % Syringe plunger movement
    %smSeringa.RPM = floor(wRPM);
    smInc.RPM = round(wRPM);
    move(smInc,floor(nPassosIncEns)); %Motor movement
end
release(smInc); %Motor release
fprintf('Continuous increase of stage tilt completed! \n');
erroIncEns = passo2inc(nPassosIncEns - floor(nPassosIncEns),passoMotor,nEntradasSemFim,nDentesCoroa); %Error in tilting during the experiment
erroIncTotal = erroIncTotal + erroIncEns; %Update of total error in tilting 
incMesaTotal = incMesaTotal - sentido*incVarEns; %Update of total stage tilt
    
%- Finishing video recording
wait(vidobj) %Wait until all frames are captured
stop(vidobj);
%-- Wait for all frames to be written to disk
while vidobj.FramesAcquired ~= vidobj.DiskLoggerFrameCount
    pause(.1);
end
fprintf('Video recording completed! \n');
clockVid = toc; %Time for creating, recording and exporting video on disk
    
%- Showing an example image
tic
close all %Close all open windows
frame = getsnapshot(vidobj); %Capture an image
imshow(frame); %Show image
movegui(gcf,'center') %Centralize image window
clockEnsaio2 = toc;
    
%- Converting video to images
tic
fprintf('Please wait, converting video to images... \n');
if totalFrame <= limTotalFrames %Continuous frame capture during video recording
    fullfilevidPath = fullfile(filePath,fullvidName); %Full path of the video
    vidFile = VideoReader(fullfilevidPath); %Read Video File
    i = 1;
    for k = 1: round(frameInterval): totalFrame %Get an image from every video second
        frame = read(vidFile,k);
        framegray = rgb2gray(frame);
        if i < 10 %Modifying the file name prefix to arrange files in ascending tilt order
            prefixo = "00";
        elseif (i >= 10) && (i < 100)
            prefixo = "0";
        elseif (i >= 100) && (i < 1000)
            prefixo = "";
        end
        %fullimName = strcat(prefixo,num2str(i),'_',fileName,'_V(uL)_',num2str((volGotaIni + k*volCapIntervalframe),'%.2f'),imFmt);
        fullimName = strcat(prefixo,num2str(i),'_',fileName,'_alpha(°)_',num2str((incIni + i*incCaptInterval),'%.2f'),imFmt);
        fullfileimPath = fullfile(filePath,fullimName);
        imwrite(framegray,fullfileimPath);
        i = i + 1; %Interval of image name
    end
else %Capturing spaced frames during video recording
    fullfilevidPath = fullfile(filePath,fullvidName); %Full path of the video
    vidFile = VideoReader(fullfilevidPath); %Read Video File
    %i = 1;
    for k = 1: nImages %Get an image from every video frame
        frame = read(vidFile,k);
        framegray = rgb2gray(frame);
        if k < 10 %Modifying the file name prefix to arrange files in ascending tilt order
            prefixo = "00";
        elseif (k >= 10) && (k < 100)
            prefixo = "0";
        elseif (k >= 100) && (k < 1000)
            prefixo = "";
        end
        %fullimName = strcat(fileName,'_',num2str(i),'_V(uL)_',num2str((volGotaIni + k*volCaptInterval),'%.2f'),imFmt);
        fullimName = strcat(prefixo,num2str(k),'_',fileName,'_alpha(°)_',num2str((incIni + k*incCaptInterval),'%.2f'),imFmt);
        fullfileimPath = fullfile(filePath,fullimName);
        imwrite(framegray,fullfileimPath);
        %i = i + 1; %Interval of image name
    end
end
delete(fullfilevidPath) %Delete video
fprintf('Conversion of video to images completed! \n');
clockImages = toc; %Time for converting video to images and export images to disk

%- Deleting video object
tic
delete(vidobj)
clear vidobj vidsrc
clockEnsaio3 = toc;

%- Return to the starting stage tilt 
tic
fprintf('Please wait, return stage to the starting tilt... \n');
moverIncMesa(smInc,passoMotor,nEntradasSemFim,nDentesCoroa,incVarRet,sentidoRet,"high");
fprintf('Return of stage completed! \n');
incMesaTotal = incMesaTotal - incVarRet; %Update of total stage tilt
clockIncRet = toc;

clockEnsaio = clockEnsaio1 + clockIncIni + clockEspGota + clockPrepVid + clockVid + clockEnsaio2 + clockImages + clockEnsaio3 + clockIncRet; %Update experiment total time 

%- Export/save experiment description
filenametxt = "Description.txt";
fullfilenametxt = fullfile(filePath,filenametxt);
varnames = {'Datetime','incIni[°]','taxaIncIni[°/s]','incMax[°]','taxaIncEns[°/s]','sentidoInc','tEspera[s]','nImages','tEnsaio[s]','frameInterval','tCaptInterval[s]','incCaptInterval[°]','erroIncIni[°]','erroIncEns[°]','errroIncTotal[°]','clockIncIni[s]','clockEspGota[s]','clockPrepVid[s]','clockVid[s]','clockImages[s]','clockIncRet[s]','clockEnsaio[s]','frameRate','VidCompression','VidQuality'}; %Variables name in txt file
T = table(datetime,incIni,taxaIncIni,incMax,round(taxaIncEns,4),sentidoInc,tEspera,nImages,round(tEnsaio,2),frameInterval,round(tCaptInterval,2),round(incCaptInterval,4),round(erroIncIni,4),round(erroIncEns,4),round(erroIncTotal,4),round(clockIncIni,2),round(clockEspGota,2),round(clockPrepVid,2),round(clockVid,2),round(clockImages,2),round(clockIncRet,2),round(clockEnsaio,2),framerate,string(vidComp),vidQual,'VariableNames',varnames);
writetable(rows2vars(T),fullfilenametxt,'Delimiter',' ','WriteVariableNames',false) %Write table to .txt file

pause(5);
end

