function [volSerTotal] = routineContPendDrop(vidAdap,vidID,vidFmt,imFmt,smSeringa,passoMot,passoSer,wMin,framerate,limTotalFrames,volSerTotal)
%ROUTINECONTPENDDROP Execution of a continuous pendant drop experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Execution of a continuous pendant drop experiment. Performs a
%   continuous increase in the drop volume until reaching a pre-established
%   maximum volume. During the volume increase, drop images are captured
%   allowing further analysis. The routine requires the user to define the 
%   initial drop volume (volGotaIni), the maximum drop volume (volGotaMax),
%   the drop flow rate (taxaDepGota), the waiting time for drop equilibrium
%   (tEsp) and the number of captured images (nImages). For initial drop 
%   formation, the minimum flow rate is used (0.3158 uL/s). The routine 
%   allows choosing the directory where the images will be exported, as 
%   well as defining the name of the images. At the end, a text file is 
%   exported with all relevant parameters of the experiment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT
% vidAdap - Character vector that specifies the name of the adaptor used to communicate with the device
% vidID - Numeric scalar value that identifies a particle device available through the specified adaptor 
% vidFmt - Character vector that specifies a particular video format supported by the device
% imFmt - Format of the captured image 
% smSerigna - Syringe motor control
% passoMot - Number of steps required for one revolution in the motor [steps/rev]
% passoSer - Amount of volume displaced in the syringe corresponding to one revolution of its plunger [uL/rev]
% wMin - Minimum stable rotation speed (RPM) in which the motor can reach. Corresponds to a sm.RPM = 2 (Microstep)
% framerate - Actual camera frame rate
% volSerTotal - Total liquid volume in the syringe in uL
%   OUTPUT
% volSerTotal - Total liquid volume in the syringe in uL

fprintf('--------------------- PENDANT DROP - CONTINUOUS EXPERIMENT ---------------------\n');
% Experiment parameters (entered by the user)
fprintf('Experiment parameters: \n');
volGotaIni = input("- Initial drop volume (uL): "); %Initial drop volume before start the measurement series
volGotaMax = input("- Maximum drop volume (uL): "); %Maximum drop volume reached by the pendant drop
while 1
    taxaDepGota = input("- Flow rate (uL/s) (max. 1 uL/s): "); %Flow rate in uL/s
    if taxaDepGota > 1 %Check that the drop flow rate does not exceed the maximum motor rotation speed limit
        fprintf("Unable to proceed with the routine. Please lower the flow rate. \n")
    else
        break;
    end
end
tEspera = input("- Waiting time for drop equilibrium (s): "); %Waiting time for drop to reach equilibrium
nImages = input("- Number of captured images: "); %Number of captured images during the continuous drop increase

% Experiment parameters (calculated)
volGotaVar = volGotaMax - volGotaIni; %Volume variation during the experiment in uL
tEnsaio = volGotaVar/taxaDepGota; %Experiment total time (increase of drop volume) in s
%tCaptInterval = tEnsaio/nImages; %Acquisition time between images in s
%volCaptInterval = taxaDepGota*tCaptInterval; %Volume variation between images in uL
totalFrame = ceil(framerate*tEnsaio); %Total number of frames acquired during the experiment
frameInterval = round(totalFrame/nImages); %Interval between frames to acquire the desired number of images
tCaptInterval = frameInterval/framerate; %Acquisition time between images in s
volCaptInterval = taxaDepGota*tCaptInterval; %Volume variation between images in uL
erroVolTotal = 0;

if volGotaMax >= volSerTotal  %Check that there is enough volume in the syringe to perfrom the procedure
    fprintf("Unable to complete the routine. Please increase the liquid volume in the syringe. \n")
else
    %- Drive motor settings
    nPassosvolDepEns = vol2passo(volGotaVar,passoMot,passoSer); %No of steps corresponding to the released volume during the experiment

    %- Image capture settings
    fprintf('Path selection to export video and images... \n');
    filePath = uigetdir('C:\'); %Path selection to export video and images
    name = input('Enter a name for the video and image series: ','s');
    vidDiskFmt = '.avi'; %Video format ('.avi','.mp4','.mj2')
    %imFmt = '.tif'; %Images format ('.jpeg','.png','.tif')
    vidComp = 'Motion JPEG AVI'; %Video compression ('Archival','Motion JPEG AVI','Motion JPEG 2000','MPEG-4','Uncompressed AVI','Indexed AVI' e 'Grayscale AVI')
    vidQual = 90; %Video quality (related with the compression ratio)

    %- Begin of experiment
    fprintf('--------------------- EXECUTION OF PENDANT DROP ROUTINE ----------------------\n');
    fprintf('Press any key to proceed with the pendant drop routine... \n');
    pause;
    tic
    close all %Close all open windows
    %delete(vidobj) %Delete any video objects that may have been created
    clear vidobj vidsrc %Clear video related objects
    clockEnsaio1 = toc;

    %- Initital drop formation
    tic
    nPassosvolGotaIni = vol2passo(volGotaIni,passoMot,passoSer); %No of steps corresponding to the initial drop volume formed at the needle tip
    smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
    fprintf('Please wait, initial drop formation... \n');
    move(smSeringa,floor(nPassosvolGotaIni)); %Drive motor connected to the syringe
    release(smSeringa); %Motor release
    fprintf('Initial drop formation completed! \n');
    erroVolGotaIni = passo2vol(nPassosvolGotaIni - floor(nPassosvolGotaIni),passoMot,passoSer); %Error in the volume during the initial drop formation (uL)
    erroVolTotal = erroVolTotal + erroVolGotaIni;
    volSerTotal = volSerTotal - volGotaIni; %Update of total liquid volume in the syringe
    clockGotaIni = toc; %Time for the initial drop formation in s

    %- Waiting time for drop equilibrium
    tic
    fprintf('Please wait, waiting time for drop reach equilibrium.... \n');
    pause(tEspera); %Waiting time for drop equilibrium
    fprintf('Waiting time completed! \n');
    clockEspGotaIni = toc; %Waiting time for drop equilibrium
    
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
    clockPrepVid = toc; %Time spent for creating and preparing the video object

    %- Start video recording
    tic
    fprintf('Start video recording. \n');
    start(vidobj) %Start video recording
    
    %- Continuous increase of drop volume
    fprintf('Please wait, continuous increase in drop volume... \n');
    if taxaDepGota <= 0.32 %Flow rate less than minimum motor rotation speed. Motor will work in pulses
        smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
        wMotorEns = uLs2RPM(taxaDepGota,passoSer); %Motor rotation speed (RPM) to achieve the desired fow rate (uL/s)
        tadpasso = tadicionalpasso(passoMot,wMotorEns,wMin); %Time added to each step (in seconds) to reach the desired motor rotational speed
        for i = 1:floor(nPassosvolDepEns) %Continuous drop volume release
            move(smSeringa,1); %Executes 1 rotation step on the motor
            pause(tadpasso);
        end
    elseif (0.32 < taxaDepGota) && (taxaDepGota <= 1)
        wRPM = fuLs2RPM(taxaDepGota); %Calculation of the necessary rotational speed to be sent to the motor [RPM] in order to reach the desired flow rate [uL/s]
        taxaDepGota = fRPM2uLs(ceil(wRPM)); %Update flow rate
        % Update variables
        tEnsaio = volGotaVar/taxaDepGota;
        %tCaptInterval = tEnsaio/nImages;
        %volCaptInterval = taxaDepGota*tCaptInterval;
        totalFrame = ceil(framerate*tEnsaio);
        frameInterval = round(totalFrame/nImages);
        tCaptInterval = frameInterval/framerate;
        volCaptInterval = taxaDepGota*tCaptInterval;
        %vidobj.FramesPerTrigger = totalFrame;
        % Syringe plunger movement
        %smSeringa.RPM = floor(wRPM);
        smSeringa.RPM = round(wRPM);
        move(smSeringa,floor(nPassosvolDepEns)); %Motor movement
    end
    release(smSeringa); %Motor release
    fprintf('Continuous increase in drop volume completed! \n');
    erroVolDepEns = passo2vol(nPassosvolDepEns - floor(nPassosvolDepEns),passoMot,passoSer); %Error in the volume released during experiment
    erroVolTotal = erroVolTotal + erroVolDepEns; %Update error in the total drop volume
    volSerTotal = volSerTotal - volGotaVar; %Update of total liquid volume in the syringe
    
    %- Finishing video recording
    wait(vidobj) %Wait until all frames are captured
    stop(vidobj);
    %-- Wait for all frames to be written to disk
    while vidobj.FramesAcquired ~= vidobj.DiskLoggerFrameCount
        pause(.1);
    end
    fprintf('Video recording completed! \n');
    clockVid = toc; %Time for the main drop formation and criation/store video on disk
    
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
            if i < 10 %Modifying the file name prefix to arrange files in ascending volume order
                prefixo = "00";
            elseif (i >= 10) && (i < 100)
                prefixo = "0";
            elseif (i >= 100) && (i < 1000)
                prefixo = "";
            end
            %fullimName = strcat(prefixo,num2str(i),'_',fileName,'_V(uL)_',num2str((volGotaIni + k*volCapIntervalframe),'%.2f'),imFmt);
            fullimName = strcat(prefixo,num2str(i),'_',fileName,'_V(uL)_',num2str((volGotaIni + i*volCaptInterval),'%.2f'),imFmt);
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
            if k < 10 %Modifying the file name prefix to arrange files in ascending volume order
                prefixo = "00";
            elseif (k >= 10) && (k < 100)
                prefixo = "0";
            elseif (k >= 100) && (k < 1000)
                prefixo = "";
            end
            %fullimName = strcat(fileName,'_',num2str(i),'_V(uL)_',num2str((volGotaIni + k*volCaptInterval),'%.2f'),imFmt);
            fullimName = strcat(prefixo,num2str(k),'_',fileName,'_V(uL)_',num2str((volGotaIni + k*volCaptInterval),'%.2f'),imFmt);
            fullfileimPath = fullfile(filePath,fullimName);
            imwrite(framegray,fullfileimPath);
            %i = i + 1; %Interval of image name
        end
    end
    delete(fullfilevidPath) %Delete video
    fprintf('Conversion of video to images completed! \n');
    clockImages = toc; %Time for converitng video to images and store images on disk
    %- Deleting video object
    tic
    delete(vidobj)
    clear vidobj vidsrc
    clockEnsaio3 = toc;

    clockEnsaio = clockEnsaio1 + clockGotaIni + clockEspGotaIni + clockPrepVid + clockVid + clockEnsaio2 + clockImages + clockEnsaio3; %Update experiment total time

    %- Export/save experiment description
    filenametxt = "Description.txt";
    fullfilenametxt = fullfile(filePath,filenametxt);
    varnames = {'Datetime','volGotaIni[uL]','volGotaMax[uL]','taxaDepGota[uL/s]','tEspera[s]','nImages','tEnsaio[s]','totalFrames','FrameGrabInterval','tCaptInterval[s]','volCaptInterval[uL]','ErroVolumeGotaIni[uL]','ErroVolumeGotaEns[uL]','ErroVolumeGotaTotal[uL]','clockGotaIni[s]','clockEspGotaIni[s]','clockPrepVid[s]','clockVid[s]','clockImages[s]','clockEnsaio[s]','frameRate','VidCompression','VidQuality'}; %Variables name in txt file
    T = table(datetime,volGotaIni,volGotaMax,round(taxaDepGota,4),tEspera,nImages,round(tEnsaio,2),totalFrame,frameInterval,round(tCaptInterval,2),round(volCaptInterval,4),round(erroVolGotaIni,4),round(erroVolDepEns,4),round(erroVolTotal,4),round(clockGotaIni,2),round(clockEspGotaIni,2),round(clockPrepVid,2),round(clockVid,2),round(clockImages,2),round(clockEnsaio,2),framerate,string(vidComp),vidQual,'VariableNames',varnames);
    writetable(rows2vars(T),fullfilenametxt,'Delimiter',' ','WriteVariableNames',false) %Write table to .txt file

    pause(5);
end
end

