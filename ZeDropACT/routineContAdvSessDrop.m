function [volSerTotal] = routineContAdvSessDrop(vidAdap,vidID,vidFmt,imFmt,smSeringa,passoMot,passoSer,wMin,framerate,limTotalFrames,volSerTotal)
%ROUTINECONTADVSESSDROP Execution of a continuous sessile drop experiment 
% (needle in). Includes only the advancing routine.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Execution of continuous sessile drop experiment (needle in). Includes
%   only the advancing routine. Images are captured concomitantly with
%   increasing drop volume. Initially, the drop volume is increased until
%   the indicated volume is reached (volGotaAdvIni) to carry out the
%   advancing routine. THe routine requires the user to define the drop
%   initial volume (volGotaStart) the drop maixmum volume after advancing
%   (volGotaAdvMax), the advancing flow rate (taxaDepAdv), the waiting time
%   for drop equilibrium during advancing (tEsperaAdv) and the number of
%   images captured during advancing (nImagesAdv). The routine allows the
%   choosing the directory where the images will be exported, as well as
%   defining the name of the images. At the end, a text file is exported
%   with all relevant parametes of the experiment. Prior to start the
%   experiment, the needle is required to be inserted into the sessile
%   drop.
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

fprintf('------------- SESSILE DROP (NEEDLE IN) - CONTINUOUS EXPERIMENT (ONLY ADVANCING) ---------------\n');
% Experiment parameters (entered by the user)
fprintf('Experiment parameters (advancing): \n');
volGotaStart = input("- Initial drop volume (uL): "); %Initial drop volume
while 1
    volGotaAdvIni = input("- Initial drop volume before advancing (uL): "); %Initial drop volume for the execution of the advancing contact angle measurement protocol
    if volGotaAdvIni < volGotaStart
        fprintf('Invalid initial drop volume before advancing. Please enter a larger volume. \n');
    else
        break;
    end
end
taxaDepAdvIni = 0.1; %Flow rate in uL/s during the initial drop advance
while 1
    volGotaAdvMax = input("- Maximum drop volume after advancing (uL): "); %Maximum drop volume reached during the execution of the advancing contact angle measurement protocol
    if volGotaAdvMax <= volGotaAdvIni
        fprintf('Invalid maximum drop volume after advancing. Please enter a larger volume. \n');
    else
        break;
    end
end
while 1
    taxaDepAdv = input("- Flow rate during advancing (uL/s) (max. 1 uL/s): "); %Drop flow rate in uL/s during advancing
    if taxaDepAdv > 1 %Check that the flow rate does not exceed the maximum rotation speed 
        fprintf("Unable to proceed with the routine. Please lower the flow rate. \n")
    else
        break;
    end
end
tEsperaAdv = input("- Waiting time for drop equilibrium during advancing (s): "); %Waiting time for the drop to reach equilibrium during advancing protocol
nImagesAdv = input("- Number of images captured during advancing: "); %Number of images captured during advancing

% Experimental parameters (calculated)
volGotaVarAdvIni = volGotaAdvIni - volGotaStart; %Volume variation during intial drop increase before advancing
volGotaVarAdv = volGotaAdvMax - volGotaAdvIni; %Volume variation during advancing

tEnsaioAdv = volGotaVarAdv/taxaDepAdv;  %Total experiment time (drop volume increase) during advancing in s
totalFrameAdv = ceil(framerate*tEnsaioAdv); %Number of total frames acquired during advancing
frameIntervalAdv = round(totalFrameAdv/nImagesAdv); %Frame interval during advancing
tCaptIntervalAdv = frameIntervalAdv/framerate; %Time interval between captured frames during advancing in s 
volCaptIntervalAdv = taxaDepAdv*tCaptIntervalAdv; %Volume interval between captured frames during advancing in uL 

erroVolTotalAdv = 0;
clockAvanco = 0;
%volCaptInterval = taxaDepGota*tCaptInterval; %Volume change between images in uL

if volGotaAdvMax >= volSerTotal  %Check that there is enough volume in the syringe to perfrom the procedure
    fprintf("Unable to complete the routine. Please increase the liquid volume in the syringe. \n")
else
    %- Drive motor settings
    nPassosvolDepAdvIni = vol2passo(volGotaVarAdvIni,passoMot,passoSer); %No of steps corresponding to the released volume during preparation for advancing
    nPassosvolDepAdv = vol2passo(volGotaVarAdv,passoMot,passoSer); %No of steps corresponding to the released volume between each measurement during advancing
     
    %- Image capture settings
    fprintf('Path selection to export video and images... \n');
    filePath = uigetdir('C:\'); %Path selection to export video and images
    name = input('Enter a name for the video and image series: ','s');
    vidDiskFmt = '.avi'; %Video format ('.avi','.mp4','.mj2')
    %imFmt = '.tif'; %Images format ('.jpeg','.png','.tif')
    vidComp = 'Motion JPEG AVI'; %Video compression ('Archival','Motion JPEG AVI','Motion JPEG 2000','MPEG-4','Uncompressed AVI','Indexed AVI' e 'Grayscale AVI')
    vidQual = 90; %Video quality (related with the compression ratio)
    
    %1) ADVANCING ROUTINE
    %- Begin of experiment
    tic
    fprintf('--------------------- EXECUTION OF DROP ADVANCING ROUTINE ----------------------\n');
    fprintf('Press any key to proceed with the advancing routine... \n');
    pause;
    close all %Close all open windows 
    %delete(vidobj) %Delete any video objects that may have been created
    clear vidobj vidsrc %Clear video related objects
    clockAvanco1= toc;
    clockAvanco = clockAvanco + clockAvanco1; %Update total advancing time
    
    %- Initial drop increase
    tic
    smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
    fprintf('Please wait, inital drop increase... \n');
    wMotorEns = uLs2RPM(taxaDepAdvIni,passoSer); %Motor rotation speed (RPM) to achieve the desired flow rate (uL/s)
    tadpasso = tadicionalpasso(passoMot,wMotorEns,wMin); %Time added to each step (in seconds) to reach the desired motor rotational speed
    for i = 1:floor(nPassosvolDepAdvIni) %Initial drop increase
        move(smSeringa,1); %Executes 1 rotation step on the motor
        pause(tadpasso);
    end
    release(smSeringa); %Motor release
    fprintf('Initial drop increase completed! \n');
    erroVolDepAdvIni = -passo2vol(nPassosvolDepAdvIni - floor(nPassosvolDepAdvIni),passoMot,passoSer); %Error in the initial drop volume formed before advancing (uL)
    erroVolTotalAdv = erroVolTotalAdv + erroVolDepAdvIni;
    volSerTotal = volSerTotal - volGotaVarAdvIni; %Update total liquid volume in the syringe
    fprintf('Please wait, waiting time for drop reach equilibrium... \n');
    pause(tEsperaAdv); %Waiting time for drop reach equilibrium
    fprintf('Waiting time completed! \n');
    clockAdvGotaIni = toc; %Time for initial drop increase in s before advancing
    clockAvanco = clockAvanco + clockAdvGotaIni; %Update total advancing time
    
    %- Continuous increase of drop volume
    tic
    %-- Create video object
    [vidobj, vidsrc] = createVidObj(vidAdap,vidID,vidFmt);

    %- Update of image names
    %fileName = strcat(name,'_Avanco_',num2str(k),'_V(uL)_',num2str((volGotaAdvIni+k*volGotaVarAdv),'%.2f'));
    fileName = name;
    fullvidName = strcat(fileName,vidDiskFmt);
    fullfilevidPath = fullfile(filePath,fullvidName);
    diskLogger = VideoWriter(fullfilevidPath,vidComp);
    diskLogger.FrameRate = framerate; %Matching the video frame rate to the actual camera frame rate
    diskLogger.Quality = vidQual; %Definition of the video quality
    vidobj.DiskLogger = diskLogger;
    if totalFrameAdv <= limTotalFrames %Continuous frame capture during video recording
        vidobj.FramesPerTrigger = totalFrameAdv; %Total frames to acquire
        vidobj.FrameGrabInterval = 1; %Frame interval
    else %Capturing spaced frames during video recording
        vidobj.FramesPerTrigger = nImagesAdv; %Total frames to acquire
        vidobj.FrameGrabInterval = frameIntervalAdv; %Frame interval
    end
    clockVidCreateAdv = toc;
    clockAvanco = clockAvanco + clockVidCreateAdv; %Update total advancing time

    %-- Start video recording
    tic
    fprintf('Start video recording. \n');
    start(vidobj) %Start video recording
    %-- Continuous increase in drop volume
    fprintf('Please wait, continuous increase in drop volume... \n');
    if taxaDepAdv <= 0.32 %Flow rate less than minimum motor rotation speed. Motor will work in pulses
        smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
        wMotorEns = uLs2RPM(taxaDepAdv,passoSer); %Motor rotation speed (RPM) to achieve the desired fow rate (uL/s)
        tadpasso = tadicionalpasso(passoMot,wMotorEns,wMin); %Time added to each step (in seconds) to reach the desired motor rotational speed
        for i = 1:floor(nPassosvolDepAdv) %Continuous drop volume release
            move(smSeringa,1); %Executes 1 rotation step on the motor
            pause(tadpasso);
        end
    elseif (0.32 < taxaDepAdv) && (taxaDepAdv <= 1)
        wRPM = fuLs2RPM(taxaDepAdv); %Calculation of the necessary rotational speed to be sent to the motor in RPM
        taxaDepAdv = fRPM2uLs(ceil(wRPM)); %Update flow rate
        % Update variables
        tEnsaioAdv = volGotaVarAdv/taxaDepAdv;
        totalFrameAdv = ceil(framerate*tEnsaioAdv);
        frameIntervalAdv = round(totalFrameAdv/nImagesAdv);
        tCaptIntervalAdv = frameIntervalAdv/framerate;
        volCaptIntervalAdv = taxaDepAdv*tCaptIntervalAdv;
        % Syringe plunger movement
        smSeringa.RPM = round(wRPM);
        move(smSeringa,floor(nPassosvolDepAdv)); %Motor movement
    end
    release(smSeringa); %Motor release
    fprintf('Continuous increase in drop volume completed! \n');
    erroVolDepAdv = passo2vol(nPassosvolDepAdv - floor(nPassosvolDepAdv),passoMot,passoSer); %Error in the volume released during advancing (uL)
    erroVolTotalAdv = erroVolTotalAdv + erroVolDepAdv; %Update error on total drop volume
    volSerTotal = volSerTotal - volGotaVarAdv; %Update total liquid volume in the syringe
    %-- Finishing video recording
    wait(vidobj) %Wait until all frames are captured
    stop(vidobj);
    %--- Wait for all frames to be written to disk
    while vidobj.FramesAcquired ~= vidobj.DiskLoggerFrameCount
        pause(.1);
    end
    fprintf('Video recording completed! \n');
    clockAdv = toc;
    clockAvanco = clockAvanco + clockAdv; %Update total advancing time
      
    %-- Showing an example image
    tic
    close all %Close all open windows
    frame = getsnapshot(vidobj); %Capture an image
    imshow(frame); %Show image
    movegui(gcf,'center') %Centralize image window
    clockExAdv = toc;
    clockAvanco = clockAvanco + clockExAdv; %Update total advancing time
   
    %- Converting video to images
    tic
    fprintf('Please wait, converting video to images... \n');
    if totalFrameAdv <= limTotalFrames %Continuous frame capture during video recording
        fullfilevidPath = fullfile(filePath,fullvidName); %Full path of the video
        vidFile = VideoReader(fullfilevidPath); %Read Video File
        iAdv = 1;
        for k = 1: round(frameIntervalAdv): totalFrameAdv %Getting an image from every video second
            frame = read(vidFile,k);
            framegray = rgb2gray(frame);
            if iAdv < 10 %Modifying the file name prefix to arrange files in ascending volume order
                prefixo = "00";
            elseif (iAdv >= 10) && (iAdv < 100)
                prefixo = "0";
            elseif (iAdv >= 100) && (iAdv < 1000)
                prefixo = "";
            end
            %fullimName = strcat(prefixo,num2str(i),'_',fileName,'_V(uL)_',num2str((volGotaIni + k*volCapIntervalframe),'%.2f'),imFmt);
            fullimName = strcat(prefixo,num2str(iAdv),'_',fileName,'_V(uL)_',num2str((volGotaAdvIni + iAdv*volCaptIntervalAdv),'%.2f'),imFmt);
            fullfileimPath = fullfile(filePath,fullimName);
            imwrite(framegray,fullfileimPath);
            iAdv = iAdv + 1; %Interval of image name
        end
    else %Capturing spaced frames during video recording
        fullfilevidPath = fullfile(filePath,fullvidName); %Full path of the video
        vidFile = VideoReader(fullfilevidPath); %Read Video File
        for k = 1: nImagesAdv %Getting a image from every video frame
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
            fullimName = strcat(prefixo,num2str(k),'_',fileName,'_V(uL)_',num2str((volGotaAdvIni + k*volCaptIntervalAdv),'%.2f'),imFmt);
            fullfileimPath = fullfile(filePath,fullimName);
            imwrite(framegray,fullfileimPath);
            %i = i + 1; %Interval of image name
        end
    end
    delete(fullfilevidPath) %Delete video
    fprintf('Conversion of video to images completed! \n');
    clockImagesAdv = toc;
    clockAvanco = clockAvanco + clockImagesAdv; %Update total advancing time

    %- Deleting video object
    tic
    delete(vidobj)
    clear vidobj vidsrc
    clockExcVideoAdv = toc;
    clockAvanco = clockAvanco + clockExcVideoAdv; %Update total advancing time

    %- Export/save experiment description
    filenametxt = "Description.txt";
    fullfilenametxt = fullfile(filePath,filenametxt);
    varnames = {'Datetime','volGotaStart[uL]','volGotaAdvIni[uL]','taxaDepAdvIni[uL]','volGotaAdvMax[uL]','taxaDepAdv[uL/s]','tEsperaAdv[s]','tEnsaioAdv[s]','nImagesAdv','totalFrameAdv','frameIntervalAdv','tCaptIntervalAdv[s]','volCaptIntervalAdv[uL]','volGotaVarAdv[uL]','erroVolDepAdvIni[uL]','erroVolDepAdv[uL]','erroVolTotalAdv[uL]','clockAdvGotaIni[s]','clockAdv[s]','clockImagesAdv[s]','clockTotalAvanco[s]','frameRate','VidCompression','VidQuality'}; %Variables name in txt file
    T = table(datetime,volGotaStart,volGotaAdvIni,taxaDepAdvIni,volGotaAdvMax,taxaDepAdv,tEsperaAdv,tEnsaioAdv,nImagesAdv,totalFrameAdv,frameIntervalAdv,tCaptIntervalAdv,volCaptIntervalAdv,volGotaVarAdv,round(erroVolDepAdvIni,4),round(erroVolDepAdv,4),round(erroVolTotalAdv,4),round(clockAdvGotaIni,2),round(clockAdv,2),round(clockImagesAdv,2),round(clockAvanco,2),framerate,string(vidComp),vidQual,'VariableNames',varnames);
    writetable(rows2vars(T),fullfilenametxt,'Delimiter',' ','WriteVariableNames',false) %Write table to .txt file
    pause(5);
end
end

