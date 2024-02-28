function [volSerTotal] = routineContRecSessDrop(vidAdap,vidID,vidFmt,imFmt,smSeringa,passoMot,passoSer,wMin,framerate,limTotalFrames,volSerTotal)
%ROUTINECONTRECSESSDROP Execution of a continuous sessile drop experiment 
% (needle in). Includes only the receding routine.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Execution of continuous sessile drop experiment (needle in). Includes
%   only the receding routine. Images are captured concomitantly with
%   decreasing drop volume. Initially, the drop volume is increased until
%   the preparation volume for receding is reached (volGotaRecPrep). Then,
%   reduction in the drop volume, at high and low flow rates, are performed
%   until the volume indicated for receding is reached (volGotaRecIni). The
%   routine requires the user to define the initial drop volume
%   (volGotaStart), the preparation drop volume for receding
%   (volGotaRecPrep), the initial drop volume before receding
%   (volGotaRecIni), the flow rate during receding (taxaDepRec), the
%   waiting time for drop equilibrium during receding (tEsperaRec) and the
%   number of images capruted during receding (nImagesRec). The routine
%   allows choocing the directory where the images will be exported, as 
%   well as defining the name of the images. At the end, a text file is
%   exported with all relevant parameters of the experiment. Prior to start
%   the experiment, the needle is required to be inserted into the sessile
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

fprintf('-------------- SESSILE DROP (NEEDLE IN) - CONTINUOUS EXPERIMENT (ONLY RECEDING) ----------------\n');
% Experiment parameters (entered by the user)
fprintf('Experiment parameters (receding): \n');
volGotaStart = input("- Initial drop volume (uL): "); %Initial drop volume
while 1
    volGotaRecPrep = input("- Preparation drop volume for receding (uL): "); %Preparation drop volume for receding
    if volGotaRecPrep < volGotaStart
        fprintf('Invalid preparation drop volume for receding. Please enter a larger volume. \n');
    else
        break;
    end
end
while 1
    volGotaRecIni = input("- Initial drop volume before receding (uL): "); %Initial drop volume before receding
    if volGotaRecIni > volGotaRecPrep
        fprintf('Invalid initial drop volume before receding. Please enter a smaller volume. \n');
    elseif volGotaRecIni < volGotaStart
        fprintf('Invalid initial drop volume before receding. Please enter a larger volume. \n');
    else
        break;
    end
end
taxaDepRecIni = 0.1; %Initial flow rate before receding in uL/s
while 1
    volGotaRecEnd = input("- Final drop volume after receding (uL): "); %Final drop volume after receding protocol
    if volGotaRecEnd > volGotaRecIni
        fprintf('Invalid final drop volume after receding. Please enter a smaller volume. \n');
    else
        break;
    end
end
while 1
    taxaDepRec = input("- Flow rate during receding (uL/s) (max. 1 uL/s): "); %Flow rate during receding in uL/s
    if taxaDepRec > 1 %Check that the flow rate does not exceed the maximum motor rotational speed limit
        fprintf("Unable to proceed with the routine. Please lower the flow rate. \n")
    else
        break;
    end
end
tEsperaRec = input("- Waiting time for drop equilibrium during receding (s): "); %Waiting time for drop equilibrium during receding
nImagesRec = input("- Number of captured images during receding: "); %Number of captured images during receding

% Experiment parameters (calculated)
volGotaVarRecPrep = volGotaRecPrep - volGotaStart; %Volume variation during preparation (drop volume increase) before receding
volGotaVarRecIni = abs(volGotaRecPrep - volGotaRecIni); %Volume variation during the initial stage of receding (drop volume decrease)
volGotaVarRecIniAlta = floor(volGotaVarRecIni*0.8); %Volume variation during the initial receding stage (drop reduction) at high flow rate
volGotaVarRecIniBaixa = volGotaVarRecIni - volGotaVarRecIniAlta; %Volume variation during the initial receding stage (drop reduction) at low flow rate
volGotaVarRec = abs(volGotaRecIni - volGotaRecEnd); %Volume variation during receding

tEnsaioRec = volGotaVarRec/taxaDepRec; %Total experiment time (drop volume decrease) during receding in s
totalFrameRec = ceil(framerate*tEnsaioRec); %Total number of frames acquired during receding
frameIntervalRec = round(totalFrameRec/nImagesRec); %Frame interval during receding
tCaptIntervalRec = frameIntervalRec/framerate; %Time interval between captured frames during receding in s
volCaptIntervalRec = taxaDepRec*tCaptIntervalRec; %Volume variation between captured frames during receding in uL

erroVolTotalRec = 0;
clockRecuo = 0;
%volCaptInterval = taxaDepGota*tCaptInterval; %Volume change between images in uL

if volGotaRecPrep >= volSerTotal  %Check that there is enough volume in the syringe to perform the procedure
    fprintf("Unable to complete the routine. Please increase the liquid volume in the syringe. \n")
else
    %- Drive motor settings
    nPassosvolDepRecPrep = vol2passo(volGotaVarRecPrep,passoMot,passoSer); %No of steps corresponding to the volume released during preparation for receding
    nPassosvolAdqRecIniAlta = vol2passo(volGotaVarRecIniAlta,passoMot,passoSer); %No of steps corresponding to the volume withdrawal before the start of receding in high flow rate
    nPassosvolAdqRecIniBaixa = vol2passo(volGotaVarRecIniBaixa,passoMot,passoSer); %No of steps corresponding to the volume withdrawal before the start of receding in low flow rate
    nPassosvolAdqRec = vol2passo(volGotaVarRec,passoMot,passoSer); %No of steps corresponding to the volume withdrawal during receding
    
    %- Image capture settings
    fprintf('Path selection to export video and images... \n');
    filePath = uigetdir('C:\'); %Path selection to export video and images
    name = input('Enter a name for the video and image series: ','s');
    vidDiskFmt = '.avi'; %Video format ('.avi','.mp4','.mj2')
    %imFmt = '.tif'; %Images format ('.jpeg','.png','.tif')
    vidComp = 'Motion JPEG AVI'; %Video compression ('Archival','Motion JPEG AVI','Motion JPEG 2000','MPEG-4','Uncompressed AVI','Indexed AVI' e 'Grayscale AVI')
    vidQual = 90; %Video quality (related with the compression ratio)
    
    %2) RECEDING ROUTINE
    %- Begin experiment
    tic
    fprintf('--------------------- EXECUTION OF DROP RECEDING ROUTINE ----------------------\n');
    fprintf('Press any key to proceed with the receding routine... \n');
    pause;
    close all %Close all open windows
    %delete(vidobj) %Delete any video objects that may have been created
    clear vidobj vidsrc %Clear video related objects
    clockRecuo1= toc;
    clockRecuo = clockRecuo + clockRecuo1; %Update total receding time

    %- Increase of drop volume
    tic
    smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
    fprintf('Please wait, increasing drop volume... \n');
    move(smSeringa,floor(nPassosvolDepRecPrep)); %Drive motor connected to syringe
    release(smSeringa); %Motor release
    fprintf('Increase of drop volume completed! \n');
    erroVolDepRecPrep = -passo2vol(nPassosvolDepRecPrep - floor(nPassosvolDepRecPrep),passoMot,passoSer); %Error in the preparation drop volume before receding (uL)
    erroVolTotalRec = erroVolTotalRec + erroVolDepRecPrep;
    volSerTotal = volSerTotal - volGotaVarRecPrep; %Update total liquid volume in the syringe
    pause(2); %Wait to start liquid withdrawal from the drop 
    clockRecGotaPrep = toc; %Time of preparaing drop volume before receding in s
    clockRecuo = clockRecuo + clockRecGotaPrep; %Update total receding time

    %- Initial withdrawal of drop volume
    tic
    smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
    fprintf('Please wait, initial withdrawal of drop volume at high flow rate... \n');
    move(smSeringa,-floor(nPassosvolAdqRecIniAlta)); %Drive motor connected to the syringe
    release(smSeringa); %Motor release
    pause(1);
    fprintf('Please wait, initial withdrawal of drop volume at low flow rate... \n');
    smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
    wMotorEns = uLs2RPM(taxaDepRecIni,passoSer); %Motor rotation speed (RPM) to achieve the desired flow rate (uL/s)
    tadpasso = tadicionalpasso(passoMot,wMotorEns,wMin); %Time added to each step (in seconds) to reach the desired motor rotational speed
    for i = 1:floor(nPassosvolAdqRecIniBaixa) %Initial withdrawal of drop volume
        move(smSeringa,-1); %Executes 1 rotation step on the motor
        pause(tadpasso);
    end
    release(smSeringa); %Motor release
    fprintf('Initial withdrawal of drop volume completed! \n');
    erroVolAdqRecIni = passo2vol((nPassosvolAdqRecIniAlta + nPassosvolAdqRecIniBaixa) - (floor(nPassosvolAdqRecIniAlta) + floor(nPassosvolAdqRecIniBaixa)),passoMot,passoSer); %Error in the initial drop volume before receding (uL)
    erroVolTotalRec = erroVolTotalRec + erroVolAdqRecIni;
    volSerTotal = volSerTotal + volGotaVarRecIni; %Update total liquid volume in the syringe
    fprintf('Please wait, waiting time for drop reach equilibrium... \n');
    pause(tEsperaRec); %Waiting time for drop reach equilibrium
    fprintf('Waiting time completed! \n');
    clockRecGotaIni = toc; %Time for initial drop receding in s before receding 
    clockRecuo = clockRecuo + clockRecGotaIni; %Update total receding time
    
    %- Continuous decrease of drop volume
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
    if totalFrameRec <= limTotalFrames %Continuous frame capture during video recording
        vidobj.FramesPerTrigger = totalFrameRec; %Total frames to acquire
        vidobj.FrameGrabInterval = 1; %Frame interval
    else %Capturing spaced frames during video recording
        vidobj.FramesPerTrigger = nImagesRec; %Total frames to acquire
        vidobj.FrameGrabInterval = frameIntervalRec; %Frame interval
    end
    clockVidCreateRec = toc;
    clockRecuo = clockRecuo + clockVidCreateRec; %Update total receding time

    %-- Start video recording
    tic
    fprintf('Start video recording. \n');
    start(vidobj) %Start video recording
    %-- Continuous decrease in drop volume
    fprintf('Please wait, continuous decrease in drop volume... \n');
    if taxaDepRec <= 0.32 %Flow rate less than minimum motor rotation speed. Motor will work in pulses
        smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
        wMotorEns = uLs2RPM(taxaDepRec,passoSer); %Motor rotation speed (RPM) to achieve the desired fow rate (uL/s)
        tadpasso = tadicionalpasso(passoMot,wMotorEns,wMin); %Time added to each step (in seconds) to reach the desired motor rotational speed
        for i = 1:floor(nPassosvolAdqRec) %Continuous drop volume withdrawal 
            move(smSeringa,-1); %Executes 1 rotation step on the motor
            pause(tadpasso);
        end
    elseif (0.32 < taxaDepRec) && (taxaDepRec <= 1)
        wRPM = fuLs2RPM(taxaDepRec); %Calculation of the necessary rotational speed to be sent to the motor in RPM
        taxaDepRec = fRPM2uLs(ceil(wRPM)); %Update flow rate
        % Update variables
        tEnsaioRec = volGotaVarRec/taxaDepRec;
        totalFrameRec = ceil(framerate*tEnsaioRec);
        frameIntervalRec = round(totalFrameRec/nImagesRec);
        tCaptIntervalRec = frameIntervalRec/framerate;
        volCaptIntervalRec = taxaDepRec*tCaptIntervalRec;
        % Syringe plunger movement
        smSeringa.RPM = round(wRPM);
        move(smSeringa,-floor(nPassosvolAdqRec)); %Motor movement
    end
    release(smSeringa); %Motor release
    fprintf('Continuous decrease in drop volume completed! \n');
    erroVolAdqRec = passo2vol(nPassosvolAdqRec - floor(nPassosvolAdqRec),passoMot,passoSer); %Error in the volume withdrawal during the experiment
    erroVolTotalRec = erroVolTotalRec + erroVolAdqRec; %Update error on total drop volume
    volSerTotal = volSerTotal + volGotaVarRec; %Update total liquid volume in the syringe
    %-- Finishing video recording
    wait(vidobj) %Wait until all frames are captured
    stop(vidobj);
    %--- Wait for all frames to be written to disk
    while vidobj.FramesAcquired ~= vidobj.DiskLoggerFrameCount
        pause(.1);
    end
    fprintf('Video recording completed! \n');
    clockRec = toc;
    clockRecuo = clockRecuo + clockRec; %Update total receding time
      
    %-- Showing an example image
    tic
    close all %Close all open windows
    frame = getsnapshot(vidobj); %Capture an image
    imshow(frame); %Show image
    movegui(gcf,'center') %Centralize image window
    clockExRec = toc;
    clockRecuo = clockRecuo + clockExRec; %Update total receding time
   
    %- Converting video to images
    tic
    fprintf('Please wait, converting video to images... \n');
    if totalFrameRec <= limTotalFrames %Continuous frame capture during video recording
        fullfilevidPath = fullfile(filePath,fullvidName); %Full path of the video
        vidFile = VideoReader(fullfilevidPath); %Read Video File
        iRec = 1;
        for k = 1: round(frameIntervalRec): totalFrameRec %Getting a image from every video second
            frame = read(vidFile,k);
            framegray = rgb2gray(frame);
            if iRec < 10 %Modifying the file name prefix to arrange files in ascending volume order
                prefixo = "00";
            elseif (iRec >= 10) && (iRec < 100)
                prefixo = "0";
            elseif (iRec >= 100) && (iRec < 1000)
                prefixo = "";
            end
            %fullimName = strcat(prefixo,num2str(i),'_',fileName,'_V(uL)_',num2str((volGotaIni + k*volCapIntervalframe),'%.2f'),imFmt);
            fullimName = strcat(prefixo,num2str(iRec),'_',fileName,'_V(uL)_',num2str((volGotaRecIni - iRec*volCaptIntervalRec),'%.2f'),imFmt);
            fullfileimPath = fullfile(filePath,fullimName);
            imwrite(framegray,fullfileimPath);
            iRec = iRec + 1; %Interval of image name
        end
    else %Capturing spaced frames during video recording
        fullfilevidPath = fullfile(filePath,fullvidName); %Full path of the video
        vidFile = VideoReader(fullfilevidPath); %Read Video File
        for k = 1: nImagesRec %Getting a image from every video frame
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
            fullimName = strcat(prefixo,num2str(k),'_',fileName,'_V(uL)_',num2str((volGotaRecIni - k*volCaptIntervalRec),'%.2f'),imFmt);
            fullfileimPath = fullfile(filePath,fullimName);
            imwrite(framegray,fullfileimPath);
            %i = i + 1; %Interval of image name
        end
    end
    delete(fullfilevidPath) %Delete video
    fprintf('Conversion of video to images completed! \n');
    clockImagesRec = toc;
    clockRecuo = clockRecuo + clockImagesRec; %Update total receding time

    %- Deleting video object
    tic
    delete(vidobj)
    clear vidobj vidsrc
    clockExcVideoRec = toc;
    clockRecuo = clockRecuo + clockExcVideoRec; %Update total receding time

    %- Export/save experiment description
    filenametxt = "Description.txt";
    fullfilenametxt = fullfile(filePath,filenametxt);
    varnames = {'Datetime','volGotaStart[uL]','volGotaRecPrep[uL]','volGotaRecIni[uL]','volGotaRecEnd[uL]','taxaDepRecIni[uL/s]','taxaDepRec[uL/s]','tEsperaRec[s]','tEnsaioRec[s]','nImagesRec','totalFrameRec','frameIntervalRec','tCaptIntervalRec[s]','volCaptIntervalRec[uL]','volGotaVarRec[uL]','erroVolDepRecPrep[uL]','erroVolAdqRecIni[uL]','erroVolAdqRec[uL]','erroVolTotalRec[uL]','clockRecGotaPrep[s]','clockRecGotaIni[s]','clockRec[s]','clockImagesRec[s]','clockTotalRecuo[s]','frameRate','VidCompression','VidQuality'}; %Variables name in txt file
    T = table(datetime,volGotaStart,volGotaRecPrep,volGotaRecIni,volGotaRecEnd,taxaDepRecIni,taxaDepRec,tEsperaRec,tEnsaioRec,nImagesRec,totalFrameRec,frameIntervalRec,tCaptIntervalRec,volCaptIntervalRec,volGotaVarRec,round(erroVolDepRecPrep,4),round(erroVolAdqRecIni,4),round(erroVolAdqRec,4),round(erroVolTotalRec,4),round(clockRecGotaPrep,2),round(clockRecGotaIni,2),round(clockRec,2),round(clockImagesRec,2),round(clockRecuo,2),framerate,string(vidComp),vidQual,'VariableNames',varnames);
    writetable(rows2vars(T),fullfilenametxt,'Delimiter',' ','WriteVariableNames',false) %Write table to .txt file
    pause(5);
end
end

