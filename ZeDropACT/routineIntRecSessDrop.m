function [volSerTotal] = routineIntRecSessDrop(vidAdap,vidID,vidFmt,imFmt,smSeringa,smSeringaZ,passoMot,passoSer,passoFuso,wMin,framerate,volSerTotal)
%ROUTINEINTRECSESSDROP Execution of a intermittent sessile drop experiment 
% (needle out). Includes only the receding routine.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Execution of a intermittent sessile drop experiment (needle out). 
%   Includes only the receding routine. Stops during withdrawal are 
%   performed to remove the needle from the drop and capture images. 
%   Initially, the drop volume is increased until the preparation volume 
%   for receding is reached (volGotaRecPrep). Then, decrease of the drop 
%   volume, at high and low speed, are carried out before drop receding is 
%   executed. The routine requires the user to define the initial drop 
%   volume (volGotaStart), volume of drop preparation for receding 
%   (volGotaRecPrep), initial drop volume before receding (volGotaRecIni), 
%   flow rate during receding (taxaDepRec), syringe displacement 
%   (dispZSer), waiting time for drop equilibrium during receding 
%   (tEsperaRec), capture time during receding (tCaptRec), number of 
%   captured images during receding (nImagesRec) and number of stops during
%   receding (nStopsRec). The routine allows choosing the directory where 
%   the images will be saved, as well as defining the name of the images. 
%   At the end, a text file is exported with all relevant parametes of the 
%   experiment. Prior to start the experiment, the needle is required to be
%   inserted into the sessile drop. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT
% vidAdap - Character vector that specifies the name of the adaptor used to communicate with the device
% vidID - Numeric scalar value that identifies a particle device available through the specified adaptor 
% vidFmt - Character vector that specifies a particular video format supported by the device
% imFmt - Format of the captured image 
% smSerigna - Syringe motor control
% smSeringaZ - Syringe displacement motor control
% passoMot - Number of steps required for one revolution in the motor [steps/rev]
% passoSer - Amount of volume displaced in the syringe corresponding to one revolution of its plunger [uL/rev]
% passoFuso - Spindle pitch [mm/rev]
% wMin - Minimum stable rotation speed (RPM) in which the motor can reach. Corresponds to a sm.RPM = 2 (Microstep)
% framerate - Actual camera frame rate
% volSerTotal - Total liquid volume in the syringe in uL
%   OUTPUT
% volSerTotal - Total liquid volume in the syringe in uL

fprintf('----------- SESSILE DROP (NEEDLE OUT) - INTERMITTENT EXPERIMENT (ONLY RECEDING) -------------\n');
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
dispZSer = input("- Syringe displacement (mm): "); %Syringe displacement during the execution of advancing and receding contact angle measurement protocol. Displacement required for the needle to move out of the camera's field of view
tEsperaRec = input("- Waiting time for drop equilibrium during receding (s): "); %Waiting time for drop equilibrium during receding
tCaptRec = input("- Capture time during a stop in receding (s): "); %Capture time during a stop in receding
nImagesRec = input("- Number of captured images during a stop in receding: "); %Number of captured images during a stop in receding
nStopsRec = input("- Number of stops during receding: "); %Number of stops during receding

% Experimental parameters (calculated)
volGotaVarRecPrep = volGotaRecPrep - volGotaStart; %Volume variation during preparation (drop increase) before receding
volGotaVarRecIni = abs(volGotaRecIni - volGotaRecPrep); %Volume variation during the initial step of receding (drop decrease)
volGotaVarRecIniAlta = floor(volGotaVarRecIni*0.8); %Volume variation during the initial step of receding (drop decrease) at high flow rate
volGotaVarRecIniBaixa= volGotaVarRecIni - volGotaVarRecIniAlta; %Volume variation during the initial step of receding (drop decrease) at low flow rate
volGotaVarRec = abs(volGotaRecIni - volGotaRecEnd)/nStopsRec; %Volume variation in each stop during receding
tCaptIntervalRec = tCaptRec/nImagesRec; %Acquired time between images during receding in s
totalFrameRec = ceil(framerate*tCaptRec); %Number of total frames acquired during receding
frameIntervalRec = round(totalFrameRec/nImagesRec); %Frame interval during receding
erroVolTotalRec = 0;
clockRecuo = 0;

if volGotaRecPrep >= volSerTotal  %Check that there is enough volume in the syringe to perfrom the procedure
    fprintf("Unable to proceed with the routine. Please increase the liquid volume in the syringe. \n")
else
    %- Drive motor settings
    nPassosvolDepRecPrep = vol2passo(volGotaVarRecPrep,passoMot,passoSer); %No of steps corresponding to the volume released during the preparation for receding
    nPassosvolAdqRecIniAlta = vol2passo(volGotaVarRecIniAlta,passoMot,passoSer); %No of steps corresponding to the volume withdrawal before receding in high flow rate
    nPassosvolAdqRecIniBaixa = vol2passo(volGotaVarRecIniBaixa,passoMot,passoSer); %No of steps corresponding to the volume withdrawal before receding in low flow rate
    nPassosvolAdqRec = vol2passo(volGotaVarRec,passoMot,passoSer); %No of steps corresponding to the volume withdrawal in each stop during receding
    
    %- Image capture settings
    fprintf('Path selection to export video and images... \n');
    filePath = uigetdir('C:\'); %Path selection to export video and images
    name = input('Enter a name for the video and image series: ','s');
    vidDiskFmt = '.avi'; %Video format ('.avi','.mp4','.mj2')
    %imFmt = '.tif'; %Images format ('.jpeg','.png','.tif')
    vidComp = 'Motion JPEG AVI'; %Video compression ('Archival','Motion JPEG AVI','Motion JPEG 2000','MPEG-4','Uncompressed AVI','Indexed AVI' e 'Grayscale AVI')
    vidQual = 90; %Video quality (related with the compression ratio)
    
    %2) RECEDING ROUTINE
    %- Begin of experiment
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
    move(smSeringa,floor(nPassosvolDepRecPrep)); %Drive motor connected to the syringe
    release(smSeringa); %Motor release
    fprintf('Increase of drop volume completed! \n');
    erroVolDepRecPrep = -passo2vol(nPassosvolDepRecPrep - floor(nPassosvolDepRecPrep),passoMot,passoSer); %Error in the drop preparation volume before receding (uL)
    erroVolTotalRec = erroVolTotalRec + erroVolDepRecPrep;
    volSerTotal = volSerTotal - volGotaVarRecPrep; %Update total liquid volume in the syringe
    pause(1); %Wait to start liquid withdrawal from the drop 
    clockRecGotaPrep = toc; %Time of preparaing drop volume before receding in s
    clockRecuo = clockRecuo + clockRecGotaPrep; %Update total receding time

    %- Initial withdrawal of drop volume
    tic
    smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
    fprintf('Please wait, initial withdrawal of drop volume at high flow rate...\n');
    move(smSeringa,-floor(nPassosvolAdqRecIniAlta)); %Drive motor connected to the syringe
    release(smSeringa); %Motor release
    pause(1);
    fprintf('Please wait, initial withdrawal of drop volume at low flow rate...\n');
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
    fprintf('Please wait, waiting time for drop reach equilibrium...\n');
    pause(tEsperaRec); %Waiting time for drop reach equilibrium
    fprintf('Waiting time completed! \n');
    clockRecGotaIni = toc; %Time for initial drop receding in s before receding
    clockRecuo = clockRecuo + clockRecGotaIni; %Update total receding time

    %- Receding measurement series
    iRec = 1; %Counter for all captured images during receding
    for k = 1:nStopsRec
        tic
        %- Create video object
        [vidobj, vidsrc] = createVidObj(vidAdap,vidID,vidFmt);

        %- Update of image names
        fileName = strcat(name,'_Recuo_',num2str(k),'_V(uL)_',num2str((volGotaRecIni-k*volGotaVarRec),'%.2f'));
        fullvidName = strcat(fileName,vidDiskFmt);
        fullfilevidPath = fullfile(filePath,fullvidName);
        diskLogger = VideoWriter(fullfilevidPath,vidComp);
        diskLogger.FrameRate = framerate; %Matching the video frame rate to the actual camera frame rate
        diskLogger.Quality = vidQual; %Definition of the video quality
        vidobj.DiskLogger = diskLogger;
        %vidobj.FramesPerTrigger = framerate*tCapt; %Total frames to acquire
        vidobj.FramesPerTrigger = nImagesRec; %Total frames to acquire
        vidobj.FrameGrabInterval = frameIntervalRec; %Frame interval

        %- Execution of measurement
        fprintf("------------------------------------------------------------------- \n");
        fprintf("Execution of receding measurement %d of %d \n",k,nStopsRec);
        clockVidCreateRec = toc;
        clockRecuo = clockRecuo + clockVidCreateRec; %Update total receding time

        %- Additional decrease of drop volume
        tic
        fprintf('Please wait, additional decrease of drop volume... \n');
        if taxaDepRec <= 0.32 %Flow rate less than minimum motor rotation speed. Motor will work in pulses
            smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
            wMotorEns = uLs2RPM(taxaDepRec,passoSer); %Motor rotation speed (RPM) to achieve the desired fow rate (uL/s)
            tadpasso = tadicionalpasso(passoMot,wMotorEns,wMin); %Time added to each step (in seconds) to reach the desired motor rotational speed
                for i = 1:floor(nPassosvolAdqRec) %Additional decrease of drop volume
                    move(smSeringa,-1); %Executes 1 rotation step on the motor
                    pause(tadpasso);
                end
        elseif (0.32 < taxaDepRec) && (taxaDepRec <= 1)
            wRPM = fuLs2RPM(taxaDepRec); %Calculation of the necessary rotational speed to be sent to the motor in RPM
            smSeringa.RPM = floor(wRPM);
            move(smSeringa,-floor(nPassosvolAdqRec)); %Motor movement
        end
        release(smSeringa); %Motor release
        fprintf('Additional decrease of drop volume completed! \n');
        erroVolAdqRec = passo2vol(nPassosvolAdqRec - floor(nPassosvolAdqRec),passoMot,passoSer); %Error in the volume acquired during receding (uL)
        erroVolTotalRec = erroVolTotalRec + erroVolAdqRec;
        volSerTotal = volSerTotal + volGotaVarRec; %Update total liquid volume in the syringe
        pause(1);
        clockRec(k) = toc;
        clockRecuo = clockRecuo + clockRec(k); %Update total receding time

        %- Syringe displacement (syringe departure from the drop)
        tic
        fprintf('Please wait, syringe displacement... \n');
        moverSeringa(passoMot,passoFuso,smSeringaZ,"up",dispZSer,"high"); %Syringe movement
        release(smSeringaZ); %Motor release
        fprintf('Syringe displacement completed! \n');
        clockDispSerRec1(k) = toc;
        clockRecuo = clockRecuo + clockDispSerRec1(k); %Update total receding time

        %- Waiting time for drop reach equilibrium
        tic
        fprintf('Please wait, waiting time for drop reach equilibrium...\n');
        if clockDispSerRec1(k) >= tEsperaRec %Syringe displacement time is longer than waiting time for drop equilibrium
            tEsperaRecReal(k) = clockDispSerRec1(k); %Calculation of actual waiting time
        else %Syringe displacement time is less than waiting time for drop equilibrium
            tEsperaEspReal(k) = tEsperaRec;
            deltat = tEsperaRec - clockDispSerRec1(k);
            pause(deltat)
        end
        %pause(tEsperaRec);
        fprintf('Waiting time completed! \n');
        clockEspRec = toc;
        clockRecuo = clockRecuo + clockEspRec;

        %- Start video recording
        tic
        fprintf('Start video recording. \n');
        start(vidobj) %Start video recording
        wait(vidobj) %Wait until all frames are captured
        stop(vidobj);
        % Wait for all frames to be written to disk
        while vidobj.FramesAcquired ~= vidobj.DiskLoggerFrameCount
            pause(.1);
        end
        fprintf('Video recording completed! \n');
        clockVidRec(k) = toc;
        clockRecuo = clockRecuo + clockVidRec(k); %Update total receding time
        
        %- Showing an example image
        tic
        close all %Close all open windows
        frame = getsnapshot(vidobj); %Capture an image
        imshow(frame); %Show image
        movegui(gcf,'center') %Centralize image window
        set(gcf,'Name',strcat(fileName,'_Example'))
        clockExRec = toc;
        clockRecuo = clockRecuo + clockExRec; %Update total receding time

        %- Converting video to images
        tic
        fprintf('Please wait, converting video to images... \n');
        %vid2frame(filePath,fullvidName,fileName,imFmt,tCaptInterval);
        fullfilevidPath = fullfile(filePath,fullvidName); %Full path of the video
        vidFile = VideoReader(fullfilevidPath); %Read Video File
        %i = 1;
        for m = 1: nImagesRec %Getting a image from every video frame
            frame = read(vidFile,m);
            framegray = rgb2gray(frame);
            if iRec < 10 %Modifying the file name prefix to arrange files in ascending volume order
                prefixo = "00";
            elseif (iRec >= 10) && (iRec < 100)
                prefixo = "0";
            elseif (iRec >= 100) && (iRec < 1000)
                prefixo = "";
            end
            fullimName = strcat(prefixo,num2str(iRec),'_',fileName,imFmt);
            fullfileimPath = fullfile(filePath,fullimName);
            imwrite(framegray,fullfileimPath);
            iRec = iRec + 1; %Interval of image name
        end
        delete(fullfilevidPath) %Delete video
        fprintf('Conversion of video to images completed! \n');
        clockImagesRec(k) = toc;
        clockRecuo = clockRecuo + clockImagesRec(k); %Update total receding time
     
        %- Syringe displacement (syringe return to the drop)
        tic
        fprintf('Please wait, syringe displacement... \n');
        moverSeringa(passoMot,passoFuso,smSeringaZ,"down",dispZSer,"high"); %Syringe movement
        release(smSeringaZ); %Motor release
        fprintf('Syringe displacement completed! \n');
        clockDispSerRec2(k) = toc;
        clockRecuo = clockRecuo + clockDispSerRec2(k); %Update total receding time
        
        %- Syringe position monitoring (showing an example image)
        tic
        close all %Close all open windows
        frame = getsnapshot(vidobj); %Capture an image
        imshow(frame); %Show image
        movegui(gcf,'center') %Centralize image window
        set(gcf,'Name',strcat(fileName,'_CheckNeedlePosition'))
        %%{
        hold on;
        vidResValues = vidobj.VideoResolution; %Get resolution of the video object
        stepSizeY = vidResValues(2)/8; % Distance between lines in Y direction
        stepSizeX = vidResValues(1)/4; %Distance between lines in X direction
        for row = 1 : stepSizeY : vidResValues(2)
        yline(row, '-.r', 'LineWidth', 1);
        end
        for col = 1 : stepSizeX : vidResValues(1)
        xline(col, '-.r', 'LineWidth', 1);
        end
        %}
        clockCheckSerRec = toc;
        clockRecuo = clockRecuo + clockCheckSerRec; %Update total receding time

        %- Deleting video object
        tic
        delete(vidobj)
        clear vidobj vidsrc
        clockExcVidRec = toc;
        clockRecuo = clockRecuo + clockExcVidRec; %Update total receding time
    end

    %- Export/save experiment description
    filenametxt = "Description.txt";
    fullfilenametxt = fullfile(filePath,filenametxt);
    varnames = {'Datetime','volGotaStart[uL]','volGotaRecPrep[uL]','volGotaRecIni[uL]','volGotaRecEnd[uL]','taxaDepRecIni[uL/s]','taxaDepRec[uL/s]','tEsperaRec[s]','tEsperaRecReal[s]','tCaptRec[s]','nImagesRec','nStopsRec','deslZSer[mm]','volGotaVarRec[uL]','totalFrameRec','frameIntervalRec','erroVolDepRecPrep[uL]','erroVolAdqRecIni[uL]','erroVolAdqRec[uL]','erroVolTotalRec[uL]','clockRecGotaPrep[s]','clockRecGotaIni[s]','clockRec[s]','clockDispSerRec[s]','clockVidRec[s]','clockImagesRec[s]','clockTotalRecuo[s]','frameRate','VidCompression','VidQuality'}; %Variables name in txt file
    T = table(datetime,volGotaStart,volGotaRecPrep,volGotaRecIni,volGotaRecEnd,taxaDepRecIni,taxaDepRec,tEsperaRec,round(mean(tEsperaRecReal),2),tCaptRec,nImagesRec,nStopsRec,dispZSer,volGotaVarRec,totalFrameRec,frameIntervalRec,round(erroVolDepRecPrep,4),round(erroVolAdqRecIni,4),round(erroVolAdqRec,4),round(erroVolTotalRec,4),round(clockRecGotaPrep,2),round(clockRecGotaIni,2),round(mean(clockRec),2),round(mean(mean(clockDispSerRec1)+mean(clockDispSerRec2)),2),round(mean(clockVidRec),2),round(mean(clockImagesRec),2),round(clockRecuo,2),framerate,string(vidComp),vidQual,'VariableNames',varnames);
    writetable(rows2vars(T),fullfilenametxt,'Delimiter',' ','WriteVariableNames',false) %Write table to .txt file
    pause(5);
end
end

