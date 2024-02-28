function [volSerTotal] = routineIntAdvSessDrop(vidAdap,vidID,vidFmt,imFmt,smSeringa,smSeringaZ,passoMot,passoSer,passoFuso,wMin,framerate,volSerTotal)
%ROUTINEINTADVSESSDROP Execution of a intermittent sessile drop experiment 
% (needle out). Includes only the advancing routine.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Execution of a intermittent sessile drop experiment (needle out).
%   Includes only the advancing routine. Stops during advancing are made to
%   remove the needle from the drop and capture images. Initially, the drop
%   volume is increased until the indicated volume is reached 
%   (volGotaAdvIni) to carry out the advancing routine. The routine 
%   requires the user to define the initial drop volume (volGotaStart), 
%   initial drop volume before advancing (volGotaAdvIni), maximum drop 
%   volume after advancing (volGotaAdvMax), the flow rate during advancing 
%   (taxaDepAdv) , syringe displacement (dispZSer), waiting time for drop 
%   equilibrium during advancing (tEsperaAdv), capture time during stop in 
%   advancing (tCaptAdv), number of images captured during advance 
%   (nImagesAdv) and number of stops during advancing (nStopsAdv). The 
%   routine allows choosing the directory where the images will be 
%   exported, as well as defining the name of the images. At the end, a 
%   text file is exported with all relevant parametes of the experiment. 
%   Prior to start the experiment, the needle is required to be inserted 
%   into the sessile drop. 
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

fprintf('------------ SESSILE DROP (NEEDLE OUT) - INTERMITTENT EXPERIMENT (ONLY ADVANCING) --------------\n');
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
dispZSer = input("- Syringe displacement (mm): "); %Syringe displacement during the execution of advancing and receding contact angle measurement protocol. Displacement required for the needle to move out of the camera's field of view
tEsperaAdv = input("- Waiting time for drop equilibrium during advancing (s): "); %Waiting time for drop equilibrium during advancing
tCaptAdv = input("- Capture time during a stop in advancing (s): "); %Capture time during a stop in advancing
nImagesAdv = input("- Number of captured images during a stop in advancing: "); %Number of captured images during a stop in advancing
nStopsAdv = input("- Number of stops during advancing: "); %Number of stops during advancing

% Experimental parameters (calculated)
volGotaVarAdvIni = volGotaAdvIni - volGotaStart; %Volume variation during intial drop increase before advancing
volGotaVarAdv = (volGotaAdvMax - volGotaAdvIni)/nStopsAdv; %Volume variation in each stop during advancing
tCaptIntervalAdv = tCaptAdv/nImagesAdv; %Acquisition time between images during advancing in s
totalFrameAdv = ceil(framerate*tCaptAdv); %Number of total frames acquired during advancing
frameIntervalAdv = round(totalFrameAdv/nImagesAdv); %Frame interval during advancing
erroVolTotalAdv = 0;
clockAvanco = 0;
clockRecuo = 0;
clockEnsaio = 0;

if volGotaAdvMax >= volSerTotal  %Check that there is enough volume in the syringe to perfrom the procedure
    fprintf("Unable to proceed with the routine. Please increase the liquid volume in the syringe. \n")
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
    fprintf('Please wait, inital drop increase...\n');
    wMotorEns = uLs2RPM(taxaDepAdvIni,passoSer); %Motor rotation speed (RPM) to achieve the desired flow rate (uL/s)
    tadpasso = tadicionalpasso(passoMot,wMotorEns,wMin); %Time added to each step (in seconds) to reach the desired motor rotational speed
    for i = 1:floor(nPassosvolDepAdvIni) %Initial drop increase
        move(smSeringa,1); %Executes 1 rotation step on the motor
        pause(tadpasso);
    end
    release(smSeringa);  %Motor release
    fprintf('Initial drop increase completed! \n');
    erroVolDepAdvIni = -passo2vol(nPassosvolDepAdvIni - floor(nPassosvolDepAdvIni),passoMot,passoSer); %Error in the initial drop volume formed before advancing (uL)
    erroVolTotalAdv = erroVolTotalAdv + erroVolDepAdvIni;
    volSerTotal = volSerTotal - volGotaVarAdvIni; %Update total liquid volume in the syringe
    fprintf('Please wait, waiting time for drop reach equilibrium...\n');
    pause(tEsperaAdv); %Waiting time for drop equilibrium
    fprintf('Waiting time completed! \n');
    clockAdvGotaIni = toc; %Time for initial drop increase in s before advancing
    clockAvanco = clockAvanco + clockAdvGotaIni; %Update total advancing time

    %- Advancing measurement series
    iAdv = 1; %Counter for all captured images during advancing
    for k = 1:nStopsAdv
        tic
        %- Create video object
        [vidobj, vidsrc] = createVidObj(vidAdap,vidID,vidFmt);

        %- Update of image names
        fileName = strcat(name,'_Avanco_',num2str(k),'_V(uL)_',num2str((volGotaAdvIni+k*volGotaVarAdv),'%.2f'));
        fullvidName = strcat(fileName,vidDiskFmt);
        fullfilevidPath = fullfile(filePath,fullvidName);
        diskLogger = VideoWriter(fullfilevidPath,vidComp);
        diskLogger.FrameRate = framerate; %Matching the video frame rate to the actual camera frame rate
        diskLogger.Quality = vidQual; %Definition of the video quality
        vidobj.DiskLogger = diskLogger;
        %vidobj.FramesPerTrigger = framerate*tCapt; %Total frames to acquire
        vidobj.FramesPerTrigger = nImagesAdv; %Total frames to acquire
        vidobj.FrameGrabInterval = frameIntervalAdv; %Frame interval

        %- Execution of measurement
        fprintf("------------------------------------------------------------------- \n");
        fprintf("Execution of advancing measurement %d of %d \n",k,nStopsAdv);
        clockVidCreateAdv = toc;
        clockAvanco = clockAvanco + clockVidCreateAdv; %Update total advancing time

        %- Additional increase of drop volume
        tic
        fprintf('Please wait, additional increase of drop volume... \n');
        if taxaDepAdv <= 0.32 %Flow rate less than minimum motor rotation speed. Motor will work in pulses
            smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
            wMotorEns = uLs2RPM(taxaDepAdv,passoSer); %Motor rotation speed (RPM) to achieve the desired fow rate (uL/s)
            tadpasso = tadicionalpasso(passoMot,wMotorEns,wMin); %Time added to each step (in seconds) to reach the desired motor rotational speed
                for i = 1:floor(nPassosvolDepAdv) %Additional increase of drop volume
                    move(smSeringa,1); %Executes 1 rotation step on the motor
                    pause(tadpasso);
                end
        elseif (0.32 < taxaDepAdv) && (taxaDepAdv <= 1)
            wRPM = fuLs2RPM(taxaDepAdv); %Calculation of the necessary rotational speed to be sent to the motor in RPM
            smSeringa.RPM = floor(wRPM);
            move(smSeringa,floor(nPassosvolDepAdv)); %Motor movement
        end
        release(smSeringa); %Motor release
        fprintf('Additional increase of drop volume completed! \n');
        erroVolDepAdv = -passo2vol(nPassosvolDepAdv - floor(nPassosvolDepAdv),passoMot,passoSer); %Error in the volume released during advancing (uL)
        erroVolTotalAdv = erroVolTotalAdv + erroVolDepAdv;
        volSerTotal = volSerTotal - volGotaVarAdv; %Update total liquid volume in the syringe
        pause(1);
        clockAdv(k) = toc;
        clockAvanco = clockAvanco + clockAdv(k); %Update total advancing time

        %- Syringe displacement (syringe departure from the drop)
        tic
        fprintf('Please wait, syringe displacement... \n');
        moverSeringa(passoMot,passoFuso,smSeringaZ,"up",dispZSer,"high"); %Syringe movement
        release(smSeringaZ); %Motor release
        fprintf('Syringe displacement completed! \n');
        clockDispSerAdv1(k) = toc;
        clockAvanco = clockAvanco + clockDispSerAdv1(k); %Update total advancing time

        %- Waiting time for drop reach equilibrium
        tic
        fprintf('Please wait, waiting time for drop reach equilibrium... \n');
        if clockDispSerAdv1(k) >= tEsperaAdv %Syringe displacement time is longer than waiting time for drop equilibrium
            tEsperaAdvReal(k) = clockDispSerAdv1(k); %Calculation of actual waiting time
        else %Syringe displacement time is less than waiting time for drop equilibrium
            tEsperaAdvReal(k) = tEsperaAdv;
            deltat = tEsperaAdv - clockDispSerAdv1(k);
            pause(deltat)
        end
        %pause(tEsperaAdv);
        fprintf('Waiting time completed! \n');
        clockEspAdv = toc;
        clockAvanco = clockAvanco + clockEspAdv; %Update total advancing time

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
        clockVidAdv(k) = toc;
        clockAvanco = clockAvanco + clockVidAdv(k); %Update total advancing time
        
        %- Showing an example image
        tic
        close all %Close all open windows
        frame = getsnapshot(vidobj); %Capture an image
        imshow(frame); %Show image
        movegui(gcf,'center') %Centralize image window
        set(gcf,'Name',strcat(fileName,'_Example'))
        clockExAdv = toc;
        clockAvanco = clockAvanco + clockExAdv; %Update total advancing time

        %- Converting video to images
        tic
        fprintf('Please wait, converting video to images... \n');
        fullfilevidPath = fullfile(filePath,fullvidName); %Full path of the video
        vidFile = VideoReader(fullfilevidPath); %Read Video File
        for m = 1: nImagesAdv %Getting a image from every video frame
            frame = read(vidFile,m);
            framegray = rgb2gray(frame);
            if iAdv < 10 %Modifying the file name prefix to arrange files in ascending volume order
                prefixo = "00";
            elseif (iAdv >= 10) && (iAdv < 100)
                prefixo = "0";
            elseif (iAdv >= 100) && (iAdv < 1000)
                prefixo = "";
            end
            fullimName = strcat(prefixo,num2str(iAdv),'_',fileName,imFmt);
            fullfileimPath = fullfile(filePath,fullimName);
            imwrite(framegray,fullfileimPath);
            iAdv = iAdv + 1; %Interval of image name
        end
        delete(fullfilevidPath) %Delete video
        fprintf('Conversion of video to images completed! \n');
        clockImagesAdv(k) = toc;
        clockAvanco = clockAvanco + clockImagesAdv(k); %Update total advancing time
        
        %- Syringe displacement (syringe return to the drop)
        tic
        fprintf('Please wait, syringe displacement... \n');
        moverSeringa(passoMot,passoFuso,smSeringaZ,"down",dispZSer,"high"); %Syringe movement
        release(smSeringaZ); %Motor release
        fprintf('Syringe displacement completed! \n');
        clockDispSerAdv2(k) = toc;
        clockAvanco = clockAvanco + clockDispSerAdv2(k); %Update total advancing time

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
        clockCheckSerAdv = toc;
        clockAvanco = clockAvanco + clockCheckSerAdv; %Update total advancing time

        %- Deleting video object
        tic
        delete(vidobj)
        clear vidobj vidsrc
        clockExcVidAdv = toc;
        clockAvanco = clockAvanco + clockExcVidAdv; %Update total advancing time
    end

    %- Export/save experiment description
    filenametxt = "Description.txt";
    fullfilenametxt = fullfile(filePath,filenametxt);
    varnames = {'Datetime','volGotaStart[uL]','volGotaAdvIni[uL]','taxaDepAdvIni[uL]','volGotaAdvMax[uL]','taxaDepAdv[uL/s]','deslZSer[mm]','tEsperaAdv[s]','tEsperaAdvReal[s]','tCaptAdv[s]','nImagesAdv','nStopsAdv','volGotaVarAdv[uL]','totalFrameAdv','frameIntervalAdv','erroVolDepAdvIni[uL]','erroVolDepAdv[uL]','erroVolTotalAdv[uL]','clockAdvGotaIni[s]','clockAdv[s]','clockDispSerAdv[s]','clockVidAdv[s]','clockImagesAdv[s]','clockTotalAvanco[s]','frameRate','VidCompression','VidQuality'}; %Variables name in txt file
    T = table(datetime,volGotaStart,volGotaAdvIni,taxaDepAdvIni,volGotaAdvMax,taxaDepAdv,dispZSer,tEsperaAdv,round(mean(tEsperaAdvReal),2),tCaptAdv,nImagesAdv,nStopsAdv,volGotaVarAdv,totalFrameAdv,frameIntervalAdv,round(erroVolDepAdvIni,4),round(erroVolDepAdv,4),round(erroVolTotalAdv,4),round(clockAdvGotaIni,2),round(mean(clockAdv),2),round(mean(mean(clockDispSerAdv1)+mean(clockDispSerAdv2)),2),round(mean(clockVidAdv),2),round(mean(clockImagesAdv),2),round(clockAvanco,2),framerate,string(vidComp),vidQual,'VariableNames',varnames);
    writetable(rows2vars(T),fullfilenametxt,'Delimiter',' ','WriteVariableNames',false) %Write table to .txt file
    pause(5);
end
end

