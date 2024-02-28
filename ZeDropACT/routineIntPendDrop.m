function [volSerTotal] = routineIntPendDrop(vidAdap,vidID,vidFmt,imFmt,smSeringa,passoMot,passoSer,wMin,framerate,volSerTotal)
%ROUTINEINTPENDDROP Execution of an intermittent pendant drop experiment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Execução de ensaio intermitente de gota pendente. Executa aumentos 
%   sucessivos e intermitentes no volume da gota até o alcance de um volume 
%   máximo pré-estalecido. Entre os aumentos, imagens da gota são 
%   capturadas possibilitando posterior análise. A rotina requer que o 
%   usuário defina o volume iniical da gota (volGotaIni), o volume máximo 
%   da gota (volGotaMax), a taxa de deposição da gota (taxaDepGota), o 
%   tempo de espera para a gota entrar em equilíbrio (tEsp), o tempo de 
%   duração da captura de imagens durante uma parada (tCapt), o número de 
%   imagens que se deseja capturar durante uma parada (nImages) e o número 
%   de paradas ao longo do aumento intermitente da gota (nStops). Para a 
%   formação da gota inicial, a mínima taxa de deposição é usada 
%   (0,3158 uL/s). A rotina permite a escolha do diretório em que as 
%   imagens serão salvas, assim como a definição do nome das imagens. Ao 
%   final, um arquivo de texto é exportado com todos os parâmetros 
%   relevantes do ensaio.
%   Execution of an intermittent pendant drop experiment. Performs
%   successive and intermittent increases in the drop volume until reaching
%   a pre-established maximum volume. During drop increase, drop images are
%   captured for further analysis. The routine requires the user to define 
%   the initial drop volume (volGotaIni), the maximum drop volume 
%   (volGotaMax), the drop flow rate (taxaDepGota), the waiting time for 
%   drop equilibrium (tEsp), the time of image capture during a stop 
%   (tCapt), the number of captured images during a stop (nImages) and the 
%   number of stops along the intermittent drop increase (nStops). For 
%   initial drop formation, the minimum flow rate is used (0.3158 uL/s). 
%   The routine allows choosing the directory where the images will be 
%   exported, as well as defining the name of the images. At the end, a 
%   text file is exported with all relevant parameters of the experiment.
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

fprintf('----------------- PENDANT DROP - INTERMITTENT EXPERIMENT -------------------\n');
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
tCapt = input("- Capture time during a stop (s): "); %Capture time during a stop
%tCaptInterval = input("- Capture time between images (s): ");
nImages = input("- Number of captured images during a stop: "); %Number of captured images during a stop
nStops = input("- Number of stops: "); %Number of stops

volGotaVar = (volGotaMax - volGotaIni)/nStops; %Volume variation during each stop
tCaptInterval = tCapt/nImages; %Acquissition time between images in s
%volCaptInterval = taxaDepGota*tCaptInterval; %Volume variation between images in uL
totalFrame = framerate*tCapt; %Total number of frames acquired during a stop
frameInterval = totalFrame/nImages; %Interval between frames to acquire the desired number of images
erroVolTotal = 0; %Error in the total drop volume [uL]
clockPrepVid = zeros(1,nStops); %Time spent in the video preparation in each stop
clockDepGota = zeros(1,nStops); %Time spent to volume addition in each stop
clockEspGota = zeros(1,nStops); %Time spent to wait for drop equilibrium before video recording in each stop
clockVid = zeros(1,nStops); %Time spent for video recording in each stop
clockImages = zeros(1,nStops); %Time spent for converitng video to images in each stop
clockStop = zeros(1,nStops); %Total time spent in each stop
clockEnsaio = 0; %Total time to execute experiment

if volGotaMax >= volSerTotal %Check that there is enough volume in the syringe to perfrom the procedure
    fprintf("Unable to complete the routine. Please increase the liquid volume in the syringe. \n")
else
    %- Drive motor settings
    nPassosvolDepEns = vol2passo(volGotaVar,passoMot,passoSer); %%No of steps corresponding to the released volume in each stop

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
    clockEnsaio = clockEnsaio + clockEnsaio1; %Update total experiment time

    %- Initital drop formation
    tic
    nPassosvolGotaIni = vol2passo(volGotaIni,passoMot,passoSer); %No of steps corresponding to the initial drop volume formed at the needle tip
    smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
    fprintf('Please wait, initial drop formation...\n');
    move(smSeringa,floor(nPassosvolGotaIni)); %Drive motor connected to the syringe
    release(smSeringa); %Motor release
    fprintf('Initial drop formation completed! \n');
    erroVolGotaIni = passo2vol(nPassosvolGotaIni - floor(nPassosvolGotaIni),passoMot,passoSer); %Error in the volume during the initial drop formation (uL)
    erroVolTotal = erroVolTotal + erroVolGotaIni;
    volSerTotal = volSerTotal - volGotaIni; %Update of total liquid volume in the syringe
    clockGotaIni = toc; %Time for the initial drop formation in s
    clockEnsaio = clockEnsaio + clockGotaIni; %Update total experiment time

    %- Waiting time for drop equilibrium
    tic
    fprintf('Please wait, waiting time for drop reach equilibrium... \n');
    pause(tEspera); %Waiting time for drop equilibrium
    fprintf('Waiting time completed! \n');
    clockEspGotaIni = toc; %Waiting time for initial drop equilibrium
    clockEnsaio = clockEnsaio + clockEspGotaIni; %Update total experiment time

    %- Start of measurement series
    for k = 1:nStops
        %clockStop(k) = 0; %Time for the execution of each stop
        tic
        %- Create video object
        [vidobj, vidsrc] = createVidObj(vidAdap,vidID,vidFmt);

        %- Update of image names
        fileName = strcat(name,'_',num2str(k),'_V(uL)_',num2str((volGotaIni+k*volGotaVar),'%.2f'));
        fullvidName = strcat(fileName,vidDiskFmt);
        fullfilevidPath = fullfile(filePath,fullvidName);

        %- Preparation of video object
        diskLogger = VideoWriter(fullfilevidPath,vidComp);
        diskLogger.FrameRate = framerate; %Matching the video frame rate to the actual camera frame rate
        diskLogger.Quality = vidQual; %Definition of the video quality
        vidobj.DiskLogger = diskLogger;
        %vidobj.FramesPerTrigger = framerate*tCapt; %Total frames to acquire
        vidobj.FramesPerTrigger = nImages; %Total frames to acquire
        vidobj.FrameGrabInterval = round(frameInterval); %Frame interval
        clockPrepVid(k) = toc; %Time spent for creating and preparing the video object
        clockStop(k) = clockStop(k) + clockPrepVid(k); %Update time to execute each stop
        clockEnsaio = clockEnsaio + clockPrepVid(k); %Update total experiment time

        %- Execution of measurement/stop
        fprintf("----------------------------------------------- \n");
        fprintf("Execution of measurement %d of %d \n",k,nStops);
        
        %- Additional increase of drop volume
        tic
        fprintf('Please wait, additional increase of drop volume... \n');
        if taxaDepGota <= 0.32 %Flow rate less than minimum motor rotation speed. Motor will work in pulses
            smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s
            wMotorEns = uLs2RPM(taxaDepGota,passoSer); %Motor rotation speed (RPM) to achieve the desired fow rate (uL/s)
            tadpasso = tadicionalpasso(passoMot,wMotorEns,wMin); %Time added to each step (in seconds) to reach the desired motor rotational speed
                for i = 1:floor(nPassosvolDepEns) %Additional increase of drop volume
                    move(smSeringa,1); %Executes 1 rotation step on the motor
                    pause(tadpasso);
                end
        elseif (0.32 < taxaDepGota) && (taxaDepGota <= 1)
            wRPM = fuLs2RPM(taxaDepGota); %Calculation of the necessary rotational speed to be sent to the motor [RPM] in order to reach the desired flow rate [uL/s]     
            smSeringa.RPM = floor(wRPM);
            move(smSeringa,floor(nPassosvolDepEns)); %Motor movement
        end
        release(smSeringa); %Motor release
        fprintf('Additional increase of drop volume completed! \n');
        erroVolDepEns = passo2vol(nPassosvolDepEns - floor(nPassosvolDepEns),passoMot,passoSer); %Error in the additional increase of drop volume (uL)
        erroVolTotal = erroVolTotal + erroVolDepEns;
        volSerTotal = volSerTotal - volGotaVar; %Update of total liquid volume in the syringe
        clockDepGota(k) = toc; %Time spent to the additional increase of drop volume
        clockStop(k) = clockStop(k) + clockDepGota(k); %Update time to execute each stop
        clockEnsaio = clockEnsaio + clockDepGota(k); %Update total experiment time

        %- Waiting time for drop equilibrium
        tic
        fprintf('Please wait, waiting time for drop reach equilibrium... \n');
        pause(tEspera);
        fprintf('Waiting time completed! \n');
        clockEspGota(k) = toc; %%Waiting time for drop equilibrium before image capture
        clockStop(k) = clockStop(k) + clockEspGota(k); %Update time to execute each stop
        clockEnsaio = clockEnsaio + clockEspGota(k); %Update total experiment time

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
        %toc
        clockVid(k) = toc;
        clockStop(k) = clockStop(k) + clockVid(k); %Update time to execute each stop
        clockEnsaio = clockEnsaio + clockVid(k); %Update total experiment time
        
        %- Showing an example image
        tic
        close all %Close all open windows
        frame = getsnapshot(vidobj); %Capture an image
        imshow(frame); %Show image
        movegui(gcf,'center') %Centralize image window
        set(gcf,'Name',strcat(fileName,'_Example'))
        clockEnsaio2 = toc;
        clockStop(k) = clockStop(k) + clockEnsaio2; %Update time to execute each stop
        clockEnsaio = clockEnsaio + clockEnsaio2; %Update total experiment time

        %- Converting video to images
        tic
        fprintf('Please wait, converting video to images... \n');
        %vid2frame(filePath,fullvidName,fileName,imFmt,tCaptInterval);
        fullfilevidPath = fullfile(filePath,fullvidName); %Full path of the video
        vidFile = VideoReader(fullfilevidPath); %Read Video File
        %i = 1;
        for m = 1: nImages %Get an image from every video frame
            frame = read(vidFile,m);
            framegray = rgb2gray(frame);
            fullimName = strcat(fileName,'_',num2str(m),imFmt);
            fullfileimPath = fullfile(filePath,fullimName);
            imwrite(framegray,fullfileimPath);
            %i = i + 1; %Interval of image name
        end
        delete(fullfilevidPath) %Delete video
        fprintf('Conversion of video to images completed! \n');
        clockImages(k) = toc;
        clockStop(k) = clockStop(k) + clockImages(k); %Update time to execute each stop
        clockEnsaio = clockEnsaio + clockImages(k); %Update total experiment time

        %- Deleting video object
        tic
        delete(vidobj)
        clear vidobj vidsrc
        clockEnsaio3 = toc;
        clockStop(k) = clockStop(k) + clockEnsaio3; %Update time to execute each stop
        clockEnsaio = clockEnsaio + clockEnsaio3; %Update total experiment time
    end

    %- Export/save experiment description
    filenametxt = "Description.txt";
    fullfilenametxt = fullfile(filePath,filenametxt);
    varnames = {'Datetime','volGotaIni[uL]','volGotaMax[uL]','nStops','taxaDepGota[uL/s]','tEspera[s]','tCaptura[s]','tIntervaloCaptura[s]','erroVolGotaIni[uL]','erroVolDepStop[uL]','erroVolTotal[uL]','clockGotaIni[s]','clockEspGotaIni[s]','clockAvgPrepVid[s]','clockTotPrepVid[s]','clockAvgDepGota[s]','clockTotDepGota[s]','clockAvgEspGota[s]','clockTotEspGota[s]','clockAvgVid[s]','clockTotVid[s]','clockAvgImages[s]','clockTotImages[s]','clockAvgStop[s]','clockTotStop[s]','clockEnsaio[s]','frameRate','VidCompression','VidQuality'}; %Variables name in txt file
    T = table(datetime,volGotaIni,volGotaMax,nStops,taxaDepGota,tEspera,tCapt,tCaptInterval,round(erroVolGotaIni,4),round(erroVolDepEns,4),round(erroVolTotal,4),round(clockGotaIni,2),round(clockEspGotaIni,2),round(mean(clockPrepVid),2),round(sum(clockPrepVid),2),round(mean(clockDepGota),2),round(sum(clockDepGota),2),round(mean(clockEspGota),2),round(sum(clockEspGota),2),round(mean(clockVid),2),round(sum(clockVid),2),round(mean(clockImages),2),round(sum(clockImages),2),round(mean(clockStop),2),round(sum(clockStop),2),round(clockEnsaio,2),framerate,string(vidComp),vidQual,'VariableNames',varnames);
    writetable(rows2vars(T),fullfilenametxt,'Delimiter',' ','WriteVariableNames',false) %Write table to .txt file

    pause(5);
end
end

