clc
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% COMMUNICATION WITH ARDUINO AND MOTORS %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1)Initializing communication
fprintf('------------------- COMMUNICATION WITH ARDUINO AND MOTORS -------------------------\n');
%porta = 'COM5'; %Connection port with Aruino
porta = input('Enter the communication port with the Arduino board [COM6]: ','s');
if isempty(porta)
    porta = 'COM6';
end
fprintf('Initializing communication with the Arduino board... \n');
a = arduino(porta,'UNO','Libraries','Adafruit\MotorShieldV2'); %Initializing communication with the Arduino board
shieldbot = addon(a,'Adafruit\MotorShieldV2'); %Initializing communication with the lower Adafruit MotorShield (microsyringe control)
shieldtop = addon(a,'Adafruit\MotorshieldV2','I2CAddress','0x61'); %Initializing communication with the upper Adafruit MotorShield (sample stage control)
%2)Initializing communication with motors
fprintf('Initializing motors... \n');
%   2.1) Microsyringe motor
smprevSer = 200; %Motor steps/rev. (read-only property)
smnumbSer = 1; %Number of the connection in which the stepper motor is plugged into (M1&M2 = 1; M3&M4 = 2)(read-only property)
smtypeSer = 'Microstep'; %Motor step type ('Single','Double','Interleave' e 'MicroStep')(read-only property)
smRPMSer = 10; %Motor rotation speed (RPM)
smSeringa = stepper(shieldbot,smnumbSer,smprevSer,'stepType',smtypeSer,'RPM',smRPMSer); %Initializing microsyringe motor control
%   2.2) Microsyringe positioning motor - Z Axis
smprevSerZ = 200; %Motor steps/rev.
smnumbSerZ = 2; %Number of the connection in which the stepper motor is plugged into (M1&M2 = 1; M3&M4 = 2)(read-only property)
smtypeSerZ = 'Microstep'; %Motor step type ('Single','Double','Interleave' e 'MicroStep')(read-only property)
smRPMSerZ = 10; %Motor rotation speed (RPM)
smSeringaZ = stepper(shieldbot,smnumbSerZ,smprevSerZ,'stepType',smtypeSerZ,'RPM',smRPMSerZ); %Initializing the control of the microsyringe positioning motor (Z Axis)
%%{
modulo = input('Do you like to activate the positioning stage (P) or the tilting stage (T) [P]? ','s');
if isempty(modulo)
    modulo = 'P';
end
if modulo == 'P'
    %   2.3) Positioning stage motor - X Axis
    smprevMesaX = 200; %Motor steps/rev. (read-only property)
    smnumbMesaX = 2; %Number of the connection in which the stepper motor is plugged into (M1&M2 = 1; M3&M4 = 2)(read-only property)
    smtypeMesaX = 'Microstep'; %Motor step type ('Single','Double','Interleave' e 'MicroStep')(read-only property)
    smRPMMesaX = 10; %Motor rotation speed (RPM)
    smMesaX = stepper(shieldtop,smnumbMesaX,smprevMesaX,'stepType',smtypeMesaX,'RPM',smRPMMesaX); %Initializing the control of the positioning stage motor (X Axis)
    %   2.3) Positioning stage motor - Y Axis
    smprevMesaY = 200; %Motor steps/rev. (read-only property)
    smnumbMesaY = 1; %Number of the connection in which the stepper motor is plugged into (M1&M2 = 1; M3&M4 = 2)(read-only property)
    smtypeMesaY = 'Microstep'; %Motor step type ('Single','Double','Interleave' e 'MicroStep')(read-only property)
    smRPMMesaY = 10; %Motor rotation speed (RPM)
    smMesaY = stepper(shieldtop,smnumbMesaY,smprevMesaY,'stepType',smtypeMesaY,'RPM',smRPMMesaY); %Initializing the control of the positioning stage motor (Y Axis)
    mode = 1; %Defines which system is connected (1 - Positioning stage system/ 2- Tilting stage system)
    %}
    %%{
elseif modulo == 'T'
    %   2.4) Tilting stage
    smprevInc = 200; %Motor steps/rev. (read-only property)
    smnumbInc = 1; %Number of the connection in which the stepper motor is plugged into (M1&M2 = 1; M3&M4 = 2)(read-only property)
    smtypeInc = 'Microstep'; %Motor step type ('Single','Double','Interleave' e 'MicroStep')(read-only property)
    smRPMInc = 10; %Motor rotation speed (RPM)
    smInc = stepper(shieldtop,smnumbInc,smprevInc,'stepType',smtypeInc,'RPM',smRPMInc); %Initializing the control of the tilting stage motor
    mode = 2; %Defines which system is connected (1 - Positioning stage system/ 2- Tilting stage system)
    %}
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN ROUTINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
% Goniometer configuration
passoSer = 13.23; %Microsyringe step 13.23 uL/rev.
passoMot = 200; %Motor step 200 steps/rev.
passoFuso = 2; %Screw spindle pitch 2 mm/rev.
nEntradasSemFim = 1; %Worm number of starts
nDentesCoroa = 104; %Number of worm gear teeth
wMin = 1.4171; %Minimum stable rotation (RPM) that motor can reach. Corresponds to sm.RPM = 2 (Microstep) (determined empirically)
wMax = 4.8526; %Maximum stable rotation (RPM) that motor can reach. Corresponds to sm.RPM = 2 (Microstep) (determined empirically)

% Video configuration
vidFmt = 'MJPG_1920x1080'; %Video format
vidAdap = 'winvideo';
vidID = 1;
limTotalFrames = 300; %Total limit of frames that the camera is able to store during a continuous capture
framerate = 30; %Real frame rate (suggestion: run a video preview and identify the frame rate)
imFmtPD = '.pnm'; %Image format for pendant drops ('.jpeg','.png','.tif','.pnm' and others)
imFmtSD = '.tif'; %Image format for sessile drops ('.jpeg','.png','.tif','.pnm' and others)
imFmtID = '.tif'; %Image format for inclined drops ('.jpeg','.png','.tif','.pnm' and others)

% Routine settings
%{
load VolumeSeringa %Load variable volSerTot (Monitorar volume total de líquido na seringa)
%volSerTot = 0; %Total volume available in syringe
%volSerTot = 200; %Total volume available in syringe
volAutLib = 20; %Released volume in uL after the automatic liquid acquisition routine
if mode == 1
    load PosMesa %Load variables posMesaX and posMesaY (actual position of the sample stage)
elseif mode == 2
    load IncMesa %Load variable incMesaTot (actual tilting of the sample stage)
end
%}
volAutLib = 20; %Released volume in uL after the automatic liquid acquisition routine

% Main
while 1
    %clc
    close all
    % Variables update
    load VolumeSeringa %Load variable volSerTot (monitor total volume of liquid in the syringe)
    displayVolSer(volSerTot); %Displays how much liquid volume is in the syringe
    if mode == 1
        load PosMesa %Load variables posMesaX and posMesaY (actual position of the sample stage)
        displayPosMesa(posMesaX,posMesaY); %Displays the current position of the sample stage
    elseif mode == 2
        load IncMesa %Load variable incMesaTot (actual tilting of the sample stage)
        displayIncMesa(incMesaTot); %Displays the current tilt of the sample stage
    end
    
    % Display live video
    [vidobj, vidsrc] = createVidObj(vidAdap,vidID,vidFmt); %Create vidobject
    gridRow = 8; %Number of grid lines placed over the video
    gridCol = 4; %Number of grid lines placed over the video
    liveVid(vidobj,gridRow,gridCol)
    
    % Menu
    fprintf('------------------------------- ZEDROP_ACT SOFTWARE ---------------------------------\n');
    %fprintf('Volume total de líquido [uL]: %.2f \n',volSerTot);
    fprintf('MENU: \n');
    fprintf('1. Adjust syringe liquid volume; \n');
    fprintf('2. Reset monitoring of syringe volume; \n');
    fprintf('3. Adjust syringe vertical position ; \n');
    if mode == 1
        fprintf('4. Adjust stage position (solid substrate); \n');
        fprintf('5. Center stage (solid substrate); \n');
        fprintf('6. Reset monitoring of stage position; \n');
    elseif mode == 2
        fprintf('4. Adjust stage tilt (solid substrate); \n');
        fprintf('5. Level stage (solid susbtrate); \n');
        fprintf('6. Reset monitoring of stage tilt; \n');
    end
    fprintf('7. Initial drop formation; \n');
    fprintf('8. Capture image; \n');
    fprintf('9. Pendant drop - intermittent experiment; \n');
    fprintf('10. Pendant drop - continuous experiment; \n');
    fprintf('11. Sessile drop (needle out) - intermittent experiment (advancing and receding); \n');
    fprintf('111. Sessile drop (needle out) - intermittent experiment (only advancing); \n');
    fprintf('112. Sessile drop (needle out) - intermittent experiment (only receding); \n');
    fprintf('12. Sessile drop (needle in) - continuous experiment (advancing and receding); \n');
    fprintf('121. Sessile drop (needle in) - cycles of continuous experiments (advancing and receding); \n');
    fprintf('122. Sessile drop (needle in) - continuous experiment (only advancing); \n');
    fprintf('123. Sessile drop (needle in) - continuous experiment (only receding); \n');
    if mode == 2
        fprintf('13. Inclined drop - intermittent experiment; \n');
        fprintf('14. Inclined drop - continuous experiment; \n');
    end
    
    opcao = input('\n Enter one of the above options [0 to exit]: ');
    
    if opcao == 1 %Adjust syringe liquid volume
        clc
        volSerTot = ajusteVolSeringa(smSeringa,passoMot,passoSer,volSerTot,volAutLib);
    elseif opcao == 2 %Reset monitoring of syringe volume
        clc
        volSerTot = 0;
    elseif opcao == 3 %Adjust syringe position (z axis)
        clc
        ajustePosSeringa(smSeringaZ,passoMot,passoFuso);
    elseif opcao == 4 %Adjust stage position (solid substrate)/ Level stage (solid susbtrate)
        clc
        if mode == 1
            [posMesaX,posMesaY] = ajustePosMesa(smMesaX,smMesaY,posMesaX,posMesaY,passoMot,passoFuso);
        elseif mode == 2
            [incMesaTot] = ajusteIncMesa(smInc,incMesaTot,passoMot,nEntradasSemFim,nDentesCoroa);
        end
    elseif opcao == 5 %Center stage (solid substrate)/ Level stage (solid susbtrate)
        clc
        if mode == 1
            [posMesaX,posMesaY] = centralPosMesa(smMesaX,smMesaY,posMesaX,posMesaY,passoMot,passoFuso);
        elseif mode == 2
            [incMesaTot] = nivelIncMesa(smInc,incMesaTot,passoMot,nEntradasSemFim,nDentesCoroa);
        end
    elseif opcao == 6 %Reset monitoring of stage position/ Reset monitoring of stage tilt
        clc
        if mode == 1
            posMesaX = 0;
            posMesaY = 0;
        elseif mode == 2
            incMesaTot = 0;
        end
    elseif opcao == 7 %Initial drop formation
        clc
        volSerTot = formGotaIni(smSeringa,passoMot,passoSer,volSerTot);
    elseif opcao == 8 %Capture image
        clc
        captImage(vidAdap,vidID,vidFmt,imFmtSD);
    elseif opcao == 9 %Pendant drop - intermittent experiment
        clc
        volSerTot = routineIntPendDrop(vidAdap,vidID,vidFmt,imFmtPD,smSeringa,passoMot,passoSer,wMin,framerate,volSerTot);
    elseif opcao == 10 %Pendant drop - continuous experiment
        clc
        volSerTot = routineContPendDrop(vidAdap,vidID,vidFmt,imFmtPD,smSeringa,passoMot,passoSer,wMin,framerate,limTotalFrames,volSerTot);
    elseif opcao == 11 %Sessile drop (needle out) - intermittent experiment (advancing and receding)
        clc
        volSerTot = routineIntSessDrop(vidAdap,vidID,vidFmt,imFmtSD,smSeringa,smSeringaZ,passoMot,passoSer,passoFuso,wMin,framerate,volSerTot);
    elseif opcao == 111 %Sessile drop (needle out) - intermittent experiment (only advancing)
        clc
        volSerTot = routineIntAdvSessDrop(vidAdap,vidID,vidFmt,imFmtSD,smSeringa,smSeringaZ,passoMot,passoSer,passoFuso,wMin,framerate,volSerTot);
    elseif opcao == 112 %Sessile drop (needle out) - intermittent experiment (only receding)
        clc
        volSerTot = routineIntRecSessDrop(vidAdap,vidID,vidFmt,imFmtSD,smSeringa,smSeringaZ,passoMot,passoSer,passoFuso,wMin,framerate,volSerTot);
    elseif opcao == 12 %Sessile drop (needle in) - continuous experiment (advancing and receding)
        clc
        volSerTot = routineContSessDrop(vidAdap,vidID,vidFmt,imFmtSD,smSeringa,passoMot,passoSer,wMin,framerate,limTotalFrames,volSerTot);
    elseif opcao == 121 %Sessile drop (needle in) - cycles of continuous experiments (advancing and receding)
        clc
        volSerTot = routineContCycleSessDrop(vidAdap,vidID,vidFmt,imFmtSD,smSeringa,passoMot,passoSer,wMin,framerate,limTotalFrames,volSerTot);
    elseif opcao == 122 %Sessile drop (needle in) - continuous experiment (only advancing)
        clc
        volSerTot = routineContAdvSessDrop(vidAdap,vidID,vidFmt,imFmtSD,smSeringa,passoMot,passoSer,wMin,framerate,limTotalFrames,volSerTot);
    elseif opcao == 123 %Sessile drop (needle in) - continuous experiment (only receding)
        clc
        volSerTot = routineContRecSessDrop(vidAdap,vidID,vidFmt,imFmtSD,smSeringa,passoMot,passoSer,wMin,framerate,limTotalFrames,volSerTot);
    elseif opcao == 13 && mode == 2 %Inclined drop - intermittent experiment
        clc
        incMesaTot = routineIntIncDrop(vidAdap,vidID,vidFmt,imFmtID,smInc,passoMot,nEntradasSemFim,nDentesCoroa,wMin,wMax,framerate,incMesaTot);
    elseif opcao == 14 && mode == 2 %Inclined drop - continuous experiment
        clc
        incMesaTot = routineContIncDrop(vidAdap,vidID,vidFmt,imFmtID,smInc,passoMot,nEntradasSemFim,nDentesCoroa,wMin,wMax,framerate,limTotalFrames,incMesaTot);
    elseif opcao == 0 %Exit
        break;
    end
    clc

    %Save variables
    save VolumeSeringa volSerTot %Save variable volSerTot (Monitor total liquid volume in the syringe)
    save PosMesa posMesaX posMesaY %Save variable volSerTot (Monitor stage positioning)
    %save IncMesa incMesaTot %Save variable incMesaTot (Monitor stage tilt)
end
clc
close all