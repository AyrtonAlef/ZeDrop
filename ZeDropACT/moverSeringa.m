function  moverSeringa(passoMotor,passoFuso,motor,sentido,deslZmm,velocidade)
%MOVERSERINGA Move the syringe to a predetermined displacement (delZmm), 
% direction (sentido) and speed (velocidade).
%   INPUT
% passoMotor - Number of steps required for one revolution in the motor [steps/rev]
% passoFuso - Linear displacement of the spindle corresponding to one revolution of the spindle [mm/rev]
% motor - Syringe displacement motor control
% sentido - Syringe displacement direction (up/down)
% deslZmm - Displacement [mm]
% velocidade - Displacement speed (low/medium/high)

nPassos = mm2passo(passoMotor,passoFuso,deslZmm); %Conversion of displacement in mm to the required number of motor steps
if strcmp(sentido,'up') %Determining the direction of movement on the desired axis
    i = 1;
elseif strcmp(sentido,'down')
    i = -1;
end
if strcmp(velocidade,'low') %Determination of displacement speed
    motor.RPM = 2; %RPM = 2 -> wreal = 1,432 RPM -> vbaixa = 2,864 mm/min
elseif strcmp(velocidade,'medium') 
    motor.RPM = 8; %RPM = 8 -> wreal = 3,045 RPM -> vmedia = 6,09 mm/min
elseif strcmp(velocidade,'high')
    motor.RPM = 113; %RPM = 113 -> wreal = 4,75 RPM -> vmedia = 9,5 mm/min
end
move(motor,i*nPassos); %Move stepper motor
end



