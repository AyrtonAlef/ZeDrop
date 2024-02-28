function  moverMesa(passoMotor,passoFuso,eixo,sentido,deslmm,velocidade)
%MOVERMESA Move the stage (solid substrate) along an axis and according to 
% a predetermined displacement (deslmm), direction (sentido) and speed 
% (velocidade).
%   INPUT
% passoMotor - Number of steps required for one revolution in the motor [steps/rev]
% passoFuso - Linear displacement of the spindle corresponding to one revolution of the spindle [mm/rev]
% eixo - Motor control of the axis on which movement will proceed
% sentido - Displacement direction (positive/negative)
% deslmm - Displacement [mm]
% velocidade - Displacement speed (high/medium/low)

nPassos = (passoMotor/passoFuso)*deslmm; %Conversion of displacement in mm to the required number of motor steps
if strcmp(sentido,'positivo') %Determining the direction of movement on the desired axis
    i = 1;
elseif strcmp(sentido,'negativo')
    i = -1;
end
if strcmp(velocidade,'low') %Determination of displacement speed
    eixo.RPM = 2; %RPM = 2 -> wreal = 1,432 RPM -> vbaixa = 2,864 mm/min
elseif strcmp(velocidade,'medium') 
    eixo.RPM = 8; %RPM = 8 -> wreal = 3,045 RPM -> vmedia = 6,09 mm/min
elseif strcmp(velocidade,'high')
    eixo.RPM = 113; %RPM = 113 -> wreal = 4,75 RPM -> vmedia = 9,5 mm/min
end
move(eixo,i*round(nPassos)); %Move stepper motor
end

