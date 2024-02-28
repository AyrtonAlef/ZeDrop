function [incMesaTot] = ajusteIncMesa(smInc,incMesaTot,passoMotor,nEntradasSemFim,nDentesCoroa)
%AJUSTEINCMESA Adjust stage tilt (solid substrate).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Allows the adjustment of the stage tilting (solid substrate). Allows
%   the user to choose the angle [°] direction (clockwise/counterclockwise)
%   and tilt rate [°/s]. Three tilt rates are available: high,
%   corresponding to a tilt rate of 0.2750 °/s, medium, corresponding to a
%   tilt rate of 0.1760 °/s, and low, corresponding to a rate of 0.0826
%   °/s.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT
% smInc - Control of the motor connected to the tilting stage
% incMesaTot - Current stage tilt [°]
% passoMotor - Number of steps required for one revolution in the motor [steps/rev]
% nEntradasSemFim -Worm number of starts
% nDentesCoroa - Number of worm gear teeth
%   OUTPUT
% incMesaTot - Current stage tilt [°]

fprintf('-------------- ADJUST STAGE TILT (SOLID SUBSTRATE) -------------- \n');
%Routine of stage tilt adjustment
while 1
    displayIncMesa(incMesaTot);
    resp = input("Do you want to proceed with stage tilt adjustment (Y/N) [N]? ","s");
    if resp == "Y"
        angle = input("Tilt angle (°): ");
        sentido = input("Tilt direction (clockwise (CW) or counterclockwise (CCW)) [C]: ", "s");
        if isempty(sentido)
            sentido = "CW";
        end
        velocidade = input("Tilt speed (low/medium/high) [high]: ", "s");
        if isempty(velocidade)
            velocidade = "high";
        end
        fprintf('Please wait, tilting stage... \n');
        moverIncMesa(smInc,passoMotor,nEntradasSemFim,nDentesCoroa,angle,sentido,velocidade); %Stage tilting
        release(smInc); %Motor release
        fprintf('Stage tilting completed! \n');
        if sentido == "CW"
            i = 1;
        elseif sentido == "CCW"
            i = -1;
        end
        incMesaTot = incMesaTot + i*angle; %Update of total stage tilt
    else
        break;
    end
end
end

