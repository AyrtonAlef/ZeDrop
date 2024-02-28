function [posMesaX,posMesaY] = ajustePosMesa(smMesaX,smMesaY,posMesaX,posMesaY,passoMot,passoFuso)
%AJUSTEPOSMESA Adjust sample stage position (solid substrate).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Allows an adjustment of the sample stage position (solid substrate). It
%   allows the user to choose the axis of the movement (X/Y) as well as the
%   value [mm], direction and speed [mm/min].  The program automatically
%   identifies the direction of the movement according to the informed
%   displacement value (positive or negative). Three speeds are available:
%   high, corresponding to a speed of 9.5 mm/min, medium corresponding to a
%   speed of 6.09 mm/min, and low, corresponding to a speed of 2,864
%   mm/min. The movement of the stage is restricted via software by
%   previously established maximum displacement limits of each axis
%   (posXlimPos, posXlimNeg, posYlimPos, posYlimNeg). *Requires the
%   connection of the positioning stage system.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT
% smMesaX - Control of the X axis motor of the positioning stage
% smMesaY - Control of the Y axis motor of the positioning stage
% posMesaX - Current X axis position of stage in mm  
% posMesaY - Current Y axis position of stage in mm
% passoMotor - Number of steps required for one revolution in the motor [steps/rev]
% passoSeringa - Amount of volume displaced in the syringe corresponding to one revolution of its plunger [uL/rev]
%   OUTPUT
% posMesaX - Current X axis position of stage in mm  
% posMesaY - Current Y axis position of stage in mm  

fprintf('--------------- ADJUST STAGE POSITION (SOLID SUBSTRATE) --------------- \n');
%Routine of stage position adjustment
posXlimPos = 8; %Positive position limit on the X axis
posXlimNeg = -8; %Negative position limit on the X axis
posYlimPos = 24; %Positive position limit on the Y axis
posYlimNeg = -24; %Negative position limit on the Y axis

while 1
    displayPosMesa(posMesaX,posMesaY);
    resp = input("Do you like to proceed with positioning stage adjustment (Y/N) [N]? ","s");
    if resp == "Y"
        eixo = input("Do you like to move the stage in the X-axis (X) or the Y-axis (Y) [Y]? ","s");
        if isempty(eixo)
            eixo = "Y";
        end
        if eixo == "X"
            eixoMotor = smMesaX;
        elseif eixo == "Y"
            eixoMotor = smMesaY;
        end
        while 1
            deslmm = input("Displacement (mm): ");
            if eixo == "X"
                if ((posMesaX + deslmm) > posXlimPos) || ((posMesaX + deslmm) < posXlimNeg)
                    fprintf("Invalid displacement. Positioning stage limit reached. \n")
                else
                    break;
                end
            elseif eixo =="Y"
                if ((posMesaY + deslmm) > posYlimPos) || ((posMesaY + deslmm) < posYlimNeg)
                    fprintf("Invalid displacement. Positioning stage limit reached. \n")
                else
                    break;
                end
            end
        end
        if deslmm >= 0
            sentido = "positivo";
        elseif deslmm < 0
            sentido = "negativo";
        end
        velDesl = input("Speed (low/medium/high) [low]: ", "s");
        if isempty(velDesl)
            velDesl = "low";
        end
        fprintf('Please wait, moving stage... \n');
        moverMesa(passoMot,passoFuso,eixoMotor,sentido,abs(deslmm),velDesl); %Movendo mesa
        if eixo == "X"
            release(smMesaX); %Liberação do motor
            posMesaX = posMesaX + deslmm; %Atualização da posição atual da mesa;
        elseif eixo =="Y"
            release(smMesaY); %Liberação do motor
            posMesaY = posMesaY + deslmm; %Atualização da posição atual da mesa;
        end
        fprintf('Stage positioning completed! \n');
    else
        break;
    end
end
end

