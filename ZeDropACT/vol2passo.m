function [nPassos] = vol2passo(vDep,passoMotor,passoSeringa)
%VOL2PASSO Returns how many motor steps are required to displace a given 
% volume [uL] of liquid in the syringe.
%   INPUT
% vDep - Volume of liquid displaced [uL]
% passoMotor - Number of steps required for one revolution in the motor [steps/rev]
% passoSeringa - Amount of volume displaced in the syringe corresponding to one revolution of its plunger [uL/rev]
%   OUPUT
% nPassos - Number of motor steps

nPassos = vDep * passoMotor/passoSeringa;
end

