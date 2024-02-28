function [vDep] = passo2vol(nPassos,passoMotor,passoSeringa)
%PASSOS2VOL Returns how much liquid volume [uL] id displaced in the syringe
% given a number of motor steps.
%   INPUT
% nPassos - NUmber of motor steps
% passoMotor - Number of steps required for one revolution in the motor [steps/rev]
% passoSeringa - Amount of volume displaced in the syringe corresponding to one revolution of its plunger [uL/rev]
%   OUPUT
% vDep - Volume displaced in the syringe [uL]

vDep = nPassos * passoSeringa/passoMotor;
end

