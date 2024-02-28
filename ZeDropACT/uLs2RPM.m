function [wMotor] = uLs2RPM(taxaDeposito,passoSeringa)
%ULS2RPM Returns the motor rotational speed [RPM] corresponding to a given flow rate (taxaDeposito)[uL/s].
%   INPUT
% taxaDeposito - flow rate [uL/s]
% passoSeringa - Amount of volume displaced in the syringe corresponding to one revolution of its plunger [uL/rev]
%   OUTPUT
% wMotor - motor rotational speed [RPM]

wMotor = 60*taxaDeposito/passoSeringa;
end