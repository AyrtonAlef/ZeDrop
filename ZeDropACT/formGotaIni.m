function [volSerTotal] = formGotaIni(smSeringa,passoMot,passoSer,volSerTotal)
%FORMGOTAINI Initial drop formation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Routine for forming an initial drop connected to the needle. Allows the
%   user to choose the drop volume [uL] as well as the drop formation flow 
%   rate [uL/s]. Three flow rates are available: high, corresponding to 
%   flow rate of 1.0607 uL/s, medium, corresponding to a flow rate of
%   0.6715 uL/s, and low, corresponding to a flow rate of 0.3158 uL/s.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT
% smSerigna - Syringe motor control
% passoMotor - Number of steps required for one revolution in the motor [steps/rev]
% passoSeringa - Amount of volume displaced in the syringe corresponding to one revolution of its plunger [uL/rev]
% volSerTotal - Total liquid volume in the syringe in uL
%   OUTPUT
% volSerTotal - Total liquid volume in the syringe in uL

fprintf('----------------------- INITIAL DROP FORMATION ------------------------\n');
volGotaIni = input("Informe o volume da gota inicial (uL): "); %Initial volume of the drop formed at the needle tip (uL)
if volGotaIni >= volSerTotal  %Check that there is enough volume in the syringe to perform the procedure
    fprintf("Unable to complete the procedure. Please increase the liquid volume of the syringe.")
else
    nPassosvolGotaIni = vol2passo(volGotaIni,passoMot,passoSer); %No of steps corresponding to the initial volume of the drop formed at the tip of the needle
    taxaCol = input("Liquid flow rate (low/medium/high) [low]: ", "s");
        if taxaCol == "high"
            smSeringa.RPM = 300; %Maximum possible RPM. Corresponds to an actual RPM = 4.8526 and a flow rate = 1.0607 uL/s 
        elseif taxaCol == "medium"
            smSeringa.RPM = 8; %Medium RPM. Corresponds to an actual RPM = 3.045 and a flow rate = 0.6715 uL/s 
        else
            smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPM = 1.4171 and a flow rate = 0.3158 uL/s 
        end
    fprintf('Please wait, drop formation...\n');
    move(smSeringa,floor(nPassosvolGotaIni)); %Drive motor connected to the syringe
    release(smSeringa); %MOtor release
    fprintf('Initial drop formation completed!\n');
    errovolGotaIni = passo2vol(nPassosvolGotaIni - floor(nPassosvolGotaIni),passoMot,passoSer); %Error in the initial volume of the drop formed at the needle tip (uL)
    volSerTotal = volSerTotal - volGotaIni; %Update total volume of liquid in the syringe
end
pause(2);
end

