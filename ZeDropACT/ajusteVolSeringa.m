function [volSerTotal] = ajusteVolSeringa(smSeringa,passoMot,passoSer,volSerTotal,volAutLib)
%AJUSTEVOLSERINGA Adjust syringe liquid volume.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Allows the control of liquid volume in the syringe. Allows the user to
%   choose a manual (M) or an automatic (A) adjustment. In manual
%   adjustment (M), the user can define the volume [uL] to be released or
%   withdrawn as well as the rate of release or withdrawal of this volume
%   [uL/s]. Two speeds are available: high, corresponding to 1.0607 uL/s,
%   and low, corresponding to 0.3158 uL/s. In automatic adjustment (A), a
%   routine intake of probing liquid is followed by a release of a portion of
%   the acquired volume. In automatic adjustment, the user only informs the
%   volume to be collected, using an average flow rate (0.6715 uL/s) during
%   withdrawal and release. The volume of liquid released corresponds to
%   volAutLib. The release aims to remove possible air bubbles and
%   impurities.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT
%smSerigna - Syringe motor control
%passoMot - Number of steps required for one revolution in the motor [steps/rev]
%passoSer - Amount of volume displaced in the syringe corresponding to one revolution of its plunger [uL/rev]
%volSerTotal - Total liquid volume in the syringe in uL
%volAutLib - Released volume of liquid in uL during automatic adjustment in order to remove bubbles and contaminations
%   OUTPUT
%volSerTotal - Total liquid volume in the syringe in uL

fprintf('------------------------- ADJUST SYRINGE LIQUID VOLUME  ---------------------------\n');
%Enchimento da seringa
%resp = "Y";
while 1
    displayVolSer(volSerTotal); %Display how much probing liquid volume is in the syringe
    resp = input("Do you like to proceed with adjusting the volume of liquid in the syringe (Y/N) [N]? ","s");
    if resp == "Y"
        %Choosing the type of volume adjustment
        tipoAjuste = input("Do you like to proceed with a manual (M) or automatic (A) adjustment [A]? ","s");
        if tipoAjuste == "M" %Manual adjustment
            tipoColeta = input("Do you like to withdrawn (W) or release (R) probing liquid [R]? ","s");
            if isempty(tipoColeta)
                tipoColeta = 'R';
            end
            if tipoColeta == "W"
                volCol = input("Withdrawn volume (uL): ");
                taxaCol = input("Liquid withdrawal rate (low/high) [low]: ", "s");
                if taxaCol == "high"
                    smSeringa.RPM = 300; %Maximum possible RPM. Corresponds to an actual RPPM = 4.8526 and an withdrawal rate = 1.0607 uL/s
                else
                    smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPPM = 1.4171 and an withdrawal rate = 0.3158 uL/s
                end
                nPassosCol = vol2passo(volCol,passoMot,passoSer); %No of steps corresponding to the volume withdrawn during preparation
                fprintf('Please wait, withdraw of probing liquid... \n');
                move(smSeringa,-floor(nPassosCol)); %Drive motor connected to syringe
                release(smSeringa); %Motor release
                errovolCol = passo2vol(nPassosCol - floor(nPassosCol),passoMot,passoSer); %Error in volume withdrawn during preparation (uL)
                volSerTotal = volSerTotal + volCol; %Update total volume of liquid in the syringe
                fprintf('Withdrawn of probing liquid completed! \n')
            elseif tipoColeta == "R"
                volLib = input("Released volume (uL): ");
                if volLib <= volSerTotal %Checking if there will be enough volume to release
                    taxaLib = input("Liquid released rate (low/high) [low]: ", "s");
                    if taxaLib == "alta"
                        smSeringa.RPM = 300; %Maximum possible RPM. Corresponds to an actual RPPM = 4.8526 and an withdrawal rate = 1.0607 uL/s
                    else
                        smSeringa.RPM = 2; %Minimum possible RPM. Corresponds to an actual RPPM = 1.4171 and an withdrawal rate = 0.3158 uL/s
                    end
                    nPassosLib = vol2passo(volLib,passoMot,passoSer); %No of steps corresponding to the volume released during preparation
                    fprintf('Please wait, release of probing liquid... \n');
                    move(smSeringa,floor(nPassosLib)); %Drive motor connected to syringe
                    release(smSeringa); %Motor release
                    errovolLib = passo2vol(nPassosLib - floor(nPassosLib),passoMot,passoSer); %Error in volume released during preparation (uL)
                    volSerTotal = volSerTotal - volLib; %Update total volume of liquid in the syringe
                    fprintf('Release of probing liquid completed! \n')
                else
                    fprintf("Unable to complete the procedure. Please decrease the released volume. \n")
                end
            end
        elseif tipoAjuste == "A" %Automatic adjustment
            %Withdrawal of probing liquid
            volAutCol = input("Withdrawal volume (uL): ");
            if (volSerTotal + volAutCol) <= volAutLib %Checking if there will be enough volume for the liquid release step
                fprintf("Unable to complete the procedure. Please increase the withdrawal volume. \n")
            else
                %Syringe filling
                nPassosAutCol = vol2passo(volAutCol,passoMot,passoSer); %No of steps corresponding to the volume withdrawal during preparation
                smSeringa.RPM = 8; %Medium possible RPM. Corresponds to an actual RPM = 3.045 and an withdrawal rate = 0.6715 uL/s
                fprintf('Please wait, withdraw of probing liquid... \n');
                move(smSeringa,-floor(nPassosAutCol)); %Drive motor connedted to syringe
                release(smSeringa); %Motor release
                errovolCol = passo2vol(nPassosAutCol - floor(nPassosAutCol),passoMot,passoSer); %Error in volume withdrawal during preparation (uL)
                volSerTotal = volSerTotal + volAutCol; %Update total volume of liquid in the syringe
            
                %Release of probing liquid (Removal of possible bubbles and impurities)
                nPassosAutLib = vol2passo(volAutLib,passoMot,passoSer); %No of steps corresponding to the volume released during preparation
                smSeringa.RPM = 8; %Medium possible RPM. Corresponds to an actual RPM = 3.045 and an withdrawal rate = 0.6715 uL/s
                fprintf('PLease wait, release of probing liquid to remove possible bubbles and impurities... \n')
                move(smSeringa,floor(nPassosAutLib)); %Drive motor connedted to syringe
                release(smSeringa); %Motor release
                errovolLib = passo2vol(nPassosAutLib - floor(nPassosAutLib),passoMot,passoSer); %Error in volume released during preparation (uL)
                volSerTotal = volSerTotal - volAutLib; %Updating total volume of liquid in the syringe
                fprintf('Withdraw of probing liquid completed! \n')
            end  
        end
    else
        break;
    end
end
    fprintf('------------------------------------------------------------------------------ \n')
end

