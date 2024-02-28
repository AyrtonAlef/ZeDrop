function displayVolSer(volSerTot)
%DISPLAYVOLSER Displays the current volume of probing liquid in the syringe
%   INPUT
% volSerTor - Liquid volume in the syringe [uL]

volTotal = 1000; %Total liquid volume in the syringe in uL
percVol = (volSerTot/volTotal); %Percetage of syringe volume occupied by probing liquid
fprintf("Liquid volume in the syringe: \n");
if percVol == 0
    fprintf('    --------------------------------------------- \n');
    fprintf('---|                                             |-----|   %.1f/%d uL \n', volSerTot,volTotal);
    fprintf('    --------------------------------------------- \n\n');
elseif (percVol > 0) && (percVol <=0.05)
    fprintf('    --------------------------------------------- \n');
    fprintf('---||                                            |-----|   %.1f/%d uL \n', volSerTot,volTotal);
    fprintf('    --------------------------------------------- \n\n');
elseif (percVol > 0.02) && (percVol <=0.05)
    fprintf('    --------------------------------------------- \n');
    fprintf('---||||                                          |-----|   %.1f/%d uL \n', volSerTot,volTotal);
    fprintf('    --------------------------------------------- \n\n');
elseif (percVol > 0.05) && (percVol <=0.1)
    fprintf('    --------------------------------------------- \n');
    fprintf('---||||||                                        |-----|   %.1f/%d uL \n', volSerTot,volTotal);
    fprintf('    --------------------------------------------- \n\n');
elseif (percVol > 0.1) && (percVol <=0.15)
    fprintf('    --------------------------------------------- \n');
    fprintf('---||||||||                                      |-----|   %.1f/%d uL \n', volSerTot,volTotal);
    fprintf('    --------------------------------------------- \n\n');
elseif (percVol > 0.15) && (percVol <=0.2)
    fprintf('    --------------------------------------------- \n');
    fprintf('---||||||||||                                    |-----|   %.1f/%d uL \n', volSerTot,volTotal);
    fprintf('    --------------------------------------------- \n\n');
elseif (percVol > 0.2) && (percVol <=0.25)
    fprintf('    --------------------------------------------- \n');
    fprintf('---|||||||||||||                                 |-----|   %.1f/%d uL \n', volSerTot,volTotal);
    fprintf('    --------------------------------------------- \n\n');
elseif (percVol > 0.25) && (percVol <=0.3)
    fprintf('    --------------------------------------------- \n');
    fprintf('---|||||||||||||||                               |-----|   %.1f/%d uL \n', volSerTot,volTotal);
    fprintf('    --------------------------------------------- \n\n');
elseif (percVol > 0.3) && (percVol <=0.35)
    fprintf('    --------------------------------------------- \n');
    fprintf('---|||||||||||||||||                             |-----|   %.1f/%d uL \n', volSerTot,volTotal);
    fprintf('    --------------------------------------------- \n\n');
elseif (percVol > 0.35) && (percVol <=0.4)
    fprintf('    --------------------------------------------- \n');
    fprintf('---|||||||||||||||||||                           |-----|   %.1f/%d uL \n', volSerTot,volTotal);
    fprintf('    --------------------------------------------- \n\n');
elseif (percVol > 0.4) && (percVol <=0.5)
    fprintf('    --------------------------------------------- \n');
    fprintf('---||||||||||||||||||||||||                      |-----|   %.1f/%d uL \n', volSerTot,volTotal);
    fprintf('    --------------------------------------------- \n\n');
elseif (percVol > 0.5) && (percVol <=0.6)
    fprintf('    --------------------------------------------- \n');
    fprintf('---||||||||||||||||||||||||||||                  |-----|   %.1f/%d uL \n', volSerTot,volTotal);
    fprintf('    --------------------------------------------- \n\n');
elseif (percVol > 0.6) && (percVol <=0.7)
    fprintf('    --------------------------------------------- \n');
    fprintf('---|||||||||||||||||||||||||||||||||             |-----|   %.1f/%d uL \n', volSerTot,volTotal);
    fprintf('    --------------------------------------------- \n\n');
elseif (percVol > 0.7) && (percVol <=0.8)
    fprintf('    --------------------------------------------- \n');
    fprintf('---|||||||||||||||||||||||||||||||||||||         |-----|   %.1f/%d uL \n', volSerTot,volTotal);
    fprintf('    --------------------------------------------- \n\n');
elseif (percVol > 0.8) && (percVol <=0.9)
    fprintf('    --------------------------------------------- \n');
    fprintf('---||||||||||||||||||||||||||||||||||||||||||    |-----|   %.1f/%d uL \n', volSerTot,volTotal);
    fprintf('    --------------------------------------------- \n\n');
elseif (percVol > 0.9) && (percVol <=0.95)
    fprintf('    --------------------------------------------- \n');
    fprintf('---||||||||||||||||||||||||||||||||||||||||||||| |-----|   %.1f/%d uL \n', volSerTot,volTotal);
    fprintf('    --------------------------------------------- \n\n');
elseif (percVol > 0.95) && (percVol <=1)
    fprintf('    --------------------------------------------- \n');
    fprintf('---|||||||||||||||||||||||||||||||||||||||||||||||-----|   %.1f/%d uL \n', volSerTot,volTotal);
    fprintf('    --------------------------------------------- \n\n');
end

end