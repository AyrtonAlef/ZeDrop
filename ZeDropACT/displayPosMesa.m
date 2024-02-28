function displayPosMesa(posMesaX,posMesaY)
%DISPLAYVOLMESA Displays the current stage position (solid substrate).
%   INPUT
% posMesaX - Current X axis position of stage in mm
% posMesaY -  Current Y axis position of stage in mm

% Current stage position
fprintf("Current stage position: X = %.2f mm Y = %.2f mm \n\n",posMesaX,posMesaY);
% Limits of stage positioning
fprintf("                                ^ Xmax = +8.0 mm \n");
fprintf("                                | \n");
fprintf("            Ymin = -24.0 mm <- Pos -> Ymax = +24.0 mm \n");
fprintf("                                | \n");
fprintf("                                v Xmin = -8 0 mm \n\n")
end