function [b,alphac,alphamL,alphamR,x_CPL,y_CPL,x_CPR,y_CPR,indexCPL,indexCPR,Contact_angle_left,Contact_angle_right,Contact_angle,TiltAngle_S] ...
    = CircleMask2(ID,I_New,maskSize,AXxL,AXyL,AXxR,AXyR,Apex_x,Apex_y,Left_edgeTx,Left_edgeDx,Right_edgeTx,Right_edgeDx,...
Num_Down_L, Num_Up_L,Num_Down_R,Num_Up_R,End_x_L,End_x_R)
%CIRCLEMASK Automatic identification of contact points and calculation of
%contact angles based on a mask method (circle mask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function extracted from Dropen_V01.m code 
% (https://board.unimib.it/datasets/wzchzbm58p/3) and modified to analyze
% left and right drop profile edges separately. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT:
% ID: image
% b: maskSize
% AXx and AXy:  x coordinates and y coordinates of the drop profile + solid subatrate (only one boundary line in each side of drop) 
% Apex_x and Apex_y: z and y coordinates of the drop apex
% I_New: binary image generated using the border data
% Left_edgeTx/Left_edgeDx/Right_edgeTx/Right_edgeDx: coordinates of left and right extreme edge points (substrate)
% Num_Up_L and Num_Down_L: number of data upper and down of left extreme edge point
% Num_Up_R and Num_Down_R: number of data upper and down of right extreme edge point
% End_x_l: x coordinates of the last drop point on the left
% End_x_R: x coordinates of the last drop point on the right
%   OUTPUT:
% b: matrix size
% alphac: calculated angle in 2D system (for all position of the 2D image)
% alphamL and alphamR: gathered data from the mask method (L - from left drop profile and R - from right drop profile)
%   alphamL(1,:) and alphamR(1,:): x coordinates of the detected edge profile (drop profile + solid substrate)
%   alphamL(2,:) and alphamR(2,:): y coordinates of the detected edge profile (drop profile + solid substrate)
%   alphamL(3,:) and alphamR(3,:): length of arc (from apex)
%   alphamL(4,:) and alphamR(4,:): convolution area
%   alphamL(5,:) and alphamR(5,:): local slope
% x_CPL,y_CPL,x_CPR and y_CPR: coordinates of the left and right contact points, respectively.
% indexCPL: index of alphamL related to the data on the contact point
% indexCPR: index of alphamR related to the data on the contact point
% Contact_angle_left: left CA
% Contact_angle_right: right CA
% Contact_angle: mean CA between left and right CAs
% TiltAngle_S: baseline tilt angle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = maskSize; % Neighbourhood index. b = 20 -> optimum value found for CA evaluation in the original work of Akbari and Antonini 2021 (doi: 10.17632/wzchzbm58p.3)
m = 2*b + 1; % Mask size = 2b + 1
R = b + 1;
Mask = zeros(m,m);% Mask, ones up
Ar = zeros(size(ID,1),size(ID,2)); % Convolution area
alphac = zeros(size(ID,1),size(ID,2)); % alphac = Calculated angle in 2D system
IA = zeros(m,m); % Locale area matrix
lambda = 1; % Sign(slope)
MP = zeros(m,m); % Mask, ones up
MM = zeros(m,m); % Mask, ones down
MS = 0; % Mask size
v = 0;
vv = 0;
vx = 0;
%alpham = Calculated angle data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Generation of masks (characterization mask)(circle mask - cMask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:m
    for jj = 1:m
        condition = ((ii - R)^2+ (jj - R)^2)/(b^2);
        if condition <=1  && ii<=b
            Mask(ii,jj) = 1;
        end
        if condition <=1  &&  ii>(b+1)
            Mask(ii,jj) = -1;
        end
    end
end
ArTotal1 = nnz(Mask);
 
%- Generation of masks (angle calculation masks)(MP = TMask)(MM = BMask)
for ii = 1:m
    for jj = 1:m
        condition = ((ii - (b+1))^2+ (jj - (b+1))^2)/(b^2);
        if condition <=1  && jj <=(b+1)
            MP(jj,ii) = 1;
            MS=MS+1;
        end
        if condition <=1  && jj >=(b+1)
            MM(jj,ii) = 1;
        end
    end
end
ArTotal = 2*MS-2*b-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Calculation of the convolution area and tangent angle in each point of border (using the circle mask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- Left side
up = 0;
down = 0;
m_up = 0;
m_down = 0;
for kk = 1:max(size(AXxL))
    rr = AXxL(kk);
    cc = AXyL(kk);
    if (rr > b) && (rr < size(ID,1)-b) && (cc > b) && (cc < size(ID,2)-b)
        IA = I_New(rr-b:rr+b,cc-b+1:cc+b+1);
        IA1 = IA.*Mask;
        IA3 = sum(IA1,'all'); %baseline identification factor
        if (rr <= Left_edgeTx) %upper side
            up = up+1;
            IA2 = IA.*MP;
            Ar(rr,cc) = sum(IA2,'all');
        end
        if (rr > Left_edgeDx) && (rr > Left_edgeTx) %down side
            down = down+1;
            IA2 = IA.*MM;
            Ar(rr,cc) = sum(IA2,'all');
        end
        if (IA3 >= 0) && ((rr > Left_edgeTx) && (rr <= Left_edgeDx)) %down part of middle
            IA2 = IA.*MM;
            Ar(rr,cc) = (sum(IA2,'all'));
            m_down = m_down+1;
        end
        if  (IA3 < 0) && ((rr > Left_edgeTx) && (rr <= Left_edgeDx) && (Num_Down_L > Num_Up_L*0.5 )) %top part of middle
            IA2 = IA.*MP;
            Ar(rr,cc) = (sum(IA2,'all'));
            m_up = m_up+1;
        end
        alphac(rr,cc) = 180-(180/pi)*(2*pi* Ar(rr,cc)/(ArTotal)); %Calculation of local slope
        if  alphac(rr,cc)> 171
            alphac(rr,cc) = alphac(rr-5,cc);
        end
    end
end
%-- Right side
up = 0;
down = 0;
m_up = 0;
m_down = 0;
for kk = 1:max(size(AXxR))
    rr = AXxR(kk);
    cc = AXyR(kk);
    if (rr > b) && (rr < size(ID,1)-b) && (cc > b) && (cc < size(ID,2)-b)
        IA = I_New(rr-b:rr+b,cc-b-1:cc+b-1);
        IA1 = IA.*Mask;
        IA3 = sum(IA1,'all'); %baseline identification factor
        if (rr <= Right_edgeTx) %upper side
            up = up+1;
            IA2 = IA.*MP;
            Ar(rr,cc) = sum(IA2,'all');
        end
        if (rr > Right_edgeDx) && (rr > Right_edgeTx)  %Down side
            down = down+1;
            IA2 = IA.*MM;
            Ar(rr,cc) = sum(IA2,'all');
        end
        if (IA3 >= 0) && ((rr> Right_edgeTx) && (rr<= Right_edgeDx))  %down part of middle
            IA2 = IA.*MM;
            Ar(rr,cc) = (sum(IA2,'all'));
            m_down = m_down+1;
        end
        if  (IA3 < 0) && ((rr > Right_edgeTx) &&  (rr <= Right_edgeDx) && (Num_Down_R > Num_Up_R*0.5)) %top part of middle
            IA2 = IA.*MP;
            Ar(rr,cc) = (sum(IA2,'all'));
            m_up = m_up+1;
        end
        alphac(rr,cc) = 180-(180/pi)*(2*pi* Ar(rr,cc)/(ArTotal)); %Calculation of local slope
        if  alphac(rr,cc)> 171
            alphac(rr,cc) = alphac(rr-5,cc);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Gathering of the area and angle data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- Left side
ss = 0;
for kk = 1:max(size(AXxL))
    rr = AXxL(kk);
    cc = AXyL(kk);
    ss = ss + 1;
    lambda = -1;
    if (rr <= size(ID,1)) && (cc <= size(ID,2)) && (Ar(rr,cc) ~= 0)
        alphamL(1,ss) = rr;
        alphamL(2,ss) = cc;
        alphamL(3,ss) = lambda*sqrt((rr-Apex_x)^2+(cc-Apex_y)^2);% length of arc (from the apex)
        alphamL(4,ss) = Ar(rr,cc);
        alphamL(5,ss) = alphac(rr,cc);
    end
end
%-- Right side
ss = 0;
for kk = 1:max(size(AXxR))
    rr = AXxR(kk);
    cc = AXyR(kk);
    ss = ss + 1;
    lambda = 1;
    if (rr <= size(ID,1)) && (cc <= size(ID,2)) && (Ar(rr,cc) ~= 0)
        alphamR(1,ss) = rr;
        alphamR(2,ss) = cc;
        alphamR(3,ss) = lambda*sqrt((rr-Apex_x)^2+(cc-Apex_y)^2);% length of arc (from the apex)
        alphamR(4,ss) = Ar(rr,cc);
        alphamR(5,ss) = alphac(rr,cc);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Identifiction of baseline and left and right contact points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_MioL = max(max(alphac(1:End_x_L-3,1:Apex_y)));
M_MioR = max(max(alphac(1:End_x_R-3,Apex_y:size(ID,2))));

colXR = 0;
[rowXL,colXL] = find(round(alphac(:,1:Apex_y),1) == round(M_MioL,1));
[rowXR,colXR] = find(round(alphac(:,Apex_y:size(ID,2)),1) == round(M_MioR,1));
%[rowXL,colXL]=find(floor(alphac(:,1:Apex_y))==floor(M_MioL));
%[rowXR,colXR]=find(floor(alphac(:,Apex_y:size(ID,2)))==floor(M_MioR));
colXR = colXR + Apex_y-1;

x_CPL = min(rowXL);
if M_MioL<89
    M_MioL = max(max(alphac(1:Left_edgeTx,1:Apex_y)));
    [rowXL,colXL] = find(round(alphac(:,1:Apex_y),1)==round(M_MioL,1));
    for ii = 1:size(rowXL)
        drL(ii) = rowXL(ii)-Left_edgeTx;
    end
    mindrL = min(abs(drL(:)));
    for ii = 1:size(rowXL)
        if drL(ii) == -mindrL
            x_CPL = rowXL(ii);
        end
    end
end
y_CPL = colXL(find(rowXL==x_CPL));

x_CPR = min(rowXR);
if M_MioR<89
    colXR = 0;
    M_MioR = max(max(alphac(1:Right_edgeTx,Apex_y:size(ID,2))));
    [rowXR,colXR] = find(round(alphac(:,Apex_y:size(ID,2)),1)==round(M_MioR,1));
    colXR = colXR+Apex_y-1;
    for ii = 1:size(rowXR)
        drR(ii) = rowXR(ii)-Right_edgeTx;
    end
    mindrR = min(abs(drR(:)));
    for ii = 1:size(rowXR)
        if drR(ii) == -mindrR
            x_CPR = rowXR(ii);
        end
    end
end
y_CPR = colXR(find(rowXR==x_CPR));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Calculation of the contact point index in alphamL and alphamR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
distL = sqrt((alphamL(1,:)-x_CPL).^2 + (alphamL(2,:)-y_CPL).^2);
indexCPL = find(distL == min(distL));
distR = sqrt((alphamR(1,:)-x_CPR).^2 + (alphamR(2,:)-y_CPR).^2);
indexCPR = find(distR == min(distR));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Calculation of baseline tilt angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TiltAngle_S = 0;
if x_CPL ~= x_CPR
   TiltAngle_S = atan2(x_CPR-x_CPL, y_CPR-y_CPL) * 180/pi;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Calculation of contact angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Contact_angle_left = round(alphac(x_CPL,y_CPL)+TiltAngle_S,2);
Contact_angle_right = round(alphac(x_CPR,y_CPR)-TiltAngle_S,2);
Contact_angle = (Contact_angle_left+Contact_angle_right)/2;
end