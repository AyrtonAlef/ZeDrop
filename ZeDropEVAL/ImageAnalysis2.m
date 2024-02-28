function [I_New,AXxL,AXyL,AXxR,AXyR,Apex_x,Apex_y,Left_edgeTx,Left_edgeDx,Right_edgeTx,Right_edgeDx,Num_Down_L, Num_Up_L,Num_Down_R,...
    Num_Up_R,End_x_L,End_x_R] = ImageAnalysis2(ID,edgeL,edgeR,fillBlankHorLines,marginMaskBaselineDet)
% IMAGEANALYSIS Automatic detection of the drop profile and identification of left and right contact points.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function extracted from Dropen_V01.m code 
% (https://board.unimib.it/datasets/wzchzbm58p/3) and modified to receive 
% edge information (longestedge) from a subpixel edge detection routine.
% Determination of left and right edge properties separately. Addition of
% filling and filtering steps during the determination of edge profile
% informations. Modified by Ayrton Pereira.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT:
% ID: image
% edgeL.x and edgeL.y: coordinates of the left side of the edge
% edgeR.x and edgeR.y: coordinates of the right side of the edge
% fillBlankHorLines: Option for filling blank horizontal lines (1 = enable and 0 = disable)
% marginMaskBaselineDet: minimum distance between edge points to do not consider an abrupt change
%   OUTPUT:
% AXxL and AXyL:  x coordinates and y coordinates of the left side of the drop profile
% AXxR and AXyR:  x coordinates and y coordinates of the right side of the drop profile
% Apex_x and Apex_y: z and y coordinates of the drop apex
% I_New: binary image generated using the border data
% Left_edgeTx/Lef_edgeDx/Right_edgeTx/Right_edgeDx: coordinates of left and right extreme edge points
% Num_Up_L and Num_Down_L: number of data upper and down of left extreme edge point
% Num_Up_R and Num_Down_R: number of data upper and down of right extreme edge point
% End_x_L and End_y_L: coordinates of the last drop point on the left
% End_x_R and End_y_R: coordinates of the last drop point on the right

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Load coordinates of the left and right detected edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AXx1L = round(edgeL.y);
AXy1L = round(edgeL.x);
AXx1R = round(edgeR.y);
AXy1R = round(edgeR.x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Find drop apex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Apex_x = round(mean(AXx1L(1),AXx1R(1)));
Apex_y = round(mean(AXy1L(1),AXy1R(1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Find the lowest point in the detected edge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Low_xL = max(AXx1L);
Low_yL = round(mean(AXy1L(find(AXx1L == Low_xL))));
Low_xR = max(AXx1R);
Low_yR = round(mean(AXy1R(find(AXx1R == Low_xR))));
Low_x = round(mean(Low_xL,Low_xR));
Low_y = round(mean(Low_yL,Low_yR));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Alteration of boundary to have only one boundary line in each side of drop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- Left Side
ll = 0;
AXxL = 0;
AXyL = 0;
for kk = 1: max(size(AXx1L))
    rr = AXx1L(kk);
    cc = AXy1L(kk);
    if not(ismember(rr,AXxL))
        cc1 = max(AXy1L(find(AXx1L == rr)));
        ll = ll + 1;
        AXxL(ll,1) = rr;
        AXyL(ll,1) = cc1;
    end
end
%-- Right Side
ll = 0;
AXxR = 0;
AXyR = 0;
for kk = 1: max(size(AXx1R))
    rr = AXx1R(kk);
    cc = AXy1R(kk);
    if not(ismember(rr,AXxR))
        cc1 = min(AXy1R(find(AXx1R == rr)));
        ll = ll + 1;
        AXxR(ll,1) = rr;
        AXyR(ll,1) = cc1;
    end
end
%%{
if fillBlankHorLines
    %-- Filling blank horizontal lines along the edge
    %--- Left Side
    AXxLFill = AXxL(1);
    AXyLFill = AXyL(1);
    jj = 1;
    for ii = 2:size(AXxL)
        if (AXxL(ii)-AXxL(ii-1)) == 1
            jj = jj + 1;
            AXxLFill(jj) = AXxL(ii);
            AXyLFill(jj) = AXyL(ii);
        else
            %diff = (AXx(ii)-AXx(ii-1));
            while (AXxL(ii)-AXxLFill(jj)) > 1
                jj = jj + 1;
                AXxLFill(jj) = AXxLFill(jj-1) + 1;
                AXyLFill(jj) = round((((AXyL(ii)-AXyL(ii-1))*(AXxLFill(jj)-AXxL(ii-1)))/(AXxL(ii)-AXxL(ii-1))) + AXyL(ii-1)); %Estimation of y coordinate based on line regression
            end
            jj = jj + 1;
            AXxLFill(jj) = AXxL(ii);
            AXyLFill(jj) = AXyL(ii);
        end
    end
    AXxL = AXxLFill';
    AXyL = AXyLFill';
    %--- Right Side
    AXxRFill = AXxR(1);
    AXyRFill = AXyR(1);
    jj = 1;
    for ii = 2:size(AXxR)
        if (AXxR(ii)-AXxR(ii-1)) == 1
            jj = jj + 1;
            AXxRFill(jj) = AXxR(ii);
            AXyRFill(jj) = AXyR(ii);
        else
            %diff = (AXx(ii)-AXx(ii-1));
            while (AXxR(ii)-AXxRFill(jj)) > 1
                jj = jj + 1;
                AXxRFill(jj) = AXxRFill(jj-1) + 1;
                AXyRFill(jj) = round((((AXyR(ii)-AXyR(ii-1))*(AXxRFill(jj)-AXxR(ii-1)))/(AXxR(ii)-AXxR(ii-1))) + AXyR(ii-1)); %Estimation of y coordinate based on line regression
            end
            jj = jj + 1;
            AXxRFill(jj) = AXxR(ii);
            AXyRFill(jj) = AXyR(ii);
        end
    end
    AXxR = AXxRFill';
    AXyR = AXyRFill';
end
%}
%%{
%- Filtering abrupt changes in the edge profile
%margin = 50; %Minimum distance between edge points to do not consider an abrupt change
%-- Left Side
for ii = 5:size(AXxL)
    if abs(AXyL(ii)-AXyL(ii-1)) > marginMaskBaselineDet
        AXyL(ii) = AXyL(ii-1);
    end
end
%-- Right side
for ii = 5:size(AXxR)
    if abs(AXyR(ii)-AXyR(ii-1)) > marginMaskBaselineDet
        AXyR(ii) = AXyR(ii-1);
    end
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Search the final data in vertical axis (=x) the left and right side of drop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
End_x_L = max(AXxL(find(AXyL < Apex_y)));
End_y_L = min(AXyL(find(AXxL == End_x_L)));
End_x_R = max(AXxR(find(AXyR > Apex_y)));
End_y_R = max(AXyR(find(AXxR == End_x_R)));
End_Data = max(End_x_R, End_x_L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Generation of new image in black and white using the border data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- Before the drop apex
I_New = zeros(size(ID,1),size(ID,2)); 
for ii = 1:Apex_x
    for jj = 1:size(ID,2)
        I_New(ii,jj) = 1;
    end
end
%-- After the final border data
for ii = End_Data+1:size(ID,1)
    for jj = 1:size(ID,2)
        I_New(ii,jj) = 1;
    end
end
%-- In drop area
%{
for kk = 1:max(size(AXx))
    rr = AXx(kk);
    cc = AXy(kk);
    if cc < Apex_y
        for jj = 1:cc
            I_New(rr,jj) = 1;
        end
    end
    if cc >= Apex_y
        for jj = cc:size(ID,2)
            I_New(rr,jj) = 1;
        end
    end
end
%}
for kk = 1:max(size(AXxL))
    rr = AXxL(kk);
    cc = AXyL(kk);
    %if cc < Apex_y
    for jj = 1:cc
        I_New(rr,jj) = 1;
    end
end
for kk = 1:max(size(AXxR))
    rr = AXxR(kk);
    cc = AXyR(kk);
    %if cc < Apex_y
    for jj = cc:size(ID,2)
        I_New(rr,jj) = 1;
    end
end
%-- Alteration in new intensity matrix to avoid unexpected mistakes in the convolution matrix
if End_x_L < End_Data
    for ii = End_x_L+1:End_Data
        for jj = 1:End_y_L-1
            I_New(ii,jj) = 1;
        end
    end
end
if End_x_R < End_Data
    for ii = End_x_R+1:End_Data
        for jj = End_y_R+1:size(ID,2)
            I_New(ii,jj) = 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Search for the left edges of drop(= minimum of horizontal(=y) data of drop border)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- Top side
Left_edgeTy = min(AXyL);
L_X = AXxL(find(AXyL == Left_edgeTy));
Left_edgeTx = min(L_X);
Right_edgeTy = max(AXyR);
R_X = AXxR(find(AXyR == Right_edgeTy));
Right_edgeTx = min(R_X);
while I_New(Left_edgeTx,Left_edgeTy) == 0
    if I_New(Left_edgeTx,Left_edgeTy) == 0
        Left_edgeTx = Left_edgeTx+1;
    end
end
rr = Left_edgeTx;
cc = Left_edgeTy;
ia = 0;
%-- Finding the centeral data of the edge part
while ismember(rr,L_X)
    if I_New(rr,cc) ~= 0
        ia = ia+1;
        rr = rr+1;
    end
end
Left_edgeTx = abs(round((rr+Left_edgeTx)/2))-1;
%-- Reflection part
Left_edgeDx = max(L_X);
Left_edgeDy = min(AXyL(find(AXxL == Left_edgeDx)));
if Left_edgeDx > Left_edgeTx+20
    while I_New(Left_edgeDx,Left_edgeDy) == 0
        if I_New(Left_edgeDx,Left_edgeDy) == 0
            Left_edgeDx = Left_edgeDx-1;
        end
    end
end
%- Finding the central data of the edge part
rr = Left_edgeDx;
cc = Left_edgeDy;
ia = 0;
while ismember(rr,L_X)
    if I_New(rr,cc) ~= 0
        ia = ia+1;
        rr = rr-1;
    end
end
Left_edgeDx = abs(round((rr+Left_edgeDx)/2))-1;
%this alteration is needed in some hydrophilic drops
if Left_edgeDx < Left_edgeTx
    Left_edgeDx = Left_edgeTx;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Search for the right edges of drop(= maximum of horizontal(=y) data of drop border)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while I_New(Right_edgeTx,Right_edgeTy) == 0
    if I_New(Right_edgeTx,Right_edgeTy) == 0
        Right_edgeTx = Right_edgeTx+1;
    end
end
rr = Right_edgeTx;
cc = Right_edgeTy;
ia = 0;
%-- Finding the centeral data of the edge part
while ismember(rr,R_X)
    if I_New(rr,cc) ~= 0
        ia = ia+1;
        rr = rr+1;
    end
end
Right_edgeTx = abs(round((rr+Right_edgeTx)/2))-1;
%-- Reflection part
Right_edgeDx = max(R_X);
Right_edgeDy = max(AXyR(find(AXxR == Right_edgeDx)));

if Right_edgeDx > Right_edgeTx+20
    while I_New(Right_edgeDx,Right_edgeDy) == 0
        if I_New(Right_edgeDx,Right_edgeDy) == 0
            Right_edgeDx = Right_edgeDx-1;
        end
    end
end
%-- Finding the central data of the edge part
rr = Right_edgeDx;
cc = Right_edgeDy;
ia = 0;
while ismember(rr,R_X)
    if I_New(rr,cc) ~= 0
        ia = ia+1;
        rr = rr-1;
    end
end
Right_edgeDx = abs(round((rr+Right_edgeDx)/2))-1;
%this alteration is needed in some hydrophilic drops
if Right_edgeDx < Right_edgeTx
    Right_edgeDx = Right_edgeTx;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Counting the number of data upper than the top edge=Num_Up and after that=Num_Down
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- Left Side
Num_Up_L = 0;
Num_Down_L = 0;
for kk = 1:max(size(AXxL))
    if AXxL(kk) <= Left_edgeTx
        Num_Up_L = Num_Up_L+1;
    end
    if AXxL(kk) > Left_edgeTx
        Num_Down_L = Num_Down_L+1;
    end
end
%-- Right Side
Num_Up_R = 0;
Num_Down_R = 0;
for kk = 1:max(size(AXxR))
    if AXxR(kk) <= Right_edgeTx
        Num_Up_R = Num_Up_R+1;
    end
    if AXxR(kk) > Right_edgeTx
        Num_Down_R = Num_Down_R+1;
    end
end

SE_L = 0;
if (Left_edgeDx == Left_edgeTx) && (Num_Down_L >= Num_Up_L*0.5) && (Low_yL ~= Apex_y) && ((End_x_L-Apex_x) > 0.6*(Right_edgeTy-Left_edgeTy))
    %&& (max(End_x_L,End_x_R)-Apex_x)<2*(Left_edgeTx-Apex_x)&& (max(End_x_L,End_x_R)-Apex_x)>(Right_edgeTy-Left_edgeTy)
    Left_edgeDx = End_x_L;
    SE_L = SE_L+1;
end
SE_R = 0;
if (Right_edgeDx == Right_edgeTx) && (Num_Down_R >= Num_Up_R*0.5) && (Low_yR ~= Apex_y) && ((End_x_R-Apex_x) > 0.6*(Right_edgeTy-Left_edgeTy))
    %&& (max(End_x_L,End_x_R)-Apex_x)<2*(Left_edgeTx-Apex_x)&& (max(End_x_L,End_x_R)-Apex_x)>(Right_edgeTy-Left_edgeTy)
    Right_edgeDx = End_x_R;
    SE_R = SE_R+1;
end