function Data = circlefit(edgeL,edgeR,baseL,baseR,n_fit)
%CIRCLEFIT Calculation of contact angle using circle fitting in both drop sides
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function extracted from Dropen_V01.m code 
% (https://board.unimib.it/datasets/wzchzbm58p/3) and modified to summarize
% the process of circle fitting. The rotation of the image step before 
% the fitting was removed. The edge detection step was removed. The 
% coordinates of the detected edges and the contact points are given as 
% inputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT:
% edgeL.x and edgeL.y: coordinates of the left side of the edge
% edgeR.x and edgeR.y: coordinates of the right side of the edge
% baseL: coordinates of the left contact point
% baseR: coordinates of the right contact point
% n_fit: number of fitting points (in the original code the value was set to 100)
%   OUTPUT:
% Data.CAL: left contact angle (CA)
% Data.CAR: right contact angle (CA)
% Data.TLL: left contact point (CP)
% Data.TLR: right contact point (CP)
% Data.CoordsL: coordinates of the data in the left side used to fit circle
% Data.CoordsR: coordinates of the data in the right side used to fit circle
% Data.EvalPolyL: estimation of polynomial function for the left side
% Data.EvalPolyR: estimation of polynomial function for the right side

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Load vairables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drop_left = [edgeL(:,1),edgeL(:,2)];
drop_right = [edgeR(:,1),edgeR(:,2)];
CPleft = [baseL(2),baseL(1)];
CPright = [baseR(2),baseR(1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Edge identification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%EDGE_canny = edge(I0gray, 'canny');
EDGE_zerocross = edge(I0gray, 'zerocross');
% alternatively, log filter can be used
% EDGE_log = edge(CI,'log',0);
EDGE = bwareaopen(EDGE_zerocross, round(size(I0gray,1)/2)); %removes small objects
EDGEa = bwmorph(EDGE,  'bridge', Inf);
EDGEb = bwmorph(EDGEa, 'close', Inf);
EDGEc = bwmorph(EDGEb, 'thin', Inf);
% figure (5)
% clf(5)
 %imshow(EDGEc);
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Identification of the longest edge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%boundaries = bwboundaries(EDGE);
[mat_boundaries, label_boundaries] = bwboundaries(EDGE);
for ii=1:size(mat_boundaries,1)
    mmx(ii)=size (mat_boundaries{ii,1},1);
end
mmx=mmx.';
BL = find (mmx==max(mmx));
boundaries = mat_boundaries;
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Identification of left and right drop boundary among all boundaries
% drop boundaries are identified as the lines, closest to CPleft and
% CPright
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Identification of the boundaries that cross the baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
xmin = size(EDGE,2);
xmax = 0;
j = BL;
element = boundaries{j};
pixel_baseline = find(element(:,1) == baseline);
if isempty(pixel_baseline) == 0 % if the line has points on the baseline
    x_baselinel = min(element(pixel_baseline,2));
    if x_baselinel < xmin
        drop_left_index = j;
        xmin = x_baselinel;
    end
    x_baseliner = max(element(pixel_baseline,2));
    if x_baseliner > xmax
        drop_right_index = j;
        xmax = x_baseliner;
    end
end
%}  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Drop left analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%drop_left = [edgeL.x,edgeL.y];
if isempty(drop_left) == 0
    posCPleft = find_index(drop_left,CPright,CPleft); %find coordinate closest to baseline
    %[x0L,y0L]
    %posCPleft = find_index(drop_left(:,1),[y0R,x0R],[y0L,x0L]); %find coordinate closest to baseline
    %CPleft = drop_left(posCPleft,:);
    %CPleft = [y0L,x0L];
    %-- Increased vector drop_left
    drop_left_inc = [drop_left; drop_left; drop_left; drop_left; drop_left];
    posCPleft_inc = posCPleft + 2*size(drop_left,1);
    %--- Check drop reflection
    %{
    if drop_left_inc(posCPleft_inc + 3,1) < drop_left_inc(posCPleft_inc,1) %For drops without reflection
        drop_left_fit  = drop_left_inc(posCPleft_inc:posCPleft_inc+n_fit-1,:);
    else
        drop_left_fit  = drop_left_inc(posCPleft_inc:-1:posCPleft_inc-n_fit+1,:);
    end
    %}
    drop_left_fit  = drop_left_inc(posCPleft_inc:-1:posCPleft_inc-n_fit+1,:);
    %-- A coordinate change is performed, using CP as center of local coordinate systems
    x_fitl = drop_left_fit(:,2) - drop_left_fit(1,2);
    y_fitl = - drop_left_fit(:,1) + drop_left_fit(1,1);

    %-- Circle fitting
    [xcl,ycl,Rl] = circlefitcon(x_fitl, y_fitl);
    dydx_l = - xcl / ycl;
    atan_circ_l = atan(dydx_l) * 180/pi;
    %-- Calculation of left CA
    if atan_circ_l > 0
        CA_circ_l =  atan_circ_l;
    else
        CA_circ_l = 180 + atan_circ_l;
    end
    
    %-- Estimation of the data using the fitted circle
    %circle_l = [CPleft(2)+xcl,CPleft(1)-ycl,Rl]; %circle fitting parameters for the left drop profile 
    circle_l = [xcl,ycl,Rl];
    %--- Evaluation only in the region where the circle was fitted
    %{
    %xEvalCircL = drop_left_fit(:,2);
    yEvalCircL = sqrt(circle_l(3)^2-(x_fitl-circle_l(1)).^2)+circle_l(2);
    %EvalCircL = [xEvalCircL,yEvalCircL];
    EvalCircL = [drop_left_fit(:,2),(yEvalCircL- drop_left_fit(1,1))*-1];
    %xEvalCircL = x_fitl;
    %yEvalCircL = sqrt(circle_l(3)^2-(xEvalCircL-circle_l(1)).^2)+circle_l(2);
    %EvalCircL = [drop_left_fit(:,2),(yEvalCircL- drop_left_fit(1,1))*-1];
    %}
    %--- Evaluation in all the region of the left edge
    %%{
    x_fitlc = drop_left(posCPleft:-1:1,2) - drop_left(posCPleft,2);
    yEvalCircLc = sqrt(circle_l(3)^2-(x_fitlc-circle_l(1)).^2)+circle_l(2);
    EvalCircL = [drop_left(posCPleft:-1:1,2),(yEvalCircLc- drop_left(posCPleft,1))*-1];
    %}
end
if  CPleft(1) == 0
    drop_left_fit = [];
    CA_circ_l = NaN;
    %CA_pol_l = NaN;
    CPleft = NaN;
    EvalCircL = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Drop right analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%drop_right = [edgeR.x,edgeR.y];
if isempty(drop_right) == 0
    posCPright = find_index(drop_right,CPright,CPleft); %find coordinate closest to baseline
    %CPright = drop_right(posCPright,:);
    %CPright = [y0L,x0L];
    %-- Increased vector drop_right (determine points to be used in polynomial fit)
    drop_right_inc = [drop_right; drop_right; drop_right; drop_right; drop_right];
    posCPright_inc = posCPright + 2* size(drop_right,1);
    %{
    %--- Check drop reflection
    if drop_right_inc(posCPright_inc + 3,1) < drop_right_inc(posCPright_inc,1)
        drop_right_fit  = drop_right_inc(posCPright_inc:posCPright_inc+n_fit-1,:);
    else
        drop_right_fit  = drop_right_inc(posCPright_inc:-1:posCPright_inc-n_fit+1,:);
    end
    %}
    drop_right_fit  = drop_right_inc(posCPright_inc:-1:posCPright_inc-n_fit+1,:);
    %-- A coordinate change is performed, using CP as center of local coordinate systems
    x_fitr = drop_right_fit(:,2) - drop_right_fit(1,2);
    y_fitr = - drop_right_fit(:,1) + drop_right_fit(1,1);
    
    %--- Circle fitting
    [xcr,ycr,Rr] = circlefitcon(x_fitr, y_fitr);
    dydx_r = - xcr / ycr;
    atan_circ_r = atan(dydx_r) * 180/pi;
    %-- Calculation of left CA
    if atan_circ_r > 0
        CA_circ_r = 180 - atan_circ_r;
    else
        CA_circ_r = - atan_circ_r;
    end
    
    %-- Estimation of the data using the fitted circle
    %circle_r = [CPright(2)+xcr,CPright(1)-ycr,Rr]; %circle fitting parameters for the right drop profile
    circle_r = [xcr,ycr,Rr];
    %--- Evaluation only in the region where the circle was fitted
    %{
    %xEvalCircR = drop_right_fit(:,2);
    %yEvalCircR = sqrt(circle_r(3)^2-(xEvalCircR-circle_r(1)).^2)+circle_r(2);
    yEvalCircR = sqrt(circle_r(3)^2-(x_fitr-circle_r(1)).^2)+circle_r(2);
    %EvalCircR = [xEvalCircR,yEvalCircR];
    EvalCircR = [drop_right_fit(:,2),(yEvalCircR- drop_right_fit(1,1))*-1];
    %}
    %--- Evaluation in all the region of the right edge
    %%{
    x_fitrc = drop_right(posCPright:-1:1,2) - drop_right(posCPright,2);
    yEvalCircRc = sqrt(circle_r(3)^2-(x_fitrc-circle_r(1)).^2)+circle_r(2);
    EvalCircR = [drop_right(posCPright:-1:1,2),(yEvalCircRc- drop_right(posCPright,1))*-1];
    %}
end
if  CPright(1) == 0
    drop_right_fit = [];
    CA_circ_r = NaN;
    %CA_pol_r = NaN;
    CPright = NaN;
    EvalCircR = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Gathering data and results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data.CAL = CA_circ_l;
Data.CAR = CA_circ_r;
Data.TLL = CPleft;
Data.TLR = CPright;
Data.CoordsL = drop_left_fit;
Data.CoordsR = drop_right_fit;
Data.EvalCircL = EvalCircL;
Data.EvalCircR = EvalCircR;
end

