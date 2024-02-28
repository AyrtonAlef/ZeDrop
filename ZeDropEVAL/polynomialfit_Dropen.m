function Data = polynomialfit_Dropen (edgeL,edgeR,baseL,baseR,n_fit,poly_ord)
%POLYNOMIALFIT_DROPEN Calculation of contact angle using polynomial fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function extracted from Dropen_V01.m code 
% (https://board.unimib.it/datasets/wzchzbm58p/3) and modified to summarize
% the process of polynomial fitting. The rotation of the image step before 
% the fitting was removed. The edge detection step was removed. The 
% coordinates of the detected edges and the contact points are given as 
% inputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT:
% edgeL.x and edgeL.y: coordinates of the left side of the edge
% edgeR.x and edgeR.y: coordinates of the right side of the edge
% poly_ord: polynomial order (the authors recommend: (i) for CA < 60° -> poly_ord = 2 and (ii) for CA > 60° -> poly_ord = 3)
% n_fit: number of fitting points (in the original code the value was set to 100)
% baseL: coordinates of the left contact point
% baseR: coordinates of the right contact point
%   OUTPUT:
% Data.CAL: left contact angle (CA)
% Data.CAR: right contact angle (CA)
% Data.TLL: left contact point (CP)
% Data.TLR: right contact point (CP)
% Data.CoordsL: coordinates of the data in the left side used to fit polynomial
% Data.CoordsR: coordinates of the data in the right side used to fit polynomial
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
    %CPleft = drop_left(posCPleft,:);
    %CPleft = [y0L,x0L];
    %-- Increased vector drop_left
    drop_left_inc = [drop_left; drop_left; drop_left; drop_left; drop_left];
    posCPleft_inc = posCPleft + 2*size(drop_left,1);
    %{
    %--- Check drop reflection
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

    %-- Polynomial fitting
    coeff = polyfitcon(y_fitl, x_fitl, poly_ord, 0);
    coeff_der = coeff(1:poly_ord) .* (poly_ord:-1:1)';
    %-- Calculation of left CA
    %--- Evaluation of the CAl at CPl, with local coordinates (0,0)
    CA_pol_l = (pi/2 - atan(coeff_der(end))) *180/pi;
    
    %-- Estimation of the data using the fitted polynomial function
    %--- Evaluation only in the region where the polynomial was fitted
    %{
    EvalPolyL = [(polyval(coeff,y_fitl) + drop_left_fit(1,2)),drop_left_fit(:,1)];
    %}
    %--- Evaluation in all the region of the left edge
    %%{
    y_fitlc = -drop_left(posCPleft:-1:1,1) + drop_left(posCPleft,1);
    EvalPolyL = [(polyval(coeff,y_fitlc) + drop_left(posCPleft,2)),drop_left(posCPleft:-1:1,1)];
    %}
end
if  CPleft(1) == 0
    drop_left_fit = [];
    %CA_circ_l = NaN;
    CA_pol_l = NaN;
    CPleft = NaN;
    EvalPolyL = [];
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

    %-- Fit the data to polynomial
    %coeff = polyfitcon(x_fitr, y_fitr, poly_ord, 0);
    coeff = polyfitcon(y_fitr, x_fitr, poly_ord, 0);
    coeff_der = coeff(1:poly_ord) .* (poly_ord:-1:1)';
    %-- Calculation of right CA
    %--- Evaluation of the CAr at CPr, with local coordinates (0,0)
    CA_pol_r = (pi/2 + atan(coeff_der(end))) *180/pi;
    
    %-- Estimation of the data using the fitted polynomial function
    %--- Evaluation only in the region where the polynomial was fitted
    %{
    EvalPolyR = [(polyval(coeff,y_fitr) + drop_right_fit(1,2)),drop_right_fit(:,1)];
    %}
    %--- Evaluation in all the region of the right edge
    %%{
    y_fitrc = -drop_right(posCPright:-1:1,1) + drop_right(posCPright,1);
    EvalPolyR = [(polyval(coeff,y_fitrc) + drop_right(posCPright,2)),drop_right(posCPright:-1:1,1)];
    %}
end
if  CPright(1) == 0
    drop_right_fit = [];
    %CA_circ_r = NaN;
    CA_pol_r = NaN;
    CPright = NaN;
    EvalPolyR = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Gathering data and results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data.CAL = CA_pol_l;
Data.CAR = CA_pol_r;
Data.TLL = CPleft;
Data.TLR = CPright;
Data.CoordsL = drop_left_fit;
Data.CoordsR = drop_right_fit;
Data.EvalPolyL = EvalPolyL;
Data.EvalPolyR = EvalPolyR;
end

