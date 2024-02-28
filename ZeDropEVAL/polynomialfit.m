function Data = polynomialfit(edgeL,edgeR,baseL,baseR,nv,poly_ord)
% POLYNOMIALFIT Fit contact angle by double ellipses or polynomials
% (https://www.mathworks.com/matlabcentral/fileexchange/57919-drop-shape-analysis-fit-contact-angle-by-double-ellipses-or-polynomials).
% More details can be found in Andersen and Taboryski 2017 work (doi:
% 10.1088/1361-6501/aa5dcf).
% Modified to enable the user to set different numbers of polynomial orders
% for the polynomial fitting function. Set polynomial order as as input
% parameter. Enable run, even fitting problems appear. 
%   INPUT:
% edgeL.x and edgeL.y: coordinates of the left side of the edge
% edgeR.x and edgeR.y: coordinates of the right side of the edge
% edgeL.ny and edgeL.nx: normal vector (normalized) of the left side of the edge
% edgeR.ny and edgeR.nx: normal vector (normalized) of the right side of the edge
% baseL: coordinates of the left contact point
% baseR: coordinates of the right contact point
% nv: number of fitting points
% poly_ord: polynomial order
%   OUTPUT:
% Data.CAL: left contact angle (CA)
% Data.CAR: right contact angle (CA)
% Data.TLL: left contact point (CP)
% Data.TLR: right contact point (CP)
% Data.CoordsL: coordinates of the data in the left side used to fit polynomial
% Data.CoordsR: coordinates of the data in the right side used to fit polynomial
% Data.EvalPolyL: estimation of polynomial function for the left side
% Data.EvalPolyR: estimation of polynomial function for the right side

%---------------------------------------------------------------------------
%                           Define variables from input
%---------------------------------------------------------------------------

left=[edgeL.x,edgeL.y];
right=[edgeR.x,edgeR.y];

x0L=baseL(1);
y0L=baseL(2);

x0R=baseR(1);
y0R=baseR(2);


%---------------------------------------------------------------------------
%                 find coordinate closest to baseline
%---------------------------------------------------------------------------

indexL=find_index(left,baseR,baseL);
indexR=find_index(right,baseR,baseL);

%nv=10; % number of pixels above baseline used to calculate approximate contact angle
if indexL-nv < 0 || indexR-nv < 0
    Data.CAL = [];
    Data.CAR = [];
    Data.TLL = [];
    Data.TLR = [];
    Data.CoordsR = [];
    Data.CoordsL = [];
    Data.EvalPolyL = [];
    Data.EvalPolyR = [];
else
    %---------------------------------------------------------------------------
    %               Determine points to be used in polyfit
    %---------------------------------------------------------------------------

    CAL=270-median(atan2(edgeL.ny(indexL-nv:indexL-1),edgeL.nx(indexL-nv:indexL-1)))*180/pi;
    if CAL>180
        %CAL=CAL-360;
        CAL = abs(CAL-360);
    end
    CAR=90+median(atan2(edgeR.ny(indexR-nv:indexR-1),edgeR.nx(indexR-nv:indexR-1)))*180/pi;

    arclength=60;
    nvL=min([round(indexL*arclength/CAL),indexL-1]);
    nvR=min([round(indexR*arclength/CAR),indexR-1]);


    %---------------------------------------------------------------------------
    %              Rotate and center datapoints to improve fit
    %---------------------------------------------------------------------------

    AnglesL=atan2(edgeL.ny(indexL-nvL:indexL-1),edgeL.nx(indexL-nvL:indexL-1));
    AnglesL(AnglesL<0)=AnglesL(AnglesL<0)+2*pi;
    AL=median(AnglesL)*180/pi;

    AR=median(atan2(edgeR.ny(indexR-nvR:indexR-1),edgeR.nx(indexR-nvR:indexR-1)))*180/pi;

    shiftcoords=@(coordinates,shift) coordinates-ones(size(coordinates,1),1)*shift;
    rotatcoords=@(coordinates,angle)([cosd(angle),-sind(angle);sind(angle),cosd(angle)]*coordinates')';

    dataL=left(indexL-nvL:indexL-1,:);
    meanL=mean(dataL);
    dataR=right(indexR-nvR:indexR-1,:);
    meanR=mean(dataR);

    RotL=-AL-90;
    RotR=-AR-90;

    RotatedL=rotatcoords(shiftcoords(dataL,meanL),RotL);
    RotatedR=rotatcoords(shiftcoords(dataR,meanR),RotR);

    %---------------------------------------------------------------------------
    % Rotate and shift baselines to the same coordinate system as datapoints to improve fit
    %---------------------------------------------------------------------------

    baseline=[x0L,y0L;x0R,y0R];
    baseL=shiftcoords(baseline,meanL);
    baseL=rotatcoords(baseL,RotL);
    baseR=shiftcoords(baseline,meanR);
    baseR=rotatcoords(baseR,RotR);

    %---------------------------------------------------------------------------
    %                    fit the data to polynomial
    %---------------------------------------------------------------------------
    pL = polyfit(RotatedL(:,1),RotatedL(:,2),poly_ord);
    pR = polyfit(RotatedR(:,1),RotatedR(:,2),poly_ord);

    %---------------------------------------------------------------------------
    %        Find triplepoint by intersection between polynomial and baseline
    %---------------------------------------------------------------------------    
    if isempty(RotatedL) || isempty(RotatedR)
        Data.CAL = [];
        Data.CAR = [];
        Data.TLL = [];
        Data.TLR = [];
        Data.CoordsR = [];
        Data.CoordsL = [];
        Data.EvalPolyL = [];
        Data.EvalPolyR = [];
    else
        [x0L,y0L]=linefitinteraction(pL,RotatedL(end,:),[baseL(1,1),baseL(1,2)],[baseL(2,1),baseL(2,2)]);
        [x0R,y0R]=linefitinteraction(pR,RotatedR(end,:),[baseR(2,1),baseR(2,2)],[baseR(1,1),baseR(1,2)]);

        slopeL=atan2(polyval(polyder(pL),x0L),1)*180/pi;
        slopeR=atan2(polyval(polyder(pR),x0R),1)*180/pi;

        %---------------------------------------------------------------------------
        %        Evaluate polynomial at triplepoint
        %---------------------------------------------------------------------------

        BaseAngleL=atan2(baseL(2,2)-baseL(1,2),baseL(2,1)-baseL(1,1))*180/pi;
        BaseAngleR=atan2(baseR(1,2)-baseR(2,2),baseR(1,1)-baseR(2,1))*180/pi;

        CAL=BaseAngleL-slopeL;
        CAR=180-(BaseAngleR-slopeR);

        Data.CAL=CAL;
        Data.CAR=CAR;

        Data.TLL=shiftcoords(rotatcoords([x0L,y0L],-RotL),-meanL);
        Data.TLR=shiftcoords(rotatcoords([x0R,y0R],-RotR),-meanR);
        Data.CoordsL=dataL;
        Data.CoordsR=dataR;

        xL=linspace(RotatedL(1,1),RotatedL(end,1))*2;
        EvalPolyL=[xL',polyval(pL,xL)'];
        EvalPolyL=rotatcoords(EvalPolyL,-RotL);
        EvalPolyL=shiftcoords(EvalPolyL,-meanL);



        xR=linspace(RotatedR(1,1),RotatedR(end,1))*2;
        EvalPolyR=[xR',polyval(pR,xR)'];
        EvalPolyR=rotatcoords(EvalPolyR,-RotR);
        EvalPolyR=shiftcoords(EvalPolyR,-meanR);


        Data.EvalPolyL=EvalPolyL;
        Data.EvalPolyR=EvalPolyR;
    end
end
