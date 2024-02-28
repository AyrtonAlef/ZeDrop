function B = AverageBaseline(xBase,yBase)
% AVERAGEBASELINE Function that averages the baseline found in several 
% images. This is done by fitting the x,y coordinates in xBase yBase to a 
% straight line. Then all outliers, defined as having a distance to the 
% found baseline of more than two standard deviations, are removed. The 
% outliers are considered to be failed baseline determinations. The 
% remaining datapoints are then fitted to a straight line and the data for 
% the line (B) is returned by the function. The averaged baseline is found 
% as y= B(1)+x*B(2);

xBase=xBase(~isnan(xBase));
yBase=yBase(~isnan(yBase));

B = [ones(length(xBase),1) xBase] \ yBase;

a=-B(2);
c=-B(1);

disttoline=abs(a*xBase+yBase+c)/sqrt(a^2+1^2);

outliers=disttoline>2*std(disttoline);
xBase=xBase(~outliers);
yBase=yBase(~outliers);
B = [ones(length(xBase),1) xBase] \ yBase;
end

