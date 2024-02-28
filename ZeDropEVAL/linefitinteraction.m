function [x0,y0]=linefitinteraction(p,P0,P1,P2)


P1x = P1(1);
P1y = P1(2);
P2x = P2(1);
P2y = P2(2);

xt=@(t)(P2x-P1x)*t+P1x;
yt=@(t)(P2y-P1y)*t+P1y;
fun1 = @(t) polyval(p,xt(t))-yt(t);
fun2 = @(t) abs(polyval(p,xt(t))-yt(t));

if P1x==P2x
    xr=P1x;
    yr=polyval(p,xr);
else 
a=(P2y-P1y)/(P2x-P1x);
b=P1y-a*P1x;

p2=p;
p2(end)=p2(end)-b;
p2(end-1)=p2(end-1)-a;

xr=roots(p2);
xr=xr(real(xr)==xr);
yr=polyval(p,xr);
end
dist=sum(sqrt(bsxfun(@minus,[xr,yr],P0).^2),2);
x0=xr(min(dist)==dist);
y0=yr(min(dist)==dist);
