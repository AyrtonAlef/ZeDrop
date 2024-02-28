function [row0,col0,index] = findreflection(trace,n,cut)
% FINDREFLECTION Algorithm for finding the reflection of the drop and 
% calculating tripple point on each side of the drop. This is done by 
% searching for the index ii where the ii-n:ii and ii+n:n points fitted to 
% two lines result in lines with a slope of equal magnitude but opposite 
% signs.

%   INPUT:
% Trace in [row,col] coordinates
% n integer number of points that are fitted 
% cut is the number of pixels from the top where the script should not look
% for reflections. This can be used to avoid problems where the script
% finds the center of the drop for contact angles >90 deg.
%   OUTPUT: 
% [row0,col0] is the coordinate for the tripple line
% index is the index of the point in trace that is closest to the triple
% line.
traceRef = trace; %Modified by Ayrton;
ncut=sum(trace(:,2)<=cut);
trace=trace(trace(:,2)>cut,:);

% predefine vectors before loop
slopevec=nan(size(trace));  % stores the slopes of the linear fits
n2=2; % number of points from ii where the fit uses points from

for ii=n:length(trace)-n
    x1=trace(ii-n+1:ii-n2,2); % exctract the coordinates for fit
    y1=trace(ii-n+1:ii-n2,1);
    x2=trace(ii+n2:ii+n-1,2);
    y2=trace(ii+n2:ii+n-1,1);
    B1 = [ones(length(x1),1) x1] \ y1; % fast linear fit
    B2 = [ones(length(x2),1) x2] \ y2;
    slopevec(ii,1)=B1(2); % store slope in vector
    slopevec(ii,2)=B2(2);
end
% find indices in slopevec with different signs
sign_slope=sign(slopevec);
clean=sum(sign_slope,2)~=0;
slopevec(clean,1)=NaN;
slopevec(clean,2)=NaN;
t=sum(slopevec,2);
% find the point where slopes are of equal magnitude but opposite signs
% (sum is close to 0)

index=find(abs(t)==min(abs(t)));

% Find the triple line by fitting lines to points above and below the
% reflection point. This is repeated twice to improve accuracy.
for ii=1:2
    if isempty(trace) %Modified by Ayrton Pereira. Check if trace is empty
        fprintf('Trace array empty! \n')
        fprintf('Probably cut value is invalid or edge detection failed. \n');
        fprintf('Suggestion1: Check that the entire edge is being identified. Review the edge detection method used. \n')
        row0 = 1; %Return the apex coordinates
        col0 = 1; %Return the apex coordinates
        break;
    elseif isempty(index) %Modified by Ayrton Pereira. Check if index is empty
        fprintf('Index position empty! \n')
        fprintf('Probably number of fitted points (n) or cut value is invalid. Please, enter new values. \n');
        fprintf('Suggestion1: Try to change cut value first to avoid problems where the script finds the center of the drop for contact angles. \n');
        fprintf('Suggestion2: Then, lower n value. \n')
        fprintf('Suggestion3: Check that the entire edge is being identified. Review the edge detection method used. \n')
        %index2 = index + ncut;
        row0 = traceRef(1,1);
        col0 = traceRef(1,2);
        break;
    elseif max(index-n+1 < 1) || max(index+n-1 > length(trace)) %Modified by Ayrton Pereira. Check if index exceeds trace length
        fprintf('Index is negative or exceeds the length of trace! \n')
        fprintf('Number of points that are fitted (n) or number of pixels from the top (cut) invalid. Please, enter new values. \n');
        fprintf('Suggestion1: Try to change cut value first to avoid problems where the script finds the center of the drop for contact angles. \n');
        fprintf('Suggestion2: Then, lower n value. \n')
        fprintf('Suggestion2: Check that the entire edge is being identified. Review the edge detection method used. \n')
        row0 = trace(round(mean(index)),1);
        col0 = trace(round(mean(index)),2);
        break;
    else
        x1=trace(index-n+1:index-n2,2);
        y1=trace(index-n+1:index-n2,1);
        x2=trace(index+n2:index+n-1,2);
        y2=trace(index+n2:index+n-1,1);
        B1 = [ones(length(x1),1) x1] \ y1;
        B2 = [ones(length(x2),1) x2] \ y2;

        a=B1(2);
        b=B2(2);
        c=B1(1);
        d=B2(1);

        % find crossing between lines
        col0=(d-c)/(a-b);
        row0=a*col0+c;
        [~,index]=min(sum(bsxfun(@minus,trace,[row0,col0]).^2,2)); %the new index is the point closest to the new intersection
    end
end

index=index+ncut;

if row0 < 0 || col0 <0 %Modified by Ayrton Pereira. Check if returned CPs are negative
    row0 = traceRef(1,1);
    col0 = traceRef(1,2);
end
% row0=trace(index,1);
% col0=trace(index,2);

