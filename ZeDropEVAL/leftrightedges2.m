function [edgeL,edgeR]=leftrightedges2(longestedge1,longestedge2,x_AxisSym,y_Int,x_IntLeft,x_IntRight,jump_dist,filterMode)
% LEFTRIGHTEDGES2 Modified leftrightedges function that identifies left and
% right edges and filter coordinates according to needle position
% Returns left and right edges in edgeL.x edgeL.y....
% Output points are sorted so e.g. edgeL.y(1) is near the apex and edgeL.y(end) is
% at the bottom of the drop.

%%{
%Identify left and right edges
if mean(longestedge1.x) > x_AxisSym
    traceR = [longestedge1.y,longestedge1.x]; %converting into row,col format
    traceL = [longestedge2.y,longestedge2.x]; %converting into row,col format
else
    traceR = [longestedge2.y,longestedge2.x]; %converting into row,col format
    traceL = [longestedge1.y,longestedge1.x]; %converting into row,col format
end

%Identify starting point and filter left and right edges
[~,indAscR]= sort(traceR(:,1),'ascend'); %Sorting traceR in ascending order according to y coordinate
if filterMode == 'H' %Horizontal filter for drop profile
    indCondR = find(traceR(:,1) > y_Int); %Finding position in traceR that are below y_Int in the image
else %Vertical filter for drop profile
    indCondR = find(traceR(:,2) > x_IntRight); %Finding position in traceR thar are above x_IntRight in the image
end
indexR = ismember(indAscR,indCondR);
divideR = indAscR(indexR); 
right = traceR(divideR,:); %Coordinates of the filtered right edge in ascending order according to y coordinate
indexStartR = divideR(1); %Index of the starting point in the right edge
xStartR = traceR(indexStartR,2);
yStartR = traceR(indexStartR,1);
    
[~,indAscL]= sort(traceL(:,1),'ascend'); %Sorting traceL in ascending order according to y coordinate
if filterMode == 'H' %Horizontal filter for drop profile
    indCondL = find(traceL(:,1) > y_Int); %Finding position in traceL that are below y_Int in the image
else %Vertical filter for drop profile
    indCondL = find(traceL(:,2) < x_IntLeft); %Finding position in traceL that are below y_Int in the image
end
indexL = ismember(indAscL,indCondL);
divideL = indAscL(indexL); 
left = traceL(divideL,:); %Coordinates of the filtered left edge in ascending order according to y coordinate
indexStartL = divideL(1); %Index of the starting point in the left edge
xStartL = traceL(indexStartL,2);
yStartL = traceL(indexStartL,1);

% Sort coordinates to start at apex
[~,li,lcut] = sortsnake(left,[xStartL,yStartL],jump_dist);
[~,ri,rcut] = sortsnake(right,[xStartR,yStartR],jump_dist);

% assign all parameters in edges to edgeL and edgeR
if mean(longestedge1.x) > x_AxisSym
    fields = fieldnames(longestedge1);
    edgeL = longestedge2;
    edgeR = longestedge1;
    for i = 1:numel(fields)
    edgeL.(fields{i}) = longestedge2.(fields{i})(divideL);
    edgeL.(fields{i}) = edgeL.(fields{i})(li);
    edgeL.(fields{i}) = edgeL.(fields{i})(lcut);

    edgeR.(fields{i})= longestedge1.(fields{i})(divideR);
    edgeR.(fields{i})= edgeR.(fields{i})(ri);
    edgeR.(fields{i})= edgeR.(fields{i})(rcut);
    %edgeR.(fields{i}) = flip(edgeR.(fields{i}));
    end
else
    fields = fieldnames(longestedge1);
    edgeL = longestedge1;
    edgeR = longestedge2;
    for i = 1:numel(fields)
    edgeL.(fields{i}) = longestedge1.(fields{i})(divideL);
    edgeL.(fields{i}) = edgeL.(fields{i})(li);
    edgeL.(fields{i}) = edgeL.(fields{i})(lcut);

    edgeR.(fields{i}) = longestedge2.(fields{i})(divideR);
    edgeR.(fields{i}) = edgeR.(fields{i})(ri);
    edgeR.(fields{i}) = edgeR.(fields{i})(rcut);
    %edgeR.(fields{i}) = flip(edgeR.(fields{i}));
    end
end
%}
%{
%Identify and filter left edge and right edge
if mean(longestedge1.x) > x_AxisSym
    divide_R = (longestedge1.y > y_Int);
    divide_L = (longestedge2.y > y_Int);
    % Assign all parameters in edges to edgeL and edgeR
    fields = fieldnames(longestedge1);
    edgeL = longestedge2;
    edgeR = longestedge1;
    for i = 1:numel(fields)
    edgeL.(fields{i})=flip(longestedge2.(fields{i})(divide_L));
    edgeR.(fields{i})=longestedge1.(fields{i})(divide_R);
    end
else
    divide_R = (longestedge2.y > y_Int);
    divide_L = (longestedge1.y > y_Int);
    % Assign all parameters in edges to edgeL and edgeR
    fields = fieldnames(longestedge1);
    edgeL = longestedge1;
    edgeR = longestedge2;
    for i = 1:numel(fields)
    edgeL.(fields{i})=flip(longestedge1.(fields{i})(divide_L));
    edgeR.(fields{i})=longestedge2.(fields{i})(divide_R);
    end
%}
%Returning array in ascending or descending order according to x position
%{
%Identify and filter left edge and right edge
if mean(longestedge1.x) > x_AxisSym
    [~,indAscR]= sort(longestedge1.x,'ascend');
    indCondR = find(longestedge1.y > y_Int);
    index_R = ismember(indAscR,indCondR);
    divide_R = indAscR(index_R); 
    [~,indAscL]= sort(longestedge2.x,'descend');
    indCondL = find(longestedge2.y > y_Int);
    index_L = ismember(indAscL,indCondL);
    divide_L = indAscL(index_L);

    % Assign all parameters in edges to edgeL and edgeR
    fields = fieldnames(longestedge1);
    edgeL = longestedge2;
    edgeR = longestedge1;
    for i = 1:numel(fields)
    edgeL.(fields{i})=longestedge2.(fields{i})(divide_L);
    edgeR.(fields{i})=longestedge1.(fields{i})(divide_R);
    end
else
    [~,indAscR]= sort(longestedge2.x,'ascend');
    indCondR = find(longestedge2.y > y_Int);
    index_R = ismember(indAscR,indCondR);
    divide_R = indAscR(index_R); 
    [~,indAscL]= sort(longestedge1.x,'descend');
    indCondL = find(longestedge1.y > y_Int);
    index_L = ismember(indAscL,indCondL);
    divide_L = indAscL(index_L);

    % Assign all parameters in edges to edgeL and edgeR
    fields = fieldnames(longestedge1);
    edgeL = longestedge1;
    edgeR = longestedge2;
    for i = 1:numel(fields)
    edgeL.(fields{i})=longestedge1.(fields{i})(divide_L);
    edgeR.(fields{i})=longestedge2.(fields{i})(divide_R);
    end
%}
%end

function [sorted,index_trace,indexcut] = sortsnake(trace,P0,jump_dist)
%function that sorts all coordinates to by starting at at specified point
%P0. Output sorted is [P0;P1;P2;...Pi] where P0 is specified starting point
%P1 is the point in trace closest to P0, P2 is the point in trace exluding
%P0, Pi is the point in trace closest to Pi-1 excluding P1:Pi-2. 
% This also means that if the data points forks the trace will choose one
% branch and avoid the other. When reaching the end of the chosen branch
% the boundary will jump to the other branch. If this happens the script will stop the
% boundary if the jump is more than jump_dist value (in pixels).
% Input
%     trace:        a list of points to be sorted in [row,col] coordinates
%     P0:           the starting coordinate in [col,row]
%     jump_dist:     number of pixels that stops the search for new branches (original = 5). According to the perspective of the author (Ayrton), this value shoud be equal to margin. 
% Output
%     sorted:       Sorted version of trace
%     index_trace:  logic vector that ensures sorted=trace(index_trace,:)
%     indexcut:     logic vector that cuts of secondary branches.
%
data=[flip(P0);trace];
dist = pdist2(data,data);

N = size(data,1);

indexcut=ones(1,N-1);
result = NaN(1,N);
result(1) = 1; % first point is first row in data matrix

for ii=2:N
    dist(:,result(ii-1)) = Inf;
    [mindist, closest_idx] = min(dist(result(ii-1),:));
        result(ii) = closest_idx;
        if mindist>jump_dist
            indexcut(ii-1:end)=0;
        end
end
result(1)=[];
result=result-1;
sorted=trace(result',:);
indexcut=indexcut==1;
index_trace=result;