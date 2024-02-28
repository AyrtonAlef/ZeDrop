function [longestedge1,longestedge2] = findlongestedgeLR(edges,imsize,margin)
% FINDLONGESTEDGELR Modified findlongestedge function to return the two 
% longest edge in the image. For needle in images, this will correspond to 
% left and right edges of the drop.
%   INPUT:
% edges in format
%       edges.x         xcoordinates, subpixel
%       edges.y         ycoordinates, subpixel
%       edges.position   linear index of x,y position in pixel
% imsize in format
%       [row,col]
% margin, pixel distance between edges that should be considered part of same contour.


% First we determine the amount of separate edges found spaced by more than
% 1 pixel. This is done in pixel resolution with bwlabel matlab function


logical=zeros(imsize);
logical(edges.position)=1; % create binary image with white pixels at edges and black where there is no edges.
[L,NumberOfboundaries]=bwlabel(logical,8); % find and label seperate edges

%Load the seperate boundaries into a cell structure that have dimension
%{number of boundaries, 1}

boundary=cell(NumberOfboundaries,1); 
for ii=1:NumberOfboundaries
    [row,col]=find(L==ii);
    boundary{ii}=[row,col];
end


% Calculate the minimum distance between all points in one boudnary and all
% other points in another boundary. This is to be used for joining
% boundaries that are closer than the threshold set in margin input
% variable. The following loop creates a matrix similar to the output from
% pdist2 but wuith the minimum distance between points boundary ii and ss
% as the value in distmatrix(ii,ss) and distmatrix(ss,ii)

distmatrix=NaN(NumberOfboundaries);

for ii=1:NumberOfboundaries
    if ii+1<=NumberOfboundaries
    for ss=ii+1:NumberOfboundaries
    dist = pdist2(boundary{ii},boundary{ss});
    distmatrix(ii,ss)=min(min(dist));
    distmatrix(ss,ii)=distmatrix(ii,ss);
    end 
    end
end

% Now we have calculated all relevant distances between boundaries found and
% can start joining boundaries that are closely spaced. We start by giving
% each boundary a seperate index and the giving changing the index of
% closely space boundaries to the index of their neighbour

BoundaryIndex=1:NumberOfboundaries;
for ii=1:NumberOfboundaries
    for ss=1:NumberOfboundaries
        if distmatrix(ii,ss)<margin
            BoundaryIndex(BoundaryIndex==BoundaryIndex(ss))=BoundaryIndex(ii);
        end
   end
end

%Now there is several boundaries with same index and we combine the
%boundaries with same index into same cell

NewIndex=unique(BoundaryIndex); % the remaining boundary indexes after renaming
NumberOfNewIndexes=length(NewIndex); %number of remaining boundaries
CombinedBoundary=cell(1,NumberOfNewIndexes);


for ii=1:NumberOfboundaries
    index=find(NewIndex==BoundaryIndex(ii));
    CombinedBoundary{index}=[CombinedBoundary{index};boundary{ii}]; %joining boundaries of same index
end


% Calculating length the joined boundaries to find the longest one 
 boundarylength=zeros(1,NumberOfNewIndexes);
for ii=1:NumberOfNewIndexes
     boundarylength(ii)=length(CombinedBoundary{ii});
end

%[~,index]=max(boundarylength);
maxBoundaries = maxk(boundarylength,2); %Length the 2 longest edges
index1 = find(boundarylength==maxBoundaries(1),1); %Find the index of the longest edge
index2 = find(boundarylength==maxBoundaries(2),1,'last'); %Find the index of the second longest edge

%mainedge=CombinedBoundary{index}; %select the longest of boundaries
mainedge1 = CombinedBoundary{index1}; %select the longest of boundaries
mainedge2 = CombinedBoundary{index2}; %Select the second longest of boundaries
[I,J] = ind2sub(imsize,edges.position); %calculate the in

pts=[I,J]; 
%[~,loc] = ismember(mainedge,pts,'rows'); %calculate the location of the longest edge in the original edge structure
[~,loc1] = ismember(mainedge1,pts,'rows'); %calculate the location of the longest edge in the original edge structure
[~,loc2] = ismember(mainedge2,pts,'rows'); %calculate the location of the second longest edge in the original edge structure

fields = fieldnames(edges);
for i = 1:numel(fields)
    longestedge1.(fields{i})=edges.(fields{i})(loc1); %asign all properties of the input edge structure to the new longestedge structure
    longestedge2.(fields{i})=edges.(fields{i})(loc2); %asign all properties of the input edge structure to the new longestedge structure
end


