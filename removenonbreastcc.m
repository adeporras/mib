%% Remove the non-breast areas (NOT using regiongrowing)
% Script written by Mariela A. Porras-Chaverri April 2012
% University of Wisconsin-Madison

% This method should be faster than the region growing method I have used
% so far. (Apr11th,2012)

% Let's mark the point where the top row and the furthest (~nipple) point
% meet:
% Let's find the countour of the breast(This may be useful to get rid of
% the use of 'region growing' in MIB):
Contour = edge(Image);
%%figure, imshow(Contour), impixelinfo
% let's fill the contour:
SE = strel('diamond',2);
Contour = imdilate(Contour, SE);
%%figure, imagesc(Contour), impixelinfo
k_Contour= find(Contour==1); % IMPORTANT TO KEEP THE k_Contour VARIABLE NAME


% Now that we have the contour points, we have to get rid of all the points
% to the right of that contour, since those will correspond to the
% 'non-breast areas', to do this, I will use the 'imfill' function.

SegmBreast = Contour;
% figure,imagesc(SegmBreast), impixelinfo

% the 'imfill' function needs the locations. To make sure , I will use
% points in the first column. 
k_loc = [450]; % this point is around the middle of the image left side
%k_loc = [m_Image*n_Image];
SegmBreast= imfill(SegmBreast, k_loc',4);
%%figure,imagesc(SegmBreast), impixelinfo
[L, numLabel] = bwlabel(SegmBreast, 4);
RGB = label2rgb(L);
