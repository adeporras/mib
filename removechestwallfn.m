function [ SegmBreastNoWall ] = removechestwallfn( m_Image,startCropChest,SegmBreastNoBott )
%% Remove chest wall structures (left side crop)
% Script written by Mariela A. Porras-Chaverri April 2012
% University of Wisconsin-Madison
iMLOB=1;
x_chest=iMLOB+startCropChest;

i=1:1:m_Image;
j = 1:1:(x_chest);

SegmBreastNoWall = SegmBreastNoBott;
SegmBreastNoWall(i,j)=0;
%figure, imshow(SegmBreastNoWall), impixelinfo


end

