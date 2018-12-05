%% Remove bottom part skin folds
% Script written by Mariela A. Porras-Chaverri April 2012
% University of Wisconsin-Madison

%% METHOD USING THE 'SHORTEST DISTANCES CHEST-WALL EDGE TO BREAST EDGE'
% to segment out the non-breast parts of the bottom of the image (skin folds):

% Get the distance between the chest-wall side of the SegmBreastNoPec and the 
% breast edge. The narrowest parts (from the bottom) should
% correspond to the points where the breast joins the skin fold.
SegmBreastNoBott = SegmBreast;
Distance = sum(SegmBreastNoBott,2);

% Now I find the minimum in the second part of 'Distance', this is where
% the curve 'under the breast' will be:
HalfImage = roundoff(m_Image/2,0);
DistanceB = Distance;

i = 1:HalfImage;
DistanceB(i)=800; % artificially increase the value here to make sure the 
%min will be found at the bottom curve of the breast. 

j= (m_Image-10):m_Image;
DistanceB(j)=800;

[min_row min_col] = min(DistanceB); % this one gives me the X and Y values for the lowest min from the bottom of the image

x_cut = min_row;
y_cut = min_col;

%Image((m_Image-5):m_Image,:)=0;

i=y_cut:1:m_Image; 
j=1:1:n_Image;

SegmBreastNoBott(i,j)=0; % '..NoBott' refers to the BOTTom skin fold
%figure, imshow(SegmBreastNoBott), impixelinfo
