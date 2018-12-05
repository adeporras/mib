%% Remove pect muscle in MLO (not present in CC)
% Script written by Mariela A. Porras-Chaverri
% Some 'major' changes on 09jan14
% University of Wisconsin-Madison

% First we remove the triangular area from the upper-left corner that contains the
% pectoral muscle. 

% Let's find the points of interest: the point furthest from the left side
% (chest) that contains breast (FURTHEST BREAST POINT), and the point where
% the pectoral muscle intersects the left side of the image (PECT MUSCLE
% POINT).
%
%% FURTHEST BREAST POINT 
% Here we change the triangle little by little to match
% the %gland for the whole breast calculated using the CC view. 

% The first triangles will have a corner where the breast intersects the
% upper side of the image: 

if strcmp(anatomTriangle,'yes')
    Furth_Point = [1 kPointA];
    Pect_Point = [kPointB, 1];
else
    stophere 
    Furth_Point = [1 topInt+iMLO];
    %% PECT MUSCLE POINT - Test 09jan14
    Pect_Point = [sideInt+iMLOC, 1];
    
end

%% BUILD TRIANGLE MASK
% Now we construct the triangle mask that gets rid of the pect muscle

m_line = (Pect_Point(2)-Furth_Point(2))/(Pect_Point(1)-Furth_Point(1));
b_line = Furth_Point(2);

[X,Y]=meshgrid(1:m_Image, 1:n_Image);

Triangle = Y - m_line * X - b_line > 0;
Triangle = Triangle';
%figure, imagesc(Triangle), impixelinfo





