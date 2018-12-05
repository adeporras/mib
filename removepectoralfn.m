function [ Triangle ] = removepectoralfn( anatomTriangle,kPointA,kPointB,m_Image,n_Image )
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
end

