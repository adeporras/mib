%% MIB algorithm for Raw density maps - May2015 update.
% Mariela A. Porras-Chaverri  February 28th, 2011
% Script was updated for ease of code maintenance -May2015

% This script divides the breast image into three sections. The sections
% will be used to determine the relative glandular composition for the
% ´perpendicular´ (or almost...) view. For example, the three sections in
% the CC view will be assumed to represent the three sections for the MLO
% view, and viceversa
% Input is a matrix file with the heights of the fibroglandular tissue 
% for each pixel (Image), CBT, and the curve P(h_int) that relates pixel value to
% column-height of glandular tissue.
% Breast edge corresponds to the rounded part of the breast that is not
% directly in contact with the compressor plate. It's assumed to be only
% glandular tissue. 
% Breast core corresponds to the part of the breast that is in direct
% contact with the compressor plates. It's assumed to have constant height
% equal to CBT. 


display(' inStruct(1).mTranfs debe ser igual a 2.03')
 inStruct(1).mTranfs

 
 mTrasnf = inStruct(1).mTranfs,
 
 
 
TH=floor(bTransf);% TH separates breast core from breast edge

%CBTmm is the compressed breast thickness in mm that Volpara states for the breast,
% it is NOT corrected for the presence of skin.

A_pix = pixSize*pixSize; % (mm^2)this is the size of the pixels. I'm using 300 microns (0.3 mm).
N_pix=numel(Image);% number of pixels in the image file

% The input matrix 'Image' should have the Volpara pixel values. NOT the h_int values.
% Here I will use the relationship between pixel value (P) and h_int (mm)
% this relationship was given by Volpara.

Image = double(Image);% Image is in 'pixel value' units
I = (Image-bTransf)./mTransf;% I is in 'mm' of h_int tissue
auxImage = I; % This matrix will be used to get the Core/Edge masks.

%% Find the glandular tissue information from the h_int image.
% Negative or zero pixel values in the 'I' matrix indicate that there is no glandular tissue for that pixel.(Negative
% values belong to the breast edge).
k=find(I<0);
I(k)=0;
I=double(I);

%% Division of density map into sections:
% Let's find the masks that will define each of the three sections
[m,n]=size(SegmBreast);
MA=zeros(size(SegmBreast)); % output matrices
MB=zeros(size(SegmBreast));
MC=zeros(size(SegmBreast));

PixNum=sum(sum(SegmBreast)); % number of pixels with value=1 for the J#neg.

Asec=(PixNum)/3;
Asum1=0;
Asum2=0;

% Mask for first section (layer1 in HLB)
for i=1:m,
    for j=1:n,
        if Asum1 <= Asec;
            Asum1 = Asum1 + SegmBreast(i,j);
            MA(i,j)= SegmBreast(i,j);
        elseif MA(i,j)==0;
        end
    end
end
SegmBreastB = SegmBreast-MA; 

% Mask for second section (layer2 in HLB)
for ii=1:m,
    for jj=1:n,
        if Asum2 <= Asec;
            Asum2 = Asum2 + SegmBreastB(ii,jj);
            MB(ii,jj)= SegmBreastB(ii,jj);
        elseif MB(ii,jj)==0;
        end
    end
end
% Mask for third section (layer3 in HLB)
MC=SegmBreast-MA-MB;

% Apply these masks to the image I:
IA=I.*MA; 
IB=I.*MB;
IC=I.*MC;

%% Division of breast into Core and Edge.
% Calculate the masks for the Core and Edge, for each of the HLB layers:
preCoreMask1 =(auxImage>0); 
preEdgeMask1 =(auxImage<=0); 
% Remove any tiny holes that may be left in the masks:
preCoreMask = imfill(preCoreMask1, 'holes');
preEdgeMask = imfill(preEdgeMask1, 'holes');
se = strel('disk',9);
CoreMask=imclose(preCoreMask, se);
EdgeMask=imclose(preEdgeMask, se);

% Breast Core, and Core for each breast section:
Icore = double(Image).* CoreMask.*SegmBreast; 
IcoreA = double(Image).* CoreMask.*SegmBreast.*MA; 
IcoreB = double(Image).* CoreMask.*SegmBreast.*MB;
IcoreC = double(Image).* CoreMask.*SegmBreast.*MC;
% Breast Edge, and Edge for each breast section:
Iedge = double(Image).* EdgeMask.*SegmBreast;
IedgeA = double(Image).* EdgeMask.*SegmBreast.*MA;
IedgeB = double(Image).* EdgeMask.*SegmBreast.*MB;
IedgeC = double(Image).* EdgeMask.*SegmBreast.*MC;
figure, imagesc(Icore), impixelinfo


%% Calculation of glandular and total volumes for breast and breast sections.

% Volume of glandular tissue for each section, whole breast:
V_glandtot = A_pix*sum(sum(I)); % in mm^3
V_glandA = A_pix*sum(sum(IA*1)); % layer1
V_glandB = A_pix*sum(sum(IB*1)); % layer2
V_glandC = A_pix*sum(sum(IC*1)); % layer3

%% Volume of glandular tissue for each section, Core only:
% For the first section (layer1)
IcoreABW=IcoreA*1;
k=find(IcoreABW>0);
IcoreABW(k)=1;
V_coreA=A_pix*CBTmm*sum(sum(IcoreABW));
clear k
% For the second section (layer2)
IcoreBBW=IcoreB*1;
k=find(IcoreBBW>0);
IcoreBBW(k)=1;
V_coreB=A_pix*CBTmm*sum(sum(IcoreBBW));
clear k
% For the third section (layer3)
IcoreCBW=IcoreC*1;
k=find(IcoreCBW>0);
IcoreCBW(k)=1;
V_coreC=A_pix*CBTmm*sum(sum(IcoreCBW));
clear k

%% Volume of glandular tissue for each section, Edge only.
%  The edge volume is calculated assuming a spherical shape fo the breast. 
% Total area of the projected breast (mm^2)
A_breast=A_pix*sum(sum(SegmBreast)); % in mm^2.

IedgeBW=Iedge*1;
jk=find(IedgeBW>0);
IedgeBW(jk)=1; % this makes a matrix with ones where the edge is, and we
% can use it to get the edge area and volume.
A_edge=A_pix*sum(sum(IedgeBW)); % in mm^2.

% First we will assume that the distance (R) from nipple to chest wall edge of
% the mammogram is half as much as the breast width (D, top to bottom in the
% images.
D = pixSize*sum(SegmBreast(:,1)); % in mm
R = D/2; % in mm.
Rcm = R/10; % in cm.

%Then we will use the area of the breast edge to calculate the 'approximate
%internal breast edge radius' (r), that is, the approx distance between the
%chest wall edge of the image and the internal border of the breast edge.

% Internal radius squared is:
r2=(R^2 - (2*A_edge)/pi );% in mm^2

% now let's get the relevant volumes, first the volume of the spherical
% zone that corresponds to the "full" compressed breast:
Vol1 = (1/12) * pi * CBTmm * ((CBTmm^2/4)+ 3*R^2 + 3*r2); % in mm^3

% now we get the volume of the cylinder that encloses the "core" part of
% the compressed breast:
Vol2= (1/2)*pi*r2*CBTmm; % in mm^3

% and now we can get the edge volume:
V_edge = Vol1 - Vol2; % in mm^3

% and I will divide the edge volume equally into the three sections (layers):
V_edgeA = V_edge/3;
V_edgeB = V_edgeA*1;
V_edgeC = V_edgeB*1;

% Now, the total volume is going to be the sum of the volume from the edge
% and from the breast core/glandular disk(without the breast edge). Note
% how it is different for each layer. The edge is also symmetric about the
% midbreast plane.

V_tot = V_edge + V_coreA + V_coreB + V_coreC;

V_totA = V_edgeA + V_coreA;
V_totB = V_edgeB + V_coreB;
V_totC = V_edgeC + V_coreC;

%% Glandularity for each breast section:
glandA = (V_glandA/V_totA)*100; % layer1
glandB = (V_glandB/V_totB)*100; % layer2
glandC = (V_glandC/V_totC)*100; % layer3
glandtot = (glandA + glandB+glandC)/3; % Whole breast glandularity.

% Glandularity values for the HLB layers of the other-view case:
HLB_layers=[glandtot, glandA, glandB, glandC];
Idist = (glandA-glandC)/glandB;

% Comparison to the value estimated by other methods, such as Volpara or
% from phantoms. 
wholeGlanDiffperc = 100*(theorGland-glandtot)/((theorGland+glandtot)/2);
wholeGlanDiff = theorGland-glandtot;