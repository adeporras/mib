%% Remove the non-breast areas (NOT regiongrowing)
% Script written by Mariela A. Porras-Chaverri May 2012
% University of Wisconsin-Madison

clear Contour
clear SegmBreast

Contour  = edge(Image, 'sobel');
%figure, imshow(Contour), impixelinfo
% let's fill the contour:
SE = strel('diamond',3);
Contour = imdilate(Contour, SE);

%figure, imshow(Contour), impixelinfo
k_Contour= find(Contour==1); % IMPORTANT TO KEEP THE k_Contour VARIABLE NAME

% Now that we have the contour points, we have to get rid of all the points
% to the right of that contour, since those will correspond to the
% 'non-breast areas'.
SegmBreast = Contour;
%figure,imshow(SegmBreast), impixelinfo

clear m_Image
clear n_Image

[m_Image n_Image]=size(Image);

k_loc = [m_Image*n_Image];

SegmBreast= imfill(SegmBreast, k_loc',8);
% locations = [450 200];
% SegmBreast= imfill(SegmBreast, locations, 8);
%figure,imshow(SegmBreast), impixelinfo

% here I find the 'islands' in the SegmBreast image, some of these
% are microcalcifications, but others are not. 
[L, numLabel] = bwlabel(SegmBreast, 8);
RGB = label2rgb(L);

% now let's find the area of each of the regions in the SegmBreast
% (including all kinds of islands)
%figure, imshow(SegmBreast)
for g =1:numLabel;
    k_L = find(L==g);
    areaL(g) = (sum(sum(L(k_L)))/g);
        
    % Now I can decide what I want to do which each kind of island.
       
    if areaL(g)==max(areaL);
%         display(strcat('I am the largest area:','__', num2str(areaL(g))))
        SegmBreast(k_L)=1;
    elseif areaL(g)~=max(areaL);
%         if areaL(g)<(min(areaL)+5);
        SegmBreast(k_L)=0; % when it was '0' in some cases it got rid of all the breast :/
      % so I changed it to '1' to see how it goes... 
        SegmBreast(k_L)=1;
              
        
        % display(strcat('This is this:','__', num2str(areaL(g))))
%         elseif areaL(g)>(min(areaL)+5);
%         SegmBreast(k_L)=0;
%         display(strcat('This is something else:','__', num2str(areaL(g))))
%         else
%         SegmBreast(k_L)=0;
%         display(strcat('That is that:','__', num2str(areaL(g))))
%         end
    end
end

%display(strcat('I am the smallest area:','__', num2str(min(areaL)))) 
%figure, imshow(SegmBreast)
    
[L, numLabel] = bwlabel(SegmBreast, 8);
RGB = label2rgb(L);
%     PerimL= bwperim(L);
%     AreaL = bwarea(L);
% figure,imshow(PerimL)
% figure,imshow(AreaL)

% % fullscreen = get(0,'ScreenSize');
% % fig = figure('Visible','off','Position',[0 -50 fullscreen(3) fullscreen(4)])
% % subaxis(1,1,1,  'Spacing', 0, 'Padding', 0, 'Margin', 0), subimage(RGB)  
% % % figure, imagesc(RGB)
% % axis off
% %     name = strcat(pathStart,'\Dropbox\From h_int to DgN-HLB\Images - Highnam100Full\\MIB figures\', char(patientID),'__Objects');
% %     saveas(fig,name,'png')


% STATS = regionprops(SegmBreast, 'BoundingBox');
% STATS = regionprops(SegmBreast);
 %[ul_corner width]= STATS.BoundingBox; % not sure what this output is...

% but now we have to take the negative of the SegmBreast matrix
% and use that one in the calculations:

SegmBreast = im2bw(ait_imneg(double(SegmBreast)));

%figure,imshow(SegmBreast), impixelinfo

%  nonBreastPixVal=Image(m_Image,1) % this is the pixel value of the non-breast part.
%  ImageEdgeMask = Image-nonBreastPixVal;
%  figure,imshow(ImageEdgeMask), impixelinfo
% SegmBreast = SegmBreast.*im2bw(ImageEdgeMask);
%ImageEdgeMask = im2bw(ait_imneg(double(ImageEdgeMask)));
% %figure,imshow(ImageEdgeMask), impixelinfo
%SegmBreast = imdilate(SegmBreast, SE);
%figure,imshow(SegmBreast), impixelinfo

testPix = SegmBreast(200,100); % for non-Pect(clinical cases)


count =1;

% % % while testPix==0;
% % %     count = count+1;
% % %     display('Let`s do this thing!')
% % % if count<=5;
% % %     removenonbreastagain;
% % %     SegmBreast = imdilate(SegmBreast, SE);
% % %         
% % % else
% % %     %display(strcat('I tried_',num2str(count),'_times and it did not work :( '))
% % %     break
% % % end
% % % 
% % % %SegmBreast = im2bw(ait_imneg(double(SegmBreast)));
% % % end

SegmBreast = imdilate(SegmBreast, SE);
%figure,imshow(SegmBreast), impixelinfo    
