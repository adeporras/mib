% Mammography-Image Based algorithm 
% maporras A. Porras-Chaverri
% Master script (21dec13, based on earlier scripts)
% Rework on May2015 - Updated and optimized for speed and code maintenance.
% ********************************************
clear all,clc, close all, warning off, fclose all; tic
display('WELCOME TO THE MIB ALGORITHM - REFRESHED IN MAY2015')
pathStart = 'C:\Users\Mariela Porras C'; % Run the scripts from the laptop.

studyID='testingnewcode';
ss1='99-None-2017'; % date for HLBlayers spreadsheet
% This corresponds to the *PATIENT BREAST* ID, meaning, there is one
% for each set of CC-MLO views.(Each breast is considered as a 'separate patient')
% Therefore, if only one patientID is read, both the CC and the MLO are
% read and used by the script.
startPat = 1;
endPat   = 1;

partialPath = '\Dropbox\aaa_UCR\Proyectos Investigación\2017\915-b7-081_estimacion tejido glandular (hlb2)\Software and codes - playground\HLB research scripts\';
StudyBaseCamp=strcat(pathStart, partialPath, studyID, '\');
% and this is the file I to have in that folder with the info:
InfoFileName= strcat(studyID,'_info.xlsx'); % File where the case parameters (CBT, glands, etc) are stored.




SheetInfoFileName1='CC'; % Name of the sheet where the parameters are stored.
SheetInfoFileName2='MLO'; % Name of the sheet where the parameters are stored.
RangeSheetInfoFileName='A2:K405'; % Range in the sheet where the parameters are located.
DensityMapsFolder =strcat(StudyBaseCamp,'Density Maps\');

[caseInfoCC caseIDsCC] = xlsread(strcat(StudyBaseCamp,InfoFileName),SheetInfoFileName1,RangeSheetInfoFileName);
[caseInfoMLO caseIDsMLO] = xlsread(strcat(StudyBaseCamp,InfoFileName),SheetInfoFileName2,RangeSheetInfoFileName);
% caseIDs(pat,1) is CaseID
% caseInfo(pat,1) is whole breast gland Perc (Volpara value)
% caseInfo(pat,2) is CBT in mm
% caseInfo(pat,3) is mTransf
% caseInfo(pat,4) is bTransf
% Now we concatenate, First matrix is the CCs, Second is the MLOs:
caseInfo(:,:,1) = caseInfoCC;
caseInfo(:,:,2) = caseInfoMLO;
caseIDs(:,:,1) = caseIDsCC;
caseIDs(:,:,2) = caseIDsMLO;

numFiles = 2*size(caseIDs,1); % The number '2' comes from it being a pair of CC-MLO images.
N_patients = numFiles; % N_patients is actually the number of JPEGS I want to run.

pixSize = 0.3; % (mm)Size of the pixels in the density maps. Volpara maps' are 300 microns.
% Tolerance of difference between CC and segmented MLO (%, not pp)
diffCCMLOtol=5;
% When would a case be considered 'homogeneous'?
homogThreshold = 0.25 % from Highnam100 study

%% Param for pectoral removal
anatomTriangle = 'yes'; % controls the test values for FurthPoint, etc. (see removepectoral)
%control iMLO and iMLOC
topTriangleStep=5; % step to increase the upper cathetus
leftTriangleStep=5; % step to increase the left side cathetus

smallerTriangleSide = 40; % controls the dimensions of the 'starting triangle', linked to size of Image
lastTriangleTop =1;  % controls the dimensions of the 'ending triangle', linked to size of Image

%%  Param for chest wall removal
stepChest = 2; %control iMLOB
startCropChest=20; % first run had value 20
maxCropChest = 21; % first run had value 21
maxAttemptsChest=maxCropChest-startCropChest;

folder=strcat(pathStart,partialPath,'1_MIBscripts\');
cd(folder);

for pat=1 %startPat:endPat;    
    % First let's work with the CC case:
    for CCorMLO = 1%:2; % i=1 for CC, i=2, MLO.
        patientID = char(caseIDs(pat,1,CCorMLO))
        densityMapFile = char(caseIDs(pat,6,CCorMLO))
        
        theorGland =caseInfo(pat,1, CCorMLO); % Whole breast percentage gland value given by Volpara
        CBTmm=caseInfo(pat,2, CCorMLO);
        mTransf=caseInfo(pat,3, CCorMLO); % slope P(t)
        bTransf=caseInfo(pat,4, CCorMLO); % intercept P(t)

        
        s_im=strcat(DensityMapsFolder,densityMapFile);
        jpeg = fopen(s_im,'r');
        if jpeg>0; % Call the Density Map image
            Image = imread(strcat(DensityMapsFolder,densityMapFile));
            % figure, imagesc(Image), impixelinfo
            Image(:,1:5)=[]; % We get rid of the first 5 columns in the Image. In some cases, there is an strip of empty pixels in the density maps from Volpara that messes up with the imfill.
            [m_Image n_Image]=size(Image);
            
            
            % DEFINIR AQUI EL STRUCTURE InStruct, usar valor temporal para
            % SegmBreast (ceros, mismo tamano que Image)
            cd(strcat(pathStart,'\Dropbox\HLB research scripts\1_MIBscripts'))
            if CCorMLO==1;% This is a CC case
                
                
                removenonbreastcc % Remove the non-breast areas of the Image; Output mask is 'SegmBreast'.
                % figure,imshow(SegmBreast),impixelinfo
                removefolds % Remove the bottom part of the breast (skin folds); Output mask is 'SegmBreastNoBott'.
                % figure,imshow(SegmBreastNoBott),impixelinfo
                % And finally, we combine all the masks together:
                SegmBreast = SegmBreast.*SegmBreastNoBott;
                % figure,imagesc(SegmBreast),impixelinfo
                AreaCC= sum(sum(SegmBreast));
                
                % ACTUALIZAR EL VALOR DEL CAMPO 'SegmBreast' en InStruct a
                % el valor real de SegmBreast que acabamos de encontrar
                
                
                % Now we call the MIB function that will get the values of the HLB
                % layers that correspond to the CC Image.Output is the SegmBreast mask and the HLB values
                
                
                stophere
                
                MIBraw;
                
                
                HLB_layersCC = HLB_layers;
                IdistCC = Idist;
                glandCC=glandtot;
                wholeGlanDiffpercCC = wholeGlanDiffperc;
                wholeGlanDiffCC = wholeGlanDiff;
                
                % And finally, save the HLB layer info to the HLBlayers-spreadsheet
                ResultsCC = horzcat(HLB_layersCC, IdistCC, theorGland, wholeGlanDiffpercCC, wholeGlanDiffCC);
                
                s_resultscc=strcat('B',num2str(pat+1),':I',num2str(pat+1));
                xlswrite(strcat(StudyBaseCamp,'MIB Outputs\',ss1,'_1-HLBlayers_',studyID,'.xls'),ResultsCC,SheetInfoFileName1,s_resultscc);
                
                s_pat = strcat('A',num2str(pat+1)); % Cell in the spreadsheet
                xlswrite(strcat(StudyBaseCamp,'MIB Outputs\',ss1,'_1-HLBlayers_',studyID,'.xls'),{patientID},SheetInfoFileName1,s_pat);
                
                s_labels='A1:I1';
                Labels = {'CaseID' 'HLB whole gland' 'lay1' 'lay2' 'lay3' 'Idist' 'theorGland' 'glandDiff(%)' 'glandDiff(pp)'};
                xlswrite(strcat(StudyBaseCamp,'MIB Outputs\',ss1,'_1-HLBlayers_',studyID,'.xls'),Labels,SheetInfoFileName1,s_labels);
               
                
            elseif CCorMLO==2; % This is a MLO case
                
                %  Let's remove most of all the non-breast structures
                % in the MLO image.
                % First the 'empty image receptor' part of the Image:
                removenonbreastmlo % Output mask is 'SegmBreast'
                %figure,imagesc(SegmBreast),impixelinfo
                
                %Then we remove the bottom part of the breast (skin folds):
                removefolds % Output mask is 'SegmBreastNoBott'.
                %figure,imagesc(SegmBreastNoBott),impixelinfo
                
                SegmBreastInitial = im2bw(SegmBreast.*SegmBreastNoBott);
                
                % Getting Point A (BB8, p74-75  - BB8 refers to lab notebook#8)
                forPointA= SegmBreastInitial(1:20,:);
                forPointA = sum(forPointA,1);
                kPointA = min(find(forPointA==0));
                
                % and Point B (BB8, p74-75)
                forPointB = SegmBreastInitial(:,:);
                forPointB = sum(forPointB,2);
                offSetPointB = 400; % look at the bottom half
                kPointB = find(forPointB(offSetPointB:size(forPointB,1),:)==0);
                kPointB = min(kPointB)+offSetPointB;
                
                lastTriangleSide = m_Image/kPointB;
                smallerTriangleTop=n_Image/kPointA;
                
                % Then we also get rid of the chest wall
                removechestwall;
                % figure,imshow(SegmBreastNoWall),impixelinfo
                
                % Now we remove the pectoral muscle
                display('USING ONLY ANATOMICAL TRIANGLE FOR NOW')
                removepectoral
                %figure, imagesc(Triangle), impixelinfo
                
                % and combine all the masks:
                SegmBreast = SegmBreastInitial.*SegmBreastNoWall.*Triangle;
                % figure, imagesc(SegmBreast), impixelinfo
                
                
                % ACTUALIZAR EL VALOR DEL CAMPO 'SegmBreast' en InStruct a
                % el valor real de SegmBreast que acabamos de encontrar
                
                
                
                % Now we call the MIB function that will get the values of the HLB
                % layers that correspond to the MLO Image. Output is the SegmBreast mask and the HLB values
                MIBraw;
                glandMLO=glandtot;
                HLB_layersMLO = HLB_layers;
                IdistMLO = Idist;
                wholeGlanDiffpercMLO = wholeGlanDiffperc;
                wholeGlanDiffMLO = wholeGlanDiff;
                
                % And finally, save the HLB layer info to the HLBlayers-spreadsheet
                ResultsMLO = horzcat(HLB_layersMLO, IdistMLO, theorGland, wholeGlanDiffpercMLO, wholeGlanDiffMLO);
                
                s_resultsmlo=strcat('B',num2str(pat+1),':I',num2str(pat+1));
                xlswrite(strcat(StudyBaseCamp,'MIB Outputs\',ss1,'_1-HLBlayers_',studyID,'.xls'),ResultsMLO,SheetInfoFileName2,s_resultsmlo);
                
                s_pat = strcat('A',num2str(pat+1)); % Cell in the spreadsheet
                xlswrite(strcat(StudyBaseCamp,'MIB Outputs\',ss1,'_1-HLBlayers_',studyID,'.xls'),{patientID},SheetInfoFileName2,s_pat);
                
                s_labels='A1:I1';
                Labels = {'CaseID' 'HLB whole gland' 'lay1' 'lay2' 'lay3' 'Idist' 'theorGland' 'glandDiff(%)' 'glandDiff(pp)'};
                xlswrite(strcat(StudyBaseCamp,'MIB Outputs\',ss1,'_1-HLBlayers_',studyID,'.xls'),Labels,SheetInfoFileName2,s_labels);
                
                
                %% Test if tolerance is met:
                diffCCMLOperc = 100*abs((glandCC-glandMLO)/(((glandCC+glandMLO)/2)));
                diffCCMLO = glandCC-glandMLO;
                
                if diffCCMLOperc<diffCCMLOtol
                    display('Anatomical triangle was OK')
                    
                    s_labels='J1:L1';
                    Labels = {'diffCCMLO(%)' 'diffCCMLO(pp)' 'Anatomical triangle used?'};
                    xlswrite(strcat(StudyBaseCamp,'MIB Outputs\',ss1,'_1-HLBlayers_',studyID,'.xls'),Labels,SheetInfoFileName2,s_labels);
                    s_resultstol=strcat('J',num2str(pat+1),':K',num2str(pat+1));
                    xlswrite(strcat(StudyBaseCamp,'MIB Outputs\',ss1,'_1-HLBlayers_',studyID,'.xls'), [diffCCMLOperc diffCCMLO] , SheetInfoFileName2, s_resultstol);
                    s_r=strcat('L',num2str(pat+1));
                    xlswrite(strcat(StudyBaseCamp,'MIB Outputs\',ss1,'_1-HLBlayers_',studyID,'.xls'), {'yes'} , SheetInfoFileName2, s_r);
                 
                    continue
                else
                    display('Anatomical triangle was not good enough')
                    % MUST ADD THE BIT WHERE THE TRIANGLE MASK IS
                    % CALCULATED ITERATIVELY AND MIBRAW IS APPLIED FOR EACH
                    % TRIANGLE
                    stophere
                end
            end
            
                     
        else % statement from jpeg if-loop
            display('There is no JPEG for that case')
            continue
        end % End of jpeg if-loop
        fclose('all');
        
    end % End of CC or MLO loop
end
toc
fclose('all');

pepTalk % Just for fun :)s