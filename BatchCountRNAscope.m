%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: lam94@cam.ac.uk, 2024
% Detect bright spots. Count them inside the mask of cells obtained via stardist
% Output: 
%   files+Summary - mean counts per cell
%   spot tables
% ! The label files always start with Lbl+"name of the image file.tif"
% ! DAPI is always the last channel
% Version for 3 or less channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath ..\bfmatlab\
% Parameters:

nrch = 2;               % numbers of channels. Last - nuclei
halfsize_est = 2;       % gaussian blur width
thInt = 10;             % threshold of spot (max) intensity 
enlarge_param = 2+1;    % Parameter for mask dilation. +1 due to rangefilter. In pixels.
thresh_ch = [20 20 20]; % Spot number per channel for the cell to be counted as positive.
skip = 0;               % skip the first slices

%% Find all label images (tif stack) in this folder. The original images should be here as well.

folder = "H:\Projects\Shel\New\CCN_AMYG_DMN_2"
%% Result folder
%resfolder = fullfile(folder,strcat('Res20x20_z0Th',num2str(thInt)));
resfolder = fullfile(folder,strcat('R_',num2str(thInt)));
mkdir(resfolder)
%% detect all label images to be analysed
files=dir(fullfile(folder,'Lbl*.tif'));

ct = 1;
allSummary = [];
clear nm
for q=1:length(files)
    disp(strcat(num2str(q), '/', num2str(length(files))))
    fn = fullfile(files(q).folder,  files(q).name); % label image
    flif = fullfile(files(q).folder, strrep(files(q).name, '.tif', '.lif'));
    tic
    data = OpenImage(strrep(files(q).name(8:end), '.tif', '.czi'), files(q).folder);
    fnres = fullfile(resfolder,files(q).name);
    if exist('data')
        T = cell(nrch-1,1);
        Lbl = tiffreadVolume(fn);
        DAPI = double(max(Lbl(:)));
        
        % dilate the nuclei masks
        for z = 1:size(data.A{1,1},3)
            % dilate cell in 2d
            aux = Lbl(:,:,z); %aux = -1*(Lbl(:,:,z)>0);
            bwper = (rangefilt(Lbl(:,:,z))>0);
            aux(bwper) =0;   
            d = bwdist(aux>0);
            L1=watershed(d);            
            enlarge = bwdist(Lbl(:,:,z)>0);
            L2 = bwlabel((enlarge<=enlarge_param).*double(L1>0));
            stat = regionprops(L2,Lbl(:,:,z), 'PixelIdxList','PixelValues');
            statOrig = regionprops(Lbl(:,:,z), 'PixelIdxList');
            aux = Lbl(:,:,z);
            Res = zeros(size(aux));
            for x = 1:size(stat,1)
                if median(stat(x).PixelValues)== 0
                    s = stat(x).PixelValues;
                    s(s==0) =[];
                    Res(stat(x).PixelIdxList)= median(s);
                else
                    Res(stat(x).PixelIdxList)= median(floor(stat(x).PixelValues));
                end
            end 
            Lbl(:,:,z) = Res;
            if (z==1)
                imwrite(uint16(Lbl(:,:,z)),strrep(fn, 'Lbl', 'DilLbl'));
            else    
                imwrite(uint16(Lbl(:,:,z)),strrep(fn, 'Lbl', 'DilLbl'),'WriteMode','append')
            end
        end
        clear data.D data.Det
        for ch = 1:nrch-1
            allspots = [];
            % skip first two slices
            for z = skip+1:size(data.A{ch,1},3)
                [data.D{ch}(:,:,z) data.Det{ch}(:,:,z) B1 n1 w1] = FindPeakWav(data.A{ch,1}(:,:,z), 0.01, 0, 2, 0, 1,0);                
                aux = Lbl(:,:,z);
                data.mx{ch}(:,:,z) = imregionalmax(data.A{ch}(:,:,z)).*data.Det{ch}(:,:,z);
                % subpixel fitting not necessary. To be removed
                spotlist = spot2dCurvefit(data.A{ch,1}(:,:,z), data.mx{ch}(:,:,z), 1);    
                if ~isempty(spotlist)
                    id = sub2ind(size(Lbl(:,:,1)),min(size(Lbl,2),max(1,round(spotlist(:,2)))), min(size(Lbl,1),max(1,round(spotlist(:,1)))));
                    % z x y,
                    allspots = [allspots; z*ones(size(spotlist,1),1) spotlist aux(id)];
                end
            end

            stat1 = regionprops3(Lbl,data.D{ch}-B1,'MeanIntensity', 'Volume');
            stat2 = regionprops3(Lbl,double(data.mx{ch}),'MeanIntensity', 'Volume');          
            ii = find(allspots(:,5)>thInt);            
            Tspot = array2table(allspots(ii,[1 2 3 5 12 13]), 'VariableNames', {'z', 'x', 'y','Ampl','exitflag', 'L'});
            % write all detected spot information
            writetable(Tspot,  strrep(fnres, '.tif', strcat('_spots_ch',num2str(ch), '.txt')));            
            statL = grpstats(Tspot, 'L', 'mean');
            st = [(1:size(stat1,1))' 0*(1:size(stat1,1))'];
            st(round(statL.L(2:end)),2) = statL.GroupCount(2:end);            
            T{ch} = table(st(1:end,1), st(1:end,2), stat2.Volume.*stat2.MeanIntensity, stat1.Volume.*stat1.MeanIntensity, DAPI*ones(size(st,1),1),'VariableNames', {'LblId', strcat('Counts','_',num2str(ch)), strcat('CountsFromMean', '_',num2str(ch)), strcat('IntegratedIntensity','_',num2str(ch)), 'NrNuclei'});
         
        end
        % summarise spots per nucleus. I assume 4 channels (modify for
        % nrch)
        if nrch==2
            allT = T{1};
        else
        if nrch == 4
            writetable(join(join(T{1},T{2}),T{3}),  strrep(fnres, '.tif', strcat('_counts.csv')));  
            allT = join(join(T{1},T{2}),T{3});
        else
            allT = join(T{1},T{2});
        end
        end
        id1 = find(allT.Counts_1>=thresh_ch(1));
        meanid1 = mean(allT.Counts_1(id1));
        if nrch>2
            meanintIid1 = mean(allT.IntegratedIntensity_1(id1));
            id2 = find(allT.Counts_2>=thresh_ch(2));
            meanid2 = mean(allT.Counts_2(id2));
            meanintIid2 = mean(allT.IntegratedIntensity_2(id2));        
            if nrch == 4
                id3 = find(allT.Counts_3>=thresh_ch(3));
                meanid3 = mean(allT.Counts_3(id3));
                meanintIid3 = mean(allT.IntegratedIntensity_3(id3));
                id13 = find((allT.Counts_1>=thresh_ch(1))&(allT.Counts_3>=thresh_ch(3)));
                meanid13 = mean(allT.Counts_1(id13));
                id23 = find((allT.Counts_2>=thresh_ch(2))&(allT.Counts_3>=thresh_ch(3)));
                meanid23 = mean(allT.Counts_2(id23));        
                summaryPos = [DAPI length(id1) length(id2) length(id3)...
                    length(find((allT.Counts_1>=thresh_ch(1))&(allT.Counts_3>=thresh_ch(3)))) length(find((allT.Counts_2>=thresh_ch(2))&(allT.Counts_3>=thresh_ch(3)))) length(find((allT.Counts_3>=thresh_ch(3))&(allT.Counts_1>=thresh_ch(1))&(allT.Counts_2>=thresh_ch(2))))...
                    meanid1, meanid2, meanid3 meanintIid1 meanintIid2 meanintIid3  meanid13  meanid23];
            else
                id12 = find((allT.Counts_1>=thresh_ch(1))&(allT.Counts_2>=thresh_ch(2)));
                meanid12 = mean(allT.Counts_1(id12));
                summaryPos = [DAPI length(id1) length(id2) ...
                    length(find((allT.Counts_1>=thresh_ch(1))&(allT.Counts_2>=thresh_ch(2)))) ...
                    meanid1, meanid2, meanintIid1 meanintIid2 meanid12];

            end
        else
             summaryPos = [DAPI length(id1)  meanid1];
   
        end
        allSummary = [allSummary; summaryPos];    
        nm{ct} = files(q).name;
        ct = ct+1;
    end
    toc
    
end
% create overall summary
atab = array2table(allSummary);
if nrch == 2
    atab.Properties.VariableNames ={'TotalNuclei', 'PosC1', 'MeanPos1'}

else
if nrch == 4
    atab.Properties.VariableNames ={'TotalNuclei', 'PosC1','PosC2','PosC3',...
        'PosC12','PosC23','PosC123','MeanPos1','MeanPos2', 'MeanPos3', 'MeanIntPos1', 'MeanIntPos2', 'MeanIntPos3', 'MeanPos1wPos3','MeanPos2wPos3'}
else
    atab.Properties.VariableNames ={'TotalNuclei', 'PosC1','PosC2',...
        'PosC12', 'MeanPos1','MeanPos2', 'MeanIntPos1', 'MeanIntPos2', 'MeanPos1wPos2'}

end
end
atab = addvars(atab, nm', 'NewVariableNames', 'File');
writetable( atab, fullfile(resfolder, 'Summary.csv'));
%%

