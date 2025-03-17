%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single Molecule detection via a trous wavelet method. Thresholding based on FDR. 
%
% Developped for:
% Leila Muresan, Jaroslaw Jacak, Erich Peter Klement, Jan Hesse and
% Gerhard J. Schutz - Microarray analysis at single molecule resolution,
% IEEE Trans. on Nanobioscience, accepted
%
% Based on:
% J.-C. Olivo-Marin. Extraction of spots in biological images using multiscale
% products. Pattern Recognition, 35:1989?1996, 2002.
%
% F. Abramovich and Y. Benjamini. Adaptive thresholding of wavelet coe?-
% cients. Computational Statistics and Data Analysis, 22:351?361, 1996.
%
% F. Abramovich, Y. Benjamini, D. Donoho, and I.M. Johnstone. Adapting
% to unknown sparsity by controlling the false discovery rate. The Annals of
% Statistics, 34(2):584?653, 2006.
%
% J. Boulanger, J.-B. Sibarita, Ch. Kervrann, and P. Bouthemy. Non-
% parametric regression for patch-based ?uorescence microscopy image se-
% quence denoising. In Proc. IEEE Int. Symp. on Biomedical Imaging: from
% nano to macro (ISBI), 2008.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [Denoised Res Backgr] = FindPeakWav(image, q, imageon, levelNo,ignoreFirst, gatP, BHpower, rSize);
%
% Beta (and sloppy) version 
% Created by L. Muresan 10.2008
% 
% image - input image
% q - rate of false discovery
% imageon - if intermediate figures displayed =1, otherwise 0
% levelNo - number of scales computed in the wavlet transform (no. of detail levels)
% ignoreFirst - the scales used for detection range from: ignoreFirst+1 to levelNo
% gatP - parameter that indicates if gat is applied in advance (= 0), or
% for each wavelet separately (= 1) or not applied (= -1, default value)
% Results: 
% Denoised - denoised image with removed background
% Res - binary image, support of detected molecules
% Backgr - estimation of the background
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Denoised Res Backgr noiseStd wave] = FindPeakWav(image, q, imageon, levelNo, ignoreFirst, gatP, BHpower, rSize);
openth = 0;
if nargin < 5
    ignoreFirst = 1;
end

A = double(image);
if(~exist('gatP'))
    gatP = -1;
end
if(~exist('BHpower'))
    BHpower = 1;
end
if(~exist('rSize'))
    rSize = 3;
end
if gatP > 0
    %[p1 p2 A] = gatB(A);
    A = 2*sqrt(A+3/8);
    %A = gat(A);
end
i0=[];
Kernel = [  1/256 1/64 3/128 1/64 1/256; ...
            1/64 1/16 3/32 1/16 1/64;...
            3/128 3/32 9/64 3/32 3/128; 
            1/64 1/16 3/32 1/16 1/64;...
            1/256 1/64 3/128 1/64 1/256 ];
kerSq=sqrt(sum(sum(Kernel.^2)));
normKernel = ones(size(Kernel));
m = 50;

%% Wavelet transform computation
wave = {};
CA = {};
CA{1} = double(padarray(A, [m m], 'symmetric'));

denoiseK = 0;

for k = 1:levelNo
    OldK = Kernel;
    CA{k+1} = conv2(CA{k}, Kernel, 'same');
    % Non Fisz
    if gatP ==0 %VS a la Zhang
        c = 0.0177;
       % wave{k} = sqrt(CA{k}+c) - sqrt(CA{k+1}+c);
        wave{k} = (CA{k}-CA{k+1})./sqrt(CA{k}+CA{k+1});
    else
        wave{k} = CA{k} -CA{k+1};
    end
    
    z = zeros(size(Kernel,1),1);
    
    NKernel = Kernel(:,1);
    for i = 2:size(Kernel,2)
        NKernel = [NKernel z Kernel(:,i)];
    end
    Kernel = NKernel;
    z = zeros(1,size(NKernel,2));    
    NKernel = Kernel(1,:);
    for i = 2:size(Kernel,1)
        NKernel = [NKernel; z; Kernel(i,:)];
    end
    Kernel = NKernel;
    normKernel = ones(size(Kernel));  %% ???
    if imageon
        figure;
        subplot(221)
        imagesc(CA{k}); title('Convolved images')

        subplot(222)
        imagesc(wave{k}); title('Wavelet');
        
        subplot(223)
        hist(CA{k}(:), 100); 

        subplot(224)
        hist(wave{k}(:), 100); 
        
        pause        
    end

end

Backgr = CA{levelNo+1}(m+1:size(CA{levelNo+1},1)-m, m+1:size(CA{levelNo+1},2)-m);

%% Spot detection
maxC = 0;
% scale-noise estimation
Denoised = zeros(size(image));
Res = ones(size(image));
for k=1:levelNo 
    wave{k} = wave{k}(m+1:size(wave{k},1)-m, m+1:size(wave{k},2)-m);
    maxC = max( maxC, max (max( abs(wave{k}) )));
    
%     if gatP>0
%         [p1 p2 wave{k} ] = gatB(wave{k} );
%     end

    % test gaussianity assumption
    vals = wave{k}(:);
   
    if k == 1
        noiseStd(k) = median(abs(vals-median(vals))/0.6745);
    else
        if noiseStd(k-1) == 0
             noiseStd(k) = median(abs(vals-median(vals))/0.6745);
        else
            noiseStd(k) = kerSq*noiseStd(k-1);
        end
    end    
    % for gaussian noise
    p =2*(1 - cdf('Normal',abs(vals),0,noiseStd(k)));
     
    n = length(vals);
    [ps idx] = sort(p);
    if BHpower
        s = (1-ps)./(n+1-(1:n)');
        id  = find(diff(s)<0);
        if ~isempty(id)
            mhat = min(floor( 1/(s(id(1)+1))), n);
            i0 =  max(find(ps <(1:n)'./mhat*q));  
            i0comp =  max(find(ps < (1:n)'./n*q));  
        end
    else                
        i0 =  max(find(ps < (1:n)'./n*q));  
    end
    if imageon
        figure; 
        plot(ps, 'r.-'); hold on; title (num2str(i0))
        plot((1:n)'./n*q, 'b*-'); hold on
        pause
    end
    
    % Modify iff found i0!!!!
    if ~isempty(i0)
        wave{k}(idx(i0+1:end)) = 0;        
    end
    
    % for poisson noise
    % p = cdf('Poisson',vals,1:6)
    
    %%%% !!!  
    wave{k}( wave{k}<0) = 0;
    if k >= ignoreFirst+1;  % To be able to ignore (noise) scales
        Denoised = Denoised +wave{k};            
        if ~isempty(i0)
            Res = Res.* max(0,wave{k}); 
        else
            Res = 0;
        end
    end
  

end

% To add the background!!! I removed it

% put back in comment, just for TIE
Denoised = Denoised +CA{levelNo+1}(m+1:size(CA{levelNo+1},1)-m, m+1:size(CA{levelNo+1},2)-m);
%Denoised(find(Denoised<0))= 0;


Res(find(Res > 0)) = 1;
Denoised(find(Denoised<0))= 0;

%
Res  = bwareaopen(gather(Res),openth);
% inverse Anscombe transform
if gatP > 0
    Denoised = 1/4*(Denoised.*Denoised-8/3);
end
clear CA
%clear wave
clear A

% %
% figure;
% l = 5
% allV = Denoised(1:37, l*37+1:l*37+37);
% for k = 1:16
%    allV =  allV+Denoised(k*37+1:k*37+37, l*17+1:l*17+37);
%    clf
%    imagesc(Denoised(k*37+1:k*37+37, l*17+1:l*17+37)-Denoised((k-1)*37+1:(k-1)*37+37, l*17+1:l*17+37))
%    pause
% end

if imageon    
    figure; 
    subplot(1,3,1); 
    imagesc(image); colormap gray; title('Original'); axis equal tight
    %imagesc(image); title('Original')    
    subplot(1,3,2); 
    %imagesc(Denoised); colormap gray; title('Denoised'); axis equal tight
    imagesc(Backgr); colormap gray; title('Denoised'); axis equal tight
    %imagesc(Denoised); title('Denoised')
    subplot(1,3,3); 
    imshow(Res, []); title('Binary')

    figure;
    hist(Backgr(:), 100); title('Background histogram')
    pause    
end