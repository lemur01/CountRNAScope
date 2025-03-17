 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LM, June 2013
% 2d gaussian fit
%
% Im - fluorescent image
% mask - detected single molecules (wuth FindPeakWav)
% nrgauss - how many gaussian to fit (1 or 2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function spotlist = spot2dCurvefit(Im, mask, nrgauss)
    %if isempty(spotlist); spotlist2 = []; return; end
    MdataSize = 6;
    %s = regionprops(mask, Im, {'Centroid','WeightedCentroid'});
    warning('off')
    
%     if nrgauss == 1
%   %      mask = bwareaopen(mask,2);
%         w = 3;
%        % L = bwlabel(mask);
%     else
%         
%   %      mask = bwareaopen(mask,9);
%       %  L = bwlabel(mask);
%         w = 9;
%     end
    L = bwlabel(mask);
    s = regionprops(L, Im, {'MinIntensity', 'MaxIntensity'});
    s1 = regionprops(L, 'Centroid');
    if isempty(s); spotlist = []; return; end
    centroids = cat(1, s1.Centroid);
    spotlist = zeros(size(centroids, 1), 11);
    spotlist(:,1:2) =centroids;
    spotlist(:,3) = 0;
    spotlist(:,4) = cat(1, s.MaxIntensity);
    spotlist(:,5) = cat(1, s.MinIntensity);
    spotlist(:,9:10) = spotlist(:,1:2);
    spotlist(:,11) = -25;
    
    wz = 3;
    % pos3 z!!!
    % tolerance
    tolerance = 1e-3;
    
    % maximum number of iterations
    max_n_iterations = 20;
    
    % estimator id
    estimator_id = EstimatorID.MLE;
    
    % model ID
    model_id = ModelID.GAUSS_2D;
    Ly = size(spotlist,1);
    for sp = 1:Ly
        
        d = Im(max(1,round(spotlist(sp,2))-wz):min(round(spotlist(sp,2))+wz,size(Im,1)), max(1,round(spotlist(sp,1))-wz):min(round(spotlist(sp,1))+wz,size(Im,2)));
        initial_parameters = single([spotlist(sp,4)  spotlist(sp,1)-round(spotlist(sp,1))+wz spotlist(sp,2)-round(spotlist(sp,2))+wz 2 spotlist(sp,5)]');
        %% run Gpufit
        [parameters, states, chi_squares, n_iterations, time] = gpufit(single(d(:)), [], ...
            model_id, initial_parameters, tolerance, max_n_iterations, [], estimator_id, []);

        %% displaying results
        %display_results('2D rotated Gaussian peak', model_id, 1, number_parameters, size_x, time, true_parameters, parameters, states, chi_squares, n_iterations);
    end



        
              %options = optimset('Display','on','MaxIter',1000, 'Diagnostics', 'on','PlotFcns', @optimplotx);
% %               if nrgauss == 1 % don't include 2 points in the window            
% %                 x0 =[1 spotlist(sp,1) 2 spotlist(sp,2)]; % 2 ];
% %                 x =x0;
% %                 [X,Y] = meshgrid(-MdataSize/2:MdataSize/2);
% %                 xm = X+floor(spotlist(sp,1))-0.5;
% %                 ym = Y+floor(spotlist(sp,2))-0.5;
% %                 rxdata = [xm(:) ym(:)];
% %                 ii = find((ym>size(Im,1))|(ym<1)|(xm>size(Im,2))|(xm<1));
% %                 exitflag = -10;
% %                 rxdata(ii,:) = [];
% %                 try
% %                     Z = Im(sub2ind(size(Im),  floor(rxdata(:,2)),  floor(rxdata(:,1))));
% %                     Z = Z(:);                     
% %                     lb = [0,spotlist(sp,1)-MdataSize/2,0,spotlist(sp,2)-MdataSize/2]; %,0]; 
% %                     ub = [realmax('double'),spotlist(sp,1)+MdataSize/2,(MdataSize/2)^2,spotlist(sp,2)+MdataSize/2]; %,(MdataSize/2)^2];
% %                     %options = optimset('Display','on','MaxIter',1000, 'Diagnostics', 'on');%,'PlotFcns', @optimplotx);
% %                     options = optimset('MaxIter',200,'MaxFunEvals', 200);
% %                    
% %                     [x,resnorm,residual,exitflag] = lsqcurvefit(@D2GaussFunction,x0,rxdata,Z,lb,ub, options); 
% %                     if ~isnan(x(2))
% %                         spotlist(sp,1) = x(2); % x
% %                         spotlist(sp,2) = x(4); % y
% %                         spotlist(sp,11) = exitflag; % exit -1 / 0 / 1
% %                         spotlist(sp,6) = x(3); % sx
% %                         %spotlist(sp,7) = x(5); % sy
% %                         spotlist(sp,3) = x(1); % A
% %                     end
% %                 catch
% %                     floor(ym(1));
% %                     floor(xm(1));
% %                 end                    
% %                
% %                 
% %                 %spotlist(sp,1:6) = dat(1:6);
% %                 % fields 5 & 6 are not filled in second pass
% % 
% %                 aux = Im;
% %                 
% %                 %plot(dat(1), dat(2), 'mo')
% %               else
% %                 MdataSize = 12;
% %                 aux = imdilate((L==sp),strel('disk',11)).*Im;
% %                 %figure; imagesc(aux); hold on
% %                 [cy cx]=find((sp==L) >0);
% %                 idx = find(imregionalmax(aux,8)>0);
% %                 [val pos] =  sort(aux(idx), 'descend');
% %                 [yy xx]= ind2sub(size(aux),idx(pos(1:2)));
% %                 %plot(xx, yy, 'k.')
% %                 % x, y, s, A, b, x, y, s, A
% %                 % initialize with the two brightest local maxima
% %                 %dat = double([spotlist(sp,1)+1, spotlist(sp,2)+1, 1.5,  spotlist(sp, 4)-spotlist(sp, 5)+0.5 spotlist(sp, 5) ...
% %                 %    spotlist(sp,1)-1, spotlist(sp,2)-1, 1.5,  spotlist(sp, 4)-spotlist(sp, 5)+0.5 ]);
% %                 x0 = double([ val(1) xx(1) 1.5 yy(1) 1.5 val(2) xx(2) 1.5 yy(2) 1.5 ]);
% %                 
% %                 % x0 =[1 0 2 0 2 ];
% %                 [X,Y] = meshgrid(-MdataSize/2:MdataSize/2);
% %                 xdata = zeros(size(X,1),size(X,2),2);
% %                 xm = X+sum(xx)/2; %+0.5;
% %                 ym = Y+sum(yy)/2; %+0.5;
% %                 
% %                 
% %                 rxdata = [xm(:) ym(:)];
% %                 ii = find((ym>size(Im,1))|(ym<1)|(xm>size(Im,2))|(xm<1));
% %                 exitflag = -10;
% %                 rxdata(ii,:) = [];
% %                 try
% %                     Z = Im(sub2ind(size(Im),  floor(rxdata(:,2)),  floor(rxdata(:,1))));
% %                     Z = Z(:);                     
% %                     %options = optimset('Display','on','MaxIter',1000, 'Diagnostics', 'on');%,'PlotFcns', @optimplotx);
% %                     options = optimset('MaxIter',2000,'MaxFunEvals', 10000);
% %                      
% %                      x =x0;
% %                      lb = [0,xx(1)-MdataSize/2,0,yy(1)-MdataSize/2,0 0,xx(2)-MdataSize/2,0,yy(1)-MdataSize/2,0]; 
% %                      ub = [realmax('double'),xx(1)+MdataSize/2,(MdataSize/2)^2,yy(1)+MdataSize/2,(MdataSize/2)^2, realmax('double'),xx(2)+MdataSize/2,(MdataSize/2)^2,yy(2)+MdataSize/2,(MdataSize/2)^2];
% %                 
% %                     [x,resnorm,residual,exitflag] = lsqcurvefit(@D2Gauss2Function,x0,rxdata,Z,lb,ub, options);    
% %                 catch
% %                     floor(ym(1))
% %                     floor(xm(1))
% %                 end
% %                 spotlist(sp,1) = x(2);
% %                 spotlist(sp,2) = x(4); 
% %                 % A, x, y, s, A, x, y, s (b?)
% %                 spotlist(sp,3) = (x(3)+x(5))/2;
% %                 spotlist(sp,6) = x(7);%+ centroids(sp,1); 
% %                 spotlist(sp,7) = x(9);%+ centroids(sp,2); 
% %                 spotlist(sp,8) = (x(8)+x(10))/2;                 
% %                 spotlist(sp,11) = exitflag; % exit -1 / 0 / 1
% %               end

    

end


function g = gaussian_2d(x, y, p)
% Generates a 2D Gaussian peak.
% http://gpufit.readthedocs.io/en/latest/api.html#gauss-2d
%
% x,y - x and y grid position values
% p - parameters (amplitude, x,y center position, width, offset)

g = p(1) * exp(-((x - p(2)).^2 + (y - p(3)).^2) / (2 * p(4)^2)) + p(5);

end
