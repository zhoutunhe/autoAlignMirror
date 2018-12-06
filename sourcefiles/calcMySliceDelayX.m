function [delay,coeff,pic3,YI] = calcMySliceDelayX(plane2, plane1,optdisp)
%[delay,coeff,pic3,YI] = calcMySliceDelayX(plane2, plane1)
%     Calculate the shift between plane2 (image to search in) and 
%      plane1 (template) using normalized 2d cross correlation. Subpixel
%      accuracy is achieved by using findpeak_max modified from findpeak from
%      MathWorks.
%
%      The margin of the template in plane1 can be set in this function. 
%      The default is edge = 27 in x direction and edgeV = 5 in y. 
%
%See also: calcMyPhaseScan, normxcorr2, findpeak_max,
if ~exist('optdisp','var') || isempty(optdisp)
    optdisp = 0; % to activate display of the peak
end
[m, n] = size(plane1);

optproc = 2;% 1 for linux 2- for windows 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% edge =27;%important for the offset
edge = 10;
edgeV =10;%
plane1t = plane1(edgeV:end-edgeV,edge:end-edge);


if optproc == 1
    [crosscorr3,pic3] = normxcorr2_mexsub(plane1t,plane2,'valid');
    xpeak = pic3(3); ypeak = pic3(4);
    delay = [xpeak ypeak];
    coeff  = pic3(5);
    crosscorr = crosscorr3;
end;


if optproc == 2
   crosscorr3 = normxcorr2(plane1t,plane2);
   [xpeak, ypeak, coeff] = findpeak_max(crosscorr3,true);
   corr_offset = [ (ypeak-size(plane1t,1)-edgeV+1)  (xpeak-size(plane1t,2)-edge+1)];
   delay = corr_offset;
   pic3 = [corr_offset(1) corr_offset(2) corr_offset(1) corr_offset(2) coeff];
   crosscorr = crosscorr3(m-2*edgeV:m,n-2*edge:n);
end;
[YI,CI] = max2(crosscorr3);

if optdisp
    figure(2)
    subplot(2,2,1)
    imagesc(plane1t)
    subplot(2,2,2)
    imagesc(plane2)
    subplot(2,2,3)
    imagesc(crosscorr3)
    subplot(2,2,4)
    
    plot(crosscorr3(CI(1),:))
    drawnow
end;






