function [delay] = calcMySliceDelay2D(plane2, plane1,a,x0v,optdisp)
%[delay,coeff,pic3,YI] = calcMySliceDelay2D(plane2, plane1,a,x0v)
%     Calculate the shift between plane2 (image to search in) and 
%      plane1 (template) using normalized 2d cross correlation. Subpixel
%      accuracy is achieved by using findpeak_max modified from findpeak from
%      MathWorks.
%      a: window width
%      x0v: cooridnates ijn x direction of tracking center
%
%      The margin of the template in plane1 can be set in this function. 
%      The default is edge = 27 in x direction and edgeV = 5 in y. 
%
%See also: calcMySliceDelayX, calcMyPhaseScan, normxcorr2, findpeak_max,

[m, n] = size(plane1);

optproc = 2;% 1 for linux 2- for windows 
if ~exist('optdisp','var') || isempty(optdisp)
    optdisp =0;% to activate display of the peak
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

edge =5;%important for the offset
edgeV =5;%

% plane1 = plane1./repmat(sum(plane1,1),[m 1])./repmat(sum(plane1,2),[1 n]);
% plane2 = plane2./repmat(sum(plane2,1),[m 1])./repmat(sum(plane2,2),[1 n]);

% plane1 = (plane1 - mean(plane1(:)))./std(plane1(:));
% plane2 = (plane2 - mean(plane2(:)))./std(plane2(:));

delay = zeros(1, length(x0v));


if optproc == 1
    [crosscorr3,pic3] = normxcorr2_mexsub(plane1t,plane2,'valid');
    xpeak = pic3(3); ypeak = pic3(4);
    delay = [xpeak ypeak];
    coeff  = pic3(5);
    crosscorr = crosscorr3;
end;


if optproc == 2
    parfor j = 1:length(x0v)
        x0 = x0v(j);
        plane1t = plane1(x0-a+1+edgeV:x0-edgeV,edge:end-edge);
        crosscorr3 = normxcorr2(plane1t,plane2(x0-a+1:x0,:));
        [xpeak, ypeak, coeff] = findpeak_max(crosscorr3,true);
        corr_offset = [ (ypeak-size(plane1t,1)-edgeV+1)  (xpeak-size(plane1t,2)-edge+1)];
        delay(1,j) = corr_offset(2);
%         pic3 = [corr_offset(1) corr_offset(2) corr_offset(1) corr_offset(2) coeff];
        
    end
end;
% [YI,CI] = max2(crosscorr3);

if optdisp
    figure(1)
    subplot(2,2,1)
    imagesc(plane1t)
    subplot(2,2,2)
    imagesc(plane2(x0-a+1:x0,:))
    subplot(2,2,3)
    imagesc(crosscorr3)
    subplot(2,2,4)
    
    plot(crosscorr3(CI(1),:))
    drawnow
end;






