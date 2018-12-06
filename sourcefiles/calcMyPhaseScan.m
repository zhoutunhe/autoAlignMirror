function [radius,elwavefront,elwgradient,ratio,delay] = calcMyPhaseScan(pixsize,dist,...
    ROI1,path2,track2D,step,maxImages,rotim,ScanDirection,a,stepSize,nstep,optdisp,format)
% [radius,elwavefront,elwgradient] = calcMyPhaseScan(pixsize,dist,ROI1,path2,track2D,step,maxImages)
%       Calculate phase map from speckle-scanning.
%       If no input of track2D, track2D  = 0 to have 1D output
%       If no input of step size, step = 0.25 um. 
%       If no input of maximum Image number, maxImages = 500;
%       
%       Reference: Wang, Opt. Express 23(2), 1605 (2015)
%
% See also: calcMySliceDelayX

%   by S. Berujon April 2012
%   modified by tunhe zhou 2018-08-17

if ~exist('track2D','var') || isempty(track2D)
    track2D = 0; %Default is tracking only 1D
end
if ~exist('step','var') || isempty(step)
    step = .25; % scan 56step in um
end
if ~exist('maxImages','var') || isempty(maxImages)
    maxImages =201;% 
end
if ~exist('rotim','var') || isempty(rotim)
    rotim =0;% 
end
if ~exist('ScanDirection','var') || isempty(ScanDirection)
    ScanDirection =-1;% 
end
if ~exist('a','var') || isempty(a)
    a =100;% 
end
if ~exist('stepSize','var') || isempty(stepSize)
    stepSize =a/2;% 
end
if ~exist('nstep','var') || isempty(nstep)
    nstep = 2;% important (3 or 2 are best)  
% The distance between the rows used for correlation analysis
% (j-i) in Wang Opt. Express 23 (2015);% 
end
if ~exist('optdisp','var') || isempty(optdisp)
    optdisp = 0; %Default is tracking only 1D
end
if ~exist('format','var') || isempty(format)
    format = 'TIF*'; %file extension
end




%----------------------------------------------------------------------
filtersz = 0; % 15;% averaging filter size
myleefilter = 10000; % It is actually porphological processing
undersampl = 0; % undersampling factor, 0 for no



% -------------------------dark field and flat field images--------------%
darkfieldImagePath = [];
darkfieldImageName =  [];

FlatfieldImagePath =  [];
FlatfieldImageName1 =  [];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Begin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ischar(path2), return;else path1 = path2;end;
display(['Path: ' path2]);
  fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    open the files from the path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n Starting loading pictures\n');


[files] = open_seq_dia(path2,ROI1,format);

if ~isempty(darkfieldImagePath) && exist(darkfieldImagePath,'dir')
    darkim1 = single(imread(fullfile(darkfieldImagePath,darkfieldImageName),'tif'))-1;
    darkim1 = darkim1(ROI1(1):ROI1(2),ROI1(3):ROI1(4));
else
    darkim1 = zeros(size(files(1).data));
end;

if ~isempty(FlatfieldImagePath) && exist(FlatfieldImagePath,'dir')
    flatim1 = single(imread(fullfile(FlatfieldImagePath,FlatfieldImageName1),'tif'));
    flatim1 = flatim1(ROI1(1):ROI1(2),ROI1(3):ROI1(4));
    flatim1 = flatim1- darkim1;
else  flatim1 = files(1).data.*0+1;
end;
    

disp('  DONE');

% error checking
if ~isstruct(files) || isempty(files(1)), error('No 1pictures in memory: path probably not correct');end;

nImages = length(files);
if nImages > maxImages, nImages = maxImages;    end
disp(['Number of pictures = ' num2str(nImages)])

disp('     ');



if (myleefilter > 0),  files(1).data = calcMyerosion(files(1).data,myleefilter);end;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                    filter
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1 : 1 : nImages, 

    if myleefilter > 0, files(k).data = calcMyerosion(files(k).data, myleefilter); end;
    if undersampl > 0,files(k).data = files(k).data(1:undersampl:end,1:undersampl:end); end;
    if filtersz > 0 ,files(k).data = imfilter(files(k).data,ones(filtersz)./(filtersz).^2,'replicate','same');end;
    files(k).data = single(files(k).data) - darkim1;
    if ~isempty(FlatfieldImagePath), files(k).data = single(files(k).data)./single(flatim1);end;
    if rotim ~= 0, files(k).data = rot90(files(k).data,rotim );end;
end;
if undersampl > 0, pixsize = pixsize *undersampl;end;
% 

[m1, n1] = size(files(1).data);
 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           build stack                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stack= zeros(m1 , n1, nImages );
for pp = 1 : 1 : nImages, 
        % make a stack of sample images
        stack(:,:,pp) = files(pp).data;
 end;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            calculate the lcoal radius of curvature        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear delay delay11 delaypix delay1 elwgradient      
if track2D

    x0v = a:stepSize:n1-a;
    delay = NaN(m1-nstep,length(x0v));  %delay in x direction
    
    
    for pp = 1 :1: m1-nstep
    
        [delay(pp,:)] = ...
            calcMySliceDelay2D(squeeze(stack(pp,:,:)), ...
            squeeze(stack(pp+nstep,:,:)),a,x0v,optdisp);
        
    end
    ratio = abs(-(delay.*step/nstep/pixsize)+ScanDirection);
    delaypix = ratio./(dist*1e-3);% [m-1], Eq. (5) in reference

    
else
% if 1D tracking
    delay = zeros(m1-nstep,2);
    
    for pp = 1 :1: m1-nstep
        [delay(pp,:),~,~,~] = ...
            calcMySliceDelayX(squeeze(stack(pp,:,:)), ...
            squeeze(stack(pp+nstep,:,:)),optdisp);
    end
    
    
    ratio = abs(-(delay(:,2).*step/nstep/pixsize)+ScanDirection);
    delaypix = ratio./(dist*1e-3);% [m-1], Eq. (5) in reference


end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration to Phase
radius = 1./delaypix;   %[m]

secderi = delaypix; 
elwgradient = cumsum(secderi).*pixsize ;%[urad]
elwavefront = -cumsum(elwgradient).*pixsize.*1e-3; %[nm]

