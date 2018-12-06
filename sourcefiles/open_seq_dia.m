function [files] = open_seq_dia(path2, ROI,format)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author S.Berujon
%OUTPUT
%files : cell structure files(k).header
%                  files(k).data 
% INPUT
%ROI=[down top left right]
%path= string directory path

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% slect a path if doesn't exist
if nargin == 0
    path2 = uigetdir('~','Folder picking');
end;
if ~exist('format','var') || isempty(format)
    format = 'TIF*'; %file extension
end
%% check if path is good
if ~ischar(path2)
    disp('Cancellation: var path is not correct...')
    return
elseif exist(path2,'dir')==0
    disp('Cancellation: path is not a correct directory...')
    return
end;

%% search for the files to store
path1 = [path2 filesep '*.' format];
file = dir(path1);
nImages = length(file);

%% limit to a stack of 50 inages

files = cell(nImages,3);
file = struct2cell(file);
file = file(1,:)';                       
file = sort_nat(file);
file_tmp = imread(fullfile(path2,file{1,1})); % open a sample picture (1 one)

%% Check ROI if exist
if nargin >1 
    if ~isequal(size(ROI),[1 4]) % bad format of ROI
        disp('Uncorrect ROI: cancellation')
        return
    elseif (ROI(1) >= ROI(2)) || (ROI(3) >= ROI(4)) % bad order 
        disp('Uncorrect ROI: cancellation')
        return;
    elseif ROI(2) >= size(file_tmp,1) || ROI(3) >= size(file_tmp,2) % ROI exeed image size
        disp('Uncorrect ROI: cancellation');
        return;
    end;
end;

%% Open and store the images
clear files
for k = 1 : nImages
    files(k).name = fullfile(path2,file{k,1});    
    file_tmp = imread(files(k).name);
    if nargin >1
        files(k).data = uint16(file_tmp(ROI(1):ROI(2),ROI(3):ROI(4)));
    else
        files(k).data = uint16(file_tmp(:,:));
    end;
    files(k).header = [];

end;