function varargout = autoAlignMirror(varargin)
%AUTOALIGNMIRROR MATLAB code file for autoAlignMirror.fig
%      AUTOALIGNMIRROR, by itself, creates a new AUTOALIGNMIRROR or raises the existing
%      singleton*.
%
%      H = AUTOALIGNMIRROR returns the handle to a new AUTOALIGNMIRROR or the handle to
%      the existing singleton*.
%
%      AUTOALIGNMIRROR('Property','Value',...) creates a new AUTOALIGNMIRROR using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to autoAlignMirror_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      AUTOALIGNMIRROR('CALLBACK') and AUTOALIGNMIRROR('CALLBACK',hObject,...) call the
%      local function named CALLBACK in AUTOALIGNMIRROR.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help autoAlignMirror

% Last Modified by GUIDE v2.5 31-Aug-2018 16:58:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @autoAlignMirror_OpeningFcn, ...
                   'gui_OutputFcn',  @autoAlignMirror_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before autoAlignMirror is made visible.
function autoAlignMirror_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

global format;

format = 'TIF';

% Choose default command line output for autoAlignMirror
handles.output = hObject;

%// Initialize flag
handles.StopNow = false;
handles.plotROI = false;
handles.mask = [];
handles.op2 = [];

% Update handles structure
guidata(hObject, handles);

% Display diamond logo on the bottom right of the UI
axes(handles.logo_axes);
imshow('Diamond_logo_col.jpg');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cumstomed functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%
% % preview image %
% %%%%%%%%%%%%%%%%%%%%%%
function previewImage(handles)
    global format
% preview the first image in a scan folder
    fileDir = get(handles.table_filePath,'Data');
    if ispc
        dirs = regexp(genpath(fileDir{1}),['[^;]*'],'match');
    else
        dirs = regexp(genpath(fileDir{1}),['[^:]*'],'match');
    end
    % output a cell with all paths levels
    path1 = [dirs{end} filesep '*.' format];
    file = dir(path1);
    file = struct2cell(file);
    file_tmp = imread(fullfile(dirs{end},file{1,1})); % open a sample picture (1 one)
    imshow(file_tmp,[]);

% %%%%%%%%%%%%%%%%%%%%%%
% % plot profile %
% %%%%%%%%%%%%%%%%%%%%%%
function plotProfile(handles)
    global format
% plot the profile of the image within the roi to adjust setting for ROI
    dirNo = str2num(get(handles.edit_profileNo,'String'));
    horizontal = get(handles.horizontal, 'Value');
    startpix = str2double(get(handles.start_pix,'String'));
    endpix = str2double(get(handles.end_pix,'String'));
    
    fileDir = get(handles.table_filePath,'Data');
    if ispc
        % windows
        dirs = regexp(genpath(fileDir{dirNo}),['[^;]*'],'match');
    else
        % others
        dirs = regexp(genpath(fileDir{dirNo}),['[^:]*'],'match');
    end
        % output a cell with all paths levels
    path1 = [dirs{end} filesep '*.' format];
    file = dir(path1);
    file = struct2cell(file);
    file_tmp = imread(fullfile(dirs{end},file{1,1})); % open a sample picture (1 one)
    
    if ~horizontal
        Imean = mean(file_tmp(:,startpix:endpix),2);
    else
        Imean = mean(file_tmp(startpix:endpix,:),1);
    end
    Imean = Imean./max(Imean(:));
    figure, plot(Imean)
    xlabel('pixel number')
    ylabel('normalized profile within roi')
    

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % define roi %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function defineROI(hObject,handles)
global format
% find all the rois from all the scan positions
% roi = [y1 y2 x1 x2];
horizontal = get(handles.horizontal, 'Value');
startpix = str2double(get(handles.start_pix,'String'));
endpix = str2double(get(handles.end_pix,'String'));
pitchN = length(get(handles.table_filePath,'Data'));
fileDir = get(handles.table_filePath,'Data');
edge = str2double(get(handles.edit_edge,'String'));
threshold = str2double(get(handles.edit_percent,'String'));

handles.roi = zeros(pitchN,4);


if handles.plotROI
    figure
end
for it  = 1:pitchN
    
    if ispc
        % windows
        dirs = regexp(genpath(fileDir{it}),['[^;]*'],'match');
    else
        % others
        dirs = regexp(genpath(fileDir{it}),['[^:]*'],'match');
    end
    path1 = [dirs{end} filesep '*.' format];
    file = dir(path1);
    file = struct2cell(file);
    I = imread(fullfile(dirs{end},file{1,1})); % open a sample picture (1 one)
    if ~horizontal
        Imean = mean(I(:,startpix:endpix),2);
        Imean = Imean./max(Imean(:));
        mirrorFOV = find(Imean>threshold);    % !!adjust according to the profile of the mirror, 0.3 for Lane3
        handles.roi(it,:) = [mirrorFOV(1)+edge,mirrorFOV(end)-edge,startpix,endpix];
    else
        Imean = mean(I(startpix:endpix,:),1);
        Imean = Imean./max(Imean(:));
        mirrorFOV = find(Imean>threshold);    % !!adjust according to the profile of the mirror, 0.3 for Lane3
        handles.roi(it,:) = [startpix,endpix,mirrorFOV(1)+edge,...
            mirrorFOV(end)-edge];
    end
    if handles.plotROI
       subplot(ceil(sqrt(pitchN)),ceil(sqrt(pitchN)),it)
       imshow(I,[])
       warning('off', 'Images:initSize:adjustingMag');
       hold on
       rectangle('Position',[handles.roi(it,3),handles.roi(it,1),...
           handles.roi(it,4)-handles.roi(it,3),...
           handles.roi(it,2)-handles.roi(it,1)],'EdgeColor','r')
       hold on
       title(['roi', num2str(it)])
    end
    hold off
end
guidata(hObject, handles); % Update handles structure
    
    
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % start calculation %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function startCalculation(hObject,handles)
global format
pixsize = str2double(get(handles.edit_pix,'String'));
DistL = str2double(get(handles.edit_MDdist,'String'));
pitchN = length(get(handles.table_filePath,'Data'));
fileDir = get(handles.table_filePath,'Data');
step = str2double(get(handles.scanStep,'String'));
rotim = get(handles.horizontal,'Value')*(-1);   %-1 means rotate clockwise
angles = get(handles.uitable_pitchAngle,'Data');
showProc = get(handles.radiobutton_showXcorr, 'Value');
scanTowards = get(handles.radiobutton_towardsFocus,'Value');
scanAway = get(handles.radiobutton_away,'Value');
ScanDirection = scanTowards*(-1)+scanAway*1;
track2D = get(handles.button_2D,'Value');
window = str2double(get(handles.edit_window,'String'));
stepSize = str2double(get(handles.edit_stepsize,'String'));
scanNo = str2double(get(handles.edit_whichScan,'String'));
parallel = get(handles.radiobutton_parallel, 'Value');
roi = handles.roi;  

% if horizontal, rotate the image

radius = cell(pitchN,1);
Rslope = zeros(pitchN,1);
meanR = zeros(pitchN,1);
if ~track2D

    if parallel

        parfor it  = 1:pitchN
            if ispc
                % windows
                dirs = regexp(genpath(fileDir{it}),['[^;]*'],'match');
            else
                % others
                dirs = regexp(genpath(fileDir{it}),['[^:]*'],'match');
            end
                
                pathName = dirs{end};
                [Radius,~,~] = calcMyPhaseScan(pixsize,...
                    DistL,roi(it,:),pathName,track2D,step,[],rotim,ScanDirection,...
                    [],[],[],showProc,format);
                radius{it} = Radius;
                meanR(it) = median(Radius(:));
                % calculate the slope of the radius
                polyR1 = polyfit([1:length(Radius)]',Radius,1);
                Rslope(it) = polyR1(1);

        end
    else
        for it  = 1:pitchN
            if handles.StopNow == false
                if ispc
                    % windows
                    dirs = regexp(genpath(fileDir{it}),['[^;]*'],'match');
                else
                    % others
                    dirs = regexp(genpath(fileDir{it}),['[^:]*'],'match');
                end
                pathName = dirs{end};
                [Radius,~,~] = calcMyPhaseScan(pixsize,...
                        DistL,roi(it,:),pathName,track2D,step,[],rotim,ScanDirection,...
                        [],[],[],showProc,format);
                radius{it} = Radius;
                meanR(it) = median(Radius(:));
                % calculate the slope of the radius
                polyR1 = polyfit([1:length(Radius)]',Radius,1);
                Rslope(it) = polyR1(1);
                pause(.5)
            else
                msgbox('Process stopped');
                return
            end
        end
    end
else
    if handles.StopNow == false
        if ispc
            % windows
            dirs = regexp(genpath(fileDir{scanNo}),['[^;]*'],'match');
        else
            % others
            dirs = regexp(genpath(fileDir{scanNo}),['[^:]*'],'match');
        end

        pathName = dirs{end};
        [radius,~,~] = calcMyPhaseScan(pixsize,...
            DistL,roi(scanNo,:),pathName,track2D,step,[],rotim,ScanDirection,...
            window,stepSize,[],showProc,format);
        
    else
        msgbox('Process stopped');
        return
    end
end
handles.Rslope = Rslope;
handles.radius = radius;
handles.meanR = meanR;
guidata(hObject, handles); % Update handles structure
% automatic plot radius in GUI
axes(handles.plotR);
cla reset
displayPlotR(handles);
set(gca,'fontsize',8)
set(handles.pushbutton_zoomR, 'Enable', 'on');
% linear fitting 
if ~track2D
    if pitchN == 1% only one scan
        set(handles.text_result,'String',...
            'Error: There is only 1 scan, fitting cannot be done.');
    else
        % linear fit the slope
        polyS = polyfit(angles,Rslope,1);
        handles.slsl = polyval(polyS,angles);
        handles.bestAngle = -polyS(2)/polyS(1);
        
        polyR = polyfit(angles,meanR,1);
        handles.slR = polyval(polyR,angles);
        handles.bestR = polyR(1)*handles.bestAngle+polyR(2);
        
        guidata(hObject, handles); % Update handles structure

        % automatic plot radius slope in GUI
        axes(handles.plotSlope);
        cla reset
        displaySlope(handles);
        set(gca,'fontsize',8)
        set(handles.pushbutton_zoomSlope, 'Enable', 'on');

        % show the best angle result in GUI
        set(handles.text_result,'String',['The best pitch angle is ', num2str(...
            handles.bestAngle), '.', ' Focus is ', num2str(...
            handles.bestR*1e3-DistL), 'mm before the diffuser.']);
    end
else
    
    calcZernike(hObject,handles);
    handles = guidata(hObject); % Update handles structure
    % automatic plot radius slope in GUI
    axes(handles.plotSlope);
    cla reset
    displaySlope(handles);
    set(gca,'fontsize',8)
    set(handles.pushbutton_zoomSlope, 'Enable', 'on');
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot radius %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function displayPlotR(handles)
track2D = get(handles.button_2D,'Value');

if ~track2D
    legendInfo = cell(length(handles.radius),1);
    for it = 1:length(handles.radius)
        plot(handles.radius{it})
        hold all
        legendInfo{it} = ['pitch ' num2str(it)];
    end
    hold off
    if length(handles.radius)<5
        legend(legendInfo)
    end
    xlabel('pixel number')
    ylabel('R (m)')
    axis tight
else
    stepSize = str2double(get(handles.edit_stepsize,'String'));
    x = 1:stepSize:size(handles.radius,2)*stepSize;
    y = 1:size(handles.radius,1);
    imagesc(x,y,handles.radius)
    colorbar;
    title('R (m)')
    xlabel('pixel')
    ylabel('pixel')
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot slope %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function displaySlope(handles)
track2D = get(handles.button_2D,'Value');

if ~track2D
    angles = get(handles.uitable_pitchAngle,'Data');
    
    
    yyaxis left
    plot(angles,handles.Rslope,'ko')
    xlabel('pitch angle')
    ylabel('slope of R')
    hold on
    yyaxis left
    plot(angles,handles.slsl,'k-')
    axis tight
    ax = gca;
    ax.YColor = 'k';

    
    yyaxis right
    plot(angles,handles.meanR,'bx')
    ylabel('average R (m)')
    yyaxis right
    hold on
    plot(angles,handles.slR,'b-')
    ax = gca;
    ax.YColor = 'b';
    
    ymax = max(handles.meanR);
    ymin = min(handles.meanR);
    ylim([ymin-0.005,ymax+0.005])
    hold off

%     axis tight
    
    legend('slope of R','fitting','average R','fitting')
else
    bar(handles.op2(2:8).*1e3)
    ylabel('zernike coeffecients')
    xlabel('zernike terms')
    str = {'Zernike terms: 1 horizontal tilt, 2 vertical tilt, 3 defocus, 4 oblique astigmatism, 5 vertical astigmatism, 6 vertical coma, 7 horizontal coma'};
    set(handles.text_result,'String',str)
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % calculation zernike coeffecient for 2D radius %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calcZernike(hObject,handles)
        handles = guidata(hObject);
        mask = ones(size(handles.radius));
        if ~isempty(handles.mask)
            mask(handles.mask(2):handles.mask(2)+handles.mask(4),...
                handles.mask(1):handles.mask(1)+handles.mask(3))=0;
        end
        mask(abs(handles.radius-median(handles.radius(:)))>0.5)=0;
        %figure,imagesc(mask)
        outMask = handles.radius;
        outMask(mask==0) = [];
        deltaR = handles.radius-median(outMask(:));
        ZernikeList = [1:8];
        ShapeParameter = size(deltaR,2)/sqrt(size(deltaR,1)^2+size(deltaR,2)^2);
        [op1,op2]=ZernikeCalc(ZernikeList, deltaR, mask, ...
            'RECTANGLE', ShapeParameter);
        handles.op2 = op2;
        guidata(hObject, handles); % Update handles structure

        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output and callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Outputs from this function are returned to the command line.
function varargout = autoAlignMirror_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%
% Callback functions %
%%%%%%%%%%%%%%%%%%%%%%
% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
delete(hObject);


% --- Executes during object creation, after setting all properties.
function plotR_CreateFcn(hObject, eventdata, handles)



% --- Executes during object deletion, before destroying properties.
function plotR_DeleteFcn(hObject, eventdata, handles)



% --- Executes on button press in pushbutton_Start.
function pushbutton_Start_Callback(hObject, eventdata, handles)
handles.StopNow = false;
handles.plotROI = false;
guidata(hObject,handles); 
set(handles.status,'String','Caculating...')
set(handles.figure1, 'pointer', 'watch')
drawnow;
handles = guidata(hObject);

defineROI(hObject,handles);
handles = guidata(hObject);
if handles.StopNow == false
    tstart=tic;
    startCalculation(hObject,handles);
    telapsed = toc(tstart);
    set(handles.figure1, 'pointer', 'arrow')
    set(handles.status,'String',sprintf('Finished! Elapse time %f s',telapsed))
else
    msgbox('Process stopped');
    set(handles.figure1, 'pointer', 'arrow')
    set(handles.status,'String','Stopped!')
end

% --- Executes on button press in pushbutton_zoomR.
function pushbutton_zoomR_Callback(hObject, eventdata, handles)
figure
displayPlotR(handles);
angles = get(handles.uitable_pitchAngle,'Data');
legendCell=strcat('angle=',strtrim(cellstr(num2str(angles(:)))));
legend(legendCell)



% --- Executes on button press in pushbutton_zoomSlope.
function pushbutton_zoomSlope_Callback(hObject, eventdata, handles)
figure
displaySlope(handles)





function edit_pix_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function edit_pix_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_MDdist_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function edit_MDdist_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbutton_browseFile_Callback(hObject, eventdata, handles)
fileDir = uigetfile_n_dir(0,'Choose the folders containing scan images');
fileDirName = sort_nat(fileDir);
set(handles.table_filePath,'Data',fileDirName');
set(handles.uitable_pitchAngle,'Data',cell(length(fileDirName),1)); 


function pushbutton_browseFile_CreateFcn(hObject, eventdata, handles)


function scanStep_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function scanStep_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbutton_debug_Callback(hObject, eventdata, handles)
keyboard;

function zoomPreview_Callback(hObject, eventdata, handles)
figure
previewImage(handles);
colorbar;

function button_preview_Callback(hObject, eventdata, handles)
axes(handles.imagePreview);
previewImage(handles);
set(handles.zoomPreview, 'Enable', 'on');

% --- Executes on button press in button_roi.
function button_roi_Callback(hObject, eventdata, handles)
horizontal = get(handles.horizontal, 'Value');
figure
previewImage(handles);
title('Draw a rectangle to choose ROI')
h = imrect;
pos = getPosition(h);   %[xstart,ystart,width,height]
title(['ROI is [', num2str(round(pos(1))),':',num2str(round(pos(1)+pos(3))),...
    '] in x, [',num2str(round(pos(2))),':',num2str(round(pos(2)+pos(4))),...
    '] in y. You can close the window now.'])
if ~horizontal
    set(handles.start_pix,'String',num2str(round(pos(1))));
    set(handles.end_pix,'String',num2str(round(pos(1)+pos(3))));
else
    set(handles.start_pix,'String',num2str(round(pos(2))));
    set(handles.end_pix,'String',num2str(round(pos(2)+pos(4))));
end




% --- Executes during object creation, after setting all properties.
function table_filePath_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function logo_axes_CreateFcn(hObject, eventdata, handles)


% --- Executes when entered data in editable cell(s) in table_filePath.
function table_filePath_CellEditCallback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function uitable_pitchAngle_CreateFcn(hObject, eventdata, handles)



% --- Executes when entered data in editable cell(s) in uitable_pitchAngle.
function uitable_pitchAngle_CellEditCallback(hObject, eventdata, handles)
data = get(handles.uitable_pitchAngle,'Data'); 
set(handles.uitable_pitchAngle,'Data',data); 


% --- Executes during object creation, after setting all properties.
function pushbutton_debug_CreateFcn(hObject, eventdata, handles)




% --- Executes on button press in vertical.
function vertical_Callback(hObject, eventdata, handles)
set(handles.horizontal, 'Value', 0);


% --- Executes on button press in horizontal.
function horizontal_Callback(hObject, eventdata, handles)
set(handles.horizontal, 'Value', 1);


% --- Executes during object creation, after setting all properties.
function imagePreview_CreateFcn(hObject, eventdata, handles)


function start_pix_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function start_pix_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function end_pix_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function end_pix_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_parallel.
function radiobutton_parallel_Callback(hObject, eventdata, handles)
parallel = get(handles.radiobutton_parallel, 'Value');

if parallel 
   
   set(handles.edit_worker,'Enable','on');
   set(handles.pushbutton_stop,'Enable','off');
   
else

   set(handles.edit_worker,'Enable','off');
   set(handles.pushbutton_stop,'Enable','on');
end


    
    


% --- Executes on button press in radiobutton_towardsFocus.
function radiobutton_towardsFocus_Callback(hObject, eventdata, handles)
set(handles.radiobutton_towardsFocus, 'Value', 1);
set(handles.radiobutton_away, 'Value', 0);

% --- Executes on button press in radiobutton_away.
function radiobutton_away_Callback(hObject, eventdata, handles)
set(handles.radiobutton_away, 'Value', 1);
set(handles.radiobutton_towardsFocus, 'Value', 0);


% --- Executes during object creation, after setting all properties.
function plotSlope_CreateFcn(hObject, eventdata, handles)




function edit_startPitch_Callback(hObject, eventdata, handles)
pitchStart = str2double(get(handles.edit_startPitch,'String'));
pitchStep = str2double(get(handles.edit_pitchStep,'String'));
filePath = get(handles.table_filePath,'Data');
pitches = [pitchStart:pitchStep:pitchStart+pitchStep*(length(filePath)-1)]';
set(handles.uitable_pitchAngle,'Data',pitches); 

function edit_startPitch_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_pitchStep_Callback(hObject, eventdata, handles)
pitchStart = str2double(get(handles.edit_startPitch,'String'));
pitchStep = str2double(get(handles.edit_pitchStep,'String'));
filePath = get(handles.table_filePath,'Data');
pitches = [pitchStart:pitchStep:pitchStart+pitchStep*(length(filePath)-1)]';
set(handles.uitable_pitchAngle,'Data',pitches); 

function edit_pitchStep_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton_stop_Callback(hObject,eventdata,handles)
handles.StopNow = true;
guidata(hObject,handles); 


% --- Executes on button press in pushbutton_plotROI.
function pushbutton_plotROI_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
handles.plotROI = true;
defineROI(hObject,handles);
handles = guidata(hObject);



function edit_edge_Callback(hObject, eventdata, handles)

function edit_edge_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_plotProfile.
function pushbutton_plotProfile_Callback(hObject, eventdata, handles)
plotProfile(handles)


function edit_profileNo_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_profileNo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_percent_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_percent_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_showXcorr.
function radiobutton_showXcorr_Callback(hObject, eventdata, handles)


function edit_format_Callback(hObject, eventdata, handles)
format = get(handles.edit_format,'String');


% --- Executes during object creation, after setting all properties.
function edit_format_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_whichScan_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_whichScan_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_stepsize_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function edit_stepsize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_mask.
function button_mask_Callback(hObject, eventdata, handles)
figure
previewImage(handles);
title('Draw a rectangle to choose mask')
h = imrect;
mask = getPosition(h);   %[xstart,ystart,width,height]
title(['mask is [', num2str(round(pos(1))),':',num2str(round(pos(1)+pos(3))),...
    '] in x, [',num2str(round(pos(2))),':',num2str(round(pos(2)+pos(4))),...
    '] in y. You can close the window now.'])
handles.mask = mask;
guidata(hObject,handles); 


function button_1D_Callback(hObject, eventdata, handles)
handles.panel_2D.Visible = 'off';

function button_2D_Callback(hObject, eventdata, handles)
handles.panel_2D.Visible = 'on';



function edit_window_Callback(hObject, eventdata, handles)

function edit_window_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_worker_Callback(hObject, eventdata, handles)
poolsize = str2num(get(handles.edit_worker,'String'));
corenum = feature('numcores');
if poolsize>corenum
    set(handles.status,'String','Error: worker cannot be more than number of cores!');
    set(handles.edit_worker,'String',num2str(corenum));
    poolsize = corenum;
end
p=gcp('nocreate');
if isempty(p)
    parpool(poolsize);
    set(handles.status,'String',sprintf('Parpool opened: connected to %d workers.',...
        poolsize));
elseif poolsize ~= p.NumWorkers
    delete(gcp('nocreate'));
    parpool(poolsize);
    set(handles.status,'String',sprintf('Parpool opened: connected to %d workers.',...
        poolsize));
end

% --- Executes during object creation, after setting all properties.
function edit_worker_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_add.
function pushbutton_add_Callback(hObject, eventdata, handles)
olddata =get(handles.table_filePath,'Data');
newdata = [olddata;['  ']];
set(handles.table_filePath,'Data',newdata);


% --- Executes on button press in pushbutton_remove.
function pushbutton_remove_Callback(hObject, eventdata, handles)
rowNo = str2num(get(handles.edit_remove,'String'));
if ~isempty(rowNo)
    
    olddata =get(handles.table_filePath,'Data');
    newdata = olddata;
    newdata(rowNo)=[];
    set(handles.table_filePath,'Data',newdata);
    
end



function edit_remove_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function edit_remove_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
