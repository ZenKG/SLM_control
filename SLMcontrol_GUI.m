function varargout = SLMcontrol_GUI(varargin)

%      SLMCONTROL_GUI MATLAB code for SLMcontrol_GUI.fig
%      SLMCONTROL_GUI, by itself, creates a new SLMCONTROL_GUI or raises the existing
%      singleton*.
%
%      H = SLMCONTROL_GUI returns the handle to a new SLMCONTROL_GUI or the handle to
%      the existing singleton*.
%
%      SLMCONTROL_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SLMCONTROL_GUI.M with the given input arguments.
%
%      SLMCONTROL_GUI('Property','Value',...) creates a new SLMCONTROL_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SLMcontrol_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SLMcontrol_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SLMcontrol_GUI

% Last Modified by GUIDE v2.5 06-Jul-2021 11:16:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SLMcontrol_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SLMcontrol_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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



% --- Executes just before SLMcontrol_GUI is made visible.
function SLMcontrol_GUI_OpeningFcn(hObject, ~, handles, varargin)

handles.output = hObject;

pos = get(0, 'MonitorPositions');
handles.pos = pos;

handles.slmPatterns.String = [];
handles.count = 1;

% Update handles structure
guidata(hObject, handles);



function varargout = SLMcontrol_GUI_OutputFcn(~, ~, handles) 
varargout{1} = handles.output;



function startStopCamera_Callback(hObject, ~, handles)
% read an image
[image, ~, ~, ~] = func_readImage;
axes(handles.camera), imshow(image,[])

%impixelinfo;
guidata(hObject, handles);
    

function initializeSLM_Callback(hObject, ~, handles)

pos = handles.pos;
load('paraV.mat')

SLMX = pos(posNum,3);
SLMY = pos(posNum,4);

%Initialize refer image data
Inpic = zeros(SLMY,SLMX);
func_displaySLM(pos, posNum,Inpic)
    
guidata(hObject, handles);




% --- Executes on button press in spotGenerator.
function spotGenerator_Callback(hObject, ~, handles)

% clear all previous spot
if  isfield(handles, 'spotTxt')
    delete(handles.spotTxt);
end

tic

% read zoom information
zoom1 = get(handles.zoom1,'Value');
zoom2 = get(handles.zoom2,'Value');
zoom3 = get(handles.zoom3,'Value');


zoom = [zoom1 zoom2 zoom3];
casenum = find(zoom == 1);
handles.casenum = casenum;

% load weight comp data
weight = importdata('weight.csv');
handles.weight = weight;


% load shift-rotation data
name = strcat('zoom',num2str(casenum),'calidata','.csv');
cali = importdata(name);
a = cali(1,1);
b = cali(1,2);
c = cali(2,1);
d = cali (2,2);
tx = round(cali (3,1));
ty = round(cali (3,2));
V = [a b ; c d];

handles.V = V;
handles.tx = tx;
handles.ty = ty;

% load spherical phase data
load('sphericalComp.mat');
handles.compendata = compendata;


% load property
pos = handles.pos;
load('paraV.mat')

sizex = pos(posNum,3);
sizey = pos(posNum,4);

[Fx,Fy] = meshgrid(-pos(posNum,3)/2:1:pos(posNum,3)/2-1,-pos(posNum,4)/2:1:pos(posNum,4)/2-1);

i = sqrt(-1);
dum = zeros(sizey,sizex);
col = rand(1,3);

% mouse click to the target position
% [clickx,clicky] = ginputc('Color','r','ShowPoints',true,'ConnectPoints',false);
% click = [clickx clicky];
% not mouse click
click = [244.8848  428.8702;239.3094  614.7142;239.3094  843.3022;514.3584  843.3022;...
    746.6633  847.0191; 741.0880  657.4583; 741.0880  466.0390; 510.6416  473.4728;512.5000  659.3167];

[m,~] = size(click);
     
x = round(click(:,1));
y = round(click(:,2));

     
% display spot numbers
for k = 1:m
    
    hold(handles.camera,'on');
    Num = sprintf('%d',k);
   
    spotTxt(k) = text(x(k),y(k),Num, ...
        'HorizontalAlignment','center', ...
        'Color', col, ...
        'FontSize',18);
    hold(handles.camera,'off');
end

     
% spotTxt is overlayed numbers on the screen.
handles.spotTxt = spotTxt;


% create SLM patterns
for k = 1:m
    
    click = [x(k);y(k)]-[tx;ty];
    shift = V\click;
    
    
    shiftx = shift(1,1);
    shifty = shift(2,1);
     
    data =sqrt(weight(y(k),x(k))).* exp(-i*2.0*pi*(shiftx*Fx/sizex+shifty*Fy/sizey));
    % Use this code when you don't have weight map data
    % data = exp(-i*2.0*pi*(shiftx*Fx/sizex+shifty*Fy/sizey));
    dum = dum+data;
    
end

data = dum.*compendata;
clear dum;


% make round shape
W = str2double(get(handles.spotSize,'String'));
W = W*(100/60);
iter = 30;

if W ~= 0  
CGHdata = func_roundCGH(W,iter,sizex,sizey);
data = CGHdata.*data;
end


% normalize phase map   
phaseMap = (angle(data)+pi)./(2*pi); 

% SLM quantitazation
SLMdata = func_quantSLM10bit(phaseMap,Gmax,sizex,sizey);

% count SLM patterns
count = handles.count;
fname = sprintf('SLMpattern_%d.bmp',count);

%add file to the list
cell_slmPatterns = get(handles.slmPatterns,'String');
length_cell_slmPatterns = length(cell_slmPatterns);

cell_slmPatterns{length_cell_slmPatterns + 1} = fname;
set(handles.slmPatterns, 'String', cell_slmPatterns);   

handles.slm{1,count} = SLMdata;
handles.memory{1,count} = m;
handles.xmemery{1,count} = x;
handles.ymemery{1,count} = y;


% display to SLM
func_displaySLM(pos,posNum,SLMdata)
        
count = count+1;
handles.count = count;

toc        
  
%impixelinfo;
guidata(hObject, handles);



function deleteSpot_Callback(hObject, eventdata, handles)
% clear all previous spot
if  isfield(handles, 'spotTxt')
    delete(handles.spotTxt);
end
guidata(hObject, handles);



function Play_Callback(hObject, eventdata, handles)

pos = handles.pos;
load('paraV.mat')

repeat = str2double(get(handles.loop,'String'));
interval = str2double(get(handles.Interv,'String'));
interval = interval/1000;

slmdata = handles.slm;
[~,l] = size(slmdata);

all_figs = findobj(0,'type','figure');
figs2keep = all_figs(1);
delete(setdiff(all_figs, figs2keep));

% full screen display
 FigH = figure('Name','SLMdata','IntegerHandle','on','MenuBar','none','ToolBar','none',...
    'OuterPosition',[pos(posNum,1) pos(posNum,2) pos(posNum,3) pos(posNum,4)],...
    'InnerPosition',[pos(posNum,1) pos(posNum,2) pos(posNum,3) pos(posNum,4)],...
    'WindowState','fullscreen');

for m = 1:repeat
    
    for n = 1:l
        figure(FigH),imshow(slmdata{n},'InitialMagnification','fit','Border','tight');
        pause(interval)
    end
end

fprintf('display OK');



    
function slmPatterns_Callback(hObject, eventdata, handles)

pos = handles.pos;
load('paraV.mat')

% clear all previous spot
if  isfield(handles, 'spotTxt')
    delete(handles.spotTxt);
end

index_selected = get(handles.slmPatterns,'Value');
slmdata = handles.slm;

% retrieving spot information
m = handles.memory{1,index_selected};
x = handles.xmemery{:,index_selected};
y = handles.ymemery{:,index_selected};

% overlay selected spots
col = rand(1,3);
for k = 1:m
    
    hold(handles.camera,'on');
    
    Num = sprintf('%d',k);
    
    spotTxt(k) = text(x(k),y(k),Num, ...
        'HorizontalAlignment','center', ...
        'Color', col, ...
        'FontSize',18);
    hold(handles.camera,'off');
           
end
handles.spotTxt = spotTxt;

imdata=slmdata{1,index_selected};
func_displaySLM(pos, posNum, imdata)    
guidata(hObject, handles);


function slmPatterns_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function spotSize_Callback(hObject, eventdata, handles)

function spotSize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Gmin_Callback(hObject, eventdata, handles)

function Gmin_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Gmax_Callback(hObject, eventdata, handles)

function Gmax_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function zoom1_Callback(hObject, eventdata, handles)
guidata(hObject, handles);


function zoom2_Callback(hObject, eventdata, handles)
guidata(hObject, handles);


function zoom3_Callback(hObject, eventdata, handles)
guidata(hObject, handles);




function loop_Callback(hObject, eventdata, handles)

function loop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Interv_Callback(hObject, eventdata, handles)

function Interv_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

