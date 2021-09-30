function varargout = calibration_GUI(varargin)
% CALIBRATION_GUI MATLAB code for calibration_GUI.fig
%      CALIBRATION_GUI, by itself, creates a new CALIBRATION_GUI or raises the existing
%      singleton*.
%
%      H = CALIBRATION_GUI returns the handle to a new CALIBRATION_GUI or the handle to
%      the existing singleton*.
%
%      CALIBRATION_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALIBRATION_GUI.M with the given input arguments.
%
%      CALIBRATION_GUI('Property','Value',...) creates a new CALIBRATION_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before calibration_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to calibration_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help calibration_GUI

% Last Modified by GUIDE v2.5 05-Sep-2021 23:33:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @calibration_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @calibration_GUI_OutputFcn, ...
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


% --- Executes just before calibration_GUI is made visible.
function calibration_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
set(handles.zScan,'Value',0);
pos = get(0, 'MonitorPositions');
handles.pos = pos;
handles.slmPatterns.String = [];
guidata(hObject, handles);


function varargout = calibration_GUI_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;



function save_Callback(~, ~, handles)
x1 = str2double(get(handles.x1,'String'));
x2 = str2double(get(handles.x2,'String'));
x3 = str2double(get(handles.x3,'String'));
x4 = str2double(get(handles.x4,'String'));
x5 = str2double(get(handles.x5,'String'));
x6 = str2double(get(handles.x6,'String'));
y1 = str2double(get(handles.y1,'String'));
y2 = str2double(get(handles.y2,'String'));
y3 = str2double(get(handles.y3,'String'));
y4 = str2double(get(handles.y4,'String'));
y5 = str2double(get(handles.y5,'String'));
y6 = str2double(get(handles.y6,'String'));

shiftx = handles.shiftx;
shifty = handles.shifty;


A1 = [shiftx(1) shifty(1) 1; shiftx(2) shifty(2) 1; shiftx(3) shifty(3) 1];
B1_x = [x1; x2; x3];
B1_y = [y1; y2; y3];

A2 = [shiftx(4) shifty(4) 1; shiftx(5) shifty(5) 1; shiftx(6) shifty(6) 1];
B2_x = [x4; x5; x6];
B2_y = [y4; y5; y6];

val1 = linsolve(A1,B1_x);
val2 = linsolve(A1,B1_y);

val3 = linsolve(A2,B2_x);
val4 = linsolve(A2,B2_y);

cal1 = (val1+val3)./2;
cal2 = (val2+val4)./2;

V = [cal1(1) cal1(2);cal2(1) cal2(2);cal1(3) cal2(3)];

name = strcat('zoom',num2str(handles.casenum),'calidata','.csv');
csvwrite(name,V);




function x1_Callback(hObject, eventdata, handles)
function x1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y1_Callback(hObject, eventdata, handles)
function y1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x2_Callback(hObject, eventdata, handles)
function x2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y2_Callback(hObject, eventdata, handles)
function y2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x3_Callback(hObject, eventdata, handles)
function x3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y3_Callback(hObject, eventdata, handles)
function y3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x4_Callback(hObject, eventdata, handles)
function x4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y4_Callback(hObject, eventdata, handles)
function y4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x5_Callback(hObject, eventdata, handles)
function x5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y5_Callback(hObject, eventdata, handles)
function y5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x6_Callback(hObject, eventdata, handles)
function x6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y6_Callback(hObject, eventdata, handles)
function y6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function spot_Callback(hObject, eventdata, handles)

% clear all previous spot
if  isfield(handles, 'spotTxt')
    delete(handles.spotTxt);
end

% clear field
set(handles.x1, 'String', ''); 
set(handles.x2, 'String', ''); 
set(handles.x3, 'String', ''); 
set(handles.x4, 'String', ''); 
set(handles.x5, 'String', ''); 
set(handles.x6, 'String', ''); 

set(handles.y1, 'String', ''); 
set(handles.y2, 'String', ''); 
set(handles.y3, 'String', ''); 
set(handles.y4, 'String', ''); 
set(handles.y5, 'String', ''); 
set(handles.y6, 'String', ''); 


% read properties

% read zoom setting
zoom1 = get(handles.zoom1,'Value');
zoom2 = get(handles.zoom2,'Value');
zoom3 = get(handles.zoom3,'Value');


zoom = [zoom1 zoom2 zoom3];
casenum = find(zoom == 1);
handles.casenum = casenum;


name = strcat('zoom',num2str(casenum),'calidata','.csv');
cali = importdata(name);
a = cali(1,1);
b = cali(1,2);
c = cali(2,1);
d = cali (2,2);
tx = cali(3,1);
ty = cali (3,2);


V = [a b ; c d];


[image, imageY, imageX, vmax] = func_readImage;
axes(handles.image), imshow(image,[-vmax/2 vmax])

% tx = imageX/2;
% ty = imageY/2;

% mouse click to the target position

[clickx,clicky] = ginput(6);
click = [clickx clicky];
     
x = round(click(:,1));
y = round(click(:,2));
     
if  isfield(handles, 'spotTxt')
    delete(handles.spotTxt);
end


for k = 1:6
         
           hold(handles.image,'on');

           Num = sprintf('%d',k);

              spotTxt(k) = text(x(k),y(k),Num, ...
                'HorizontalAlignment','center', ...
                'FontSize',10);
           hold(handles.image,'off');
end
 
handles.spotTxt = spotTxt;

pos = handles.pos;
posNum = str2double(get(handles.posN,'String'));
sizex = pos(posNum,3);
sizey = pos(posNum,4);

dum = zeros(sizey,sizex);

[Fx,Fy] = meshgrid(-pos(posNum,3)/2:1:pos(posNum,3)/2-1,-pos(posNum,4)/2:1:pos(posNum,4)/2-1);

Gmax = str2double(get(handles.Gmax,'String'));
load('sphericalComp.mat');

i = sqrt(-1);

for k = 1:6
    
    click = [x(k);y(k)]-[tx;ty];
    shift = V\click;
    
    shiftx(k) = shift(1,1);
    shifty(k) = shift(2,1);
    
    
    data = exp(-i*2.0*pi*(shiftx(k)*Fx/sizex+shifty(k)*Fy/sizey));
    dum = dum+data;
    
  
end

data = dum.*compendata;
% data = dum;
clear dum;

phase = (angle(data)+pi)./(2*pi); 
imdata = func_quantSLM10bit(phase,Gmax,sizex,sizey);

func_displaySLM(pos, posNum, imdata)
handles.shiftx=shiftx;
handles.shifty =shifty;

%impixelinfo;
 
guidata(hObject, handles);

function zoom1_Callback(hObject, eventdata, handles)
function zoom2_Callback(hObject, eventdata, handles)
function zoom3_Callback(hObject, eventdata, handles)


function slmPatterns_Callback(hObject, eventdata, handles)
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
    
    hold(handles.image,'on');
    
    Num = sprintf('%d',k);
    axis = sprintf('%d,%d',x(k),y(k));
    
    spotTxt(k) = text(x(k),y(k),Num, ...
        'HorizontalAlignment','center', ...
        'Color', col, ...
        'FontSize',18);
    hold(handles.image,'off');
           
end
handles.spotTxt = spotTxt;

imdata=slmdata{1,index_selected};
pos = handles.pos;
posNum = str2double(get(handles.posN,'String'));

func_displaySLM(pos, posNum, imdata)    
guidata(hObject, handles);


function slmPatterns_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in scan.
function scan_Callback(hObject, eventdata, handles)

set(handles.slmPatterns, 'String','');

% clear all previous spot
if  isfield(handles, 'spotTxt')
    delete(handles.spotTxt);
end

% count SLM patterns
count = 1;

% read zoom information
zoom1 = get(handles.zoom1,'Value');
zoom2 = get(handles.zoom2,'Value');
zoom3 = get(handles.zoom3,'Value');

zoom = [zoom1 zoom2 zoom3];
casenum = find(zoom == 1);

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

% load spherical phase data
load('sphericalComp.mat');
handles.compendata = compendata;

% load other properties
Gmax = str2double(get(handles.Gmax,'String'));

pos = handles.pos;
posNum = str2double(get(handles.posN,'String'));

sizex = pos(posNum,3);
sizey = pos(posNum,4);

[Fx,Fy] = meshgrid(-pos(posNum,3)/2:1:pos(posNum,3)/2-1,-pos(posNum,4)/2:1:pos(posNum,4)/2-1);


imageX = str2double(get(handles.imageX,'String'));
imageY = str2double(get(handles.imageY,'String'));

i = sqrt(-1);

grid = 21;
if imageX > imageY
    step = round((imageY-100)/(grid-1));
else
    step = round((imageX-100)/(grid-1));
end


for x = imageX/2-10*step:step:imageX/2+10*step
    for y = imageY/2-10*step:step:imageY/2+10*step

        click = [x;y]-[tx;ty];
        shift = V\click;
        
        shiftx = shift(1,1);
        shifty = shift(2,1);
        
        data = exp(-i*2.0*pi*(shiftx*Fx/sizex+shifty*Fy/sizey)).*compendata;
        phase = (angle(data)+pi)./(2*pi);
        
        SLMdata = func_quantSLM10bit(phase,Gmax,sizex,sizey);

        % add to slmlist
        fname = sprintf('Scan_%d',count);

        cell_slmPatterns = get(handles.slmPatterns,'string');
        length_cell_slmPatterns = length(cell_slmPatterns);

        cell_slmPatterns{length_cell_slmPatterns + 1} = fname;
        set(handles.slmPatterns, 'String', cell_slmPatterns);  
        
                
        handles.slm{1,count} = SLMdata;
        handles.memory{1,count} = 1;
        handles.xmemery{1,count} = x;
        handles.ymemery{1,count} = y;
        count = count+1; 

    end
end

guidata(hObject, handles);





function loop_Callback(hObject, eventdata, handles)
function loop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function interv_Callback(hObject, eventdata, handles)
function interv_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in play.
function play_Callback(hObject, eventdata, handles)
repeat = str2double(get(handles.loop,'String'));
interval = str2double(get(handles.interv,'String'));
interval = interval/1000;

slmdata = handles.slm;
[~,l] = size(slmdata);

pos = handles.pos;
posNum = str2double(get(handles.posN,'String'));

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


function slider_Callback(hObject, ~, handles)
zScandata = handles.zScandata;
sphe = handles.spheData;

num = round(get(hObject, 'value'));
imdata = zScandata{num};

% display to SLM
pos = handles.pos;
posNum = str2double(get(handles.posN,'String'));

func_displaySLM(pos, posNum, imdata)

handles.zdata = sphe{num};
set(handles.lev,'string',num);
guidata(hObject, handles);




function slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function zScan_Callback(hObject, eventdata, handles)
set(handles.slider, 'Visible','On');
zMin = str2double(get(handles.zMin,'String'));
zMax = str2double(get(handles.zMax,'String')); 
level = 21;
interval = (zMax-zMin)*(10^-3)/20;

set(handles.slider, 'SliderStep', [1/(level-1), 10/(level-1)], 'max', level , 'min', 1, 'value', (1+level)/2);

Gmax = str2double(get(handles.Gmax,'String'));

pos = handles.pos;
posNum = str2double(get(handles.posN,'String'));
[Fx,Fy] = meshgrid(-pos(posNum,3)/2:1:pos(posNum,3)/2-1,-pos(posNum,4)/2:1:pos(posNum,4)/2-1);

sizex = pos(posNum,3);
sizey = pos(posNum,4);

f = 200/str2double(get(handles.mag,'String'));
lambda = str2double(get(handles.wave,'String'))*10^-6;

dx = str2double(get(handles.pixelSize,'String'))*10^-3;
dy = str2double(get(handles.pixelSize,'String'))*10^-3;

i = sqrt(-1);
dum = zeros(sizey,sizex);

% create SLM patterns
Nsample = 20;

for k = 1:Nsample
    
    shiftx = randi([-round(sizex/10),round(sizex/10)],1);
    shifty = randi([-round(sizex/10),round(sizey/10)],1);
    % shiftx = sizex/30;
    % shifty = sizey/30;
   
    data =exp(-i*2.0*pi*((-1)^k*shiftx*Fx/sizex+(-1)^k*shifty*Fy/sizey));
    dum = dum+data; 
end

data = dum;
clear dum;

f1 =10000;
f2=-8000;

for h = 1:level
    sphe{h} = exp(pi*i*(Fx.^2*dx^2+Fy.^2*dy^2)*(h-11)*interval/(lambda*f^2)).*exp(pi*i*((Fx.^2*dx^2)/(lambda*f1)+(Fy.^2*dy^2)/(lambda*f2)));
    phase = (angle(data.*sphe{h})+pi)./(2*pi);
    Hdata = func_quantSLM10bit(phase,Gmax,sizex,sizey);
    zScandata{h} = Hdata;
end

handles.zScandata = zScandata;
handles.spheData = sphe;

func_displaySLM(pos, posNum, zScandata{11})

guidata(hObject, handles);


function zMin_Callback(hObject, eventdata, handles)
function zMin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zMax_Callback(hObject, eventdata, handles)
function zMax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function saveZ_Callback(hObject, eventdata, handles)
compendata = handles.zdata;
save('sphericalComp.mat','compendata');



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



function pixelSize_Callback(hObject, eventdata, handles)

function pixelSize_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function posN_Callback(hObject, eventdata, handles)

function posN_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wave_Callback(hObject, eventdata, handles)

function wave_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mag_Callback(hObject, eventdata, handles)

function mag_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function generateWM_Callback(hObject, eventdata, handles)

[scan,imageY, imageX, maxv] = func_readImage;
scan = scan-scan(1,1);
scan_copy =scan;


casenum = 1;
% load shift-rotation data
name = strcat('zoom',num2str(casenum),'calidata','.csv');
cali = importdata(name);
tx = round(cali (3,1));
ty = round(cali (3,2));

grid = 21;
if imageX > imageY
    step = round((imageY-100)/(grid-1));
else
    step = round((imageX-100)/(grid-1));
end
i = 0;

for x = imageX/2-10*step-round(step/2):step:imageX/2+10*step-round(step/2)
        i = i+1;
        j = 0;
    
    for y = imageY/2-10*step-round(step/2):step:imageY/2+10*step-round(step/2)
        j = j+1; 

        
        if  y+step>imageY && x+step<imageX 
            weightMap(j,i) = max(max(scan(y:2048,x:imageX)));
        end
        if  x+step>imageX && y+step<imageY 
            weightMap(j,i) = max(max(scan(y:step,x:imageX)));
        end
        
        if  x+step>imageX && y+step>imageY 
            weightMap(j,i) = max(max(scan(y:imageY,x:imageX)));
        end
        
        if  x+step<imageX && y+step<imageY    
            weightMap(j,i) = max(max(scan(y:y+step,x:x+step)));
        end

    end
end

figure(),imshow(weightMap,[]);

weightMap = weightMap/max(max(weightMap));
[sizem,sizen]=size(weightMap);
totalsize = sizem*sizen;

a = ones(7,7);
convmap = conv2(weightMap,a)./49;
weightMap = convmap(4:25,4:25);

[~, numzerox] = find(weightMap<0.1);
[sizep,~]=size(numzerox);

ave = sum(sum(weightMap))/(totalsize-sizep);

inverseMap = zeros(sizen,sizem);


for i = 1:sizen
    for j = 1:sizem
        
        if weightMap(i,j)==0
            weightMap(i,j)=0;
        else
            
            inverseMap(i,j) = ave./weightMap(i,j);
        end
    end
end

inverseMap(inverseMap>9)=0;

dum1 = zeros(imageY,imageX);

if imageX > imageY
   compenMap = imresize(inverseMap,imageY/sizem,'nearest');
   weight_zoom1 = abs(compenMap.^(1));
   [wx,wy] = size(weight_zoom1);
 
   dum1(1:wx,1:wy)= weight_zoom1;
 
else
   compenMap = imresize(inverseMap,imageX/sizem,'nearest');
   weight_zoom1 = abs(compenMap.^(1));
   [wx,wy] = size(weight_zoom1);
   dum1(1:wx,1:wy)= weight_zoom1;
end

weight_zoom1 = dum1;


figure(),imshow(weight_zoom1,[]);
name1 = strcat('weight','.csv');


csvwrite(name1,weight_zoom1);





function imageX_Callback(hObject, eventdata, handles)
function imageX_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function imageY_Callback(hObject, eventdata, handles)
function imageY_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function savePara_Callback(hObject, eventdata, handles)
savefile = 'paraV.mat';
mag = str2double(get(handles.mag,'String'));
wave = str2double(get(handles.wave,'String'));
posNum = str2double(get(handles.posN,'String'));
pixelSize = str2double(get(handles.pixelSize,'String'));
Gmin = str2double(get(handles.Gmin,'String'));
Gmax = str2double(get(handles.Gmax,'String'));
imageX = str2double(get(handles.imageX,'String'));
imageY = str2double(get(handles.imageY,'String'));
save(savefile,'mag','wave','posNum','pixelSize','Gmin','Gmax','imageX','imageY');
