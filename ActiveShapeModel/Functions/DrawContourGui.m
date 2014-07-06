function varargout = DrawContourGui(varargin)
% This GUI DrawContourGui can be used to draw a contour line on 
% a picture  (saddly, only very basic functionality). 
%
% How does it work :
% - Start DrawContourGui
% - Load a picture
% Drawing the line
% - Left click with the mouse on a landmark point for 
%   instance the top of a finger
% - Move the mouse along the structure and use some right clicks to
%   annotate approximately the contour
% - Left click again when you arrive on the next recognizible landmark point
% When done 
% - Use save in the menu
%
% The resulting .mat file will contain a structure "p"
%  p.n : Number of contour points clicked
%  p.x, p.y : The location of the contour poinst
%  p.I : The image
%  p.t : same length as the coordinates, with "0" a major landmark point
%        and "2" only a simple point on the contour.
%
% When read with LoadDataSetNiceContour.m, the function will keep the 
% major landmark points, and interpolate a number of evenly distributed
% minor contour poinst between them.
%   
% Function written by D.Kroon University of Twente (February 2010)

% Edit the above text to modify the response to help DrawContourGui

% Last Modified by GUIDE v2.5 10-Feb-2010 11:15:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DrawContourGui_OpeningFcn, ...
                   'gui_OutputFcn',  @DrawContourGui_OutputFcn, ...
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


% --- Executes just before DrawContourGui is made visible.
function DrawContourGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DrawContourGui (see VARARGIN)

% Choose default command line output for DrawContourGui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DrawContourGui wait for user response (see UIRESUME)
% uiwait(handles.figure1);
data.handles=handles;
data.mouse_position_last=[0 0];
data.mouse_position=[0 0];
data.axes_size=[100 100];
data.Image=0;
data.npoints=0;
data.handleblueline=nan;
data.pointx=0;
data.pointy=0;
data.pointt=0;
setMyData(data);

% --- Outputs from this function are returned to the command line.
function varargout = DrawContourGui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_load_Callback(hObject, eventdata, handles)
% hObject    handle to menu_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.jpg', 'JPG';'*.png', 'PNG';'*.*', 'All files (*.*)'});
if isequal(filename,0) || isequal(pathname,0)
    disp('no valid file selected');
    return;
end
data=getMyData();
data.Image=im2double(imread(fullfile(pathname, filename)));
[fx,fy]=gradient(data.Image);
data.Speedmap=mean(fx.^2+fy.^2,3);
setMyData(data);
display_photo()


function display_photo()
data=getMyData();
imshow(data.Image);   
hold on;
data.axes_size=get(data.handles.axes1,'PlotBoxAspectRatio');
set(get(data.handles.axes1,'Children'),'ButtonDownFcn','DrawContourGui(''axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))');
setMyData(data);

    

% --------------------------------------------------------------------
function menu_save_Callback(hObject, eventdata, handles)
% hObject    handle to menu_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData();
p.n=data.npoints;
p.x=data.pointx;
p.y=data.pointy;
p.t=data.pointt;
p.I=data.Image;
uisave('p');


function setMyData(data)
% Store data struct in figure
setappdata(gcf,'data2d',data);

function data=getMyData()
% Get data struct stored in figure
data=getappdata(gcf,'data2d');


% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cursor_position_in_axes(hObject,handles);
data=getMyData(); if(isempty(data)), return, end
if(data.npoints>0)
    pos=data.mouse_position.*[size(data.Image,2) size(data.Image,1)];
    pointx1=data.pointx(data.npoints); pointy1=data.pointy(data.npoints);
    pointx2=pos(2); pointy2=pos(1);
    if(ishandle(data.handleblueline))
        delete(data.handleblueline)
    end
    data.handleblueline=plot([pointy1 pointy2],[pointx1 pointx2],'b');
    set(data.handleblueline,'ButtonDownFcn','DrawContourGui(''axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))');
end
setMyData(data); 


function cursor_position_in_axes(hObject,handles)
data=getMyData(); if(isempty(data)), return, end;
data.mouse_position_last=data.mouse_position;
% Get position of the mouse in the large axes
% p = get(0, 'PointerLocation');
% pf = get(hObject, 'pos');
% p(1:2) = p(1:2)-pf(1:2);
% set(gcf, 'CurrentPoint', p(1:2));
p = get(handles.axes1, 'CurrentPoint');
data.mouse_position=[p(1, 1) p(1, 2)]./data.axes_size(1:2);
setMyData(data);


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData();
pos=data.mouse_position.*[size(data.Image,2) size(data.Image,1)];
data.mouse_button=get(handles.figure1,'SelectionType');
disp([data.mouse_button ' click']);
switch(data.mouse_button)
    case 'normal'
        plot(pos(1),pos(2),'r*');
        type=0;
    case 'extend'
        plot(pos(1),pos(2),'b.');
        type=1;
    case 'alt'
        plot(pos(1),pos(2),'g.');
        type=2;
    case 'open'
    otherwise
end
data.npoints=data.npoints+1;
data.pointx(data.npoints)=pos(2);
data.pointy(data.npoints)=pos(1);
data.pointt(data.npoints)=type;
setMyData(data);

