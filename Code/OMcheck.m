function varargout = OMcheck(varargin)
%OMCHECK M-file for OMcheck.fig
%      OMCHECK, by itself, creates a new OMCHECK or raises the existing
%      singleton*.
%
%      H = OMCHECK returns the handle to a new OMCHECK or the handle to
%      the existing singleton*.
%
%      OMCHECK('Property','Value',...) creates a new OMCHECK using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to OMcheck_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      OMCHECK('CALLBACK') and OMCHECK('CALLBACK',hObject,...) call the
%      local function named CALLBACK in OMCHECK.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help OMcheck

% Last Modified by GUIDE v2.5 15-Sep-2015 03:42:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OMcheck_OpeningFcn, ...
                   'gui_OutputFcn',  @OMcheck_OutputFcn, ...
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


% --- Executes just before OMcheck is made visible.
function OMcheck_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for OMcheck
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes OMcheck wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = OMcheck_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in butOpen.
function butOpen_Callback(hObject, eventdata, handles)
% hObject    handle to butOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname, ~] = uigetfile( '.mat', 'Select XX_check.mat File', '../data/');
filenamefull = [pathname, filename];
set(handles.txtDir, 'String', filenamefull)
setappdata( 0, 'filename', filenamefull)
% load data
load( filenamefull);
% show img
axes( handles.axeImg);
hold off
imshow( sum(datamat, 3)/size(datamat,3), []);
colormap('Hot')
hold on
maxOUF = max(OUF(:));
[x,y] = meshgrid( 1:size(v1,2), 1:size(v1, 1));
quiver(x,y,v1,u1,0.5*maxOUF, 'color', 'b', 'LineStyle', '-');
quiver(x,y,v2,u2,0.5*maxOUF, 'color', 'b', 'LineStyle', '-');

% --- Executes on button press in butPt.
function butPt_Callback(hObject, eventdata, handles)
% hObject    handle to butPt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% select pt
[xx yy] = getpts( handles.axeImg);
xx = xx(1); yy = yy(1);
xx = round(xx); yy = round(yy);
% load data
load( getappdata( 0, 'filename'));
% show pt
axes( handles.axeImg);
hold off
imshow( sum(datamat,3)/size(datamat,3), []);
colormap( 'Hot')
hold on
maxOUF = max(OUF(:));
[x,y] = meshgrid( 1:size(v1,2), 1:size(v1, 1));
quiver(x,y,v1,u1,0.5*maxOUF, 'color', 'b', 'LineStyle', '-');
quiver(x,y,v2,u2,0.5*maxOUF, 'color', 'b', 'LineStyle', '-');
scatter( xx, yy, 'fill', 'blue');
%%
axes( handles.axePlot);
hold off
scatter( Ang, squeeze(datamat( yy, xx, :)/(A(yy,xx)+B(yy,xx))), 'red');
hold on
ang = 0:1:180;
fit = A(yy,xx)*cos( 2*(ang-phy(yy,xx))/180*pi)+B(yy,xx);
fit = fit /(A(yy,xx)+B(yy,xx));
plot( ang, fit)
set( handles.axePlot, 'XTickLabel',{'0','30','60','90','120','150','180'},...
    'XTick',[0 30 60 90 120 150 180], 'XLim', [0, 180], 'YTickLabel', ...
    {'0','0.2','0.4','0.6','0.8','1.0'}, 'YTick',0:0.2:1,'YLim', [-0.1, 1.2]);
%%
set( handles.txtOri, 'String', ['Orientation: ', num2str( round(phy(yy,xx)))]);
set( handles.txtOUF, 'String', ['OUF: ', num2str( OUF(yy,xx), '%.3d')]);
set( handles.txtRmse, 'String', ['R2: ', num2str( adjR2(yy,xx), '%.3d')]);
