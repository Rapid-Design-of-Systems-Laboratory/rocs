function varargout = aigui(varargin)
% AIGUI M-file for aigui.fig
%      AIGUI, by itself, creates a new AIGUI or raises the existing
%      singleton*.
%
%      H = AIGUI returns the handle to a new AIGUI or the handle to
%      the existing singleton*.
%
%      AIGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AIGUI.M with the given input arguments.
%
%      AIGUI('Property','Value',...) creates a new AIGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before aigui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to aigui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help aigui

% Last Modified by GUIDE v2.5 21-Mar-2006 17:54:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @aigui_OpeningFcn, ...
                   'gui_OutputFcn',  @aigui_OutputFcn, ...
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


% --- Executes just before aigui is made visible.
function aigui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to aigui (see VARARGIN)

% Choose default command line output for aigui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes aigui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% My code

% Initialize
clear all; clc;

% Obtain data
load('results.mat');

% Plot Pareto front
if in.num_obj == 2
  % 2-objective optimization
  plot(od.fit_arch(1,:,od.iter),od.fit_arch(2,:,od.iter),'b.');
  grid on;
  xlabel('Objective 1');
  ylabel('Objective 2');
elseif in.num_obj == 3
  % 3-objective optimization
  plot3(od.fit_arch(1,:,od.iter),od.fit_arch(2,:,od.iter), ...
    od.fit_arch(3,:,od.iter),'b.');
  grid on;
  xlabel('Objective 1');
  ylabel('Objective 2');
  zlabel('Objective 3');
else
  error('Too many objectives');
end

% Turn on figure toolbar
set(gcf,'toolbar','figure');

% --- Outputs from this function are returned to the command line.
function varargout = aigui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function pareto_text_Callback(hObject, eventdata, handles)
% hObject    handle to pareto_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pareto_text as text
%        str2double(get(hObject,'String')) returns contents of pareto_text as a double


% --- Executes during object creation, after setting all properties.
function pareto_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pareto_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in trade_space_button.
function trade_space_button_Callback(hObject, eventdata, handles)
% hObject    handle to trade_space_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of trade_space_button


function trade_space_text_Callback(hObject, eventdata, handles)
% hObject    handle to trade_space_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trade_space_text as text
%        str2double(get(hObject,'String')) returns contents of trade_space_text as a double


% --- Executes during object creation, after setting all properties.
function trade_space_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trade_space_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in trade_space_locator_button.
function trade_space_locator_button_Callback(hObject, eventdata, handles)
% hObject    handle to trade_space_locator_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Obtain data
load('results.mat');

h = datacursormode;
info = getCursorInfo(h);

% Get location on Pareto Front
x = info.Position(1);
y = info.Position(2);
if in.num_obj == 3
  z = info.Position(3);
end

% Convert position to text for output
trade_space_string = [num2str(x),',',num2str(y)];
if in.num_obj == 3
  trade_space_string = [trade_space_string,',',num2str(z)];
end

% Output location to text box
set(handles.pareto_text,'string',trade_space_string);

fit_vector = [x;y];
if in.num_obj == 3
  fit_vector = [fit_vector;z];
end

% Locate corresponding trade space location
[TF,LOC] = ismember(fit_vector',od.fit_arch(:,:,od.iter)','rows');
trade_space_pos = od.pos_arch(:,LOC,od.iter);

% Output to text box
text_output = [];
for counter = 1 : 1 : length(trade_space_pos)
  text_output = [text_output,num2str(trade_space_pos(counter)),','];
end
text_output = text_output(1:end-1);
set(handles.trade_space_text,'string',text_output);
