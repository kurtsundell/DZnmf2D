function varargout = Example_Data_Set_1D(varargin)
% EXAMPLE_DATA_SET_1D M-file for Example_Data_Set_1D.fig
%      EXAMPLE_DATA_SET_1D, by itself, creates a new EXAMPLE_DATA_SET_1D or raises the existing
%      singleton*.
%
%      H = EXAMPLE_DATA_SET_1D returns the handle to a new EXAMPLE_DATA_SET_1D or the handle to
%      the existing singleton*.
%
%      EXAMPLE_DATA_SET_1D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EXAMPLE_DATA_SET_1D.M with the given input arguments.
%
%      EXAMPLE_DATA_SET_1D('Property','Value',...) creates a new EXAMPLE_DATA_SET_1D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Example_Data_Set_1D_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Example_Data_Set_1D_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Example_Data_Set_1D

% Last Modified by GUIDE v2.5 08-Nov-2020 13:56:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Example_Data_Set_1D_OpeningFcn, ...
                   'gui_OutputFcn',  @Example_Data_Set_1D_OutputFcn, ...
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


% --- Executes just before Example_Data_Set_1D is made visible.
function Example_Data_Set_1D_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Example_Data_Set_1D (see VARARGIN)

% Choose default command line output for Example_Data_Set_1D
handles.output = hObject;



% --- Outputs from this function are returned to the command line.
function varargout = Example_Data_Set_1D_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
