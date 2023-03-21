function varargout = parameter_visualizer(varargin)
% PARAMETER_VISUALIZER MATLAB code for parameter_visualizer.fig
%      PARAMETER_VISUALIZER, by itself, creates a new PARAMETER_VISUALIZER or raises the existing
%      singleton*.
%
%      H = PARAMETER_VISUALIZER returns the handle to a new PARAMETER_VISUALIZER or the handle to
%      the existing singleton*.
%
%      PARAMETER_VISUALIZER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PARAMETER_VISUALIZER.M with the given input arguments.
%
%      PARAMETER_VISUALIZER('Property','Value',...) creates a new PARAMETER_VISUALIZER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before parameter_visualizer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to parameter_visualizer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help parameter_visualizer

% Last Modified by GUIDE v2.5 27-Sep-2022 10:28:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @parameter_visualizer_OpeningFcn, ...
                   'gui_OutputFcn',  @parameter_visualizer_OutputFcn, ...
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


% --- Executes just before parameter_visualizer is made visible.
function parameter_visualizer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to parameter_visualizer (see VARARGIN)

% Choose default command line output for parameter_visualizer
handles.output = hObject;
handles.popupmenu1.String = {'pt1','pt2','pt3','pt4','pt5','pt6','pt7','pt8','pt9','pt10','pt11','pt12','pt13','pt14'};
load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\rescaled_Sep2022\rescaled_data.mat')
handles.ff = ff;
handles.pRescaled = pRescaled;
handles.tRescaled = tRescaled;
[handles.full_model,handles.synFun] = fittingmodel;

for i = 1:14
    load(['E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\rescaled_manual\fitted\pt' int2str(i) '_rescaled_28-Sep-2022.mat']);
    handles.pars(:,:,i) = pars;
end
for i = 1:13
    str = ['edit' int2str(i)];
    if(i<3)
        handles.(str).String = num2str(1e3*handles.pars(i,1,1));
    else
        handles.(str).String = num2str(handles.pars(i,1,1));
    end
end
plot(handles.ff,pRescaled(:,1,1),'k','LineWidth',1); hold on;
plot(handles.ff,handles.synFun(ff,handles.pars(:,1,1)),'color','b');
plot(handles.ff,handles.full_model(ff,handles.pars(:,1,1)),'color','r');
set(gca,'xscale','log');
xlim([0.5,150]);
ylim([-4,4]);


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes parameter_visualizer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = parameter_visualizer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


idx = get(hObject,'Value');
tIdx = min(floor(handles.slider1.Value*200+1)+1,200);
for i = 1:13
    str = ['edit' int2str(i)];
    if(i<3)
        handles.(str).String = num2str(1e3*handles.pars(i,tIdx,idx));
    else
        handles.(str).String = num2str(handles.pars(i,tIdx,idx));
    end
end
hold off
plot(handles.ff,handles.pRescaled(:,tIdx,idx),'k','LineWidth',1); hold on;
plot(handles.ff,handles.synFun(handles.ff,handles.pars(:,tIdx,idx)),'color','b');
plot(handles.ff,handles.full_model(handles.ff,handles.pars(:,tIdx,idx)),'color','r');
set(gca,'xscale','log');
xlim([0.5,150]);
ylim([-4,4]);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
idx = handles.popupmenu1.Value;
tIdx = min(floor(handles.slider1.Value*200+1)+1,200);
for i = 1:13
    str = ['edit' int2str(i)];
    if(i<3)
        handles.(str).String = num2str(1e3*handles.pars(i,tIdx,idx));
    else
        handles.(str).String = num2str(handles.pars(i,tIdx,idx));
    end
end
hold off;
plot(handles.ff,handles.pRescaled(:,tIdx,idx),'k','LineWidth',1);
hold on;
plot(handles.ff,handles.synFun(handles.ff,handles.pars(:,tIdx,idx)),'color','b');
plot(handles.ff,handles.full_model(handles.ff,handles.pars(:,tIdx,idx)),'color','r');
set(gca,'xscale','log');
xlim([0.5,150]);
ylim([-4,4]);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
idx = 1;
newValue = str2num(get(hObject,'String'))*1e-3;
handles = updateParams(handles,idx,newValue,hObject);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
idx = 2;
newValue = str2num(get(hObject,'String'))*1e-3;
handles = updateParams(handles,idx,newValue,hObject);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
idx = 3;
newValue = str2num(get(hObject,'String'));
handles = updateParams(handles,idx,newValue,hObject);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
idx = 4;
newValue = str2num(get(hObject,'String'));
handles = updateParams(handles,idx,newValue,hObject);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
idx = 5;
newValue = str2num(get(hObject,'String'));
handles = updateParams(handles,idx,newValue,hObject);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
idx = 6;
newValue = str2num(get(hObject,'String'));
handles = updateParams(handles,idx,newValue,hObject);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
idx = 7;
newValue = str2num(get(hObject,'String'));
handles = updateParams(handles,idx,newValue,hObject);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
idx = 8;
newValue = str2num(get(hObject,'String'));
handles = updateParams(handles,idx,newValue,hObject);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double
idx = 9;
newValue = str2num(get(hObject,'String'));
handles = updateParams(handles,idx,newValue,hObject);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double
idx = 10;
newValue = str2num(get(hObject,'String'));
handles = updateParams(handles,idx,newValue,hObject);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double
idx = 11;
newValue = str2num(get(hObject,'String'));
handles = updateParams(handles,idx,newValue,hObject);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double
idx = 12;
newValue = str2num(get(hObject,'String'));
handles = updateParams(handles,idx,newValue,hObject);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double
idx = 13;
newValue = str2num(get(hObject,'String'));
handles = updateParams(handles,idx,newValue,hObject);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function handles = updateParams(handles,idx,newValue,hObject)
    ptIdx = handles.popupmenu1.Value;
    tIdx = min(floor(handles.slider1.Value*200+1)+1,200);
    handles.pars(idx,tIdx,ptIdx) = newValue;
    updateplot(handles);

function updateplot(handles)
    ptIdx = handles.popupmenu1.Value;
    tIdx = min(floor(handles.slider1.Value*200+1)+1,200);
    hold off;
    plot(handles.ff,handles.pRescaled(:,tIdx,ptIdx),'k','LineWidth',1);
    hold on;
    plot(handles.ff,handles.synFun(handles.ff,handles.pars(:,tIdx,ptIdx)),'color','b');
    plot(handles.ff,handles.full_model(handles.ff,handles.pars(:,tIdx,ptIdx)),'color','r');
    set(gca,'xscale','log');
    xlim([0.5,150]);
    ylim([-4,4]);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    ptIdx = handles.popupmenu1.Value;
    tIdx = min(floor(handles.slider1.Value*200+1)+1,199);
    handles.pars(:,tIdx+1,ptIdx) = handles.pars(:,tIdx,ptIdx);
    handles.slider1.Value = handles.slider1.Value+1/200;
    updateplot(handles);
    guidata(hObject,handles);
    pars = handles.pars(:,:,ptIdx)';
    csvwrite(['E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\rescaled_Sep2022\pars0_pt' int2str(ptIdx) '.csv'],pars)
