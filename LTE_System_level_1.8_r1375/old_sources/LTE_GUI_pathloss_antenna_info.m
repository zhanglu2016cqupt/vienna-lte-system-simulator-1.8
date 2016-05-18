function varargout = LTE_GUI_pathloss_antenna_info(varargin)
% LTE_GUI_PATHLOSS_ANTENNA_INFO M-file for LTE_GUI_pathloss_antenna_info.fig
%      LTE_GUI_PATHLOSS_ANTENNA_INFO, by itself, creates a new LTE_GUI_PATHLOSS_ANTENNA_INFO or raises the existing
%      singleton*.
%
%      H = LTE_GUI_PATHLOSS_ANTENNA_INFO returns the handle to a new LTE_GUI_PATHLOSS_ANTENNA_INFO or the handle to
%      the existing singleton*.
%
%      LTE_GUI_PATHLOSS_ANTENNA_INFO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LTE_GUI_PATHLOSS_ANTENNA_INFO.M with the given input arguments.
%
%      LTE_GUI_PATHLOSS_ANTENNA_INFO('Property','Value',...) creates a new LTE_GUI_PATHLOSS_ANTENNA_INFO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LTE_GUI_pathloss_antenna_info_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LTE_GUI_pathloss_antenna_info_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LTE_GUI_pathloss_antenna_info

% Last Modified by GUIDE v2.5 03-Feb-2010 17:36:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LTE_GUI_pathloss_antenna_info_OpeningFcn, ...
                   'gui_OutputFcn',  @LTE_GUI_pathloss_antenna_info_OutputFcn, ...
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


% --- Executes just before LTE_GUI_pathloss_antenna_info is made visible.
function LTE_GUI_pathloss_antenna_info_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LTE_GUI_pathloss_antenna_info (see VARARGIN)

% Choose default command line output for LTE_GUI_pathloss_antenna_info
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LTE_GUI_pathloss_antenna_info wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% % Draw TU Logo
% logo = importdata('./documentation/TULogo.jpg');
% axes(handles.axes3);
% %place image onto the axes
% image(logo);
% %remove the axis tick marks
% axis(handles.axes3,'off')
% axis(handles.axes3,'equal')

% Draw first run
frequency = str2double(get(handles.edit1,'String'));
mcl       = str2double(get(handles.edit3,'String'));
y_ax_max  = str2double(get(handles.edit4,'String'));
model_choice = get(handles.popupmenu1,'Value');
plot_macroscopic_pathloss_model(handles,model_choice,frequency,mcl,y_ax_max);


% --- Outputs from this function are returned to the command line.
function varargout = LTE_GUI_pathloss_antenna_info_OutputFcn(hObject, eventdata, handles) 
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

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
frequency = str2double(get(handles.edit1,'String'));
mcl       = str2double(get(handles.edit3,'String'));
y_ax_max  = str2double(get(handles.edit4,'String'));
model_choice = get(handles.popupmenu1,'Value');
plot_macroscopic_pathloss_model(handles,model_choice,frequency,mcl,y_ax_max);


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


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
mean_antenna_gain = str2double(get(handles.edit2,'String'));
antenna_choice = get(handles.popupmenu2,'Value');
plot_antenna_gain_pattern(handles,antenna_choice,mean_antenna_gain);


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
frequency = str2double(get(handles.edit1,'String'));
mcl       = str2double(get(handles.edit3,'String'));
y_ax_max  = str2double(get(handles.edit4,'String'));
model_choice = get(handles.popupmenu1,'Value');
plot_macroscopic_pathloss_model(handles,model_choice,frequency,mcl,y_ax_max);


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
mean_antenna_gain = str2double(get(handles.edit2,'String'));
antenna_choice = get(handles.popupmenu2,'Value');
plot_antenna_gain_pattern(handles,antenna_choice,mean_antenna_gain);


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

function plot_macroscopic_pathloss_model(handles,popup_sel_index,frequency,mcl,y_ax_max)
switch popup_sel_index
    case 1
        % COST231 urban micro
        macroscopic_pathloss_model = 'cost231';
        environment = 'urban_micro';
    case 2
        % COST231 urban macro
        macroscopic_pathloss_model = 'cost231';
        environment = 'urban_macro';
    case 3
        % COST231 suburban macro
        macroscopic_pathloss_model = 'cost231';
        environment = 'suburban_macro';
    case 4
        % free space
        macroscopic_pathloss_model = 'free space';
    case 5
        % TS 36.942 recommended (urban)
        macroscopic_pathloss_model = 'TS36942';
        environment = 'urban';
    case 6
        % TS 36.942 recommended (rural)
        macroscopic_pathloss_model = 'TS36942';
        environment = 'rural';
end

switch macroscopic_pathloss_model
    case 'free space'
        macroscopic_pathloss_model = macroscopic_pathloss_models.freeSpacePathlossModel(frequency);
    case 'cost231'
        macroscopic_pathloss_model = macroscopic_pathloss_models.cost231PathlossModel(frequency,environment);
    case 'TS36942'
        macroscopic_pathloss_model = macroscopic_pathloss_models.TS36942PathlossModel(frequency,environment);
    otherwise
        error(['"' LTE_config.macroscopic_loss_model '" macroscopic pathloss model not supported']);
end

% Will set the maximum distance to 1000 m
range = 1000;
distances = 0:10:range;
pathlosses = macroscopic_pathloss_model.pathloss(distances);
pathlosses = max(pathlosses,mcl);
plot(handles.axes1,distances,pathlosses);
title(handles.axes1,['Macroscopic pathloss, using ' macroscopic_pathloss_model.name ' model']);
xlabel(handles.axes1,'Distance [m]');
ylabel(handles.axes1,'Pathloss [dB]');
box(handles.axes1,'on');
grid(handles.axes1,'on');
y_limits = ylim;
ylim(handles.axes1,[mcl-0.1*abs(mcl) max(y_ax_max,mcl)]);

function plot_antenna_gain_pattern(handles,popup_sel_index,mean_antenna_gain)
switch popup_sel_index
    case 1
        % Omnidirectional antenna
        an_antenna = antennas.omnidirectionalAntenna;
    case 2
        % Berger antenna
        an_antenna = antennas.bergerAntenna(mean_antenna_gain);
    case 3
        % TS 36.942 antenna
        an_antenna = antennas.TS36942Antenna(mean_antenna_gain);
end
angle = -180:0.1:180;

gain = zeros(1,length(angle));
for i_=1:length(angle)
    gain(i_) = an_antenna.gain(angle(i_));
end
plot(handles.axes1,angle,gain);
box(handles.axes1,'on');
grid(handles.axes1,'on');
ylim(handles.axes1,ylim(handles.axes1)*1.1);
title(handles.axes1,{['Antenna gain, ' an_antenna.antenna_type ' antenna']});
xlabel(handles.axes1,{'\theta [°]'});
ylabel(handles.axes1,{'gain [dB]'});


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
frequency = str2double(get(handles.edit1,'String'));
mcl       = str2double(get(handles.edit3,'String'));
y_ax_max  = str2double(get(handles.edit4,'String'));
model_choice = get(handles.popupmenu1,'Value');
plot_macroscopic_pathloss_model(handles,model_choice,frequency,mcl,y_ax_max);


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
