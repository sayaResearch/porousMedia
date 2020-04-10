function varargout = DEBRIS_HEAD_LOSS_POROUS_MEDIA_MODELING_GUI(varargin)
% DEBRIS_HEAD_LOSS_POROUS_MEDIA_MODELING_GUI MATLAB code for DEBRIS_HEAD_LOSS_POROUS_MEDIA_MODELING_GUI.fig
%      DEBRIS_HEAD_LOSS_POROUS_MEDIA_MODELING_GUI, by itself, creates a new DEBRIS_HEAD_LOSS_POROUS_MEDIA_MODELING_GUI or raises the existing
%      singleton*.
%
%      H = DEBRIS_HEAD_LOSS_POROUS_MEDIA_MODELING_GUI returns the handle to a new DEBRIS_HEAD_LOSS_POROUS_MEDIA_MODELING_GUI or the handle to
%      the existing singleton*.
%
%      DEBRIS_HEAD_LOSS_POROUS_MEDIA_MODELING_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DEBRIS_HEAD_LOSS_POROUS_MEDIA_MODELING_GUI.M with the given input arguments.
%
%      DEBRIS_HEAD_LOSS_POROUS_MEDIA_MODELING_GUI('Property','Value',...) creates a new DEBRIS_HEAD_LOSS_POROUS_MEDIA_MODELING_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DEBRIS_HEAD_LOSS_POROUS_MEDIA_MODELING_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DEBRIS_HEAD_LOSS_POROUS_MEDIA_MODELING_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DEBRIS_HEAD_LOSS_POROUS_MEDIA_MODELING_GUI

% Last Modified by GUIDE v2.5 26-Feb-2014 22:16:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DEBRIS_HEAD_LOSS_POROUS_MEDIA_MODELING_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @DEBRIS_HEAD_LOSS_POROUS_MEDIA_MODELING_GUI_OutputFcn, ...
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


% --- Executes just before DEBRIS_HEAD_LOSS_POROUS_MEDIA_MODELING_GUI is made visible.
function DEBRIS_HEAD_LOSS_POROUS_MEDIA_MODELING_GUI_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DEBRIS_HEAD_LOSS_POROUS_MEDIA_MODELING_GUI (see VARARGIN)

% Choose default command line output for DEBRIS_HEAD_LOSS_POROUS_MEDIA_MODELING_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DEBRIS_HEAD_LOSS_POROUS_MEDIA_MODELING_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DEBRIS_HEAD_LOSS_POROUS_MEDIA_MODELING_GUI_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Flow Setup %%%%%%%%%%%%%%%%%%%%%%%%%

function flow_rate_Callback(hObject, eventdata, handles)
% hObject    handle to flow_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of flow_rate as text
%        str2double(get(hObject,'String')) returns contents of flow_rate as a double
flow_rate = str2double(get(hObject, 'String'));
if isnan(flow_rate)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
% Save the new density value
handles.flow_setup.flow_rate = flow_rate;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function flow_rate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to flow_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function velo_Callback(hObject, eventdata, handles)
% hObject    handle to velo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of velo as text
%        str2double(get(hObject,'String')) returns contents of velo as a double
velo = str2double(get(hObject, 'String'));
if isnan(velo)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
% Save the new density value
handles.flow_setup.velo = velo;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function velo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function temp_water_Callback(hObject, eventdata, handles)
% hObject    handle to temp_water (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of temp_water as text
%        str2double(get(hObject,'String')) returns contents of temp_water as a double
temp = str2double(get(hObject, 'String'));
if isnan(temp)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
% Save the new density value
handles.flow_setup.temp = temp;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function temp_water_CreateFcn(hObject, eventdata, handles)
% hObject    handle to temp_water (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function viscosity_Callback(hObject, eventdata, handles)
% hObject    handle to viscosity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of viscosity as text
%        str2double(get(hObject,'String')) returns contents of viscosity as a double
viscosity = str2double(get(hObject, 'String'));
if isnan(viscosity)
    viscosity = 1.002*10^-3;
    
end
% Save the new density value
handles.flow_setup.viscosity = viscosity;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function viscosity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to viscosity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Debris Property %%%%%%%%%%%%%%%%%%%%

function w_deb_Callback(hObject, eventdata, handles)
% hObject    handle to w_deb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of w_deb as text
%        str2double(get(hObject,'String')) returns contents of w_deb as a double
w_deb = str2double(get(hObject, 'String'));
if isnan(w_deb)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
% Save the new density value
handles.debris_setup.w_deb = w_deb;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function w_deb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to w_deb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function den_deb_Callback(hObject, eventdata, handles)
% hObject    handle to den_deb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of den_deb as text
%        str2double(get(hObject,'String')) returns contents of den_deb as a double
den_deb = str2double(get(hObject, 'String'));
if isnan(den_deb)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
% Save the new density value
handles.debris_setup.den_deb = den_deb;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function den_deb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to den_deb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function conc_deb_Callback(hObject, eventdata, handles)
% hObject    handle to conc_deb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of conc_deb as text
%        str2double(get(hObject,'String')) returns contents of conc_deb as a double
conc_deb = str2double(get(hObject, 'String'));
if isnan(conc_deb)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
handles.debris_setup.conc_deb = conc_deb;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function conc_deb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to conc_deb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function dia_deb_Callback(hObject, eventdata, handles)
% hObject    handle to dia_deb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dia_deb as text
%        str2double(get(hObject,'String')) returns contents of dia_deb as a double
dia_deb = str2double(get(hObject, 'String'));
if isnan(dia_deb)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
handles.debris_setup.dia_deb = dia_deb;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function dia_deb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dia_deb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Debris Property %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in poro_cal.
function poro_cal_Callback(hObject, eventdata, handles)
% hObject    handle to poro_cal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
velo_cal = handles.flow_setup.flow_rate/handles.system_setup.size_strainer;
flow_rate = handles.flow_setup.velo*handles.system_setup.size_strainer;
conc_deb = handles.debris_setup.w_deb / handles.debris_setup.den_deb / handles.system_setup.vol_water;
set(handles.velo_cal, 'String', velo_cal);
set(handles.flow_rate_cal, 'String', flow_rate);
set(handles.conc_cal, 'String', conc_deb);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fluid conditions
U = handles.flow_setup.velo * 0.01; % Approach Velocity (m/s)
u = handles.flow_setup.viscosity; % Viscosity (Pa.s)
rw = 998;
% Strainer
A = handles.system_setup.size_strainer; % Strainer surface area (cm^2)

% Material porperties
r = handles.debris_setup.den_deb; % Density (g/cm^3)
D = handles.debris_setup.dia_deb*10^-6; % Fiber diameter (m)
Sv = (4/D);
Sv2 = (4/D)^2; % Specific surface area

S2uU = Sv2*u*U; %Specific surface^2 * Viscosity * Velocity

% Kozeny constant model coefficient
% Davies %%%%%%%%%%%%%%%%%%%%%%%%%%
% a = 4.0;
% b = 56;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ingmanson et al. %%%%%%%%%%%%%%%%
a = 3.5;
b = 57;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lee et al. %%%%%%%%%%%%%%%%%%%%%%
% a = 2.1;
% b = 146;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lee et al. (Lord) %%%%%%%%%%%%%%%
% a = 1.41 - 2;
% b = 0.903;

% NUREG-1862


% Compression model coefficient
e0 = 0.9998;
N = 0.0022;
% N = 0.002;
M = 0.288;
% M = 0.255;
% Bed thickness

dL = 0.001;

i = 0;
w_tot = 0;
w_inj = handles.debris_setup.w_deb;
while w_tot < w_inj
    
    i = i + 1;
    L = i*0.001;
    n = L/dL;
    x = L/n : L/n : L;

    e(1) = e0; % Porosity
    p(1) = 0; % Pressure

%     if dL <= 0.001
%         n = n + 1;
%     end
    
    for j = 2 : n
        jj = j-1;

%Davies - Lee         
%         p(j) = p(jj) + (1-e(jj))^1.5*(a*(1+b*(1-e(jj))^3))*S2uU*dL + 0.66*Sv*(1-e(jj))/e(jj)^3*rw*U^2*dL;
%Lord - Lee         
%         p(j) = p(jj) + 1/b*(1-e(jj))^a/e(jj)^2*(1-e(jj))^2/e(jj)^3*S2uU*dL + 0.66*Sv*(1-e(jj))/e(jj)^3*rw*U^2*dL;

%NUREG-1862
        X(jj) = e(jj)/(1-e(jj));
        K(jj) = -0.5+0.5*log(1+X(jj))+1/(2+2*X(jj)+X(jj)^2);
        
        p(j) = p(jj) + S2uU*X(jj)^3/(K(jj)*(1+X(jj))^2)*(1-e(jj))^2/e(jj)^3*dL + 1.95*((1-e(jj))/(rw*U*D/u))^0.071*rw*U^2*Sv*(1-e(jj))/e(jj)^3*dL;

        e(j) = e0 - N*p(j)^M;
        
    end
    w_tot = sum((1-e)*r*A*dL*100);
    
    Bed(i).porosity = e;
    Bed(i).pressure = p;
    Bed(i).thickenss = L*100;
    Bed(i).weight = w_tot;
    Bed(i).x = 100*(max(x) - x);
    
    i_final = i
    j
    n

end

% for i = 1 : i_final
%     dP(i) = max(Bed(i).pressure);
%     dx(i) = max(Bed(i).x);
% end

% popup_sel_index = get(handles.popupmenu1, 'Value');
% switch popup_sel_index
%     case 1
%         plot(rand(5));
%     case 2
%         plot(sin(1:0.01:25.99));
%     case 3
%         bar(1:.5:10);
%     case 4
%         plot(membrane);
%     case 5
%         surf(peaks);
% end

L = L*100;
P = max(Bed(i_final).pressure);

axes(handles.porosity_ax);
cla;
plot(Bed(i_final).x,Bed(i_final).porosity);

axes(handles.pressure_ax);
cla;
plot(Bed(i_final).x,Bed(i_final).pressure);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(handles.bed_thickness, 'String', L);
set(handles.head_loss, 'String', P);
set(handles.deb_trans, 'String', w_tot);

% --- Executes on button press in pre_cal.
function pre_cal_Callback(hObject, eventdata, handles)
% hObject    handle to pre_cal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in tick_cal.
function tick_cal_Callback(hObject, eventdata, handles)
% hObject    handle to tick_cal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% System Setup %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function num_to_Callback(hObject, eventdata, handles)
% hObject    handle to num_to (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_to as text
%        str2double(get(hObject,'String')) returns contents of num_to as a double


% --- Executes during object creation, after setting all properties.
function num_to_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_to (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function vol_water_Callback(hObject, eventdata, handles)
% hObject    handle to vol_water (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vol_water as text
%        str2double(get(hObject,'String')) returns contents of vol_water as a double
vol_water = str2double(get(hObject, 'String'));
if isnan(vol_water)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
handles.system_setup.vol_water = vol_water;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function vol_water_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vol_water (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function size_strainer_Callback(hObject, eventdata, handles)
% hObject    handle to size_strainer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of size_strainer as text
%        str2double(get(hObject,'String')) returns contents of size_strainer as a double
size_strainer = str2double(get(hObject, 'String'));
if isnan(size_strainer)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
% Save the new density value
handles.system_setup.size_strainer = size_strainer;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function size_strainer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to size_strainer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function flow_rate_cal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velo_cal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function velo_cal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velo_cal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function porosity_ax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to porosity_ax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate porosity_ax


% --- Executes during object creation, after setting all properties.
function pressure_ax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pressure_ax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate pressure_ax


% --- Executes during object creation, after setting all properties.
function bed_thickness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bed_thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function head_loss_CreateFcn(hObject, eventdata, handles)
% hObject    handle to head_loss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function deb_trans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to deb_trans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function uipanel16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
