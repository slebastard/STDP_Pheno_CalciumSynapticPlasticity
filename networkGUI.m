function varargout = networkGUI(varargin)
% NETWORKGUI MATLAB code for networkGUI.fig
%      NETWORKGUI, by itself, creates a new NETWORKGUI or raises the existing
%      singleton*.
%
%      H = NETWORKGUI returns the handle to a new NETWORKGUI or the handle to
%      the existing singleton*.
%
%      NETWORKGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NETWORKGUI.M with the given input arguments.
%
%      NETWORKGUI('Property','Value',...) creates a new NETWORKGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before launcher_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to launcher_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help networkGUI

% Last Modified by GUIDE v2.5 27-Nov-2018 22:26:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @launcher_OpeningFcn, ...
                   'gui_OutputFcn',  @launcher_OutputFcn, ...
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

env = getEnv();
addpath(genpath(env.functionsRoot));

% --- Executes just before networkGUI is made visible.
function launcher_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to networkGUI (see VARARGIN)

% Choose default command line output for networkGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Loading default values in table
set(handles.synParams,'Data',cell(1,11));

defSyn = getSynapse();
defSynData = [defSyn.tau_Ca, ...
    defSyn.tau_rho, ...
    defSyn.tau_x, ...
    defSyn.C_pre, ...
    defSyn.C_post, ...
    defSyn.delay_pre, ...
    defSyn.theta_pot, ...
    defSyn.theta_dep, ...
    defSyn.S_attr, ...
    defSyn.noise_lvl, ...
    defSyn.dampFactor];

set(handles.synParams,'Data',defSynData);

% UIWAIT makes networkGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = launcher_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% %%%%%% GIFS PANEL %%%%%%%%%%%%%%%%%%%%%%%%%%

function graphCheck_Callback(hObject, eventdata, handles)

function laplCheck_Callback(hObject, eventdata, handles)


function graphNedit_Callback(hObject, eventdata, handles)
    handles.gif.N = str2double(get(hObject,'String'));

function graphNedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function laplKedit_Callback(hObject, eventdata, handles)
    handles.gif.K = str2double(get(hObject,'String'));

function laplKedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5




% %%%%%% NETWORK PANEL %%%%%%%%%%%%%%%%%%%%%%%%%%

function netNedit_Callback(hObject, eventdata, handles)
    handles.net.N = str2double(get(hObject,'String'));

function netNedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function netNEedit_Callback(hObject, eventdata, handles)
    handles.net.NE = str2double(get(hObject,'String'));


function netNEedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function netConnectedit_Callback(hObject, eventdata, handles)
    handles.net.Connectivity = str2double(get(hObject,'String'));

function netConnectedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function netDedit_Callback(hObject, eventdata, handles)
    handles.net.D = str2double(get(hObject,'String'));

function netDedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function neuVredit_Callback(hObject, eventdata, handles)
    handles.neu.Vr = str2double(get(hObject,'String'));

function neuVredit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function neuVtedit_Callback(hObject, eventdata, handles)
    handles.neu.Vt = str2double(get(hObject,'String'));

function neuVtedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function synJedit_Callback(hObject, eventdata, handles)
    handles.syn.J = str2double(get(hObject,'String'));

function synJedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function netGedit_Callback(hObject, eventdata, handles)
    handles.net.g = str2double(get(hObject,'String'));

function netGedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function InhPlastCheck_Callback(hObject, eventdata, handles)
    handles.simu.inhPlast = get(hObject,'Value');

    
    
% %%%%% SIMULATION PANEL %%%%%%%%%%%%%%
    

function simuTedit_Callback(hObject, eventdata, handles)
    handles.simu.T = str2double(get(hObject,'String'));

function simuTedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function simudtedit_Callback(hObject, eventdata, handles)
    handles.simu.dt = str2double(get(hObject,'String'));

function simudtedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% %%%%%%%%%%%%%%% PLOT PANEL %%%%%%%%%%%%%%%%%%%%

% %% Stack subgroup %%%%%%%%%%%%%%%%%

% --- Executes on button press in pltNoStack.
function pltNoStack_Callback(hObject, eventdata, handles)

function pltSmallStack_Callback(hObject, eventdata, handles)

function pltFullStack_Callback(hObject, eventdata, handles)


% %% Histogram subgroup %%%%%%%%%%%%%%%%%

% --- Executes on button press in pltSplHist.
function pltSplHist_Callback(hObject, eventdata, handles)

function pltFullHist_Callback(hObject, eventdata, handles)

function pltNoHist_Callback(hObject, eventdata, handles)


% %%%%%%%% INIT STATE CALLBACKS %%%%%%%%%%%%%%%%%%%%%%
    
function netInputedit_Callback(hObject, eventdata, handles)

function netInputedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function weightInitDrop_Callback(hObject, eventdata, handles)

function weightInitDrop_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function initCedit_Callback(hObject, eventdata, handles)

function initCedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function initStrapedit_Callback(hObject, eventdata, handles)

function initStrapedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% %%%%%%%%%%%%%%% LAUCH BUTTON %%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in LaunchButton.
function LaunchButton_Callback(hObject, eventdata, handles)
% hObject    handle to LaunchButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   
    % Packing variables from table
    synData = get(handles.synParams,'Data');
    synData = synData(1,:);
    syn = getSynapse();
    syn.tau_Ca = synData(1,1);
    syn.tau_rho = synData(1,2);
    syn.tau_x = synData(1,3);
    syn.C_pre = synData(1,4);
    syn.C_post = synData(1,5);
    syn.delay_pre = synData(1,6);
    syn.theta_pot = synData(1,7);
    syn.theta_dep = synData(1,8);
    syn.S_attr = synData(1,9);
    syn.noise_lvl = synData(1,10);
    syn.dampFactor = synData(1,11);
    
    % Packing simulation parameters
    simu.dt = 1e-3*str2double(get(handles.simudtedit,'String'));
    simu.T = str2double(get(handles.simuTedit,'String'));

    % Packing network parameters
    net.N = str2double(get(handles.netNedit,'String'));
    net.NE = str2double(get(handles.netNEedit,'String'));
    net.Connectivity= str2double(get(handles.netConnectedit,'String'));
    net.D = 1e-3*str2double(get(handles.netDedit,'String'));
    
    neu.V_r = str2double(get(handles.neuVredit,'String'));
    neu.V_t = str2double(get(handles.neuVtedit,'String'));
    
    syn.J = str2double(get(handles.synJedit,'String'));
    net.g = str2double(get(handles.netGedit,'String'));
    net.rExtRel = str2double(get(handles.netInputedit,'String'));
    simu.inhPlast = get(handles.InhPlastCheck,'Value');
    
    % Preparing INIT STATE
    init.c = str2double(get(handles.initCedit,'String'));
    
    initOut = get(handles.weightInitDrop,'Value');
    if initOut == 1
            init.mode = 'rand';
    else
            init.mode = 'deter';
    end
    
    init.strap = str2double(get(handles.initStrapedit,'String'));
    
    % Preparing STACKS & HISTOS
    rasterOut = get(handles.rasterButtonGroup, 'SelectedObject');
    switch get(rasterOut, 'String')
        case 'Full'
            plt.all.raster = 1;
        case 'Samples'
            plt.all.raster = 2;
        otherwise
            plt.all.raster = 0;
    end
    
    stackOut = get(handles.stackButtonGroup, 'SelectedObject');
    switch get(stackOut, 'String')
        case 'Small'
            plt.spl.pres = 1;
        case 'Full'
            plt.spl.pres = 2;
        case 'Samples'
            plt.spl.pres = 3;
        otherwise
            plt.spl.pres = 0;
    end
 
    histOut = get(handles.histoButtonGroup, 'SelectedObject');
    switch get(histOut, 'String')
        case 'Full'
            plt.spl.hist = 1;
        case 'Time samples'
            plt.spl.hist = 2;
        otherwise
            plt.spl.hist = 0;
    end
    
    plt.spl.ca = 0;
    plt.spl.rho = 0;
    plt.spl.w = 0;
    
    phaseOut = get(handles.phaseButtonGroup, 'SelectedObject');
    switch get(phaseOut, 'String')
        case 'None'
            plt.spl.phase = 0;
        case 'Relative'
            plt.spl.phase = 2;
    end
    
    
    % Preparing GIFS
    if get(handles.graphCheck,'Value')
        gif.graph = 1;
        gif.N = str2double(get(handles.graphNedit,'String'));
    else
        gif.graph = 0;
        gif.N = 0;
    end
    
    if get(handles.laplCheck,'Value')
        gif.lapl = 1;
        gif.K = str2double(get(handles.laplKedit,'String'));
    else
        gif.lapl = 0;
        gif.K = 0;
    end

    simuNetwork(syn, simu, net, neu, init, plt, gif)
    


% --- Executes on button press in clearAllPush.
function clearAllPush_Callback(hObject, eventdata, handles)
% hObject    handle to clearAllPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    clear all

function checkbox16_Callback(hObject, eventdata, handles)

function checkbox15_Callback(hObject, eventdata, handles)


% --- Executes on button press in closeAllPush.
function closeAllPush_Callback(hObject, eventdata, handles)
% hObject    handle to closeAllPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    figs2keep = gcf;
    all_figs = findobj(0, 'type', 'figure');
    delete(setdiff(all_figs, figs2keep));


% --- Executes when entered data in editable cell(s) in synParams.
function synParams_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to synParams (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
