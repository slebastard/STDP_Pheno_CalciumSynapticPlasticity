function varargout = synapseGUI(varargin)
% SYNAPSEGUI MATLAB code for synapseGUI.fig
%      SYNAPSEGUI, by itself, creates a new SYNAPSEGUI or raises the existing
%      singleton*.
%
%      H = SYNAPSEGUI returns the handle to a new SYNAPSEGUI or the handle to
%      the existing singleton*.
%
%      SYNAPSEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SYNAPSEGUI.M with the given input arguments.
%
%      SYNAPSEGUI('Property','Value',...) creates a new SYNAPSEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before launcher_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to launcher_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help synapseGUI

% Last Modified by GUIDE v2.5 27-Nov-2018 23:13:14

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

% --- Executes just before synapseGUI is made visible.
function launcher_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to synapseGUI (see VARARGIN)

% Choose default command line output for synapseGUI
handles.output = hObject;

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

handles.env = getEnv();
addpath(genpath(handles.env.functionsRoot));
addpath(handles.env.dataRoot);
handles.data.path = strcat(handles.env.dataRoot,'Venance2016/');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes synapseGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = launcher_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% %%%%%% FREQUENCY MAP PANEL %%%%%%%%%%%%%%%%%%%%%%%%%%

function FMminDtedit_Callback(hObject, eventdata, handles)
    handles.dataFit.dt.min = str2double(get(hObject,'String'));
    

function FMminDtedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FMmaxDtedit_Callback(hObject, eventdata, handles)
    handles.dataFit.dt.max = str2double(get(hObject,'String'));


function FMmaxDtedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FMdtStepedit_Callback(hObject, eventdata, handles)
    handles.dataFit.dt.step = str2double(get(hObject,'String'));

function FMdtStepedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FMmaxFreqedit_Callback(hObject, eventdata, handles)
    handles.dataFit.freq.max = str2double(get(hObject,'String'));

function FMmaxFreqedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FMstepFreqedit_Callback(hObject, eventdata, handles)
    handles.dataFit.freq.step = str2double(get(hObject,'String'));

function FMstepFreqedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function modelistbox_Callback(hObject, eventdata, handles)

function modelistbox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in freqMapPush.
function freqMapPush_Callback(hObject, eventdata, handles)
    handles.data.freqSTDP.data = csvread(strcat(handles.data.path, 'STDP_Frequency.csv'),1,0);
    handles.data.freqSTDP.length = size(handles.data.freqSTDP.data,1);
    
    handles = getSimuVars(handles);
    
    dataFit = handles.simu;
    modeID = get(handles.modelistbox,'Value');
    dataFit.mode = handles.modelistbox.String{modeID};

    dataFit.dt.min = str2double(get(handles.FMminDtedit,'String'));
    dataFit.dt.max = str2double(get(handles.FMmaxDtedit,'String'));
    dataFit.dt.step = str2double(get(handles.FMdtStepedit,'String'));

    dataFit.freq.max = str2double(get(handles.FMmaxFreqedit,'String'));
    dataFit.freq.step = str2double(get(handles.FMstepFreqedit,'String'));

    dataFit.heat = get_freq_heatmap(dataFit, handles.params);

    figure('Name','SYN_DataFit','NumberTitle','off')
    handles.data.freqSTDP.freqs=unique(handles.data.freqSTDP.data(:,5));
    n_data_freqs=length(handles.data.freqSTDP.freqs); 

    scatter3(handles.data.freqSTDP.data(:,5), handles.data.freqSTDP.data(:,2), handles.data.freqSTDP.data(:,3)./100, 50*ones(size(handles.data.freqSTDP.data,1),1), '*r')
    hold on
    
    [freq_grid, dt_grid] = meshgrid(1:dataFit.freq.step:dataFit.freq.max, dataFit.dt.min:dataFit.dt.step:dataFit.dt.max);
    dataFit.interpol = griddata(dataFit.heat(:,1), dataFit.heat(:,2), dataFit.heat(:,3), freq_grid, dt_grid);
    ribboncoloredZ(gca,dt_grid,dataFit.interpol);
    surf(freq_grid, dt_grid, dataFit.interpol);
    colormap(bluewhitered), colorbar;
    alpha 0.3
    
    for f=1:n_data_freqs
        hold on
        ids=find(handles.data.freqSTDP.data(:,5)==handles.data.freqSTDP.freqs(f) & handles.data.freqSTDP.data(:,7)~=0);
        filtered_freq=handles.data.freqSTDP.data(ids,:);
        [a,b]=sort(filtered_freq(:,2));
        h = ribbon(filtered_freq(b,2), filtered_freq(b,3)./100, 0.15);
        set(h, 'XData', filtered_freq(b,5)-1 + get(h, 'XData'));
    end   
    
    dataFit.paramPos.x = 10.0;
    dataFit.paramPos.y = 65.0;
    dataFit.paramPos.z = 2.6;
    dataFit.paramPos.colSepY = 25.0;
    dataFit.paramPos.colSepZ = 0.15;
    dataFit.paramPos.ftSize = 8;
    stampParams(handles.params, dataFit.paramPos);
   
    


% %%%%%%%%%%%%%%% STDP PANEL %%%%%%%%%%%%%%%%%%%%


function freqSTDPedit_Callback(hObject, eventdata, handles)
    handles.simu.frequency = str2double(get(hObject,'String'));

function freqSTDPedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function nPairsSTDPedit_Callback(hObject, eventdata, handles)
    handles.simu.n_iter = str2double(get(hObject,'String'));

function nPairsSTDPedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function STDPDurationedit_Callback(hObject, eventdata, handles)
    handles.simu.T = str2double(get(hObject,'String'));


function STDPDurationedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function intStepedit_Callback(hObject, eventdata, handles)
    handles.simu.int_step = str2double(get(hObject,'String'));

function intStepedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in STDPPush.
function STDPPush_Callback(hObject, eventdata, handles)
    handles = getSimuVars(handles);

    STDP = handles.simu;
    STDP.dt.min = -75;
    STDP.dt.max = 75;
    STDP.dt.step = 1;
    STDP.mode = 'rel';

    STDP.function = get_STDP_CaProd(STDP, handles.params);

    STDP.integral = sum(STDP.function(:,2))*STDP.dt.step;
    STDP.expectation = STDP.integral/(STDP.dt.max - STDP.dt.min);
    
    axes(handles.STDPaxes);
    STDP.plot = plot(STDP.function(:,1), STDP.function(:,2), '.b');
    
    dim = [.15 .6 .3 .3];
    sidestr = strcat('STDP expectation: ', num2str(100*STDP.expectation, 3+1), '%');
    annotation('textbox',dim,'String',sidestr,'FitBoxToText','on');

    neutral_hline = refline([0 1]);
    neutral_hline.Color = 'b';
    
    % title('Plasticity as a function of pre-post spike delay')
    xlabel('Pre-post spike delay (ms)')
    ylabel('Relative change in synaptic strength')

    
    
% %% POISSON MAP PANEL %%%%%%%%%%%%%%%%%%%
    
function poissonNPairsedit_Callback(hObject, eventdata, handles)
    handles.pSTDP.nTry = str2double(get(hObject,'String'));

    
function poissonNPairsedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function poissonStepfPostedit_Callback(hObject, eventdata, handles)
    handles.pSTDP.nuPost.step = str2double(get(hObject,'String'));


function poissonStepfPostedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function poissonMaxfPostedit_Callback(hObject, eventdata, handles)
    handles.pSTDP.nuPost.max = str2double(get(hObject,'String'));


function poissonMaxfPostedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function poissonMinfPostedit_Callback(hObject, eventdata, handles)
    handles.pSTDP.nuPost.min = str2double(get(hObject,'String'));


function poissonMinfPostedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function poissonStepfPreedit_Callback(hObject, eventdata, handles)
    handles.pSTDP.nuPre.step = str2double(get(hObject,'String'));


function poissonStepfPreedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function poissonMaxfPreedit_Callback(hObject, eventdata, handles)
    handles.pSTDP.nuPre.max = str2double(get(hObject,'String'));


function poissonMaxfPreedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function poissonMinfPreedit_Callback(hObject, eventdata, handles)
    handles.pSTDP.nuPre.min = str2double(get(hObject,'String'));


function poissonMinfPreedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function poissonTedit_Callback(hObject, eventdata, handles)
    handles.pSTDP.T = str2double(get(hObject,'String'));


function poissonTedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function poissonCorsTypelistbox_Callback(hObject, eventdata, handles)
    handles.pSTDP.corr.type = str2double(get(hObject,'String'));


function poissonCorsTypelistbox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function poissonCorsC12edit_Callback(hObject, eventdata, handles)
    handles.pSTDP.corr.c12 = str2double(get(hObject,'String'));


function poissonCorsC12edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function poissonCorstcedit_Callback(hObject, eventdata, handles)
    handles.pSTDP.corr.tc = str2double(get(hObject,'String'));


function poissonCorstcedit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in poissonMapPush.
function poissonMapPush_Callback(hObject, eventdata, handles)
    env = getEnv();
    addpath(genpath(env.functionsRoot), env.dataRoot);
    data.path = strcat(env.dataRoot,'Venance2016/');

    handles = getSimuVars(handles);
    
    pSTDP = handles.simu;
    pSTDP.nuPre.min = str2double(get(handles.poissonMinfPreedit,'String'));
    pSTDP.nuPre.max = str2double(get(handles.poissonMaxfPreedit,'String'));
    pSTDP.nuPre.step = str2double(get(handles.poissonStepfPreedit,'String'));
    pSTDP.nuPost.min = str2double(get(handles.poissonMinfPostedit,'String'));
    pSTDP.nuPost.max = str2double(get(handles.poissonMaxfPostedit,'String'));
    pSTDP.nuPost.step = str2double(get(handles.poissonStepfPostedit,'String'));
    pSTDP.nTry = str2double(get(handles.poissonNPairsedit,'String'));
    
    typeID = get(handles.poissonCorsTypelistbox,'Value');
    pSTDP.corr.type = handles.poissonCorsTypelistbox.String{typeID};
    pSTDP.corr.c12 = str2double(get(handles.poissonCorsC12edit,'String'));
    pSTDP.corr.tc = str2double(get(handles.poissonCorstcedit,'String'));
    
    pSTDP.map = poissonMap(handles.params, pSTDP);

    [pre_grid_intp, post_grid_intp] = meshgrid(pSTDP.nuPre.min:0.2*pSTDP.nuPre.step:pSTDP.nuPre.max, pSTDP.nuPost.min:0.2*pSTDP.nuPost.step:pSTDP.nuPost.max);
    [pre_grid_spl, post_grid_spl] = meshgrid(pSTDP.nuPre.min:pSTDP.nuPre.step:pSTDP.nuPre.max, pSTDP.nuPost.min:pSTDP.nuPost.step:pSTDP.nuPost.max);
    pSTDP.interpol = griddata(pre_grid_spl, post_grid_spl, pSTDP.map, pre_grid_intp, post_grid_intp);
    % ribboncoloredZ(gca,dt_grid,dataFit.interpol);
    figure('Name','SYN_PoissonMap','NumberTitle','off')
    pSTDP.plot = log(pSTDP.interpol);
    surf(pre_grid_intp, post_grid_intp, pSTDP.plot);
    xlabel('Presyn rate')
    ylabel('Postsyn rate')
    zlabel('Log potentiation')
    colormap(bluewhitered), colorbar;
    alpha 0.3
    
    pSTDP.paramPos.x = pSTDP.nuPre.min + 0.85*(pSTDP.nuPre.max - pSTDP.nuPre.min);
    pSTDP.paramPos.y = pSTDP.nuPost.min + 0.85*(pSTDP.nuPost.max - pSTDP.nuPost.min);
    pSTDP.paramPos.z = max(max(pSTDP.plot)) - 0.6*(max(max(pSTDP.plot)) - min(min(pSTDP.plot)));
    pSTDP.paramPos.colSepY = 0.15*(pSTDP.nuPost.max - pSTDP.nuPost.min);
    pSTDP.paramPos.colSepZ = 0.05*(max(max(pSTDP.plot)) - min(min(pSTDP.plot)));
    pSTDP.paramPos.ftSize = 8;
    stampParams(handles.params, pSTDP.paramPos);



% %% CLEAR & CLOSE PUSH BUTTONS %%%%%%%%%%%%%%%%%%%%%%%%%%%    

% --- Executes on button press in clearAllPush.
function clearAllPush_Callback(hObject, eventdata, handles)
% hObject    handle to clearAllPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    clear all


% --- Executes on button press in closeAllPush.
function closeAllPush_Callback(hObject, eventdata, handles)
% hObject    handle to closeAllPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    figs2keep = gcf;
    all_figs = findobj(0, 'type', 'figure');
    delete(setdiff(all_figs, figs2keep));


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function handles = getSimuVars(handles)
    handles.simu.n_iter = str2double(get(handles.nPairsSTDPedit,'String'));
    handles.simu.frequency = str2double(get(handles.freqSTDPedit,'String'));
    handles.simu.int_step = str2double(get(handles.intStepedit,'String'));
    handles.simu.T = str2double(get(handles.STDPDurationedit,'String'));
    
    handles.simu.int_scheme = 'euler_expl';
    handles.simu.model = 'caProd';

    synData = get(handles.synParams,'Data');
    synData = synData(1,:);
    handles.params = getSynapse();
    handles.params.tau_Ca = synData(1,1);
    handles.params.tau_rho = synData(1,2);
    handles.params.tau_x = synData(1,3);
    handles.params.C_pre = synData(1,4);
    handles.params.C_post = synData(1,5);
    handles.params.delay_pre = synData(1,6);
    handles.params.theta_pot = synData(1,7);
    handles.params.theta_dep = synData(1,8);
    handles.params.S_attr = synData(1,9);
    handles.params.noise_lvl = synData(1,10);
    handles.params.dampFactor = synData(1,11);
    handles.params.theta_act = handles.params.theta_dep;
    handles.params.TD = 0;

    handles.params.tau_x = 1e3 .* handles.params.tau_x;
    handles.params.tau_Ca = 1e3 .* handles.params.tau_Ca;
    handles.params.tau_rho = 1e3 .* handles.params.tau_rho;
    handles.params.delay_pre = 1e3 .* handles.params.delay_pre;
    handles.params.tau_w = 50000;
