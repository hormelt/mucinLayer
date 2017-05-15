function varargout = FittingGUI(varargin)
% FITTINGGUI MATLAB code for FittingGUI.fig
%      FITTINGGUI, by itself, creates a new FITTINGGUI or raises the existing
%      singleton*.
%
%      H = FITTINGGUI returns the handle to a new FITTINGGUI or the handle to
%      the existing singleton*.
%
%      FITTINGGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FITTINGGUI.M with the given input arguments.
%
%      FITTINGGUI('Property','Value',...) creates a new FITTINGGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FittingGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FittingGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FittingGUI

% Last Modified by GUIDE v2.5 28-Sep-2015 17:08:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @FittingGUI_OpeningFcn, ...
    'gui_OutputFcn',  @FittingGUI_OutputFcn, ...
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



% --- Executes just before FittingGUI is made visible.
function FittingGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FittingGUI (see VARARGIN)

% check if data is loaded; if not, prompt user for correct file
% TO DO: check that formatting is correct

if nargin < 4 % why does this work with 4?
    
    x = inputdlg('Enter space-separated numbers:',...
        'Sample', [1 50]);
    data = str2num(x{:});
    
    %     filename = uigetfile('*.mat','Input file containing rheoarray in correct format- see makeRheologyArray.m');
    %     handles.rheoarray = load(filename,'rheoarray');
    
else
    
    handles.rheoarray = varargin{1};
    handles.datasets = varargin{2};
    
end

% initial values from GUI

handles.fitOption = 'onepoint';
handles.errorOption = 'jackknife';
handles.geometry = 'planar';
handles.frames = 100;
handles.bins = 1;
numSteps = 100;
handles.tracerOption = 'a';
handles.Rmin = 1;
handles.Rmax = 15;
handles.dedrift = 0;
set(handles.frameSlider, 'Min', 100);
set(handles.frameSlider, 'Max', handles.frameSlider.Min + numSteps);
set(handles.frameSlider, 'Value', 100);
set(handles.binSlider, 'Min', 1);
set(handles.binSlider, 'Max', handles.binSlider.Min + numSteps*2);
set(handles.binSlider, 'Value', 1);
set(handles.verifySlider, 'Value', NaN);


% Choose default command line output for FittingGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



% UIWAIT makes FittingGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FittingGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in fitMenu.
function fitMenu_Callback(hObject, eventdata, handles)
% hObject    handle to fitMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Determine the selected fit option.
val = get(hObject,'Value');
% Set current fit option to the selected one.
switch val;
    case 1 % User selects one point.
        handles.fitOption = 'onepoint';
        set(handles.binText,'String','Tracers / Bin');
        set(handles.binSlider, 'Min', 1);
        set(handles.binSlider, 'Value', 1);
        set(handles.binEdit,'String',1);
    case 2 % User selects two point.
        handles.fitOption = 'twopoint';
        set(handles.binText,'String','Number of Bins');
        set(handles.binSlider, 'Min', 2);
        set(handles.binSlider, 'Value', 2);
        set(handles.binEdit,'String',2);
end
% Save the handles structure.
guidata(hObject,handles)

% Hints: contents = cellstr(get(hObject,'String')) returns fitMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fitMenu


% --- Executes during object creation, after setting all properties.
function fitMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fitMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in errorMenu.
function errorMenu_Callback(hObject, eventdata, handles)
% hObject    handle to errorMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Determine the selected fit option.
str = get(hObject, 'String');
val = get(hObject,'Value');
% Set current error option to the selected one.
switch val;
    case 1 % User selects jackknife.
        handles.errorOption = 'jackknife';
    case 2 % User selects bootstrap.
        handles.fitOption = 'bootstrap';
    case 3
        handles.fitOption = 'FWHM';
end
% Save the handles structure.
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function errorMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to errorMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in geometryMenu.
function geometryMenu_Callback(hObject, eventdata, handles)
% hObject    handle to geometryMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Determine the selected fit option.
val = get(hObject,'Value');
% Set current error option to the selected one.
switch val;
    case 1 % User selects planar.
        handles.geometry = 'planar';
    case 2 % User selects spherical.
        handles.geometry = 'spherical';
end
% Save the handles structure.
guidata(hObject,handles)

% Hints: contents = cellstr(get(hObject,'String')) returns geometryMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from geometryMenu

% --- Executes during object creation, after setting all properties.
function geometryMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to geometryMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function minEdit_Callback(hObject, eventdata, handles)
% hObject    handle to minEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.trackButton,'Enable','off');

handles.Rmin = str2double(get(handles.minEdit,'String'));

try
    if handles.Rmin > handles.Rmax;
        set(handles.trackButton,'String','Rmin < Rmax')
    elseif le(handles.Rmin,0)
        set(handles.trackButton,'String','Rmin must be >0')
    else
        set(handles.trackButton,'String','"Lets do this."')
        set(handles.trackButton,'Enable','on')
    end
catch EM
    set(handles.trackButton,'String','Awooga Awooga')
end

guidata(hObject,handles)


% Hints: get(hObject,'String') returns contents of minEdit as text
%        str2double(get(hObject,'String')) returns contents of minEdit as a double


% --- Executes during object creation, after setting all properties.
function minEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxEdit_Callback(hObject, eventdata, handles)
% hObject    handle to maxEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.trackButton,'Enable','off');

handles.Rmax = str2double(get(handles.maxEdit,'String'));

try
    if handles.Rmax < handles.Rmin
        set(handles.trackButton,'String','Rmax < Rmin')
    elseif le(handles.Rmax,0)
        set(handles.trackButton,'String','Rmax must be >0')
    else
        set(handles.trackButton,'String','"Lets do this."')
        set(handles.trackButton,'Enable','on')
    end
catch EM
    set(handles.trackButton,'String','Awooga Awooga')
end

guidata(hObject,handles)

% Hints: get(hObject,'String') returns contents of maxEdit as text
%        str2double(get(hObject,'String')) returns contents of maxEdit as a double


% --- Executes during object creation, after setting all properties.
function maxEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function frameSlider_Callback(hObject, eventdata, handles)
% hObject    handle to frameSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.frames = get(handles.frameSlider,'Value');
set(handles.frameEdit,'String',num2str(handles.frames));

guidata(hObject,handles)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function frameSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frameSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function frameEdit_Callback(hObject, eventdata, handles)
% hObject    handle to frameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.trackButton,'Enable','off');

handles.frames = str2double(get(handles.frameEdit,'String'));
if ge(handles.frames, get(handles.frameSlider,'Min')) && le(handles.frames,get(handles.frameSlider,'Max'))
    set(handles.frameSlider,'Value',handles.frames);
elseif handles.frames < get(handles.frameSlider,'Min')
    set(handles.frameSlider,'Value', get(handles.frameSlider,'Min'));
else
    set(handles.frameSlider,'Value',get(handles.frameSlider,'Max'));
end

frames = handles.frames;

try
    if le(frames,0)
        set(handles.trackButton,'String','min frames < 0')
    elseif ~(rem(frames,1)==0)
        set(handles.trackButton,'String','min frames ~int')
    else
        set(handles.trackButton,'String','"Lets do this."')
        set(handles.trackButton,'Enable','on')
    end
catch EM
    set(handles.trackButton,'String','Awooga Awooga')
end

guidata(hObject,handles)

% Hints: get(hObject,'String') returns contents of frameEdit as text
%        str2double(get(hObject,'String')) returns contents of frameEdit as a double


% --- Executes during object creation, after setting all properties.
function frameEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function binEdit_Callback(hObject, eventdata, handles)
% hObject    handle to binEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of binEdit as text
%        str2double(get(hObject,'String')) returns contents of binEdit as a double

handles.bins = str2double(get(handles.binEdit,'String'));
% if ge(handles.bins, get(handles.binSlider,'Min')) && le(handles.frames,get(handles.binSlider,'Max'))
% set(handles.binSlider,'Value',handles.bins);
% elseif handles.bins < get(handles.binSlider,'Min')
%     set(handles.binSlider,'Value', get(handles.binSlider,'Min'));
% else
%     set(handles.binSlider,'Value',get(handles.binSlider,'Max'));
% end

set(handles.trackButton,'Enable','off');

bins = handles.bins;

try
    if le(bins,0)
        set(handles.trackButton,'String','bins must be > 0')
    elseif ~(rem(bins,1)==0)
        set(handles.trackButton,'String','bins ~int')
    else
        set(handles.trackButton,'String','"Lets do this."')
        set(handles.trackButton,'Enable','on')
    end
catch
    set(handles.trackButton,'String','Awooga Awooga')
end

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function binEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to binEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function binSlider_Callback(hObject, eventdata, handles)
% hObject    handle to binSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.bins = round(get(handles.binSlider,'Value'));
set(handles.binEdit,'String',num2str(handles.bins));

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function binSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to binSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function verifySlider_Callback(hObject, eventdata, handles)
% hObject    handle to verifySlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.verifySlider,'Enable','off');

if isreal(get(handles.verifySlider,'Value'))
    set(handles.verifySlider,'Enable','on')
end

domainind = round(get(handles.verifySlider,'Value')+1);
set(handles.verifySlider,'String',[num2str(handles.verificationIDs(2,domainind)) ' - ' num2str(handles.verificationIDs(1,domainind))]);
set(handles.verifyEdit,'String',[num2str(handles.verificationIDs(2,domainind)) ' - ' num2str(handles.verificationIDs(1,domainind))]);

data = handles.rheoarray(handles.verificationIDs(2,domainind));

if handles.dedrift
    objs = dedrift_rp(data.objs,0);
else
    objs = data.objs;
end

if strcmp(handles.geometry,'spherical')
    
    [ objs_corrected ] = arcTransform( objs, data.centerposition, data.vesicle_radius, 0);
    [ handles.objs_out, ~ ] = CalcAvgRadius( objs_corrected, [], 0, 1 ); % rescaled in last step
else
    [ handles.objs_out, ~ ] = CalcAvgRadius( objs, [], 0, data.scale );
end

trtmp = handles.objs_out(:, ismember(handles.objs_out(6,:),handles.verificationIDs(1,domainind))); % separating out track j

N = size(trtmp,2)-1; % total number of displacements in track

switch lower(handles.geometry)
    
    case 'planar'
        
        delta1 = diff(trtmp(1,:))*data.scale; % displacements in x
        
    case 'spherical'
        
        delta1 = 2*data.vesicle_radius*data.scale.*asin(sqrt(sin(diff(trtmp(1,:))./2).^2)); % polar angular displacements
end

deltaxhat = zeros(1,N);

for k = 1:N
    deltaxhat(k) = data.timestep*sum(sin(pi*k.*(1:N)/(N+1)).*delta1); % discrete sine transform of x
end

Phatexp = (2.*deltaxhat.^2)/((N+1)*data.timestep); % experimental periodogram

% (2) get theoretical periodogram

R = 1/6;

D = mean(delta1.^2)/2/data.timestep + mean(delta1(1:(end-1)).*delta1(2:end))/data.timestep; % estimate for 1st coord. from cve
sigloc = sqrt(R*mean(delta1.^2)+(2*R-1)*mean(delta1(1:(end-1)).*delta1(2:end))); % estimate for 2nd coord. from cve

Phattheory = 2*D*data.timestep^2 + 2*(sigloc^2*data.timestep-2*D*R*data.timestep^2).*(1-cos((pi.*(1:N))./(N+1)));

% (3) get histogram of the ratio of experimental to theoretical
% periodogram values.

epsilonhatexp = Phatexp./Phattheory;

shapeparam = 1/2; % true for pure diffusion
scaleparam = 2; % ditto
binwidth = .5;
binmax = 5;
binnumber = binmax/binwidth;

bincenters = binwidth/2:binwidth:(binmax+binwidth/2);
binedges = eps:binwidth:(binmax+eps); %binwidth for calculating expected form- note that the gamma function diverges at x = 0 so I've added a small offset

epsilondist = hist(epsilonhatexp, bincenters); % experimental distribution of values

linepoints = .01:.01:5.01; % finer spacing so that plotted line will look nice

epsilonhattheory = zeros(1,numel(bincenters)-1); % theoretical distribution

fun = @(x) gampdf(x,shapeparam,scaleparam)*N;

for k = 1:numel(binedges)-1
    epsilonhattheory(k) = integral(fun,binedges(k),binedges(k+1)); %theoretical distribution of values
end

chi2 = sum(((epsilondist(1:binnumber)-epsilonhattheory).^2)./epsilonhattheory)/7;

set(handles.chi2edit,'String',num2str(chi2));



pdfline = gampdf(linepoints,shapeparam,scaleparam)*N*binwidth;

axes(handles.verifyAxes);

bar(bincenters(1:binnumber),epsilondist(1:binnumber));
hold on
plot(linepoints,pdfline,'m','Linewidth',2);
axis([0 binmax+1 0 max(epsilondist+10)]);
ylabel 'count'
xlabel '\epsilon'
hold off

if get(handles.binSlider,'Value')==1
    
    axes(handles.fitAxes)
    
    hold on
    
    plot(handles.a(1,:),handles.Dt(1,:),'bo');
    plot(handles.a(1,handles.a(3,:)==handles.verificationIDs(1,domainind)),handles.Dt(1,handles.Dt(6,:)==handles.verificationIDs(1,domainind)),'ro');
    
    hold off
    
end



% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function verifySlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to verifySlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function verifyEdit_Callback(hObject, eventdata, handles)
% hObject    handle to verifyEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.domainID = str2double(get(handles.verifyEdit,'String'));
if ge(handles.domainID, get(handles.domainID,'Min')) && le(handles.domainID,get(handles.domainID,'Max'))
    set(handles.verifySlider,'Value',handles.domainID);
elseif handles.domainID < get(handles.verifySlider,'Min')
    set(handles.domainID,'Value', get(handles.verifySlider,'Min'));
else
    set(handles.domainID,'Value',get(handles.verifySlider,'Max'));
end



% Hints: get(hObject,'String') returns contents of verifyEdit as text
%        str2double(get(hObject,'String')) returns contents of verifyEdit as a double


% --- Executes during object creation, after setting all properties.
function verifyEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to verifyEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in tracerMenu.
function tracerMenu_Callback(hObject, eventdata, handles)
% hObject    handle to tracerMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Determine the selected fit option.
val = get(hObject,'Value');
% Set current error option to the selected one.
switch val;
    case 1 % User selects rotating.
        handles.tracerOption = 'rotating';
    case 2 % User selects known radius.
        handles.tracerOption = 'a';
end
% Save the handles structure.
guidata(hObject,handles)

% Hints: contents = cellstr(get(hObject,'String')) returns tracerMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tracerMenu


% --- Executes during object creation, after setting all properties.
function tracerMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tracerMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in trackButton.
function trackButton_Callback(hObject, eventdata, handles, stuff)
% hObject    handle to trackButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch lower(handles.fitOption)
    
    case 'onepoint'
        
        handles.Dt = [];
        handles.a = [];
        verificationIDstmp = [];
        datasettmp = [];
        
        for j = handles.datasets
            
            data = handles.rheoarray(j);
            
            
            if handles.dedrift
                objs = dedrift_rp(data.objs,0);
            else
                objs = data.objs;
            end
            
            if strcmp(handles.geometry,'spherical')
                
                [ objs_corrected ] = arcTransform( objs, data.centerposition, data.vesicle_radius, 0);
                [ handles.objs_out, a ] = CalcAvgRadius( objs_corrected, [], 0, 1 ); % rescaled in last step
                a(1:2,:) = a(1:2,:)*data.scale;
            else
                [ handles.objs_out, a ] = CalcAvgRadius( objs, [], 0, data.scale );
            end
            
            verificationIDs = unique(handles.objs_out(6,:));
            
            verificationIDstmp = [verificationIDstmp verificationIDs];
            datasettmp = [datasettmp j.*ones(1,numel(verificationIDs))];
            
            handles.verificationIDs = [verificationIDstmp; datasettmp];
            
            [ Dt, Dr, chi2 ] = FindD( handles.objs_out, data.timestep, data.scale, 'cve', handles.geometry, strcmp(handles.tracerOption,'rotating'), [], [], data.vesicle_radius);
            
            handles.Dt = [handles.Dt Dt];
            handles.a = [handles.a a];
            
            
            
        end
        
        mean(handles.Dt(1,:))
        
        if strcmp(handles.tracerOption,'rotating')
            [meanDr, stdDr, meanDt, stdDt] = logbinD(Dr, handles.Dt, handles.bins);
            
        else
            
            if handles.bins > 1
                
                [handles.meana, handles.stda, handles.meanDt, handles.stdDt] = numbinD(handles.a, handles.Dt, handles.bins);
                
            else
                
                handles.meana = handles.a(1,:); % stupidly named, but saves a bunch of rewritten lines
                handles.stda = handles.a(2,:);
                handles.meanDt = handles.Dt(1,:);
                handles.stdDt = handles.Dt(2,:);
                
            end
            
            [ handles.etafit, handles.chifit, handles.chi2 ] = singlePointViscosity( [handles.meanDt; handles.stdDt], [handles.meana; handles.stda], data.temperature, 1, handles.errorOption, 1, 0 );
            
            axes(handles.fitAxes);
            plot(handles.fitAxes,handles.meana,handles.meanDt,'o');
            
            hold on
            
            aval = (1e-6).*handles.meana;
            etaw = 1.053*8.9e-4; % Pa s, viscosity of .1M Sucrose
            handles.aline = logspace(log10(min(aval)),log10(max(aval)));
            kB = 1.38e-23;  % Boltzmann's constant, J/K
            gam = 0.577*ones(1,numel(handles.aline)); % Euler's constant
            e = 2*etaw*handles.aline/handles.etafit(1);
            
            % constantes from Petrov, Petrosyan and Schwille Soft Matter 2012
            p = 2.74819; q = 0.51465; v = 0.73761; w = 0.52119;
            
            handles.DPetrov = kB*data.temperature/4/pi./handles.etafit(1).*(log(2./e)- gam + 4.*e/pi - (e.^2)/2.*log(2./e)).*...
                (ones(1,length(handles.aline)) - (e.^3)/pi.*log(2./e) + (v.*e.^p)./(ones(1,length(handles.aline)) + w.*e.^q)).^(-1);
            
            plot(handles.fitAxes,handles.aline*1e6,handles.DPetrov*1e12,'r');
            
            hold off
            
            ylabel('D (microns/s)');
            
            xlabel('a (microns)');
            
        end
        
        set(handles.estimateOutput,'String',num2str(handles.etafit(1)),'Enable','On');
        set(handles.uncertaintyOutput,'String',num2str(handles.etafit(2)),'Enable','On');
        
        hold off
        
        set(handles.verifySlider,'Min',0); % array indices
        set(handles.verifySlider,'Max',numel(handles.verificationIDs(1,:))-1);
        set(handles.verifySlider,'Value',0);
        set(handles.verifySlider, 'SliderStep', [1/(numel(handles.verificationIDs(1,:))-1) , 1/(numel(handles.verificationIDs(1,:))+1) ]);
        
        domainID = 0;
        
        data = handles.rheoarray(min(handles.datasets));
        
        if strcmp(handles.geometry,'spherical')
            [ objs_corrected ] = arcTransform( data.objs, data.centerposition, data.vesicle_radius, 0);
            [ objs_out, ~ ] = CalcAvgRadius( objs_corrected, [], 0, 1 ); % rescaled in last step
        else
            [ objs_out, ~ ] = CalcAvgRadius( data.objs, [], 0, data.scale );
        end
        
        set(handles.verifyEdit,'String',[num2str(min(handles.datasets)) ' - ' num2str(handles.verificationIDs(domainID+1))]);
        
        trtmp = objs_out(:, ismember(objs_out(6,:),handles.verificationIDs(domainID+1))); % separating out track j
        
        N = size(trtmp,2)-1; % total number of displacements in track
        
        switch lower(handles.geometry)
            
            case 'planar'
                
                delta1 = diff(trtmp(1,:))*data.scale; % displacements in x
                
            case 'spherical'
             
                delta1 = diff(trtmp(1,:))*data.scale; % displacements in x
                
        end
        
        deltaxhat = zeros(1,N);
        
        for k = 1:N
            deltaxhat(k) = data.timestep*sum(sin(pi*k.*(1:N)/(N+1)).*delta1); % discrete sine transform of x
        end
        
        Phatexp = (2.*deltaxhat.^2)/((N+1)*data.timestep); % experimental periodogram
        
        % (2) get theoretical periodogram
        
        R = 1/6;
        
        D = mean(delta1.^2)/2/data.timestep + mean(delta1(1:(end-1)).*delta1(2:end))/data.timestep; % estimate for 1st coord. from cve
        sigloc = sqrt(R*mean(delta1.^2)+(2*R-1)*mean(delta1(1:(end-1)).*delta1(2:end))); % estimate for 2nd coord. from cve
        
        Phattheory = 2*D*data.timestep^2 + 2*(sigloc^2*data.timestep-2*D*R*data.timestep^2).*(1-cos((pi.*(1:N))./(N+1)));
        
        % (3) get histogram of the ratio of experimental to theoretical
        % periodogram values.
        
        epsilonhatexp = Phatexp./Phattheory;
        
        shapeparam = 1/2; % true for pure diffusion
        scaleparam = 2; % ditto
        binwidth = .5;
        binmax = 5;
        binnumber = binmax/binwidth;
        
        bincenters = binwidth/2:binwidth:(binmax+binwidth/2);
        binedges = eps:binwidth:(binmax+eps); %binwidth for calculating expected form- note that the gamma function diverges at x = 0 so I've added a small offset
        
        epsilondist = hist(epsilonhatexp, bincenters); % experimental distribution of values
        
        linepoints = .01:.01:5.01; % finer spacing so that plotted line will look nice
        
        pdfline = gampdf(linepoints,shapeparam,scaleparam)*N*binwidth;
        
        epsilonhattheory = zeros(1,numel(bincenters)-1); % theoretical distribution
        
        fun = @(x) gampdf(x,shapeparam,scaleparam)*N;
        
        for k = 1:numel(binedges)-1
            epsilonhattheory(k) = integral(fun,binedges(k),binedges(k+1)); %theoretical distribution of values
        end
        
        chi2 = sum(((epsilondist(1:binnumber)-epsilonhattheory).^2)./epsilonhattheory)/7;
        
        set(handles.chi2edit,'String',num2str(chi2));
        
        axes(handles.verifyAxes);
        
        bar(bincenters(1:binnumber),epsilondist(1:binnumber));
        hold on
        plot(linepoints,pdfline,'m','Linewidth',2);
        axis([0 binmax+1 0 max(epsilondist+10)]);
        ylabel 'count'
        xlabel '\epsilon'
        hold off
        
    case 'twopoint'
        
        corrs = [];
        
        for j = handles.datasets
            
            data = handles.rheoarray(j);
            
            if strcmp(handles.geometry,'spherical')
                objs_out = arcTransform( data.objs, data.centerposition, data.vesicle_radius);
                corrs = [corrs getTwoPointCorrs( objs_out, data.scale, handles.geometry, data.vesicle_radius)];
            else
                objs_out = data.objs;
                corrs = [corrs getTwoPointCorrs( objs_out, data.scale, handles.geometry, data.vesicle_radius)];
            end
            
            
            
        end
        
        handles.corrs_out = binTwoPointCorrs( corrs, [handles.Rmin, handles.Rmax, handles.bins], data.timestep );
        
        [ handles.etafitr, handles.chifitr, etafitt, chifitt] = fitD2pt( handles.corrs_out, data.temperature, handles.errorOption, 0);
        
        handles.Rlong = logspace(log10(1/4*min(handles.corrs_out(1,:))*1e-6),log10(1.25*max(handles.corrs_out(1,:))*1e-6),1000);
        
        %     etaw = 8.9e-4; %Pa s, viscosity of H2O
        etaw = 1.053*8.9e-4; %Pa s, viscosity of .1M Sucrose
        kB = 1.38066e-23; %m^2 kg s^-2 K^-1, Boltzmann Constant
        
        Betar = 2*etaw./handles.etafitr(1).*handles.Rlong; %scaled separation
        handles.Drrpredictionfit = (pi./Betar.*stvh1(Betar) - 2./Betar./Betar - pi/2.*(bessely(0,Betar)+bessely(2,Betar)))*kB*data.temperature/2/pi./handles.etafitr(1)*data.timestep;
        
        strr = {['\eta = ', num2str(handles.etafitr(1)), ' \pm ', num2str(handles.etafitr(2))]};
        
        axes(handles.fitAxes);
        
        plot(handles.Rlong(1,:)*1e6,handles.Drrpredictionfit*1e12,'r');
        hold on
        errorxy([handles.corrs_out(1,:)',handles.corrs_out(3,:)',handles.corrs_out(5,:)'],'ColX', 1, 'ColY', 2, 'ColYe', 3, 'EdgeColor', [0 0.3 0.7], 'FaceColor', [0.3 0.6 0.9], 'MarkSize',8,'WidthEB',3);
        
        ylabel('Drr/s (microns/s)')
        
        xlabel('R (microns)')
        
        set(handles.estimateOutput,'String',num2str(handles.etafitr(1)),'Enable','On');
        set(handles.uncertaintyOutput,'String',num2str(handles.etafitr(2)),'Enable','On');
        
        hold off
        
        params = [handles.Rmin/data.scale handles.Rmax/data.scale get(handles.binSlider,'Value')];
        
        [g] = pairCorrelate( data.objs, params, [min(data.objs(1,:)) max(data.objs(1,:)) min(data.objs(2,:)) max(data.objs(2,:))], data.scale );
        
        axes(handles.verifyAxes);
        
        %         g(1)=[];
        % edges(1)=[];
        plot(g(2,:),g(1,:));
        xlabel({'Separation'});
        ylabel({'Number Density'});
        
        
end

guidata(hObject,handles)



function estimateOutput_Callback(hObject, eventdata, handles)
% hObject    handle to estimateOutput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of estimateOutput as text
%        str2double(get(hObject,'String')) returns contents of estimateOutput as a double


% --- Executes during object creation, after setting all properties.
function estimateOutput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to estimateOutput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function uncertaintyOutput_Callback(hObject, eventdata, handles)
% hObject    handle to uncertaintyOutput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of uncertaintyOutput as text
%        str2double(get(hObject,'String')) returns contents of uncertaintyOutput as a double


% --- Executes during object creation, after setting all properties.
function uncertaintyOutput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uncertaintyOutput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function chi2edit_Callback(hObject, eventdata, handles)
% hObject    handle to chi2edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of chi2edit as text
%        str2double(get(hObject,'String')) returns contents of chi2edit as a double


% --- Executes during object creation, after setting all properties.
function chi2edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chi2edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in dedriftbutton.
function dedriftbutton_Callback(hObject, eventdata, handles)
% hObject    handle to dedriftbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (get(hObject,'Value') == get(hObject,'Max'))
    handles.dedrift = 1;
else
    handles.dedrift = 0;
end

guidata(hObject,handles)

% Hint: get(hObject,'Value') returns toggle state of dedriftbutton


% --- Executes on button press in savebutton.
function savebutton_Callback(hObject, eventdata, handles)
% hObject    handle to savebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


geometry = handles.geometry;
tracer = handles.tracerOption;
minFrames = handles.frames;
dedrift = handles.dedrift;

switch lower(handles.fitOption)
    
    case 'onepoint'

etafit1pt = handles.etafit;
Dt = [handles.meanDt; handles.stdDt];
chi21pt = handles.chi2;
numberPerBin = handles.bins;
fitline = [handles.aline*1e6,handles.DPetrov*1e12];
switch lower(handles.tracerOption)
    case 'a'
a = [handles.meana; handles.stda];
save FittingGUIOutput.mat geometry tracer minFrames dedrift etafit1pt Dt chi21pt numberPerBin fitline a
    case 'rotating'
Dr = [handles.meanDr,handles.stdDr];
save FittingGUIOutput.mat geometry tracer minFrames dedrift etafit1pt Dt chi21pt numberPerBin fitline Dr
end

    case 'twopoint'

etafit2pt = handles.etafitr;
chi22pt = handles.chifitr;
Rseparations = [handles.Rmin handles.Rmax];
numberOfBins = handles.bins;
corrs = handles.corrs_out;
fitline = [handles.Rlong(1,:)*1e6,handles.Drrpredictionfit*1e12];
save FittingGUIOutput.mat geometry tracer minFrames dedrift etafit2pt chi22pt Rseparations numberOfBins corrs fitline

end


