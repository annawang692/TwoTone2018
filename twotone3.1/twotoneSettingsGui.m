function  varargout = twotoneSettingsGui( varargin)
%  varargout = twotoneSettingsGui( varargin)
%
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0, released 110426
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only” license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%
% FUNCTION: twotoneSettingsGui

% Last Modified by GUIDE v2.5 05-Nov-2010 13:05:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @twotoneSettingsGui_OpeningFcn, ...
                   'gui_OutputFcn',  @twotoneSettingsGui_OutputFcn, ...
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


% --- Executes just before twotoneSettingsGui is made visible.
function twotoneSettingsGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to twotoneSettingsGui (see VARARGIN)

%MUST supply twotoneDataIn as the first argument - otherwise load the default
% CW settings. If the default CW settings dont exist, create them.
if (numel(varargin) > 0)
  twotoneDataIn = varargin{1};
elseif exist(twotoneInstallDir('','twotoneDefaultSettings_CW.mat'),'file')
  load(twotoneInstallDir('','twotoneDefaultSettings_CW.mat'),'twotoneData');
  twotoneDataIn = twotoneData;
  clear twotoneData;
else % create a new default file
  initializeTwotoneSetting('CWonly')
  load(twotoneInstallDir('','twotoneDefaultSettings_CW.mat'),'twotoneData');
  twotoneDataIn = twotoneData;
  clear twotoneData;
end

% load the general options
if (numel(varargin) > 1)
  twotoneGeneralSettingIn = varargin{2};
elseif exist(twotoneInstallDir('','twotoneGeneralSetting.mat'),'file')
  load(twotoneInstallDir('','twotoneGeneralSetting.mat'),'twotoneGeneralSetting');
  twotoneGeneralSettingIn = twotoneGeneralSetting;
  clear twotoneGeneralSetting;
else % create a new default file
  initializeTwotoneSetting('algSettingOnly')
  load(twotoneInstallDir('','twotoneGeneralSetting.mat'),'twotoneGeneralSetting');
  twotoneGeneralSettingIn = twotoneGeneralSetting;
  clear twotoneGeneralSetting;
end

if numel(varargin)>2
  imSize = varargin{3};
  imSizeString = [num2str(imSize(1)),' x ',num2str(imSize(2))];
  set(handles.TextBoxCurrentImSize,'String',imSizeString);
end


%settings initially all unchanged
twotoneDataOut = twotoneDataIn;
twotoneGeneralSettingOut = twotoneGeneralSettingIn;
% SET ALL THE PARAMETERS TO CURRENT VALUES
% %image properties:
%   twotoneDataIn.imageLim = [1, 256, 1, 256; ...
%		      257, 512, 1, 256];
set(handles.editDXStart,'String',num2str(twotoneDataIn.settings.imageSettings.channelImageLim(1,1)));
set(handles.editDXEnd,	'String',num2str(twotoneDataIn.settings.imageSettings.channelImageLim(1,2)));
set(handles.editDYStart,'String',num2str(twotoneDataIn.settings.imageSettings.channelImageLim(1,3)));
set(handles.editDYEnd,	'String',num2str(twotoneDataIn.settings.imageSettings.channelImageLim(1,4)));
set(handles.editAXStart,'String',num2str(twotoneDataIn.settings.imageSettings.channelImageLim(2,1)));
set(handles.editAXEnd,	'String',num2str(twotoneDataIn.settings.imageSettings.channelImageLim(2,2)));
set(handles.editAYStart,'String',num2str(twotoneDataIn.settings.imageSettings.channelImageLim(2,3)));
set(handles.editAYEnd,	'String',num2str(twotoneDataIn.settings.imageSettings.channelImageLim(2,4)));
%   twotoneDataIn.avgFirst 
set(handles.editFrameAvgMin, 'String',num2str(twotoneDataIn.settings.autoDetectSettings.averageFrameLim(1)));
%   twotoneDataIn.avgLast 
set(handles.editFrameAvgMax, 'String',num2str(twotoneDataIn.settings.autoDetectSettings.averageFrameLim(2)));
% %autodetection parameters:
%   twotoneDataIn.BPdiscDiametre 
set(handles.editBPDiameter, 'String',num2str( twotoneDataIn.settings.autoDetectSettings.bandPassKernelDiameter) );
%   twotoneDataIn.WindowSize 
set(handles.editWindowRadius, 'String',num2str( twotoneDataIn.settings.autoDetectSettings.fitSubImageRadius) );
set(handles.editFileAppendString,'String',twotoneDataIn.settings.fileAppendString);
%psf estimate
set(handles.editPSFestimate, 'String',num2str( twotoneGeneralSettingIn.PSFwidthEstimate));

% fit parameters
%
% use current fit parameters
currentAlgorithm = twotoneDataIn.settings.fitSettings.algorithm;
currentFitParamName =twotoneDataIn.settings.fitSettings.fitParamName;
currentParamVal = twotoneDataIn.settings.fitSettings.param;
set(handles.popupmenuAnalysisAlgorithm, 'String', twotoneGeneralSettingIn.analysisAlgorithm);
%work out which algorithm its currently set to
currentAlgNo = find(strcmp(twotoneGeneralSettingIn.analysisAlgorithm,currentAlgorithm));
set(handles.popupmenuAnalysisAlgorithm, 'Value', currentAlgNo);

set(handles.popupmenuFitAlgorithmParam, 'String',twotoneDataIn.settings.fitSettings.fitParamName);
set(handles.editFitAlgorithmParam,'String',num2str(currentParamVal(1)));

% Update handles structure
%update setappdata
setappdata(handles.figure1, 'twotoneDataIn', twotoneDataIn);
setappdata(handles.figure1, 'twotoneDataOut', twotoneDataOut);
setappdata(handles.figure1, 'twotoneGeneralSettingIn', twotoneGeneralSettingIn);
setappdata(handles.figure1, 'twotoneGeneralSettingOut', twotoneGeneralSettingOut);
%update the old guidata handles with the new one
% for some reason you have to reference it as handles.figure1
% not entirely sure why
guidata(handles.figure1,handles);
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

%UIWAIT pauses the execution of other programs until this one is done
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = twotoneSettingsGui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

varargout{1}=getappdata(handles.figure1, 'twotoneDataOut');
varargout{2}=getappdata(handles.figure1, 'twotoneGeneralSettingOut');
% Hint: delete(hObject) closes the figure
delete(hObject);

%-------------------------------------------
%-------------------------------------------
function figure1_CloseRequestFcn(hObject, eventdata, handles)

%settings all unchanged - close is equivalent to cancel
twotoneDataIn = getappdata(handles.figure1, 'twotoneDataIn');
twotoneDataOut = twotoneDataIn;
twotoneGeneralSettingIn = getappdata(handles.figure1, 'twotoneGeneralSettingIn');
twotoneGeneralSettingOut = twotoneGeneralSettingIn;
setappdata(handles.figure1, 'twotoneDataOut', twotoneDataOut);
setappdata(handles.figure1, 'twotoneGeneralSettingOut', twotoneGeneralSettingOut);

%update the old guidata handles with the new one
% for some reason you have to reference it as handles.figure1
% not entirely sure why
guidata(handles.figure1,handles);

% Running uiresume runs twotoneSettingsGui_OutputFcn, setting the output and closing the figure
uiresume(handles.figure1);


%-------------------------------------------
%-------------------------------------------
function figure1_CreateFcn(hObject, eventdata, handles)


%-------------------------------------------
%-------------------------------------------
function pushButtonOK_Callback(hObject, eventdata, handles)

%finished so update settings
twotoneDataOut = getappdata(handles.figure1, 'twotoneDataOut');
setappdata(handles.figure1, 'twotoneDataOut', twotoneDataOut);

%update the old guidata handles with the new one
% for some reason you have to reference it as handles.figure1
% not entirely sure why
guidata(handles.figure1,handles);

% Running uiresume runs twotoneSettingsGui_OutputFcn, setting the output and closing the figure
uiresume(handles.figure1);

%-------------------------------------------
%-------------------------------------------

function pushbuttonCancel_Callback(hObject, eventdata, handles)

%settings all unchanged -ie cancel
twotoneDataIn = getappdata(handles.figure1, 'twotoneDataIn');
twotoneDataOut = twotoneDataIn;
setappdata(handles.figure1, 'twotoneDataOut', twotoneDataOut);

%update the old guidata handles with the new one
% for some reason you have to reference it as handles.figure1
% not entirely sure why
guidata(handles.figure1,handles);

% Running uiresume runs twotoneSettingsGui_OutputFcn, setting the output and closing the figure
uiresume(handles.figure1);


%-------------------------------------------
function editDXStart_Callback(hObject, eventdata, handles)
twotoneDataOut = getappdata(handles.figure1, 'twotoneDataOut');
twotoneDataOut.settings.imageSettings.channelImageLim(1,1) = str2num(get(handles.editDXStart,'String'));
setappdata(handles.figure1, 'twotoneDataOut', twotoneDataOut);


%-------------------------------------------
function editDXStart_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%-------------------------------------------
function editDXEnd_Callback(hObject, eventdata, handles)
twotoneDataOut = getappdata(handles.figure1, 'twotoneDataOut');
twotoneDataOut.settings.imageSettings.channelImageLim(1,2) = str2num(get(handles.editDXEnd,'String'));
setappdata(handles.figure1, 'twotoneDataOut', twotoneDataOut);


%-------------------------------------------
function editDXEnd_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%-------------------------------------------
function editDYStart_Callback(hObject, eventdata, handles)
twotoneDataOut = getappdata(handles.figure1, 'twotoneDataOut');
twotoneDataOut.settings.imageSettings.channelImageLim(1,3) = str2num(get(handles.editDYStart,'String'));
setappdata(handles.figure1, 'twotoneDataOut', twotoneDataOut);


%-------------------------------------------
function editDYStart_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%-------------------------------------------
function editDYEnd_Callback(hObject, eventdata, handles)
twotoneDataOut = getappdata(handles.figure1, 'twotoneDataOut');
twotoneDataOut.settings.imageSettings.channelImageLim(1,4) = str2num(get(handles.editDYEnd,'String'));
setappdata(handles.figure1, 'twotoneDataOut', twotoneDataOut);


%-------------------------------------------
function editDYEnd_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%-------------------------------------------
function editAXStart_Callback(hObject, eventdata, handles)
twotoneDataOut = getappdata(handles.figure1, 'twotoneDataOut');
twotoneDataOut.settings.imageSettings.channelImageLim(2,1) = str2num(get(handles.editAXStart,'String'));
setappdata(handles.figure1, 'twotoneDataOut', twotoneDataOut);


%-------------------------------------------
function editAXStart_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%-------------------------------------------
function editAXEnd_Callback(hObject, eventdata, handles)
twotoneDataOut = getappdata(handles.figure1, 'twotoneDataOut');
twotoneDataOut.settings.imageSettings.channelImageLim(2,2) = str2num(get(handles.editAXEnd,'String'));
setappdata(handles.figure1, 'twotoneDataOut', twotoneDataOut);


%-------------------------------------------
function editAXEnd_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%-------------------------------------------
function editAYStart_Callback(hObject, eventdata, handles)
twotoneDataOut = getappdata(handles.figure1, 'twotoneDataOut');
twotoneDataOut.settings.imageSettings.channelImageLim(2,3) = str2num(get(handles.editAYStart,'String'));
setappdata(handles.figure1, 'twotoneDataOut', twotoneDataOut);


%-------------------------------------------
function editAYStart_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%-------------------------------------------
function editAYEnd_Callback(hObject, eventdata, handles)
twotoneDataOut = getappdata(handles.figure1, 'twotoneDataOut');
twotoneDataOut.settings.imageSettings.channelImageLim(2,4) = str2num(get(handles.editAYEnd,'String'));
setappdata(handles.figure1, 'twotoneDataOut', twotoneDataOut);


%-------------------------------------------
function editAYEnd_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%-------------------------------------------
function editBPDiameter_Callback(hObject, eventdata, handles)
twotoneDataOut = getappdata(handles.figure1, 'twotoneDataOut');
twotoneDataOut.settings.autoDetectSettings.bandPassKernelDiameter= str2num(get(handles.editBPDiameter, 'String'));
setappdata(handles.figure1, 'twotoneDataOut', twotoneDataOut);

%-------------------------------------------
function editBPDiameter_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%-------------------------------------------
function editWindowRadius_Callback(hObject, eventdata, handles)
twotoneDataOut = getappdata(handles.figure1, 'twotoneDataOut');
twotoneDataOut.settings.autoDetectSettings.fitSubImageRadius= str2num(get(handles.editWindowRadius, 'String'));
setappdata(handles.figure1, 'twotoneDataOut', twotoneDataOut);

%-------------------------------------------
function editWindowRadius_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%-------------------------------------------
function editFrameAvgMin_Callback(hObject, eventdata, handles)
twotoneDataOut = getappdata(handles.figure1, 'twotoneDataOut');
twotoneDataOut.settings.autoDetectSettings.averageFrameLim(1) = str2num(get(handles.editFrameAvgMin, 'String'));
setappdata(handles.figure1, 'twotoneDataOut', twotoneDataOut);


%-------------------------------------------
function editFrameAvgMin_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%-------------------------------------------
function editFrameAvgMax_Callback(hObject, eventdata, handles)
twotoneDataOut = getappdata(handles.figure1, 'twotoneDataOut');
twotoneDataOut.settings.autoDetectSettings.averageFrameLim(2) = str2num(get(handles.editFrameAvgMax, 'String'));
setappdata(handles.figure1, 'twotoneDataOut', twotoneDataOut);


%-------------------------------------------
function editFrameAvgMax_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1


% --- Executes on selection change in popupmenuAnalysisAlgorithm.
function popupmenuAnalysisAlgorithm_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuAnalysisAlgorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
twotoneDataOut = getappdata(handles.figure1, 'twotoneDataOut');
twotoneGeneralSettingOut= getappdata(handles.figure1, 'twotoneGeneralSettingOut');

%check which algorithm has been selected
currentAlgNo = get(handles.popupmenuAnalysisAlgorithm, 'Value');
currentAlgorithm = twotoneGeneralSettingOut.analysisAlgorithm{currentAlgNo};
twotoneDataOut.settings.fitSettings.algorithm = currentAlgorithm;
twotoneDataOut.settings.fitSettings.fitParamName = twotoneGeneralSettingOut.algorithmParam{currentAlgNo};
twotoneDataOut.settings.fitSettings.param= twotoneGeneralSettingOut.algorithmDefaultVal{currentAlgNo};

% update the algorithm options
set(handles.popupmenuFitAlgorithmParam, 'String',twotoneDataOut.settings.fitSettings.fitParamName);
set(handles.popupmenuFitAlgorithmParam, 'Value',1);
currentParamVal = twotoneDataOut.settings.fitSettings.param;
%update the current parameters so there is the correct number of elements
%set any new elements as zero
nParamName = numel(twotoneDataOut.settings.fitSettings.fitParamName);
nParamVal = numel(currentParamVal);
if nParamName > nParamVal
  currentParamVal(nParamName)	  = 0;% initialize the extra elements as 0
elseif nParamName < nParamVal
  currentParamVal(nParamName+1:end) = [];%delete the extra elements
end
%update the param values
twotoneDataOut.settings.fitSettings.param=currentParamVal;

set(handles.editFitAlgorithmParam,'String',num2str(currentParamVal(1)));

setappdata(handles.figure1, 'twotoneDataOut', twotoneDataOut);

% --- Executes during object creation, after setting all properties.
function popupmenuAnalysisAlgorithm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuAnalysisAlgorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuFitAlgorithmParam.
function popupmenuFitAlgorithmParam_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuFitAlgorithmParam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
twotoneDataOut   = getappdata(handles.figure1, 'twotoneDataOut');
twotoneGeneralSettingOut = getappdata(handles.figure1, 'twotoneGeneralSettingOut');

%check which algorithm has been selected
currentParamNo  = get(handles.popupmenuFitAlgorithmParam, 'Value');
currentParamVal = twotoneDataOut.settings.fitSettings.param;
set(handles.editFitAlgorithmParam,'String',num2str(currentParamVal(currentParamNo)));

setappdata(handles.figure1, 'twotoneDataOut', twotoneDataOut);

% --- Executes during object creation, after setting all properties.
function popupmenuFitAlgorithmParam_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuFitAlgorithmParam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editFitAlgorithmParam_Callback(hObject, eventdata, handles)
% hObject    handle to editFitAlgorithmParam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
twotoneDataOut   = getappdata(handles.figure1, 'twotoneDataOut');
twotoneGeneralSettingOut = getappdata(handles.figure1, 'twotoneGeneralSettingOut');

%update the algorithm param
currentParamNo  = get(handles.popupmenuFitAlgorithmParam, 'Value');
paramVal = str2num(get(handles.editFitAlgorithmParam,'String'));
twotoneDataOut.settings.fitSettings.param(currentParamNo) = paramVal;
%update the default alg param also
currentAlg = twotoneDataOut.settings.fitSettings.algorithm;
currentAlgNo = find(strcmp(twotoneGeneralSettingOut.analysisAlgorithm,currentAlg));
twotoneGeneralSettingOut.algorithmDefaultVal{currentAlgNo}(currentParamNo) = paramVal;

setappdata(handles.figure1, 'twotoneDataOut', twotoneDataOut);
setappdata(handles.figure1, 'twotoneGeneralSettingOut', twotoneGeneralSettingOut);

% --- Executes during object creation, after setting all properties.
function editFitAlgorithmParam_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFitAlgorithmParam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editFileAppendString_Callback(hObject, eventdata, handles)
% hObject    handle to editFileAppendString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
twotoneDataOut   = getappdata(handles.figure1, 'twotoneDataOut');
%update the algorithm param
fileAppendString= get(handles.editFileAppendString,'String');
if isempty(fileAppendString)%sanity check to make sure the data is not overwritten
  fileAppendString = twotoneDataOut.settings.fileAppendString ;
  set(handles.editFileAppendString,'String',fileAppendString);
end
twotoneDataOut.settings.fileAppendString = fileAppendString;
setappdata(handles.figure1, 'twotoneDataOut', twotoneDataOut);

% --- Executes during object creation, after setting all properties.
function editFileAppendString_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFileAppendString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editPSFestimate_Callback(hObject, eventdata, handles)


function editPSFestimate_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------
function pushbuttonAutoSetParam_Callback(hObject, eventdata, handles)
% update the parameters based on the automatic estimate
twotoneDataOut   = getappdata(handles.figure1, 'twotoneDataOut');
twotoneGeneralSettingOut = getappdata(handles.figure1, 'twotoneGeneralSettingOut');

PSFwidthEstimate = str2num(get(handles.editPSFestimate,'String'));
twotoneGeneralSettingOut.PSFwidthEstimate =PSFwidthEstimate ;
[twotoneDataOut, twotoneGeneralSettingOut] = autoGenerateFitParam(twotoneDataOut, twotoneGeneralSettingOut, PSFwidthEstimate);

%update the value of the currently selected fitting algorithm parmater
currentAlgParamNo = get(handles.popupmenuFitAlgorithmParam,'Value');
currentValue = twotoneDataOut.settings.fitSettings.param(currentAlgParamNo);
set(handles.editFitAlgorithmParam,'String',num2str(currentValue));
%update the rest
set(handles.editBPDiameter,'String',num2str(twotoneDataOut.settings.autoDetectSettings.bandPassKernelDiameter))
set(handles.editWindowRadius,'String',num2str(twotoneDataOut.settings.autoDetectSettings.fitSubImageRadius))

setappdata(handles.figure1, 'twotoneDataOut', twotoneDataOut);
setappdata(handles.figure1, 'twotoneGeneralSettingOut', twotoneGeneralSettingOut);

%--------------------------------------------------------------
%--------------------------------------------------------------
%--------------------------------------------------------------
%--------------------------------------------------------------
function [twotoneData, twotoneGeneralSetting]= autoGenerateFitParam(twotoneData, twotoneGeneralSetting,width)
% function autoParamOut = autoGenerateFitParam(width)
%  calculate default values for the paramters based on 
%  the size of the PSF

autoParam = estimateParam(width);
[twotoneData, twotoneGeneralSetting] = updateParam(twotoneData, twotoneGeneralSetting,autoParam );


%-----------------------------------------------
function autoParam = estimateParam(width)
% calculate the appropriate settings

autoParam.autoDetWindowRadius = max(2, round(2*width));
autoParam.profileFittingRadius = max(2, round(4*width));
autoParam.profileFittingMaxPosChange = autoParam.profileFittingRadius/2;
fitMinWidth = max(0,width - 0.5);
fitMaxWidth = width + 0.5;
autoParam.widthLim = [fitMinWidth fitMaxWidth];
autoParam.nnLim = 4*width;

autoParam.aperturePhotInner = floor(3*width);
autoParam.aperturePhotOuter = round(sqrt(2)*3*width); %the sqrt 2 factor gives roughly equal areas for the inner and outer radii
if autoParam.aperturePhotInner == autoParam.aperturePhotOuter
  autoParam.aperturePhotOuter = autoParam.aperturePhotOuter + 1;
end

autoParam.aperturePhotFrameAvg = 3; %default setting
autoParam.bpasskerneldiametre = max(round(2*width),1);

%-----------------------------------------------
function [twotoneData, twotoneGeneralSetting] = updateParam(twotoneData, twotoneGeneralSetting,autoParam)


%update twotoneGeneralSetting first

%update each algorithms parameter set
%update fixedGaussEllipse
currentAlg='fixedGaussEllipse';
cellNo = find(strcmp(twotoneGeneralSetting.analysisAlgorithm,currentAlg));
currentParam = 'windowRadius';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) = autoParam.profileFittingRadius;
currentParam = 'minWidth';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) = autoParam.widthLim(1);
currentParam = 'maxWidth';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) = autoParam.widthLim(2);

%update fixedGauss
currentAlg='fixedGauss';
cellNo = find(strcmp(twotoneGeneralSetting.analysisAlgorithm,currentAlg));
currentParam = 'windowRadius';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) = autoParam.profileFittingRadius;
currentParam = 'minWidth';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) = autoParam.widthLim(1);
currentParam = 'maxWidth';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) = autoParam.widthLim(2);

%update freeGaussEllipse
currentAlg='freeGaussEllipse';
cellNo = find(strcmp(twotoneGeneralSetting.analysisAlgorithm,currentAlg));
currentParam = 'windowRadius';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) = autoParam.profileFittingRadius;
currentParam = 'minWidth';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) = autoParam.widthLim(1);
currentParam = 'maxWidth';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) = autoParam.widthLim(2);
currentParam = 'maxPositionChange';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) = autoParam.profileFittingMaxPosChange;

%update fixedGauss
currentAlg='ring';
cellNo = find(strcmp(twotoneGeneralSetting.analysisAlgorithm,currentAlg));
currentParam = 'innerCircleRadius';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) =autoParam.aperturePhotInner;
currentParam = 'outerCircleRadius';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) = autoParam.aperturePhotOuter;
currentParam = 'ringTime';
paramRow = find(strcmp(twotoneGeneralSetting.algorithmParam{cellNo},currentParam));
twotoneGeneralSetting.algorithmDefaultVal{cellNo}(paramRow) =autoParam.aperturePhotFrameAvg;

% now update twotoneData fit algorithm settings
selectedAlgNo = find(strcmp(twotoneGeneralSetting.analysisAlgorithm,twotoneData.settings.fitSettings.algorithm));
twotoneData.settings.fitSettings.param=twotoneGeneralSetting.algorithmDefaultVal{selectedAlgNo};

%now update the other twotoneData settings
twotoneData.settings.autoDetectSettings.fitSubImageRadius	= autoParam.autoDetWindowRadius ;
twotoneData.settings.autoDetectSettings.bandPassKernelDiameter	= autoParam.bpasskerneldiametre ;
twotoneData.settings.autoDetectSettings.nearestNeighbor.thresh	= autoParam.nnLim;
twotoneData.settings.autoDetectSettings.PSFwidth.lim		= autoParam.widthLim;
