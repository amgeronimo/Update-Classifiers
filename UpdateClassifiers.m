function varargout = UpdateClassifiers(varargin)
% UPDATECLASSIFIERS MATLAB code for UpdateClassifiers.fig
%      UPDATECLASSIFIERS, by itself, creates a new UPDATECLASSIFIERS or raises the existing
%      singleton*.
%
%      H = UPDATECLASSIFIERS returns the handle to a new UPDATECLASSIFIERS or the handle to
%      the existing singleton*.
%
%      UPDATECLASSIFIERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UPDATECLASSIFIERS.M with the given input arguments.
%
%      UPDATECLASSIFIERS('Property','Value',...) creates a new UPDATECLASSIFIERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before UpdateClassifiers_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to UpdateClassifiers_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UpdateClassifiers

% Last Modified by GUIDE v2.5 22-Feb-2017 16:24:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @UpdateClassifiers_OpeningFcn, ...
    'gui_OutputFcn',  @UpdateClassifiers_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...ui
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


% --- Executes just before UpdateClassifiers is made visible.
function UpdateClassifiers_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to UpdateClassifiers (see VARARGIN)





% Choose default command line output for UpdateClassifiers
handles.output = hObject;
handles.pushbutton1.Enable = 'off';
handles.listbox1.Enable = 'off';
if handles.nodir==1
    handles.pushbutton3.Enable = 'off';
    set(handles.edit2,'BackgroundColor',[1 .6 .6]);
else
    set(handles.edit2,'BackgroundColor',[1 1 1]);
    set(handles.edit1,'Enable','inactive');
end
% Update handles structure
guidata(hObject, handles);


% UIWAIT makes UpdateClassifiers wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = UpdateClassifiers_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
handles.participant_id = get(hObject,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
handles.participant_id = 'Subject Code';
guidata(hObject, handles);
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fileidx=get(handles.listbox1,'value');
handles.csess = {}; handles.crun = {};
for i = 1:length(fileidx)
    handles.csess = [handles.csess handles.sessions{fileidx(i)}];
    handles.crun = [handles.crun {{handles.runs{fileidx(i)}}}];
    if i == length(fileidx)
        handles.pushbutton1.Enable = 'on';
    end
end



%     csess =
% selectedfiles =
guidata(hObject, handles);


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




% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles) %CREATE CLASSIFIER
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



BuildClassifier=1;
Cfier = 1;
Batch = 1;
RerefVal = handles.checkbox3.Value;
Art = handles.checkbox2.Value;
Figures_On = 0;
DoCSP = 0;
CSession = [];
CRun = [];

try
    if ismember(handles.ctype,1:5)
        [CM, allsums, targetstim, WORD_Acc, skiptest, TTSpell, ChosenLetterF,...
            TrainD,ClassLabel, ChanLabel, SC, SR, ttimes, prmfile, eyye, range, NotDownD] = ...
            BuildClassifier_teleBCI7(handles.participant_id,handles.csess,...
            handles.crun,BuildClassifier,Cfier,...
            Batch,RerefVal,Art,DoCSP,CSession,CRun, Figures_On, handles.dataloc, handles.ctype);
        
        if Figures_On == 1
            cpos = [100 60 150 30];
        else
            cpos = get(gcf,'Position');
        end
        cfig = figure('Name',['Classifier for ' handles.participant_id],'units','characters',...
            'Position',[cpos(1),cpos(2)-33,cpos(3),30],'menubar','none');
        tpts = size(TrainD,2)/(length(unique(ChanLabel)));
        
        lincol = [1 0 0;0 0 1];
        for i=1:8
            xi(i) = subplot(3,4,i);
            plot(ttimes/handles.params.fs,nanmean(TrainD(ClassLabel(ChanLabel==1)==2,tpts*(i-1)+1:tpts*i)),'r');
            hold on;
            plot(ttimes/handles.params.fs,nanmean(TrainD(ClassLabel(ChanLabel==1)==1,tpts*(i-1)+1:tpts*i)),'b');
            %           plot(range/handles.params.fs,nanmean(NotDownD(ClassLabel(ChanLabel==1)==2,length(range)*(i-1)+1:length(range)*i)),'Color',[1 .5 .5]);
            %         plot(range/handles.params.fs,nanmean(NotDownD(ClassLabel(ChanLabel==1)==1,length(range)*(i-1)+1:length(range)*i)),'Color',[.5 .5 1]);
            cfier = CM{end};
            
            if i>1, set(gca,'Xticklabel',[],'Yticklabel',[]); end;
            plims(i,:) = get(gca,'ylim');
            if isempty(cfier)
            else
                maxc = max(abs(cfier(:,4)));
                tmp = find(cfier(:,1) == i);
                for j = tmp'
                    if sign(cfier(j,4))==1
                        line([1; 1]*cfier(j,2)'/handles.params.fs,[-50 50],'Linewidth',2,'Color',...
                            [1-abs(cfier(j,4))/maxc 1-abs(cfier(j,4))/maxc 1]);
                    else
                        line([1; 1]*cfier(j,2)'/handles.params.fs,[-50 50],'Linewidth',2,'Color',...
                            [1 1-abs(cfier(j,4))/maxc 1-abs(cfier(j,4))/maxc]);
                    end
                end
            end
        end
        linkaxes(xi); ylim(minmax(plims(:)'))
        subplot(3,4,9); bar([2:2:length(WORD_Acc)*2],WORD_Acc); xlabel('Number of Sequences');
        text(1,120,sprintf('%s',['Parameter file written to: ' prmfile]),'FontSize',8);
        ylabel('Accuracy (%)'); xlim([1 length(WORD_Acc)*2+1]); ylim([0 100])
        if ~isempty(eyye)
            subplot(3,4,10:12); bar([eyye.gazeacc sqrt(eyye.gazevar)]);
            legend('Acc_X','Acc_y','Dev_x','Dev_y','Orientation','Horizontal','Location','Best');
            xlabel('Gaze statistics');
        end
        
    else
        
        [CM, acc, TrainD, ClassLabel, ChanLabel,freqs, prmfile, eyye] = ...
            BuildMIClassifier_teleBCI2(handles.participant_id,handles.csess,...
            handles.crun,BuildClassifier,Cfier,...
            Batch,RerefVal,Art,DoCSP,CSession,CRun, Figures_On, handles.dataloc, handles.ctype);
        
        if Figures_On == 1
            cpos = [100 60 150 30];
        else
            cpos = get(gcf,'Position');
        end
        cfig = figure('Name',['Classifier for ' handles.participant_id],'units','characters',...
            'Position',[cpos(1),cpos(2)-33,cpos(3),30],'menubar','none');
        fpts = size(TrainD,2)/(length(unique(ChanLabel)));
        
        lincol = [1 0 0;0 0 1];
        for i=1:8
            subplot(3,4,i);
            plot(freqs,nanmean(TrainD(ClassLabel==-1,fpts*(i-1)+1:fpts*i)),'r');
            hold on;
            plot(freqs,nanmean(TrainD(ClassLabel==1,fpts*(i-1)+1:fpts*i)),'b');
            cfier = CM{end};
            if i>1, set(gca,'Xticklabel',[],'Yticklabel',[]); end;
            tmp = find(cfier(:,1) == i); plims = get(gca,'ylim');
            for j = tmp'
                line([1; 1]*cfier(j,2)',plims,'Linewidth',2,'Color',lincol(sign(cfier(j,4))*.5+1.5,:))
            end
        end
        subplot(3,4,9); bar([.2 .4 .6 .8 1],acc*100); xlabel('Fraction of windows user');
        text(1,120,sprintf('%s',['Parameter file written to: ' prmfile]),'FontSize',8);
        ylabel('Accuracy (%)'); xlim([0 1.2]); ylim([0 100])
        if ~isempty(eyye)
            subplot(3,4,10:12); bar([eyye.gazeacc sqrt(eyye.gazevar)]);
            legend('Acc_X','Acc_y','Dev_x','Dev_y','Orientation','Horizontal','Location','Best');
            xlabel('Gaze statistics');
        end
    end
catch
    close(gcf)
    msgbox('Mismatched parameters or invalid files chosen');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles) %CLEAR SELECTIONS
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.listbox1, 'Value', []);
set(handles.pushbutton1,'Enable','off');

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles) %FIND FILES
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.listbox1,'Value',[]);
set(handles.pushbutton1,'Enable','off');
if isempty(get(handles.edit1,'String'))
    set(handles.listbox1,'String','<<<Enter a Subject Code','Enable','off');
else
    
    if handles.radiobutton1.Value==1, handles.ctype = 1; end;
    if handles.radiobutton7.Value==1, handles.ctype = 2; end;
    if handles.radiobutton8.Value==1, handles.ctype = 3; end;
    if handles.radiobutton9.Value==1, handles.ctype = 4; end;
    if handles.radiobutton11.Value==1, handles.ctype = 5; end;
    if handles.radiobutton13.Value==1, handles.ctype = 6; end;
    
    [handles.filestring, handles.viability, handles.errorcode,...
        handles.sessions, handles.runs, handles.date] = FindFiles(handles);
    
    goodjazzystring = cell(1,sum(handles.viability>0));
    invalidjazzystring = cell(1,sum(handles.viability==0));
    mm=1; nn=1; invind = []; goodind = [];
    [hsort, hvind] = sort(handles.viability);
    viarange = hsort(hsort~=0);
    if isempty(viarange)
    else
        viarange = diff(minmax(viarange'));
    end
    for i = flipud(hvind)'%1:length(handles.filestring)
        if handles.viability(i) == 0 %Bad files
            invalidjazzystring{mm} = ['<html><FONT color=#D3D3D3><b> S' ...
                handles.sessions{i} 'R' handles.runs{i} ' </b>' handles.date{i}...
                '<FONT color=#FF0000>' handles.errorcode{i} '</html>'];
            mm=mm+1; invind = [invind i];
        else
            clr = floor(255*([(max(hsort)-handles.viability(i))/viarange...
                1 (max(hsort)-handles.viability(i))/viarange]));
            tmp = dec2hex(clr)'; tmp = ['#' tmp(:)'];
            goodjazzystring{nn} = ['<html><body bgcolor=' tmp '><b> S' ...
                handles.sessions{i} 'R' handles.runs{i} ' </b>' handles.date{i} '</html>'];
            nn=nn+1; goodind = [goodind i];
        end
    end
    
    [recentfile,recentind] = max(goodind);
    nowjazzystring = goodjazzystring(recentind);
    goodjazzystring(recentind) = []; goodind(recentind) = [];
    jazzystring = [nowjazzystring goodjazzystring  invalidjazzystring];
    fileshuffle = [recentfile goodind invind];
    handles.filestring = handles.filestring(fileshuffle);
    handles.viability = handles.viability(fileshuffle);
    handles.errorcode = handles.errorcode(fileshuffle);
    handles.sessions = handles.sessions(fileshuffle);
    handles.runs = handles.runs(fileshuffle);
    handles.date = handles.date(fileshuffle);
    
    if isempty(jazzystring)
        set(handles.listbox1,'String',['No files found of "'...
            handles.participant_id '"!'],'Enable','off');
        
    else
        set(handles.listbox1,'String',jazzystring,'Enable','on');
    end
end
handles.output = hObject;
guidata(hObject, handles);




function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
handles.dataloc = get(hObject,'String')
if ~isdir(handles.dataloc)
    set(hObject,'BackgroundColor',[1 .6 .6]);
    set(handles.pushbutton3,'Enable','off');
    set(handles.pushbutton1,'Enable','off');
    set(handles.listbox1,'Enable','off');
    set(handles.edit1,'Enable','off','String','Subject Code');
else
    set(hObject,'BackgroundColor',[1 1 1]);
    set(handles.pushbutton3,'Enable','on');
    set(handles.edit1,'Enable','inactive');
end
guidata(hObject, handles);


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


pcusr=getenv('USERNAME');
[~, inst] = regexp(cd,'\Program Files\');
[~, inst2] = regexp(cd,'\Desktop');
if isempty(inst) && isempty(inst2) %testing
    handles.dataloc= ['C:\Users\' pcusr '\Documents\BCI2000']
else %installed
    %Identify the user (BCIcomputer01 eg.), specify associated Box folder
    if strcmp(pcusr,'ageronimo') || strcmp('BCIcomputer01',pcusr) % For testing on 'ageronimo' and 'BCIcomputer'
        busr = 'Shared teleBCI Folders\';
    else
        busr = '';
    end
    fldr = dir(['C:\Users\' pcusr '\Box Sync\' busr]);
    fldr_loc = find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'teleBCI')),{fldr.name},'UniformOutput',false)));
    if isempty(fldr_loc) %Maybe using the miniputer?
        fldr_loc = find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'BCIminiputer')),{fldr.name},'UniformOutput',false)));
    end
    handles.dataloc = ['C:\Users\' pcusr '\Box Sync\'  busr fldr(fldr_loc(1)).name '\BCI2000']
end

set(hObject,'String',handles.dataloc);
if ~isdir(handles.dataloc)
    set(hObject,'BackgroundColor',[1 .6 .6]);
    handles.nodir = 1;
else
    handles.nodir = 0;
end
guidata(hObject, handles);



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
handles.params.fs = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
handles.params.fs = str2double(get(hObject,'String'));
guidata(hObject, handles);
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
handles.params.sbs = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
handles.params.sbs = str2double(get(hObject,'String'));
guidata(hObject, handles);
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
handles.params.numch = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
handles.params.numch= str2double(get(hObject,'String'));
guidata(hObject, handles);
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end







function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
handles.params.timeon = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
handles.params.timeon = str2double(get(hObject,'String'));
guidata(hObject, handles);
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
handles.params.timeoff = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
handles.params.timeoff = str2double(get(hObject,'String'));
guidata(hObject, handles);
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit1.
function edit1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.edit1,'Enable'),'inactive')
    set(handles.edit1,'Enable','on','String','');
    uicontrol(handles.edit1);
end


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in radiobutton11.
function radiobutton11_Callback(hObject, eventdata, handles)
set(handles.uibuttongroup1,'Visible','on');
handles.edit7.String = '0.5'; handles.params.timeon = 0.5;
handles.edit8.String = '0.2'; handles.params.timeoff = 0.2;
guidata(hObject, handles);

function radiobutton1_Callback(hObject, eventdata, handles)
set(handles.uibuttongroup1,'Visible','on');
handles.edit7.String = '0.1'; handles.params.timeon = 0.1;
handles.edit8.String = '0.1'; handles.params.timeoff = 0.1;
guidata(hObject, handles);

function radiobutton7_Callback(hObject, eventdata, handles)
set(handles.uibuttongroup1,'Visible','on');
handles.edit7.String = '0.1'; handles.params.timeon = 0.1;
handles.edit8.String = '0.1'; handles.params.timeoff = 0.1;
guidata(hObject, handles);

function radiobutton8_Callback(hObject, eventdata, handles)
set(handles.uibuttongroup1,'Visible','on');
handles.edit7.String = '0.1'; handles.params.timeon = 0.1;
handles.edit8.String = '0.1'; handles.params.timeoff = 0.1;
guidata(hObject, handles);

function radiobutton9_Callback(hObject, eventdata, handles)
set(handles.uibuttongroup1,'Visible','on');
handles.edit7.String = '0.1'; handles.params.timeon = 0.1;
handles.edit8.String = '0.1'; handles.params.timeoff = 0.1;
guidata(hObject, handles);

function radiobutton13_Callback(hObject, eventdata, handles)
handles.edit7.String = '0.1'; handles.params.timeon = 0.1;
handles.edit8.String = '0.1'; handles.params.timeoff = 0.1;
set(handles.uibuttongroup1,'Visible','off');
guidata(hObject, handles);
