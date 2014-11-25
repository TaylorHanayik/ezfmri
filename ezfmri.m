function varargout = ezfmri(varargin)
% EZFMRI MATLAB code for ezfmri.fig
%      EZFMRI, by itself, creates a new EZFMRI or raises the existing
%      singleton*.
%
%      H = EZFMRI returns the handle to a new EZFMRI or the handle to
%      the existing singleton*.
%
%      EZFMRI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EZFMRI.M with the given input arguments.
%
%      EZFMRI('Property','Value',...) creates a new EZFMRI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ezfmri_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ezfmri_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ezfmri

% Last Modified by GUIDE v2.5 19-Nov-2014 10:08:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ezfmri_OpeningFcn, ...
                   'gui_OutputFcn',  @ezfmri_OutputFcn, ...
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


% --- Executes just before ezfmri is made visible.
function ezfmri_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ezfmri (see VARARGIN)

% Choose default command line output for ezfmri
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ezfmri wait for user response (see UIRESUME)
% uiwait(handles.ezfmrifig);


% --- Outputs from this function are returned to the command line.
function varargout = ezfmri_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in studyFolderBtn.
function studyFolderBtn_Callback(hObject, eventdata, handles)
%SYSTEM GUI for selecting a folder
d.studyFolder = uigetdir(pwd,'Choose your experiment folder');
h = findobj('Tag','studyFolderLabel');
set(h,'String',d.studyFolder);
%LIST CONTENTS OF THE FOLDER (UNIX AND MAC OSX ONLY)
[d.status, d.cmdout] = system(sprintf('find "%s" -iname "*.nii" ! -name ".*.nii"',d.studyFolder));
h = findobj('Tag','studyFolderText');
%SPLIT d.cmdout BASED ON NEWLINE CHARACTER
files = regexp(d.cmdout,'\n','split');
%FILL IN THE FILE LIST BOX BY NUMBER OF LINES IN files
for i = 1:size(files,2)
    file = char(files(1,i));
    if i == 1
        set(h,'String',{file});
    else
        oldString = get(h,'String');
        set(h,'String',{char(oldString);file});
    end
end
d.allFiles = get(h,'String');
%SAVE TO d IN GUIDATA
guidata(hObject,d);













