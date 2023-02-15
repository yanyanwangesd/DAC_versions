function varargout = dac_set_field(varargin)
% DAC_SET_FIELD MATLAB code for dac_set_field.fig
%      DAC_SET_FIELD, by itself, creates a new DAC_SET_FIELD or raises the existing
%      singleton*.
%
%      H = DAC_SET_FIELD returns the handle to a new DAC_SET_FIELD or the handle to
%      the existing singleton*.
%
%      DAC_SET_FIELD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DAC_SET_FIELD.M with the given input arguments.
%
%      DAC_SET_FIELD('Property','Value',...) creates a new DAC_SET_FIELD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dac_set_field_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dac_set_field_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dac_set_field

% Last Modified by GUIDE v2.5 04-Feb-2018 15:46:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dac_set_field_OpeningFcn, ...
                   'gui_OutputFcn',  @dac_set_field_OutputFcn, ...
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


% --- Executes just before dac_set_field is made visible.
function dac_set_field_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dac_set_field (see VARARGIN)

% Choose default command line output for dac_set_field
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global hdx
hdx=struct('xsize',1,...
          'ysize',1,...
          'xnumx',1,...
          'ynumy',1,...
          'numpol',0,...
          'curpol',0,...
          'polstate',0,...
          'poly',cell(1),...
          'hpoly',0,...
          'hmesh',0,...
          'hgrid',0,...
          'hxmax',0,...
          'gridnx',1,...
          'gridny',1,...
          'grid',zeros(1),...
          'filename','dac_raster1.txt',...
          'fnnode',1,...
          'fpoints',zeros(1,2),...
          'fntri',1 ...
);


% UIWAIT makes dac_set_field wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dac_set_field_OutputFcn(hObject, eventdata, handles)
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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdx
vl=num2str(get(handles.edit4,'String'));
hdx.numpol=hdx.numpol+1;
hdx.poly{hdx.numpol,1}=vl;
hdx.polstate=1;

set(handles.togglebutton1,'Value',1);
x=[];
y=[];
hloc=0;
hold on
pp=0;
while(hdx.polstate==1)
    hdx.polstate=get(handles.togglebutton1,'Value');
    [xa,ya]=ginput(1);
    pp=pp+1;
    hdx.polstate=get(handles.togglebutton1,'Value');
    if(hdx.polstate==1)
        x=[x,xa];
        y=[y,ya];

        if(xa<handles.axes1.XLim(2)*1.05)
            hloc(pp)=plot(x,y,'r');
        else
            hdx.polstate=1;
            break
        end
    else
        hdx.polstate=1;
        break;
    end
end
hold off
    delete(hloc);
set(handles.togglebutton1,'Value',0);

x(numel(x))=x(1);
y(numel(y))=y(1);

hdx.poly{hdx.numpol,2}=x;
hdx.poly{hdx.numpol,3}=y;

plot_poly(hObject, eventdata, handles, hdx.numpol);

update_polytable(hObject, eventdata, handles)



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdx
hdx.filename = get(handles.edit3,'String');
if(exist(hdx.filename,'file')>0)
    btn=questdlg('Overwrite existing file?',['File ',hdx.filename,' exists'],'Overwrite','Return','Return');
    if(btn==1)
        return
    else
        write_matrix(hObject, eventdata, handles);
    end
else
    write_matrix(hObject, eventdata, handles);
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close()



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


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


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdx
hdx.gridnx = str2num(get(handles.edit5,'String'));
hdx.gridny = str2num(get(handles.edit6,'String'));

preview_grid(hObject, eventdata, handles)



% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdx
hdx.xsize=str2num(get(handles.edit1, 'String'));
hdx.ysize=str2num(get(handles.edit2, 'String'));
hdx.xc=str2num(get(handles.edit7, 'String'));

hdx.xnumx = ceil(hdx.xsize/1000/1.5 *700/hdx.xc)-1;
hdx.ynumy = ceil(hdx.ysize/1000/1.5 *700/hdx.xc)-1;

nodestr=[num2str(hdx.xnumx),'x',num2str(hdx.ynumy)];
set(handles.text9,'string',nodestr);

make_tri(hObject, eventdata, handles)


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


function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


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
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

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

global hdx
hdx.polstate=0;


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1



% --- Executes when selected cell(s) is changed in uitable1.
function uitable1_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
global hdx
idx=eventdata.Indices;
if(numel(idx)==0)
   return
end

row=idx(1);
hdx.curpol=row;
if(row>=0)
    plot_poly(hObject, eventdata, handles, row)
end


% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
global hdx
idx=eventdata.Indices;
row=idx(1);
vl = str2num(eventdata.EditData);
hdx.poly{row,1}=vl;



% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hdx

if(hdx.curpol>hdx.numpol || hdx.polstate==1)
    return
end

%shadow structure
j=1;
for i=1:hdx.numpol
   if(i~=hdx.curpol)
       spoly{j,1}=hdx.poly{i,1};
       spoly{j,2}=hdx.poly{i,2};
       spoly{j,3}=hdx.poly{i,3};
       hploc(j)=hdx.hpoly(i);
       j=j+1;
   else
       delete(hdx.hpoly(i));
   end
end

for i=1:hdx.numpol-1
    hdx.poly{i,1}=spoly{i,1};
    hdx.poly{i,2}=spoly{i,2};
    hdx.poly{i,3}=spoly{i,3};
    hdx.hpoly(i)=hploc(i);
end
hdx.numpol=hdx.numpol-1;
update_polytable(hObject, eventdata, handles)
plot_poly(hObject, eventdata, handles, hdx.numpol)














function make_tri(hObject, eventdata, handles)
global hdx
xres= hdx.xsize/(hdx.xnumx-1);
yres= hdx.ysize/(hdx.ynumy-1);
x=[0:hdx.xnumx-1] * xres;
y=[0:hdx.ynumy-1] * yres;
[X,Y]=meshgrid(x,y);
hdx.fnnode=numel(X);
np=1;
for i =1:1:hdx.xnumx
    for j=1:1:hdx.ynumy
        xr=x(i);
        yr=y(j);
        if(i>1 && i<hdx.xnumx && j>1 && j<hdx.ynumy )
            xr=xr + xres*0.33*(rand(1)-0.5);
            yr=yr + yres*0.33*(rand(1)-0.5);
        end
        hdx.fpoints(np,1)=xr;
        hdx.fpoints(np,2)=yr;
        np=np+1;
    end
end
plot_mesh(hObject, eventdata, handles)



function plot_mesh(hObject, eventdata, handles)
global hdx

for i=1:numel(hdx.hpoly)
    h=  hdx.hpoly(i);
    if( ishandle( h ) && h>0)
        delete(hdx.hpoly(i));
    end
end


hdx.fntri=delaunay(hdx.fpoints(1:hdx.fnnode,1), hdx.fpoints(1:hdx.fnnode,2));
plot(handles.axes1,hdx.fpoints(1:hdx.fnnode,1), hdx.fpoints(1:hdx.fnnode,2),'r.')
hdx.hmesh=triplot(hdx.fntri,hdx.fpoints(1:hdx.fnnode,1), hdx.fpoints(1:hdx.fnnode,2),'Color',[1 1 1]*0.7 );
axis ij
axis(handles.axes1, [-0.1*hdx.xsize hdx.xsize*1.1 -0.1*hdx.ysize hdx.ysize*1.1])
hdx.hxmax=hdx.xsize*1.1;
axis(handles.axes1, 'equal')

if(hdx.numpol>0)
    plot_poly(hObject, eventdata, handles, hdx.numpol)
end



function plot_poly(hObject, eventdata, handles, active)
global hdx

k=min(hdx.numpol,numel(hdx.hpoly));
for i=1:k
    h=  hdx.hpoly(i);
    if( ishandle( h ) && h>0)
        delete(hdx.hpoly(i));
    end
end


hold on
for i=1:hdx.numpol
    x=hdx.poly{i,2};
    y=hdx.poly{i,3};
    %hdx.hpoly(i)=plot(x,y,'Color',[1 1 1]*0.3);
    hdx.hpoly(i)=patch('XData',x,'YData',y);
    set(hdx.hpoly(i),'FaceAlpha',0.6)
    set(hdx.hpoly(i),'FaceColor',[1 1 1]*0.4)
    if(i==active)
        set(hdx.hpoly(i),'FaceColor',[1 0 0 ]);
    end
end
hold off


function update_polytable(hObject, eventdata, handles)
global hdx
for i=1:hdx.numpol
    handles.uitable1.Data{i,1}=i;
    handles.uitable1.Data{i,2}=hdx.poly{i,1};
end
for i=hdx.numpol+1:20
    handles.uitable1.Data{i,1}=[];
    handles.uitable1.Data{i,2}=[];
end


function preview_grid(hObject, eventdata, handles)
global hdx
hdx.grid = zeros(hdx.gridnx,hdx.gridny);
xres = hdx.xsize/(hdx.gridnx-1);
x = [0:xres:hdx.xsize];
yres = hdx.ysize/(hdx.gridny-1);
y = [0:yres:hdx.ysize];

for k=1:hdx.numpol
    for j=1:hdx.gridny
        for i=1:hdx.gridny
            in = inpolygon(x(i),y(j), hdx.poly{k,2}, hdx.poly{k,3});
            if(in==1)
                hdx.grid(j,i)=hdx.poly{k,1};
            end
        end
    end
end


it = get(handles.checkbox1,'Value');
if(it==1)
    figure(9)
    hps=surf(x,y,hdx.grid);
    shading interp
    colormap(jet)
    axis ij
    view([14,56])
else
    hold on
    hp=pcolor(handles.axes1,x,y,hdx.grid);
    shading flat
    hold off
    set(hp,'FaceAlpha','0.4')
    colormap(handles.axes1,jet)
    pause(3)
    delete(hp)
    clc
    if(hdx.numpol>0)
        plot_poly(hObject, eventdata, handles, hdx.numpol)
    end
end





function write_matrix(hObject, eventdata, handles)
global hdx

hdx.gridnx = str2num(get(handles.edit5,'String'));
hdx.gridny = str2num(get(handles.edit6,'String'));

hdx.grid = zeros(hdx.gridnx,hdx.gridny);
xres = hdx.xsize/(hdx.gridnx-1);
x = [0:xres:hdx.xsize];
yres = hdx.ysize/(hdx.gridny-1);
y = [0:yres:hdx.ysize];

for k=1:hdx.numpol
    for j=1:hdx.gridny
        for i=1:hdx.gridny
            in = inpolygon(x(i),y(j), hdx.poly{k,2}, hdx.poly{k,3});
            if(in==1)
                hdx.grid(j,i)=hdx.poly{k,1};
            end
        end
    end
end

grid = hdx.grid;
save(hdx.filename,'grid','-ascii')
write_sample_snippet;



function write_sample_snippet
  global hdx

  flt_len = length(hdx.filename);
  filename = [hdx.filename(1:flt_len-4),'_sample_snippet.txt'];
  fl = fopen(filename,'w');

  fprintf(fl,'//---------------- raster block \n');
  fprintf(fl,'/f_input_type\n');
  fprintf(fl,'R\n');  
  fprintf(fl,'/raster filename\n');
  fprintf(fl,'%s\n',hdx.filename);  
  fprintf(fl,'/f_raster_format : 0 - plain text,  1 - asc, 2 - csv\n');
  fprintf(fl,'0\n');
  fprintf(fl,'/f_variable_determined  - choice:  u,v,w, p,k\n');
  fprintf(fl,'<YOUR CHOICE>\n');
  fprintf(fl,'/superposition  0 - overwrite, 1 - add, 2 - multiply\n');
  fprintf(fl,'<YOUR CHOICE>\n');
  fprintf(fl,'/time bounds [a], if any -  write 0, 0 for eternity; e.g. 0, 1.d6 makes condition only valid for first 1 Ma\n');
  fprintf(fl,'<YOUR CHOICE>,<YOUR CHOICE>\n');
  fprintf(fl,'/margin exceedance - use to set conditions where DAC is advected to. Left, right, front(bottom), back(top)\n');
  fprintf(fl,'0, 0, 0, 0\n');
  fprintf(fl,'//---------------- raster block\n');

  fclose(fl);
