function varargout = MetricTensGUI(varargin)
% MetricTens (Version: 1.0)-Created by Giovanni Esteves
% Department of Materials Science and Engineering
% North Carolina State University, Nov. 29th, 2015
% Email: gesteve@ncsu.edu

% Anytime this GUI is updated the version info at the very bottom that
% correlates to the about section, should also be upated.

%For those interested in further adding features to the GUI, you can
% send me the version of both .fig and .m files with the implemented features along with your name
% and I will integrate the file and add your name to the .m file and
% highlight the contributions made.

% METRICTENSGUI MATLAB code for MetricTensGUI.fig
%      METRICTENSGUI, by itself, creates a new METRICTENSGUI or raises the existing
%      singleton*.
%
%      H = METRICTENSGUI returns the handle to a new METRICTENSGUI or the handle to
%      the existing singleton*.
%
%      METRICTENSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in METRICTENSGUI.M with the given input arguments.
%
%      METRICTENSGUI('Property','Value',...) creates a new METRICTENSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MetricTensGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MetricTensGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MetricTensGUI

% Last Modified by GUIDE v2.5 07-Apr-2016 13:04:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MetricTensGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @MetricTensGUI_OutputFcn, ...
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

% --- Executes just before MetricTensGUI is made visible.
function MetricTensGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MetricTensGUI (see VARARGIN)

% Choose default command line output for MetricTensGUI
handles.output = hObject;
handles.n=1; % Counter for hkl list
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MetricTensGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = MetricTensGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;

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
data=get(hObject,'data');
% --------------------------------------------------------------------
function uitable1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 a=[1 1 1; 1 1 1];
 set(hObject,'data')=a;
   
    function Reset_Callback(handles)
    set(handles.edit2,'Visible','Off')
    set(handles.edit3,'Visible','Off')
    set(handles.edit4,'Visible','Off')
    set(handles.edit5,'Visible','Off')
    set(handles.edit6,'Visible','Off')
    set(handles.text5,'Visible','Off')
    set(handles.text6,'Visible','Off')
    set(handles.text7,'Visible','Off')
    set(handles.text8,'Visible','Off')
    set(handles.text9,'Visible','Off')

%% Metric Tensor    
function MetricTens_Callback(hObject, handles, a, b, c, alpha, beta, gam)
%Directions
handles.d1=[str2double(get(handles.edit13,'String')) str2double(get(handles.edit14,'String')) str2double(get(handles.edit15,'String'))];
handles.d2=[str2double(get(handles.edit16,'String')) str2double(get(handles.edit17,'String')) str2double(get(handles.edit18,'String'))];
%Planes
handles.p1=[str2double(get(handles.edit7,'String')) str2double(get(handles.edit8,'String')) str2double(get(handles.edit9,'String'))];
handles.p2=[str2double(get(handles.edit10,'String')) str2double(get(handles.edit11,'String')) str2double(get(handles.edit12,'String'))];
%Computation
g=[a^2 a*b*cosd(gam) a*c*cosd(beta); b*a*cosd(gam) b^2 b*c*cosd(alpha); c*a*cosd(beta) c*b*cosd(alpha) c^2];
handles.g=g; % real space
G=inv(g);
handles.G=G;
handles.Dir_norm_p1=handles.p1*handles.G*10; % direction normal to plane 1 (multiplicative factor added)
handles.Dir_norm_p2=handles.p2*handles.G*10; % direction normal to plane 2
handles.Plane_norm_d1=handles.d1*handles.g/10; % Plane normal to direction 1 (dividing factor added)
handles.Plane_norm_d2=handles.d2*handles.g/10; % Plane normal to direction 2
handles.Ang_Directions=acosd(((handles.d1*handles.g*handles.d2')/(sqrt(handles.d1*handles.g*handles.d1')*sqrt(handles.d2*handles.g*handles.d2')))); % real space angle between directions
handles.Ang_Planes=acosd(((handles.p1*handles.G*handles.p2')/(sqrt(handles.p1*handles.G*handles.p1')*sqrt(handles.p2*handles.G*handles.p2')))); % Angle between planes
handles.dspace1=(sqrt(handles.p1*handles.G*handles.p1'))^-1; % d-space of plane 1
handles.dspace2=(sqrt(handles.p2*handles.G*handles.p2'))^-1; % d-space of plane 2
V_Cell=sqrt(a^2*b^2*c^2*(1-cosd(alpha)^2-cosd(beta)^2-cosd(gam)^2+2*cosd(alpha)*cosd(beta)*cosd(gam)));
% Assigns d-space
set(handles.text2,'String',handles.dspace1)
set(handles.text3,'String',handles.dspace2)
% Assigns Planes normal to direction
set(handles.text20,'String',num2str(handles.Plane_norm_d1)) % round function implemented here in the past to avoid uneccessary decimals
set(handles.text21,'String',num2str(handles.Plane_norm_d2)) % round (see above)
% Assigns Direction normals to plane
set(handles.text23,'String',num2str(handles.Dir_norm_p1)) % round (see above)
set(handles.text24,'String',num2str(handles.Dir_norm_p2)) % round (see above)
% Assigns Angles between Planes and Directions
set(handles.text16,'String',num2str(handles.Ang_Planes))
set(handles.text17,'String',num2str(handles.Ang_Directions))
% Assigns Unit Cell Volume
set(handles.text31,'String',strcat('Unit Cell Volume:',num2str(V_Cell)))

% Bulk  of computation for Stereographic Projections
Mat_dir=[sqrt(g(1,1)) g(2,1)/sqrt(g(1,1)) g(3,1)/sqrt(g(1,1)); 0 V_Cell*sqrt(G(3,3)/g(1,1)) -V_Cell*G(3,2)/sqrt(G(3,3)*g(1,1));0 0 1/sqrt(G(3,3))]; % this works then above
Mat_pole=Mat_dir*handles.G; % Mulitply this by direction to get plane normal i.e. pole
handles.Mat_pole=Mat_pole;
handles.Mat_dir=Mat_dir;

if and(get(handles.checkbox1,'Value')==0,get(handles.checkbox2,'Value')==0)
    % Using GUI for just Crystallographic Computations
    set(handles.axes2,'Visible','off')
    cla
else
        set(handles.axes2,'Visible','on')

%Generates List of HKLs to plot
handles.hkld(1+(handles.n-1)*2,1:3)=handles.d1;
handles.hkld(2+(handles.n-1)*2,1:3)=handles.d2;
handles.hklp(1+(handles.n-1)*2,1:3)=handles.p1;
handles.hklp(2+(handles.n-1)*2,1:3)=handles.p2;
handles.n=handles.n+1;

 guidata(hObject, handles); % this was originally the last line
StereographicProj(hObject,handles, handles.hklp, handles.hkld)

end

function StereographicProj(hObject, handles, hklp, hkld)
    % Computing Stereographic Projections       
hklp=unique(hklp,'rows'); % Filters out only unique HKLs
hkld=unique(hkld,'rows'); % Filters out only unique HKLs

try
   hklp= cat(1,hklp,handles.hklpL);
   hkld= cat(1,hkld,handles.hkldL);
catch
end

danx=str2double(get(handles.edit19,'String'));
dany=str2double(get(handles.edit20,'String'));
danz=str2double(get(handles.edit21,'String'));
def=[danx dany danz];
if get(handles.radiobutton1, 'Value')==1
set1=handles.Mat_dir*def'; % For Setting View by directions
set1=set1/norm(set1);
else
set1=handles.Mat_pole*def'; % For Setting View by pole
set1=set1/norm(set1);
end
[azi,ele,~]=cart2sph(set1(1),set1(2),set1(3));
sanz=azi*180/pi-90;
sanx=ele*180/pi-90; sany=0;
% First Rotations to Set View [Order dependent Z-Y-X]
transx=[1 0 0; 0 cosd(sanx) -sind(sanx); 0 sind(sanx) cosd(sanx)];
transy=[cosd(sany) 0 -sind(sany); 0 1 0; sind(sany) 0 cosd(sany)];
transz=[cosd(sanz) -sind(sanz) 0; sind(sanz) cosd(sanz) 0; 0 0 1];
% Second Rotations for Additional rotations [Order dependent X-Y-Z]
sanx2=180+str2double(get(handles.edit25,'String'));
sany2=str2double(get(handles.edit26,'String'));
sanz2=-90+str2double(get(handles.edit27,'String')); % Sets 111 reflection in IV quadarant so that 010
transx2=[1 0 0; 0 cosd(sanx2) -sind(sanx2); 0 sind(sanx2) cosd(sanx2)];
transy2=[cosd(sany2) 0 -sind(sany2); 0 1 0; sind(sany2) 0 cosd(sany2)];
transz2=[cosd(sanz2) -sind(sanz2) 0; sind(sanz2) cosd(sanz2) 0; 0 0 1];

%% Compute Cartesian and Stereographic Projection Coordinates 
    
    for jj=1:1:size(hklp,1) % Computes cartesian and stereographic for poles
    % 3D plot- Cartesian Coordinates
poles_cart(jj,1:3)=handles.Mat_pole*hklp(jj,1:3)';
poles_cart(jj,1:3)=poles_cart(jj,1:3)/norm(poles_cart(jj,1:3));
poles_cart(jj,1:3)=poles_cart(jj,1:3)*transz*transy*transx*transx2*transy2*transz2; % Applies rotations for setting view and additionals
    % Wulff plot coordinates
poles_stereo_x(jj)=poles_cart(jj,1)/(1+abs(poles_cart(jj,3))); % abs forces the values of Z to the top of the hemisphere
poles_stereo_y(jj)=poles_cart(jj,2)/(1+abs(poles_cart(jj,3))); % Reasoning for this is that diffraction is equal in top and bottom of the sphere
    end
    for jj=1:1:size(hkld,1) % Computes cartesian and stereographic for directions
  % 3D plot-= Cartesian Coordinates
dir_cart(jj,1:3)=handles.Mat_dir*hkld(jj,1:3)';
dir_cart(jj,1:3)=dir_cart(jj,1:3)/norm(dir_cart(jj,1:3));
dir_cart(jj,1:3)=dir_cart(jj,1:3)*transz*transy*transx*transx2*transy2*transz2; % Applies rotations for setting view and additionals
    % Wulff plot coordinates
dir_stereox(jj)=dir_cart(jj,1)/(1+abs(dir_cart(jj,3)));
dir_stereoy(jj)=dir_cart(jj,2)/(1+abs(dir_cart(jj,3)));
    end
    
%3D Diffract Points to Plot
handles.poles_cart=poles_cart;
handles.dir_cart=dir_cart;
    
% Stereographic/Wulff Plot
handles.hklp=hklp; %stores hklp into guidata
handles.hkld=hkld; %stores hkld into guidata
handles.poles_stereo_x=poles_stereo_x;
handles.poles_stereo_y=poles_stereo_y;
handles.dir_stereox=dir_stereox;
handles.dir_stereoy=dir_stereoy;
    
guidata(hObject, handles);
cla(handles.axes2,'reset') % resets to re-plot

if get(handles.popupmenu2,'Value')==1
    Sphere(hObject, handles, hklp, poles_cart, hkld, dir_cart)   
else
    
% Begins Plottin of Stereographic Projections
 Wulff(hObject, handles, str2double(get(handles.edit28,'String')))
 hold on
% Handles plotting and checkbox options
if and(get(handles.checkbox1,'Value')==1, get(handles.checkbox2,'Value')==1)
 
 for u=1:size(hklp,1)
scatter(poles_stereo_x(u),poles_stereo_y(u),'MarkerFaceColor','black')
text(poles_stereo_x(u),poles_stereo_y(u),strcat(num2str(hklp(u,1:3))),'horizontal','left','vertical','bottom'); % in GUI you will need to specify a direction vs. a pole
 end
 
 for u=1:size(hkld,1)
scatter(dir_stereox(u),dir_stereoy(u),'square','MarkerFaceColor',[.8 .1 .1])
text(dir_stereox(u),dir_stereoy(u),strcat(num2str(hkld(u,1:3))),'horizontal','left','vertical','bottom'); % in GUI you will need to specify a direction vs. a pole
 end    
        
elseif get(handles.checkbox1,'Value')==1
 for u=1:size(hklp,1)
scatter(poles_stereo_x(u),poles_stereo_y(u),'MarkerFaceColor','black')
text(poles_stereo_x(u),poles_stereo_y(u),strcat(num2str(hklp(u,1:3))),'horizontal','left','vertical','bottom'); % in GUI you will need to specify a direction vs. a pole
 end

elseif get(handles.checkbox2,'Value')==1
 for u=1:size(hkld,1)
scatter(dir_stereox(u),dir_stereoy(u),'square','MarkerFaceColor',[.8 .1 .1])
text(dir_stereox(u),dir_stereoy(u),strcat(num2str(hkld(u,1:3))),'horizontal','left','vertical','bottom'); % in GUI you will need to specify a direction vs. a pole
 end
    
end
end
guidata(hObject, handles);
%% Create Sphere for 3D Diffract Option
function Sphere(hObject, handles, hklp, poles_cart, hkld, dir_cart)
guidata(hObject, handles)
[X,Y,Z] = sphere(50);
[A, B]=find(X>=0 & Y>=0 & Z>=0);
axes(handles.axes2)
mesh(X(A,B),Y(A,B),Z(A,B),'FaceColor','none')
colormap bone
xlabel('X')
ylabel('Y')
zlabel('Z')
set(handles.axes2,'FontName','Palatino Linotype')

hold on
% Handles plotting and checkbox options
if and(get(handles.checkbox1,'Value')==1,get(handles.checkbox2,'Value')==1)
 
 for u=1:size(hklp,1)
scatter3(poles_cart(u,1),poles_cart(u,2),abs(poles_cart(u,3))','MarkerFaceColor','black')
text(poles_cart(u,1),poles_cart(u,2),abs(poles_cart(u,3)),strcat(num2str(hklp(u,1:3))),'horizontal','left','vertical','bottom'); % abs to force values to top of hemisphere
 end
  
 for u=1:size(hkld,1)
 scatter3(dir_cart(u,1),dir_cart(u,2),abs(dir_cart(u,3))','square','MarkerFaceColor',[.8 .1 .1])
text(dir_cart(u,1),dir_cart(u,2),abs(dir_cart(u,3)),strcat(num2str(hkld(u,1:3))),'horizontal','left','vertical','bottom'); % abs to force values to top of hemisphere
 end    
        
elseif get(handles.checkbox1,'Value')==1
 for u=1:size(hklp,1)
scatter3(poles_cart(u,1),poles_cart(u,2),abs(poles_cart(u,3))','MarkerFaceColor','black')
text(poles_cart(u,1),poles_cart(u,2),abs(poles_cart(u,3)),strcat(num2str(hklp(u,1:3))),'horizontal','left','vertical','bottom'); % abs to force values to top of hemisphere
 end
elseif get(handles.checkbox2,'Value')==1
 for u=1:size(hkld,1)
scatter3(dir_cart(u,1),dir_cart(u,2),abs(dir_cart(u,3))','square','MarkerFaceColor',[.8 .1 .1])
text(dir_cart(u,1),dir_cart(u,2),abs(dir_cart(u,3)),strcat(num2str(hkld(u,1:3))),'horizontal','left','vertical','bottom'); % abs to force values to top of hemisphere
 end
    
end
%% Create Wulff Net Plot
function Wulff(hObject, handles, Ang)
% Creates a Wulf Plot and allows rotation of the Wulff plot by specified angle Ang
% Created by Giovanni Esteves, Department of Materials Science and Engineering, North Carolina State University, 11-30-2015
% This section of code was inspired by the Wulff.m script created by Gerald Middleton 
% Creates Outer Circle that corresponds to 90°
if nargin==0
    handles.Ang=0;
else
end
outter=0:.01:2*pi;
xc=cos(outter);
yc=sin(outter);
plot(xc,yc,'Color','black')
axis('square')
hold on
% Create X-Y axis centered at [0,0]
xa=[-1 1 0; 0 0 0];
ya=[0 0 0; -1 1 0];
ang=Ang; % Specify angle of rotation
trans=[cosd(ang) sind(ang) 0; -sind(ang) cosd(ang) 0; 0 0 1];
transx=trans*[1 0 0]';
transy=trans*[0 1 0]';
nxa=cat(2,transx,-transx);
nya=cat(2,-transy,transy);
if ang==0;
plot(xa(1,1:2),ya(1,1:2),xa(2,1:2),ya(2,1:2),'Color','red')
end
plot(nxa(1,1:2),nya(1,1:2),nxa(2,1:2),nya(2,1:2),'Color','black')

% Create Great Circles
phi=0:0.0001:pi;
for i=1:8;  % speifies the amount of circles between 0-90 (default is 10°)
Beta=atan(tan(i*pi/18)*sin(phi)); % Project half a circle that is tilited by angled i*pi/18.
xgc1=sin(phi).*tan(pi/4-Beta/2); % Equation comes from Priest 1985 for equal angle projections
xgc2=-sin(phi).*tan(pi/4-Beta/2); % Mentioned online as well by various sources
ygc=cos(phi).*tan(pi/4-Beta/2);
comb(1,:)=xgc1';
comb(2,:)=ygc';
comb(3,:)=zeros(1,length(xgc1));
comb2(1,:)=xgc2';
comb2(2,:)=comb(2,:);
comb2(3,:)=comb(3,:);
nrotgc=trans*comb;
nrotgc2=trans*comb2;
plot(nrotgc(1,:),nrotgc(2,:),':black',nrotgc2(1,:),nrotgc2(2,:),':black')
end
% Create Small Circles
for o=1:8
alpha=o*pi/18; % speifies the amount of circles between 0-90 (default is 10°)
rcone=sin(alpha); % radius of cone based on input angle
xres=-rcone:.001:rcone; % Creates the resolution of the base of the cone
openangle=tan(alpha); % find the r/h ratio based on opening angle of a cone
d=1/cos(alpha); % finds the the slant to the cone which is what positions in in the correct space in Z
sc1=sqrt(openangle^2-xres.^2); % creates the base of a cone based on open angle, radius, and slant
sctop=d-sc1; % used to create a mirror to plot the top small circles
scbot=-d+sc1; % used to plot the bottom small circles
combsc=[]; % resets due to changing size of xres
combsc(1,:)=xres;
combsc(2,:)=scbot;
combsc(3,:)=zeros(1,length(scbot));
nrotsc=trans*combsc;
plot(nrotsc(1,:),nrotsc(2,:),':black',-nrotsc(1,:),-nrotsc(2,:),':black')
end
axis off; % removes the axis for better visuallization
axis('square')
shg
guidata(hObject, handles)
%% Crystal Selection [Visible Options and Conditions]
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
content=get(hObject,'Value');
    switch content    
        case 1 % Cubic
Reset_Callback(handles) 

        case 2 %Tetragonal
        %Reset
Reset_Callback(handles)
        % Initial Settings
            set(handles.edit3,'Visible','On')
            set(handles.text6,'Visible','On')
            
        case 3 % Hexagonal
        % Reset    
Reset_Callback(handles)       
        % Initial Settings
            set(handles.edit3,'Visible','On')
            set(handles.text6,'Visible','On')
            set(handles.edit6,'Visible','On')
            set(handles.text9,'Visible','On')
            set(handles.edit6,'String','120','Enable','on')
            
        case 4 % Rhombohedral
        % Reset    
        Reset_Callback(handles)
        % Intial Parameters     
    set(handles.text7,'Visible','On')
    set(handles.text8,'Visible','On')
    set(handles.text9,'Visible','On')
    set(handles.edit4,'Visible','On')
    set(handles.edit5,'Visible','On','Enable','inactive')
    set(handles.edit6,'Visible','On','Enable','inactive')
    set(handles.edit4,'Visible','On','Enable','on')
    %Constraints
%alpha=beta=gamma
        set(handles.edit4,'String','90');
        set(handles.edit5,'String','90');
        set(handles.edit6,'String','90');
        try
        handles.beta=handles.alpha;
        handles.gamma=handles.alpha;
        catch
        end
        % Catches inconsitency when switching from hex to rhomb
            if or(strcmp(get(handles.edit5,'String'),get(handles.edit6,'String'))~=1,strcmp(get(handles.edit4,'String'),get(handles.edit6,'String'))~=1)
        set(handles.edit4,'String','90');
        set(handles.edit5,'String','90');
        set(handles.edit6,'String','90');
            end
            
        case 5 % Ortho
                    % Reset
                    Reset_Callback(handles)
                    
                    % Initial Parameters
    set(handles.edit2,'Visible','On')
    set(handles.edit3,'Visible','On')
    set(handles.text5,'Visible','On')
    set(handles.text6,'Visible','On')
                                     
        case 6 % Monoclinic
                        % Reset
                       Reset_Callback(handles)
                       
                       % Initial parameters
    set(handles.edit2,'Visible','On')
    set(handles.edit3,'Visible','On')
    set(handles.text5,'Visible','On')
    set(handles.text6,'Visible','On')                 
    %Angles
    set(handles.text7,'Visible','On')
    set(handles.text8,'Visible','On')
    set(handles.text9,'Visible','On')
    set(handles.edit4,'Visible','On','Enable','inactive')
    set(handles.edit5,'Visible','On')
    set(handles.edit5,'Visible','On','Enable','on')
    set(handles.edit6,'Visible','On','Enable','inactive')                       
    % Constraint 
        set(handles.edit5,'Visible','On') %beta is not equal to 90 and can change
        set(handles.edit4,'String','90');
        set(handles.edit5,'String','90');
        set(handles.edit6,'String','90');
                       
        case 7 % Triclinic
                            % Reset
                            Reset_Callback(handles)
                            
    % Initial Parameters                          
    set(handles.edit2,'Visible','On')
    set(handles.edit3,'Visible','On')
    set(handles.text5,'Visible','On')
    set(handles.text6,'Visible','On')                 
    %Angles
    set(handles.text7,'Visible','On')
    set(handles.text8,'Visible','On')
    set(handles.text9,'Visible','On')
    set(handles.edit4,'Visible','On')
    set(handles.edit4,'Visible','On','Enable','on')
    set(handles.edit5,'Visible','On')
    set(handles.edit5,'Visible','On','Enable','on')
    set(handles.edit6,'Visible','On')
    set(handles.edit6,'Visible','On','Enable','on')          
        otherwise
                                    
    end

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

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
handles.a=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
edit1_Callback(hObject,eventdata,handles)
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
handles.b=str2double(get(hObject,'String'));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
edit2_Callback(hObject,eventdata,handles)

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
handles.c=str2double(get(hObject,'String'));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
edit3_Callback(hObject,eventdata,handles)

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
try
handles.alpha=str2double(get(hObject,'String'));
if get(handles.popupmenu1,'Value')==4
set(handles.edit5,'String',handles.alpha)
set(handles.edit6,'String',handles.alpha)
end
catch
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
edit4_Callback(hObject,eventdata,handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles.beta=str2double(get(hObject,'String'));
if get(handles.popupmenu1,'Value')==6
end
catch
end

guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
edit5_Callback(hObject,eventdata,handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.gamma=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
edit6_Callback(hObject,eventdata,handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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

function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double

% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double

% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double

% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Compute Button [Crystal Constraints]
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles) % handles the Compute button of GUI
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.popupmenu1,'Value')==1 % Cubic Constraints
    
    handles.a=str2double(get(handles.edit1,'String'));
    handles.b=handles.a;
    handles.c=handles.a;
    handles.alpha=90;
    handles.beta=90;
    handles.gamma=90;
    
elseif get(handles.popupmenu1,'Value')==2 % Tetra Constraints
    
    handles.a=str2double(get(handles.edit1,'String'));
    handles.b=handles.a;
    handles.c=str2double(get(handles.edit3,'String'));
    handles.alpha=90;
    handles.beta=90;
    handles.gamma=90;
    
elseif get(handles.popupmenu1,'Value')==3 % Hexagonal Constraint
    handles.a=str2double(get(handles.edit1,'String'));

    handles.b=handles.a;
    handles.c=str2double(get(handles.edit3,'String'));

    handles.alpha=90;
    handles.beta=90;
    handles.gamma=str2double(get(handles.edit6,'String'));
        
elseif get(handles.popupmenu1,'Value')==4 % Rhombohedral Constraints
    
    handles.a=str2double(get(handles.edit1,'String'));
    handles.b=handles.a;
    handles.c=handles.a;
    handles.beta=handles.alpha;
    handles.gamma=handles.alpha;
        
elseif get(handles.popupmenu1,'Value')==5 % Orthorhombic Constraints
    handles.a=str2double(get(handles.edit1,'String'));
    handles.b=str2double(get(handles.edit2,'String'));
    handles.c=str2double(get(handles.edit3,'String'));
    handles.alpha=90;
    handles.beta=90;
    handles.gamma=90;
        
elseif get(handles.popupmenu1,'Value')==6 % Monoclinic Constraint
    handles.a=str2double(get(handles.edit1,'String'));
    handles.b=str2double(get(handles.edit2,'String'));
    handles.c=str2double(get(handles.edit3,'String'));
    handles.beta=str2double(get(handles.edit5,'String'));
    handles.alpha=90;
    handles.gamma=90;
        
elseif get(handles.popupmenu1,'Value')==7 % Triclinic Constraint
    handles.a=str2double(get(handles.edit1,'String'));
    handles.b=str2double(get(handles.edit2,'String'));
    handles.c=str2double(get(handles.edit3,'String'));
    handles.alpha=str2double(get(handles.edit4,'String'));
    handles.beta=str2double(get(handles.edit5,'String'));
    handles.gam=str2double(get(handles.edit6,'String'));
    
end



MetricTens_Callback(hObject, handles,handles.a,handles.b,handles.c, handles.alpha, handles.beta, handles.gamma)

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles) % deals with the drop down menu to select plot types
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
guidata(hObject,handles)
content=get(hObject,'Value');

    switch content  
    case 1
set(handles.text26,'String','3D Cartesian Plot')
cla(handles.axes2,'reset') % resets to re-plot

    case 2
set(handles.text26,'String','Wulff Net')
         if and(get(handles.checkbox1,'Value')==0,get(handles.checkbox2,'Value')==0) % stylastic change (turns off plot when neither "plot poles" or "plot directions" is checked
            set(handles.axes2,'Visible','off')
         else
try
Wulff(hObject, handles, handles.Ang)

catch
    Wulff(hObject, handles, 0)
end
         end
    end % end switch

if and(get(handles.checkbox1,'Value')==0,get(handles.checkbox2,'Value')==0) % stylastic change (turns off plot when neither "plot poles" or "plot directions" is checked
 set(handles.axes2,'Visible','off')
else
    
cla(handles.axes2,'reset') % resets to re-plot

if and(content==1, isempty(handles.hklp)~=1)
    Sphere(hObject, handles, handles.hklp, handles.poles_cart, handles.hkld, handles.dir_cart)   
elseif isempty(handles.hklp)~=1    
% Begins Plottin of Stereographic Projections
 Wulff(hObject, handles, 0)
 StereographicProj(hObject, handles, handles.hklp, handles.hkld)
end

end
guidata(hObject, handles)

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
handles.hklp=[];
handles.hkld=[];
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function text26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.hklp=[];
handles.hkld=[];
handles.hklpL=[];
handles.hkldL=[];
handles.n=1;
clear global handles.g
cla(handles.axes2,'reset') % resets to re-plot

if get(handles.popupmenu2,'Value')==2
 Wulff(hObject, handles, 0)
elseif get(handles.popupmenu2,'Value')==1
[X,Y,Z] = sphere(50);
[A, B]=find(X>=0 & Y>=0 & Z>=0);
mesh(X(A,B),Y(A,B),Z(A,B));
mesh(X(A,B),Y(A,B),Z(A,B),'FaceColor','none')
colormap bone
xlabel('X')
ylabel('Y')
zlabel('Z')
end
guidata(hObject, handles)

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.chkbox1=get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of checkbox1

% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.chkbox2=get(hObject,'Value');
guidata(hObject,handles)

function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function text11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if and(get(handles.radiobutton2,'Value')==0,get(handles.radiobutton1,'Value')==0)
set(handles.radiobutton1,'Value',1)
elseif get(handles.radiobutton1,'Value')==1
    set(handles.radiobutton2,'Value',0)
    
end
% Hint: get(hObject,'Value') returns toggle state of radiobutton1

% --- Executes during object creation, after setting all properties.
function radiobutton1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if and(get(handles.radiobutton1,'Value')==0,get(handles.radiobutton2,'Value')==0)
    set(handles.radiobutton2,'Value',1)
elseif get(handles.radiobutton2,'Value')==1
    set(handles.radiobutton1,'Value',0)
end
% Hint: get(hObject,'Value') returns toggle state of radiobutton2

function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit28_Callback(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function edit28_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
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

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles) % Handles the loading of data
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, datapath]=uigetfile('.xlsx','MultiSelect','off'); %will ask user to input a new file 
inFile = strcat(datapath, {filename});
HKL=xlsread(inFile{1});
hklp=HKL(:,1:3);
hkld=HKL(:,4:6);
handles.hklpL=hklp;
handles.hkldL=hkld;

guidata(hObject, handles)
% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles) % About Section, this needs to be updated everytime there is a new release of the GUI
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --------------------------------------------------------------------
function About_Callback(hObject, eventdata, handles)
% hObject    handle to About (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=msgbox({'MetricTens, version: 1.0' 'Created by Giovanni Esteves' 'North Carolina State University' 'Email: gesteve@ncsu.edu' 'Additional Contact: Jacob Jones' 'Email: jacobjones@ncsu.edu'},'About');
set(h, 'Position',[500 440 200 125])
ah=get(h,'CurrentAxes');
c=get(ah,'Children');
set(c,'FontSize',12,'FontName','Palatino Linotype');


% --------------------------------------------------------------------
function Help_Callback(hObject, eventdata, handles)
% hObject    handle to Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=msgbox('Full documentation coming soon','Help');
set(h, 'Position',[500 440 200 50])
ah=get(h,'CurrentAxes');
c=get(ah,'Children');
set(c,'FontSize',12,'FontName','Palatino Linotype');


% --------------------------------------------------------------------
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uiputfile('*.tiff;*.jpg;*.bmp;*.png','Save Window As');

if FileName==0 % this just prevents an error if user hits 'cancel' on file save
else

fileID=strcat(PathName,FileName);
img = getframe(gcf);
imwrite(img.cdata, fileID);
end
