function varargout = Three_Link_HZD(varargin)
% THREE_LINK_HZD MATLAB code for Three_Link_HZD.fig
%      THREE_LINK_HZD, by itself, creates a new THREE_LINK_HZD or raises the existing
%      singleton*.
%
%      H = THREE_LINK_HZD returns the handle to a new THREE_LINK_HZD or the handle to
%      the existing singleton*.
%
%      THREE_LINK_HZD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in THREE_LINK_HZD.M with the given input arguments.
%
%      THREE_LINK_HZD('Property','Value',...) creates a new THREE_LINK_HZD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Three_Link_HZD_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Three_Link_HZD_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Three_Link_HZD

% Last Modified by GUIDE v2.5 12-Mar-2013 18:53:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Three_Link_HZD_OpeningFcn, ...
                   'gui_OutputFcn',  @Three_Link_HZD_OutputFcn, ...
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


% --- Executes just before Three_Link_HZD is made visible.
function Three_Link_HZD_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Three_Link_HZD (see VARARGIN)

% Choose default command line output for Three_Link_HZD
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global sim_tag
sim_tag = 0;
% UIWAIT makes Three_Link_HZD wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Three_Link_HZD_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function modelParams_Callback(hObject, eventdata, handles)
% hObject    handle to modelParams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of modelParams as text
%        str2double(get(hObject,'String')) returns contents of modelParams as a double


% --- Executes during object creation, after setting all properties.
function modelParams_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modelParams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RunButton.
function RunButton_Callback(hObject, eventdata, handles)
% hObject    handle to RunButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global t_2 torque y force

torque = [];
t_2 = [];
y = [];
force = [];

tstart = 0;
tfinal = 13;

%% The optimization parameters
%
a1 = str2num(get(handles.a1,'string'));
a2 = str2num(get(handles.a2,'string'));
a  = [a1 a2];
params = str2num(get(handles.modelParams,'string'));
%= [0.512 0.073 0.035 -0.819 -2.27 3.26 3.11 1.89];

omega_1 = 1.55;
x0 = sigma_three_link(omega_1,a);
x0 = transition_three_link(x0).';
x0 = x0(1:6);

options = odeset('Events','on','Refine',4,'RelTol',10^-5,'AbsTol',10^-6);

tout = tstart;
xout = x0.';
teout = []; xeout = []; ieout = [];

disp('(impact ratio is the ratio of tangential to normal');
disp('forces of the tip of the swing leg at impact)');

for i = 1:5 % run five steps
	% Solve until the first terminal event.
	[t,x,te,xe,ie] = ode45('walker_main',[tstart tfinal],x0,options,a);

	% Accumulate output.  tout and xout are passed out as output arguments
	nt = length(t);
	tout = [tout; t(2:nt)];
	xout = [xout;x(2:nt,:)];
	teout = [teout; te]; % Events at tstart are never reported.
	xeout = [xeout; xe];
	ieout = [ieout; ie];

	% Set the new initial conditions (after impact).
	x0=transition_three_link(x(nt,:));

	% display some useful information to the user
	disp(['step: ',num2str(i),', impact ratio:  ',num2str(x0(7)/x0(8))])

	% Only positions and velocities needed as initial conditions
	x0=x0(1:6);

	tstart = t(nt);
	if tstart>=tfinal
		break
	end
end

disp('by Eric R. Westervelt, Jessy W. Grizzle,');
disp('Christine Chevallereau, Jun-Ho Choi, and Benjamin Morris');

% Run the animation
axes(handles.AnimationAxes);
hold off
% ax = gcf; 

anim(tout,xout,1/30,1,params);
%% Draw some useful graphs
%
fig_hl=figure(2);
set(fig_hl,'Position', [200 100 400 450]);
set(fig_hl,'PaperPosition',[1.25 1.5 6 8])
subplot(2,1,1)
plot(tout,xout(:,1),'-',tout,xout(:,2),'--',tout,xout(:,3),'-.')
legend('\theta_1','\theta_2','\theta_3')
grid
title('Joint Positions')
subplot(2,1,2)
plot(tout,xout(:,4),'-',tout,xout(:,5),'--',tout,xout(:,6),'-.')
legend('\theta_1 (dot)','\theta_2 (dot)','\theta_3 (dot)');
xlabel('time (sec)')
grid
title('Joint Velocities')

fig_hl=figure(3);
set(fig_hl,'Position', [220 120 400 450])
set(fig_hl,'PaperPosition',[1.25 1.5 6 8])
subplot(2,1,1)
plot(t_2,force(:,1),'-b',t_2,force(:,2),'--r')
legend('F_{tan} (N)','F_{norm} (N)')
grid
title('Forces on End of Stance Leg')
subplot(2,1,2)
plot(t_2,force(:,1)./force(:,2))
ylabel('F_{tan}/F_{norm}');
xlabel('time (sec)')
grid

fig_hl=figure(4);
set(fig_hl,'Position', [240 140 400 450]);
set(fig_hl,'PaperPosition',[1.25 1.5 6 8])
subplot(2,1,1)
plot(t_2,torque(:,1))
ylabel('\tau_1 (Nm)')
grid
title('Control Signals')
subplot(2,1,2)
plot(t_2,torque(:,2))
ylabel('\tau_2 (Nm)');
xlabel('time (sec)')
grid

fig_hl=figure(5);
set(fig_hl,'Position', [260 160 400 450]);
set(fig_hl,'PaperPosition',[1.25 1.5 6 8])
subplot(2,1,1)
plot(t_2,y(:,1))
ylabel('y_1')
title('Outputs')
grid
subplot(2,1,2)
plot(t_2,y(:,2))
ylabel('y_2');
xlabel('time (sec)')
grid




%% --------------------------------------------------------------------------
%% the animation function

function anim(t,x,ts,speed,params)

[n,m]=size(x);
[vV,vH]=hip_vel(x); % convert angles to horizontal position of hips
pH_horiz = zeros(n,1);

% Estimate hip horizontal position by estimating integral of hip velocity
for j=2:n
	pH_horiz(j)=pH_horiz(j-1)+(t(j)-t(j-1))*vH(j-1,1);
end

[te,pH_horiz]=even_sample(t,pH_horiz,1/ts);
[te,xe]=even_sample(t,x,1/ts);
[n,m]=size(xe);

k=0;

q=xe(1,1:3);
[pFoot1,pFoot2,pH,pT]=limb_position(q,pH_horiz(1));
% fig_hl = ax;
% fig_hl=figure(1);
% set(fig_hl,'Position', [280 180 400 350]);
% set(fig_hl,'PaperPosition',[0 0 6 5]);
% clf
cla
anim_axis=[-2.2 2.2 -2.2 2.2];
axis off
axis(anim_axis)
grid

% Use actual relations between masses in animation
r =params(1);
m =params(2);
Mh =params(3);
Mt =params(4);
L =params(5);
g =params(6);
%[r,m,Mh,Mt,L,g] = model_params_three_link;
scl=0.04; % factor to scale masses
mr_legs=m^(1/3)*scl; % radius of mass for legs
mr_torso=Mt^(1/3)*scl; % radius of mass for torso
leg1_color='g';
leg2_color='r';
torso_color='b';
ground_color='k'; % a.k.a. black

% Approximate circular shape of mass
param=linspace(0,2*pi+2*pi/50,50);
xmass_legs=mr_legs*cos(param);
ymass_legs=mr_legs*sin(param);

xmass_torso=mr_torso*cos(param);
ymass_torso=mr_torso*sin(param);

% Draw ground
buffer=5;
ground=line([-buffer pH_horiz(n)+buffer],[0 0]);
set(ground,'Color',ground_color,'LineWidth',2);
for k=-buffer:floor(pH_horiz(n)+buffer)
	ref_tick(k+buffer+1)=line([k k],[-0.1 0]);
	set(ref_tick(k+buffer+1),'Color',ground_color);
	ref_label(k+buffer+1)=text(-0.03+k,-0.2,num2str(k));
end

% Draw leg one
leg1=line([pFoot1(1) pH(1)],[pFoot1(2) pH(2)]);
mass1=patch(xmass_legs+(pH(1)-pFoot1(1))/2,...
	ymass_legs+(pH(2)-pFoot1(2))/2,leg1_color);
set(mass1,'EdgeColor',leg1_color)
set(leg1,'LineWidth',2,'Color',leg1_color);

% Draw leg two
leg2=line([pFoot2(1) pH(1)],[pFoot2(2) pH(2)]);
mass2=patch(xmass_legs+pH(1)-(pH(1)-pFoot2(1))/2,...
	ymass_legs+pH(2)-(pH(2)-pFoot2(2))/2,leg2_color);
set(mass2,'EdgeColor',leg2_color)
set(leg2,'LineWidth',2,'Color',leg2_color);

% Draw torso
torso=line([pH(1) pT(1)],[pH(2)*2 pT(2)*2]);
torso_mass=patch(xmass_torso+pT(1),ymass_torso+pT(2),torso_color);
set(torso_mass,'EdgeColor',torso_color)
set(torso,'LineWidth',2,'Color',torso_color);

for k=2:n
	q=xe(k,1:3);
	[pFoot1,pFoot2,pH,pT]=limb_position(q,pH_horiz(k));

	set(leg1,'XData',[pFoot1(1) pH(1)],'YData',[pFoot1(2) pH(2)]);

	set(mass1,'XData',xmass_legs+(pH(1)-pFoot1(1))/2+pH_horiz(k),...
		'YData',ymass_legs+(pH(2)-pFoot1(2))/2);

	set(leg2,'XData',[pFoot2(1) pH(1)],'YData',[pFoot2(2) pH(2)]);

	set(mass2,'XData',xmass_legs+pH(1)-(pH(1)-pFoot2(1))/2,...
		'YData',ymass_legs+pH(2)-(pH(2)-pFoot2(2))/2);

	set(torso,'XData',[pH(1) pT(1)+(pT(1)-pH(1))],...
		'YData',[pH(2) pT(2)+(pT(2)-pH(2))]);

	set(torso_mass,'XData',xmass_torso+pT(1),'YData',ymass_torso+pT(2));

	title(['T_{est} = ',num2str(te(k),'%.1f')])

	new_axis=anim_axis+[pH(1) pH(1) 0 0];
	axis(new_axis);

	for j=1:length(ref_label)
		if (j-buffer-1.05<new_axis(1)) || (j-buffer-1>new_axis(2))
			set(ref_label(j),'Visible','off')
			set(ref_tick(j),'Visible','off')
		else
			set(ref_label(j),'Visible','on');
			set(ref_tick(j),'Visible','on')
		end
	end
hold off
	drawnow;
	pause(ts*speed);
end

%% --------------------------------------------------------------------------
%% a function to calculate hip velocity

function [vV,vH] = hip_vel(x)

vV=zeros(length(x),1);
vH=cos(x(:,1)).*x(:,4); % estimate of horizontal velocity of hips


%% --------------------------------------------------------------------------
%% a function to calculate the limb position

function [pFoot1,pFoot2,pH,pT] = limb_position(q,pH_horiz)

% Use position of hips as location of stance leg foot.
pFoot1=[pH_horiz; 0];
[r,m,Mh,Mt,L,g]=model_params_three_link;
pH=[pFoot1(1)+r*sin(q(1)); pFoot1(2)+r*cos(q(1))];
pFoot2=[pH(1)-r*sin(q(2)); pH(2)-r*cos(q(2))];
pT=[pH(1)+L*sin(q(3)); pH(2)+L*cos(q(3))];


%% --------------------------------------------------------------------------
%% CONVERTS A RANDOMLY SAMPLED SIGNAL SET INTO AN EVENLY SAMPLED SIGNAL SET
%% (by interpolation)
%%
%% written by Haldun Komsuoglu, 7/23/1999

function [Et, Ex] = even_sample(t, x, Fs)

% Obtain the process related parameters
N = size(x, 2);    % number of signals to be interpolated
M = size(t, 1);    % Number of samples provided
t0 = t(1,1);       % Initial time
tf = t(M,1);       % Final time
EM = (tf-t0)*Fs;   % Number of samples in the evenly sampled case with
% the specified sampling frequency
Et = linspace(t0, tf, EM)';

% Using linear interpolation (used to be cubic spline interpolation)
% and re-sample each signal to obtain the evenly sampled forms
for s = 1:N
	Ex(:,s) = interp1(t(:,1), x(:,s), Et(:,1));
end


function a1_Callback(hObject, eventdata, handles)
% hObject    handle to a1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a1 as text
%        str2double(get(hObject,'String')) returns contents of a1 as a double


% --- Executes during object creation, after setting all properties.
function a1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a2_Callback(hObject, eventdata, handles)
% hObject    handle to a2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a2 as text
%        str2double(get(hObject,'String')) returns contents of a2 as a double


% --- Executes during object creation, after setting all properties.
function a2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in GenButton.
function GenButton_Callback(hObject, eventdata, handles)
% hObject    handle to GenButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% set(hObject,'Enable','off')
% set(handles.GenButton,'Enable','off')
filename = strrep(mfilename,'gen_model','model_save');

if ~isempty(dir([filename,'.mat']))
 delete('transition_three_link.m', 'stance_force_three_link.m', ...
        'sigma_three_link.m', 'model_save_three_link.mat' , ...
        'control_three_link.m', 'dynamics_three_link.m');
end
%  set(handles.GenButton,'Enable','off')

  gen_model_three_link;
%    set(hObject,'Enable','on')
% set(handles.GennButton,'Enable','on')

