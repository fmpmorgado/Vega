
clc; clear all

global m0 g0 T A Cd rh0 H0 Re hgr_turn tf md h_orb
Alt = 1;              %[m] Alt above sea level
n_stages = 3;         %number of stages
% VEGA Rocket
%m_stage_gross = [96243, 26300,12000];% 1st, 2nd,3d

%m_prop = [87710, 23814, 10567];       % [kg] Propellant mass
Isp    = [280,287.5, 295.9] ;        % [s]  Specific impulse
d      = 2.6;           % [m]  Diameter fairing
g0     = 9.81;        % [m/s^2] Constant at its sea-level value
m0  = 136348;         % [kg] Initial mass
Cd  = 0.5 ;             % Drag coefficient,assumed to have the constant value
rh0 = 1.225;          % [kg/m^3]
H0 = 7500;            % [m] Density scale height
Re = 6378e3;          % [m] Earth's radius
hgr_turn = 200;       % [m] Rocket starts the gravity turn when h = hgr_turn
tburn = [109.9, 77.1, 119.6];        % [s] Fuell burn time
%T   = md*(Isp*g0);
Thrust   = [3015000, 1120000, 317000];    % [N] Thrust (mean)
%t0 = 0;               % Rocket launch time

h_orb = 176e3;
niu = 3.986e14;


%v_orbit = 7506;
v_orbit = (niu/(Re+h_orb))^(0.5);
v_gravity = 1000; %m/s
v_drag = 50; %m/s

%delta_v = [1831, 1540, 2209];
%delta_v = [2831, 2540, 4209];
delta_v = v_orbit - v_gravity - v_drag;


z = 0;
while z==0
    

mass = mass_model(Isp, delta_v, Thrust); %Bad results

m_prop = [mass(1,2), mass(2,2), mass(3,2)];
m_stage_gross = [mass(1,1)+mass(1,2), mass(2,1)+mass(2,2), mass(3,1)+mass(3,2)];



%trajectory  bad results

for i=1:n_stages
    t0 = 0;
   Area(1,i)   = pi*d^2/4;       % [m^2]Frontal area
m_flow(1,i) = m_prop (1,i)/tburn(1,i); % [kg/s]Propellant mass flow rate
tf = t0 + tburn(1,i);      % The time when propellant is completely burned
%and the thrust goes to zero
mf(1,i) = m0 - m_prop(1,i); %final mass of the stage (no propellant)
t_range     = [t0,tf];  % Integration interval
md = m_flow(1,i);
A=Area(1,i);
%T   = md*(Isp(1,i)*g0);
T=Thrust(1,i);

if i == 1 % Launch initial conditions:
    
gamma0 = 89.5/180*pi;       % Initial flight path angle
v0 = 0;   % Velocity (m/s)  % Earth's Rotation considered in eq of motion.
x0 = 0;   % Downrange distance [km]
h0 = Alt; % Launch site altitude [km]
vD0 = 0;  % Loss due to drag (Velocity)[m/s]
vG0 = 0;  % Loss due to gravity (Velocity)[m/s]
   state0   = [v0, gamma0, x0, h0, vD0, vG0];
%state0 = stage1; 
% Solve initial value problem for ordinary differential equations
[t,state_1] = ode45(@RocketDynEq2,t_range,state0) ;
v1     = state_1(:,1)/1000;      % Velocity [km/s]
gamma1 = state_1(:,2)*180/pi;    % Flight path angle  [deg]
x1     = state_1(:,3)/1000;      % Downrange distance [km]
h1     = state_1(:,4)/1000;      % Altitude[km]
vD1    = -state_1(:,5)/1000;     % Loss due to drag (Velocity)[m/s]
vG1    = -state_1(:,6)/1000;     % Loss due to gravity (Velocity)[m/s]


%plot(t,h1,'b');
%hold on;
%grid on;
%plot(t,h,'.b');
%xlabel('time[s]');
%ylabel('Altitude[km]');
% text(80,5,'','Color',[0 0 1], 'VerticalAlignment','middle',...
%	'HorizontalAlignment','left','FontSize',14 );


elseif i == 2
    
 %   state1= stage
    state1 = [v1(end)*1000,gamma1(end)/180*pi,x1(end)*1000,h1(end)*1000,-vD1(end)*1000,-vG1(end)*1000];

    [t,state_2] = ode45(@RocketDynEq2,t_range,state1) ;
v2     = state_2(:,1)/1000;      % Velocity [km/s]
gamma2 = state_2(:,2)*180/pi;    % Flight path angle  [deg]
x2     = state_2(:,3)/1000;      % Downrange distance [km]
h2     = state_2(:,4)/1000;      % Altitude[km]
vD2    = -state_2(:,5)/1000;     % Loss due to drag (Velocity)[m/s]
vG2    = -state_2(:,6)/1000;     % Loss due to gravity (Velocity)[m/s]
%plot(t,h2,'r');
%hold on;
%grid on;
%plot(t,h1,'.r');
%xlabel('time[s]');
%ylabel('Altitude[km]');
% text(80,5,'','Color',[0 0 1], 'VerticalAlignment','middle',...
%	'HorizontalAlignment','left','FontSize',14 );

elseif i==3
    
        state2 = [v2(end)*1000,gamma2(end)/180*pi,x2(end)*1000,h2(end)*1000,-vD2(end)*1000,-vG2(end)*1000];
    [t,state_3] = ode45(@RocketDynEq2,t_range,state2) ;
v3     = state_3(:,1)/1000;      % Velocity [km/s]
gamma3 = state_3(:,2)*180/pi;    % Flight path angle  [deg]
x3     = state_3(:,3)/1000;      % Downrange distance [km]
h3     = state_3(:,4)/1000;      % Altitude[km]
vD3    = -state_3(:,5)/1000;     % Loss due to drag (Velocity)[m/s]
vG3    = -state_3(:,6)/1000;     % Loss due to gravity (Velocity)[m/s]
%plot(t,h3,'g');
%hold on;
%grid on;
%plot(t,h2,'.g');
%xlabel('time[s]');
%ylabel('Altitude[km]');
% text(80,5,'','Color',[0 0 1], 'VerticalAlignment','middle',...
%	'HorizontalAlignment','left','FontSize',14 );


end

t0=tf;
m0=m0-m_flow(1,i)*tburn(1,i);
end

    if h3(end) >= h_orb-5000 && h3(end) <= h_orb+5000
            z = 1;
    else
        v_reached = (niu/(Re+(h3(end)*1000)))^(0.5);
        v_diference = v_reached-v_orbit;
        delta_v = delta_v + v_diference;
        fprintf('\n %4.2f ',h3(end))
        fprintf('\n %4.2f ',v_reached)
        fprintf('\n %4.2f ',v_diference)
        fprintf('\n %4.2f ',delta_v)
        pause();
    end
end

dist('acabou')


% VEGA Rocket: First Stage P80
%fprintf('\n VEGA Rocket: First Stage P80\n')
%fprintf('\n Propellant mass           = %4.2f [kg]',m_prop)
%fprintf('\n Gross mass                = %4.2f [kg]',m_stage_gross(1))
%fprintf('\n Isp                       = %4.2f [s]',Isp)
%fprintf('\n Thrust(mean)              = %4.2f [kN]',T/1000)
%fprintf('\n Initial flight path angle = %4.2f [deg]',gamma0*180/pi)
%fprintf('\n Final speed               = %4.2f [km/s]',v(end))
%fprintf('\n Final flight path angle   = %4.2f [deg]',gamma(end))
%fprintf('\n Altitude                  = %4.2f [km]',h(end))
%fprintf('\n Downrange distance        = %4.2f [km]',x(end))
%fprintf('\n Drag loss                 = %4.2f [km/s]',vD(end))
%fprintf('\n Gravity loss              = %4.2f [km/s]',vG(end))
%fprintf('\n');