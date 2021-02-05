%% Author: Tyler Candler
%% ASEN 3200 Orbital Mechanics
%% Lab O1
%% Written: 1/20/2021



%% Housekeeping
clear; close; clc;

%% ODE call

%Sets constants
constants = calculateConstants();

% Time
T = 2*pi*sqrt(constants.a^3/constants.mu);
tfin = T*1.001;
tspan = [0:tfin];
% Initial Conditions
% rx0     = constants.a*(1-constants.e^2)/(1+constants.e*cos(constants.theta0));
rx0 = constants.rp;
vy0 = sqrt(2*constants.mu/rx0-constants.mu/constants.a);
IC  = [rx0,0,0,0,vy0,0];
options = odeset('Reltol',1e-15,'AbsTol',1e-15);



%Calling ode45 for Phase1
[t,Phase1] = ode45(@(t,input) orbitalcalc(t,input,constants),tspan,IC);


x = Phase1(1:end,1);
y = Phase1(1:end,2);
z = Phase1(1:end,3);
x_dot = Phase1(1:end,4);
y_dot = Phase1(1:end,5);
z_dot = Phase1(1:end,6);
figure(1)
plot(x,y)
hold on 
plot(1, 1, '.r', 'MarkerSize',69)
plot(x(1), y(1), '.b', 'MarkerSize',20)
% circles(0,0,500,'color','black')

title("One Period of Satellites Orbit")
xlabel("Semimajor Axis (m)")
ylabel("Semiminor Axis (m)")
legend("Satellite Orbit","Earth","Satellite Initial Location")


L = [];
m = 1000
p = [];
 for i = 1:length(x)
     r = [x(i),y(i),z(i)];
     v = [x_dot(i), y_dot(i), z_dot(i)];
     L(i) = norm(cross(r, v));
 end 
figure(2)
plot(t,L)
title("Angular Momentum over one orbital period")
xlabel("Time (s)")
ylabel("Magnitude of Angular Momentum (kg*m^2/s)")
ylim([2e11 3e11])



e_vec = [];
for i = 1:length(x)
    r = [x(i),y(i),z(i)];
    v = [x_dot(i), y_dot(i), z_dot(i)];
    temp_vec = (1/constants.mu)*((norm(v)^2*r)-(dot(r,v)*v)) - (r/norm(r));
    e_vec(i) = norm(temp_vec);
end
figure(3)
plot(t,e_vec)
title("Eccentricity Over One Orbital Period")
xlabel("Time (s)")
ylabel("Magnitude of Eccentricity Vector")
ylim([0 1])



%% Orbit Function  

function [results] = orbitalcalc(t,input,constants)

    %Calculate the current states
    rx = input(1);
    ry = input(2);
    rz = input(3);
    vx = input(4);
    vy = input(5);
    vz = input(6);
    
    %Magnitude of radius:
    r  = sqrt(rx.^2+ry.^2+rz.^2);
    
    x_vel = vx;
    y_vel = vy;
    z_vel = vz;
    x_accel = -constants.mu*rx/r.^3;
    y_accel = -constants.mu*ry/r.^3;
    z_accel = -constants.mu*rz/r.^3;

    %Equate all ODES to one vector for ODE45
    results = [x_vel; y_vel; z_vel; x_accel; y_accel; z_accel];
end



%% Constants Structure
function constants = calculateConstants(constants)

    constants.a = 37.9735 * 6378 *1000; % semimajor axis (m)
    constants.e = 0.587; %orbit eccentricity

    constants.ra = constants.a *(1+constants.e); %radius of apoapsis (m) (longest)
    constants.rp = constants.a *(1-constants.e); %radius of periapsis (m) (shortest)

    constants.me = 5.972e24 % (kg) (mass of earth)
    constants.G  = 6.67408e-11; %universal gravitational constant
    constants.mu = constants.G*constants.me; %mu of the earth

    constants.i = 0; %(degrees) orbit plane inclination
    constants.theta0 = 0;
    
end