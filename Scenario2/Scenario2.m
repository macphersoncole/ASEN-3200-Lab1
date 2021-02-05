%% Author: Tyler Candler
%% ASEN 3200 Orbital Mechanics
%% Lab O1
%% Written: 1/20/2021


%% Housekeeping
clear; close; clc;

%% Initial Conditions and constants
t = 13.5*60; %time of interest (s)
t0 = 0; %time of periapsis

rp = 1840*1000; %m
P = 6.3*3600; %seconds
r0 = 2400*1000; %m, moving away from periapsis

%% Solving for orbit parameters

G  = 6.67408e-11; %universal gravitational constant
Mm = 7.34767309e22; %mass of the moon (kg)
mu = G*Mm; %mu of the moon

n = 2*pi/P; %mean motion
a = (mu/n^2)^(1/3); % semimajor axis (m)

ra = 2*a-rp; %radius of apoapsis (m)

e = abs(ra-rp)/(ra+rp); %eccentricity

M = n*(t - t0);

%% solving for E

E = M
    for k = 1:100
        E = M + e*sin(E);
    end
    
%% Solving for location at t = 15 minutes after periapsis

r = a*(1-e*cos(E))

theta = sqrt((1+e)/(1-e))*tan(E/2);
theta = atan(theta)*2;

% theta = acos((a*(1-e^2)/(r*e))-(1/e))

theta = (theta/pi)*180;

r = r/1000;
h = r -1737.1


%% delta t = 0.5 P + 13.5 min