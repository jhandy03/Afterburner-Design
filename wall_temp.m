% ----------------------------------------------
%   Written by Jordan Handy
%   03/21/2026
%   Given the outer wall temperature of the afterburner, 
%   code calculates the inner wall temp using 1-D heat transfer
% ----------------------------------------------

clc;clear;close all;

data = readtable('h_lookup.csv');

% Time info
Nt = 2e6; % number of timesteps
time_final = 10; % end time (sec)
t = linspace(0,time_final,Nt); % time vector
dt = t(2)-t(1); % timestep

% Position info
Nx = 500; % number of spatial points
p_inner = 0.0355; % inner position within the annulus (m)
p_outer = p_inner + 0.005; % outer position of the annulus (m)
x = linspace(p_inner,p_outer,Nx); % position vector
dx = x(2)-x(1); % position step

% Gas info
T_inner = 2200; % combustion gas temp (K)
T_outer = 300; % ambient gas temp (K)
T = ones(Nt,Nx)*300; % temperature setup

% Material info
k = 15; % thermal conductivity of wall (W/mK)
rho = 8000; % density of wall (kg/m^3)
Cp = 510; % specific heat of wall (J/kgK)

% Other setup info
h_in = 1000; % convection heat transfer coef for inner gas to wall (W/m^2K)
h_out = 10; % convection heat transfer coef for outside wall to amb (W/m^2K)
CFL = k/(rho*Cp)*dt/dx^2; % Courant-Friedrichs-Lewy Number

% Run a check on the CFL number before starting calculations
if CFL >= 0.5
    warning('For a stable solution, the CFL number needs to be below 0.5. The CFL number is currently %.3f \n',CFL);
    fprintf('Check the following variables for incorrect values: \n%.4f \n%.4f \n%.4f\n',k,rho,Cp);
    CFL_check_override = input('Continue anyways? (y/n)','s');
    if strcmpi(CFL_check_override,'n')
        error('Stopping');
    end
end

