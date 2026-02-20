%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Before running this code, you will need the following:
%   - Combustion Toolbox - Alberto Cuadra Lara
%   - Optimization Toolbox
%   - Symbolic Math Toolbox (depricated)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;close all;

%For the following, nozzle exit is treated as section '2' and turbine exit is
%treated as '1'

%Given / Known
d_turbine = 71e-3; %m
d_nozzle = 51e-3; %m
gamma = 1.4;
R = 287; %J/kgK

%Find Area at the nozzle exit (2) and turbine exit (1)
A1 = (pi*d_turbine^2)/4;
A2 = (pi*d_nozzle^2)/4;

%Convert from C to K
T2 = 680;
T2 = T2 + 273.15;

%Convert from km/h to m/s
u2 = 1700;
u2 = u2*1000*(1/3600);

M2 = u2/sqrt(gamma*R*T2);

%Isentropic Area Relations
a2astar = IsenAAstar(M2,gamma);
Astar = A2/a2astar;

%Finds M1 using area relations. Outputs a vector but answer takes only
%real, subsonic answers (see M1 = M1(1))
M1 = fsolve(@(M1) A1/Astar - IsenAAstar(M1, gamma), 0.5,optimoptions('fsolve','Display','off'));

%Find the stagnation temperature for the nozzle and then finds temp at
%section '1'
T2T0 = IsenStagTemp(M2,gamma);

T0 = T2/T2T0;

T1T0 = IsenStagTemp(M1,gamma);
T1 = T1T0*T0;

%calculate velocity coming out of the turbine => goes into the afterburner
u1 = M1*sqrt(gamma*R*T1);

%%%%%%%%%%%%%%%%%%%%% Flame holder Geometry %%%%%%%%%%%%%%%%%%%%%
% theta = 15; %standard value (deg)
% D = 0.25; %got from textbook (height of flame holder)
% D = D*0.0254;%Converts from in to m
b_ratio = 0.2; % blockage ratio given in flameholder textbook - 20%
fh_height = (d_turbine/2)*b_ratio; %height of the flame holder

% H = AB_outer/2;
% V2max = 0.9516*H/tc; % H is the radius of the AB, tc is the ignition time parameter




%%%%%%%%%%%%%%%%%%%%% Combustion Calculations %%%%%%%%%%%%%%%%%%%%%

%Given
mdot_in = 300; %in L/s
mdot_in = mdot_in*10^-3; %m^3/s
mdot_fuel = 540; %g/min
mdot_fuel = mdot_fuel*(1/1000)*(1/60); %kg/s
rho_air = 1.225; %density of air at sea level. Will need to be more accurate for optimization

%Find the mass flow rate of the incoming air only
% mdot_total = rho_air*mdot_in;
% mdot_air = mdot_total - mdot_fuel;
mdot_air = 0.34; % kg/s from Radek

f_act = mdot_fuel/mdot_air;

% P3 = 2.179;%atm %Pressure needs to be in bar to pass into the combustion toolbox
% P3 = P3*1.01325;
P3 = 330e3;
P3 = P3/100000;

%%%%%%%%%%%%%%%%%%%%% Combustion Toolbox %%%%%%%%%%%%%%%%%%%%%
MW_carbon = 12.011;
MW_oxygen = 15.999;
MW_nitrogen = 14.007;
MW_hydrogen = 1.0078;
MW_fuel = 12*MW_carbon + 24*MW_hydrogen;
MW_oxidizer = 2*MW_oxygen + 3.74*2*MW_nitrogen; %Not sure if this is correct (ie supposed to not include the nitrogen)

f_st = MW_fuel/(18*MW_oxidizer);
phi = f_act/f_st;
mol_air_new = 18/phi;

% Import packages
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.equilibrium.*
import combustiontoolbox.common.Units

DB = NasaDatabase();

system = ChemicalSystem(DB);
mix = Mixture(system);

set(mix,{'Jet_AbLb'},'fuel', 1.723);
set(mix,{'O2','N2'},'oxidizer',[31.02,315.969]/31.02); %needs to normalize around 1 mol of O2 (Alberto said this)
mixArray1 = setProperties(mix,'temperature',T1,'pressure',P3,'equivalenceRatio',1); % 1=> stoiciometric combustion

mixArray2 = mixArray1.copy(); %For the products

solver = EquilibriumSolver('problemType','HP'); %adiabatic and constant P

combustion = solver.solveArray(mixArray2);

Tpost = combustion.T;
Ppost = combustion.p;
gammapost = combustion.gamma;
fprintf('Stagnation Temperature = %f K\n',Tpost);
Ppost = Units.convert(Ppost, 'bar','Pa');
fprintf('Stagnation Pressure = %f kPa\n',Ppost*10^-3);

fprintf('New gamma = %f\n',gammapost);

T0T0star = RayStagTemp(M1,gamma);
T0star = T0/T0T0star;

TpostT0star = Tpost/T0star;
Mpost = fsolve(@(Mpost) TpostT0star - ((1+gammapost)/(1+gammapost*Mpost^2))^2*Mpost^2*(1+((gammapost-1)/2)*Mpost^2)/(1+((gammapost-1)/2)), 0.8,optimoptions('fsolve','Display','off'));
AabAstar = IsenAAstar(Mpost,gammapost);
A1 = pi*((35.5)*10^-3)^2;
Astar = A1/AabAstar;
Patm = 95e3; % Pa
Pcomp = 330e3; % need to redeclare since the 330kPa value used above is converted to bar for combustion toolbox stuff
PatmP0in = Patm/Pcomp;
Mexit = fsolve(@(Mexit) PatmP0in - IsenStagPressure(Mexit, gammapost), 2,optimoptions('fsolve','Display','off'));
AeAstar = IsenAAstar(Mexit,gammapost);
Aexit = AeAstar*(AabAstar)^-1*A1;
dexit = sqrt((4*Aexit)/pi);
dthroat = sqrt((4*Astar)/pi);
fprintf('Throat diameter = %.3f mm\n',dthroat*10^3);
fprintf('Exit diameter = %.3f mm\n',dexit*10^3);
Pstatic_post = Ppost*IsenStagPressure(Mpost,gammapost);
fprintf('Static Pressure = %f Pa\n',Pstatic_post);

% Thrust Calc
Texit = Tpost*IsenStagTemp(Mexit,gammapost);
vexit = Mexit*sqrt(gammapost*R*Texit);
mdot_fuel_AB = 930.42; % g/min
mdot_fuel_AB = mdot_fuel_AB*1/1000*1/60; % convert to kg/s
mdot_exit = mdot_air + mdot_fuel + mdot_fuel_AB;
Pexit = Pcomp*IsenStagPressure(Mexit,gammapost);
Thrust = mdot_exit*vexit + (Pexit-Patm)*Aexit;
fprintf('Total Thrust = %.2f N\n',Thrust);
per_diff = 100*(Thrust/220-1);
fprintf('This is an %g%% increase over the baseline of 220N\n',per_diff);


Cp = gammapost*R/(gammapost-1);
fprintf('Cp = %f\n',combustion.cp*10^-1);
fprintf('Flame holder height = %.3f mm\n',fh_height*10^3);




function aas = IsenAAstar(M,gamma)
    aas = (1/M)*(2/(gamma+1)*(1+((gamma-1)/2)*M^2))^((gamma+1)/(2*(gamma-1)));
end

function TT0 = IsenStagTemp(M,gamma)
    TT0 = (1+((gamma-1)/2)*M^2)^-1;
end

function T0T0star = RayStagTemp(M,gamma)
    T0T0star = ((1+gamma)/(1+gamma*M^2))^2*M^2*(1+((gamma-1)/2)*M^2)/(1+((gamma-1)/2));
end

function PP0 = IsenStagPressure(M,gamma)
    PP0 = (1+((gamma-1)/2)*M^2)^-(gamma/(gamma-1));
end