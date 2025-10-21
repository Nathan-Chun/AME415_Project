%Code working Session 1
% Nathan Chun, Lilly Foord, Elijah Romero
%% --------------------- 11.1 Top Code Input ------------------------ %%

% Boundary Conditions
P01 = 185000.00; % Pa
T01 = 800.00; %K
P3 = 170000.00; %Pa
mdot = 3.00 %kg/s
Nrpm = 10000.00; %RPM

% Geometry and Design Parameters
Nb2 = 85.00 % # of Airfoils stator
Nb3 = 80.00 % # of Airfoils rotor
Phi = 0.45 % Flow parameter
Psi = 1.10 % Loading factor
BL = 0.90 % Flow blockage

% Initial Guess at Efficiency and Loss Factors
ett = 0.900 % Total to total efficiency
Kloss_N = 0.998 % Stator total pressure loss factor
Kloss_R = 0.9920 % Rotor total pressure loss factor


%% --------------------- 11.2 First Pass, Include Loss System Calculations ------------------------ %%
gamma = 1.4056;
Rgas = 4124.400
cp = 14292.300

T3is = T01*(P3/P01)^((gamma-1)/gamma);
C0 = sqrt(2*cp*(T01-T3is));
alpha2 = atan(Phi/Psi);
Rc = 1-Psi/2;
Ca = sqrt((C0^2/2)/(cot(alpha2)/(Phi*ett)+0.5));
U = Ca/Phi;
delh0is = C0^2/2 - Ca^2/2;
w = delh0is * ett;
pow = w* mdot;
rmean = 60*U/(2*pi*Nrpm); 
% Pick up at Lecture 4 Eqn 20


