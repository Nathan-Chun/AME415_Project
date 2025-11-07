%Code working Session 1
% Nathan Chun, Lilly Foord, Elijah Romero
%% --------------------- 11.1 Top Code Input ------------------------ %%
clear
% Boundary Conditions
P01 = 11617665.50; % Pa
T01 = 460.00; %K
P3 = 10804084.00; %Pa
mdot = 32.00 %kg/s
Nrpm = 15000.00; %RPM

% Geometry and Design Parameters
Nb2 = 85.00 % # of Airfoils stator
Nb3 = 80.00 % # of Airfoils rotor
Phi = 0.05 % Flow parameter
Psi = 0.99310 % Loading factor
BL = 0.90 % Flow blockage

% Initial Guess at Efficiency and Loss Factors
ett_guess = 0.900 % Total to total efficiency
Kloss_N_guess = 0.998 % Stator total pressure loss factor
Kloss_R_guess = 0.9920 % Rotor total pressure loss factor


%% --------------------- 11.2 First Pass, Include Loss System Calculations ------------------------ %%
clc;
gamma = 1.4056;
Rgas = 4124.400;
cp = 14292.300;

T3is = T01*(P3/P01)^((gamma-1)/gamma);
C0 = sqrt(2*cp*(T01-T3is));
alpha2 = atan(Phi/Psi); % REMEMBER TO CONVERT TO DEGREES FOR PRINTED OUTPUT
Rc = 1-Psi/2;
Ca = sqrt((C0^2/2)/(cot(alpha2)/(Phi*ett_guess)+0.5));
U = Ca/Phi;
delh0is = C0^2/2 - Ca^2/2;
w = delh0is * ett_guess;
pow = w* mdot;
% pow = 4800;
rmean = 60*U/(2*pi*Nrpm); 

% --------------------- Stator Exit Kinematics ------------------------ %

% Starting from Lecture 4 Eqn 20
C2 = Ca/sin(alpha2);
C2u = Ca*cot(alpha2);
W2u = C2u - U;
W2a = Ca;
C2a = Ca;
W2 = sqrt(W2u^2 + W2a^2);
alpha2p = atan(W2a/W2u);

% --------------------- Rotor Exit Kinematics ------------------------ %

% Starting from Lecture 4 Eqn 26
W3a = Ca; C3a = Ca;
W3u = -U; W3 = sqrt(W3a^2 + W3u^2);
C3 = Ca; C3a = Ca;
alpha3p = atan(W3a/W3u);

% --------------------- Thermodynamic Quantities Stator Exit ------------------------ %

% Starting from Lecture 4 Eqn 31
T2 = T01-C2^2/(2*cp);
T02=T01;
a2 = sqrt(gamma*Rgas*T2);
M2 = C2/a2; 
M2a=M2;
M2w = W2/a2;
T2w = T2*(1 + ((gamma-1)/2)*M2w^2);
T3w=T2w;
T3 = T2w - (W3^2)/(2*cp);
a3 = sqrt(gamma*Rgas*T3);
M3w = W3/a3;
M3 = C3/a3; M3a = M3;

% --------------------- Pressures/Convergence ------------------------ %
P02 = Kloss_N_guess * P01;
P2 = P02/(1+(gamma-1)/2*M2^2)^((gamma)/(gamma-1));
P2w = P2*(1+(gamma-1)/2*M2w^2)^(gamma/(gamma-1));
P3w = Kloss_R_guess * P2w;
P3ver = P3w/(1 + (gamma-1)/2 * M3w^2)^(gamma/(gamma-1));
rho2 = P2/(Rgas*T2);
rho3 = P3/(Rgas*T3);

% --------------------- Geometry ------------------------ %
A2 = mdot/(rho2*C2a); A3 = mdot/(rho3*C3a);
L2 = A2/(2*pi*rmean*BL); 
L3 = A3/(2*pi*rmean*BL);
L2D = L2/(rmean*2); 
L3D = L3/(rmean*2);

% --------------------- stator/nozzle ------------------------ %
sc0stator = 0.427 + abs(alpha2*180/pi)/58 - (abs(alpha2*180/pi)/93)^2;
sstator = 2*pi*rmean/(Nb2); % Pitch noz
cstator = sstator/sc0stator;
sbzstator = 0.8/(2*sin(alpha2)^2*(cot(alpha2 - cot(90*pi/180)))); % Assume alpha3 = 90
bzstator = sstator/sbzstator; % bz noz
betas = asin(bzstator/cstator)*180/pi;

% --------------------- Rotor ------------------------ %
sc0 = 0.427+abs(alpha3p*180/pi)/58 - (abs(alpha3p*180/pi)/93)^2;
sc1 = 0.224 + (1.575 - abs(alpha3p*180/pi)/90)*(abs(alpha3p*180/pi)/90);
xi_rot = (90 - abs(alpha2p)*180/pi)/(90-abs(alpha3p)*180/pi);
scopt = sc0 + (sc1 - sc0)*abs(xi_rot)*xi_rot; % sc rot
sbzrot = 0.8/(2*sin(alpha3p)^2*(cot(alpha2p) - cot(alpha3p)));
srot = 2*pi*rmean/Nb3; % Pitch rot
crot = srot/scopt; % Chord rot
bzrot = srot/sbzrot % bz rot
betasrot = asin(bzrot/crot)*180/pi
%%
% --------------------- Loss calculations ------------------------ %
alpha1=90; % degree..Lecture 5 Row 1 calculation...need to follow up
yp0 =Yp0(sc0stator,alpha2*180/pi)
yp1 =Yp1(sc0stator,alpha2*180/pi)
xi =(90-abs(alpha1))/(90-abs(alpha2)*180/pi) % CONFIRM WITH PROF...alpha 1 = 90?
yp_noz=abs(yp0+(xi^2)*(yp1-yp0)) % CONFIRM WITH PROF

mu=DynVisc_H2(T2)
nu=mu/rho2
[KRe,Re_c] = K_Re(C2, cstator, nu)

ys_noz= Ys(alpha1,alpha2*180/pi,cstator,L2)

dL=0; % CONFIRM WITH PROF: IS dL 0 BECAUSE NOZZLE CLEARANCE IS 0 (LEC5 PG 16)
ycl_noz=Ycl(alpha2p*180/pi, alpha3p*180/pi, cstator, L2, 0)

y_noz = KRe*yp_noz + ys_noz + ycl_noz;

Kloss_N=(((1+((gamma-1)*M2^2)/2))^(gamma/(gamma-1)))/(((1+((gamma-1)*M2^2)/2)^(gamma/(gamma-1)))*(y_noz+1)-y_noz)

% --------------------- Loss calculations Rotor ------------------------ %
alpha1=alpha2p*180/pi; % degree..Lecture 5 Row 1 calculation...need to follow up
yp0_rot =Yp0(sc0,alpha3p*180/pi)
yp1_rot =Yp1(scopt,alpha3p*180/pi)
xi_rot =(90-abs(alpha1))/(90-abs(alpha3p)*180/pi) % CONFIRM WITH PROF...alpha 1 = 90?
yp_rot=abs(yp0_rot+(xi_rot^2)*(yp1_rot-yp0_rot)) % CONFIRM WITH PROF

mu3=DynVisc_H2(T3)
rho_rotor = (rho2+rho3)/2;
nu3 = mu3/rho3;
[KRe_rotor,Re_c] = K_Re(W3, crot, nu3)

ys_rot= Ys(alpha1,alpha3p*180/pi,crot,L3)

dL=0.0075; % CONFIRM WITH PROF: Given in Lecture 5
ycl_rot=Ycl(alpha1, alpha3p*180/pi, crot, L3, dL)

y_rot = KRe*yp_rot + ys_rot + ycl_rot;
% y_rot = 0.10828

exp = gamma/(gamma-1);
exp1 =(gamma-1)/(gamma);

Kloss_R=(((1+((gamma-1)*M3w^2)/2))^(gamma/(gamma-1)))/(((1+((gamma-1)*M3w^2)/2)^(gamma/(gamma-1)))*(y_rot+1)-y_rot)
T03 = T3*(1+((gamma-1)/2)*M3^2);

P03 = P3ver*((1+((gamma-1)/2)*M3^2)^exp);

ett = (1-(T03/T01)) / (1-(P03/P01)^exp1); % Check with Prof - input values may be incorrect

%Power Calculations
% massflow = 3;
work = U * C2u;
Power = work * mdot/745.7;

%% --------------------- Set B ------------------------ %%
%% --------------------- Input ------------------------ %%

fscale=P3/P3ver
Ca_new=Ca/fscale
U_new=U/fscale



%% --------------------- Functions ------------------------ %%

% Calc Yp1
function [val,X] = Yp0(sc, alpha2)
    alpha2 = abs(alpha2);
    % Equation 4 
    if alpha2<=30
        scmin = 0.46 + alpha2/77;
    else
        scmin = 0.614+alpha2/130;
    end
    % Equation 5
    X=sc-scmin;

    % Equation 6
    if alpha2<=27
        A = 0.025 + (27-alpha2)/530;
    else
        A=0.025+(27-alpha2)/3085;
    end
    B = 0.1583-alpha2/1640;
    C = 0.08*((alpha2/30)^2-1);
    n=1+alpha2/30;

    % Equation 7
    if alpha2<=30
        val = A + B*X^2 + C*X^3;
    else
        val = A +B*(abs(X))^n;
    end
end

% Calc Yp2
function [val,X] = Yp1(sc, alpha2)
    alpha2 = abs(alpha2);
    % Equation 8
    scmin = 0.224 + 1.575*(alpha2/90) - (alpha2/90)^2;

    % Equation 9
    X=sc-scmin;

    % Equation 6
    A = 0.242 - alpha2/151 + (alpha2/127)^2;
    if alpha2 <= 30
        B = 0.3 + (30 - alpha2)/50;
    else
        B = 0.3 + (30 - alpha2)/275;
    end
    C = 0.88 - alpha2/42.4 +(alpha2/72.8)^2;

    % Equation 7
    val = A + B*X^2 - C*X^3;
end

% Calc K_Re
function [val,Re_c] = K_Re(C2, c, nu)
    % C2: stator velocity
    % c: chord value
    % nu: kinematic viscosity
    Re_c = (C2*c)/nu;

    % surface roughness
    e = 5e-6;

    % critical roughness-based Re
    Re_r=100*c/e;

    % Reynolds number correction factor
    if Re_c <= 1e5
        val = sqrt(1e5/Re_c);
    elseif Re_c > 1e5 && Re_c <=5e5
        val = 1;
    else
        val = 1+((log10(5e5)/log10(Re_r))^2.58 - 1)*(1 - (5e5/Re_c))
    end
    
end

% Ys
function val = Ys(alpha1, alpha2, c, L)
    alpha1 = deg2rad(alpha1);
    alpha2 = deg2rad(alpha2);
    alpham = abs(atan( 2/(cot(alpha2) + cot(alpha1)) ));
    % C_l/sc
    Cl_sc = abs(2*(cot(alpha1)-cot(alpha2))*sin(alpham)); 
    % Ainley loading parameter Z
    Z = abs( (Cl_sc)^2 * ( sin(alpha2)^2/sin(alpham)^3 ) );

    % cascade aspect correction F_AR
    if L/c >= 2.0
        F_AR = c/L;
    else
        F_AR = 0.5*(2*(c/L))^0.7;
    end

    % Ys
    val = 0.0334*F_AR*Z*(sin( abs(alpha2) )/sin( abs(alpha1) ));
end

% C_l;
function val = Ycl(alpha1, alpha2, c, L, dL)
    alpha1 = deg2rad(alpha1);
    alpha2 = deg2rad(alpha2);
    alpham = abs(atan( 2/(cot(alpha2) + cot(alpha1)) ));
    % C_l/sc
    Cl_sc = abs(2*(cot(alpha1)-cot(alpha2))*sin(alpham)); 
    % Ainley loading parameter Z
    Z = abs( (Cl_sc)^2 * ( sin(alpha2)^2/sin(alpham)^3 ) );
    Fb=0.36;

    val = Fb*Z*(c/L)*(dL*L/c)^0.78
    
end

% dynamic viscosity, mu
function [DynVisc] = DynVisc_H2(Temp)
    % From Lecture 5 eqns 14-16
    mu_ref=8.76e-6;
    T_ref=293.85;
    C_ref=72;
    
    % T = inlet temp
    DynVisc=mu_ref*((T_ref+C_ref)/(Temp+C_ref))*((Temp/T_ref)^(3/2));
end

