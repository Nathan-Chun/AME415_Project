clear, clc

disp('Running main script');
run('Coding_Working_Session_1.m');  % calls another script or function

fileID = fopen('ordered_variables.txt','w');

%% --------------------- PRINT VALIDATION: SET A ------------------------ %%
%% --------------------- PRINT VALIDATION INPUT SECTION ------------------------ %%
fprintf('\n=========== INPUT ==============================================\n');
fprintf(' Po1          =  %12.4f Pa\n', P01);
fprintf(' To1          =  %12.4f K\n', T01);
fprintf(' P3           =  %12.4f Pa\n', P3);
fprintf(' Massflow     =  %12.4f Kg/s\n', mdot);
fprintf(' Shaft Speed  =  %12.4f RPM\n', Nrpm);
fprintf(' Noz In angle =  %12.4f deg\n', 90.0);   % fixed inlet angle per calibration case

fprintf('\n------- fluid thermodynamic properties ------------\n');
fprintf(' gamma =  %12.4f\n', gamma);
fprintf(' Rgas  =  %12.4f J/kg-K\n', Rgas);
fprintf(' Cp    =  %12.4f J/kg-K\n', cp);

fprintf('\n------- Choice of design parameters ---------------\n');
fprintf(' Phi     =  %12.4f\n', Phi);
fprintf(' Psi     =  %12.4f\n', Psi);
fprintf(' Stator nr vanes  =  %12.4f\n', Nb2);
fprintf(' Rotor nr blades  =  %12.4f\n', Nb3);

fprintf('\n------- First guess values of parameters ----------\n');
fprintf(' Efficiency  =  %12.4f\n', ett_guess);
fprintf(' Kloss_N     =  %12.4f\n', Kloss_N_guess);
fprintf(' Kloss_R     =  %12.4f\n', Kloss_R_guess);
fprintf(' Blockage    =  %12.4f\n', BL);
fprintf('===================================================\n\n');


fprintf('\n=========== FIRST PASS Validation Output ===========\n');
fprintf('Co     = %12.4f m/s\n', C0);
fprintf('U      = %12.4f m/s\n', U);
fprintf('Ca     = %12.4f m/s\n', Ca);

fprintf('\n----- Stator Exit Kinematics -----\n');
fprintf('C2     = %12.4f m/s\n', C2);
fprintf('Cu2    = %12.4f m/s\n', C2u);
fprintf('Ca2    = %12.4f m/s\n', C2a);
fprintf('Wu2    = %12.4f m/s\n', W2u);
fprintf('Wa2    = %12.4f m/s\n', W2a);
fprintf('W2     = %12.4f m/s\n', W2);
fprintf('alfa2  = %12.4f deg\n', alpha2*180/pi);
fprintf('alfa2p = %12.4f deg\n', alpha2p*180/pi);

fprintf('\n----- Rotor Exit Kinematics -----\n');
fprintf('Wa3    = %12.4f m/s\n', W3a);
fprintf('Wu3    = %12.4f m/s\n', W3u);
fprintf('W3     = %12.4f m/s\n', W3);
fprintf('Ca3    = %12.4f m/s\n', C3a);
fprintf('Cu3    = %12.4f m/s\n', 0);
fprintf('C3     = %12.4f m/s\n', C3);
fprintf('alfa3p = %12.4f deg\n', alpha3p*180/pi);

fprintf('\n----- Thermodynamic Quantities (Stator Exit) -----\n');
fprintf('a2     = %12.4f m/s\n', a2);
fprintf('Mw2    = %12.4f\n', M2w);
fprintf('Ma2    = %12.4f\n', M2a);
fprintf('T2     = %12.4f K\n', T2);
fprintf('To2    = %12.4f K\n', T02);
fprintf('Tw2    = %12.4f K\n', T2w);

fprintf('\n----- Thermodynamic Quantities (Rotor Exit) -----\n');
fprintf('a3     = %12.4f m/s\n', a3);
fprintf('Mw3    = %12.4f\n', M3w);
fprintf('Ma3    = %12.4f\n', M3a);
fprintf('T3     = %12.4f K\n', T3);
fprintf('To3    = %12.4f K\n', T03);
fprintf('Tw3    = %12.4f K\n', T3w);

fprintf('\n----- Pressures and Densities -----\n');
fprintf('Po1    = %12.4f Pa\n', P01);
fprintf('Po2    = %12.4f Pa\n', P02);
fprintf('P2     = %12.4f Pa\n', P2);
fprintf('Pw2    = %12.4f Pa\n', P2w);
fprintf('rho2   = %12.4f kg/m^3\n', rho2);
fprintf('Pw3    = %12.4f Pa\n', P3w);
fprintf('P3    = %12.4f Pa\n', P3);
fprintf('P3ver  = %12.4f Pa\n', P3ver);
fprintf('P03    = %12.4f Pa\n', P03);
fprintf('rho3   = %12.4f kg/m^3\n', rho3);

fprintf('\n----- Geometry -----\n');
fprintf('rmean  = %12.5f m\n', rmean);
fprintf('L2     = %12.5f m\n', L2);
fprintf('L3     = %12.5f m\n', L3);
fprintf('L2/D   = %12.5f\n', L2D);
fprintf('L3/D   = %12.5f\n', L3D);

fprintf('\n----- Stator/Nozzle Geometry -----\n');
fprintf('(s/c)noz = %12.4f\n', sc0stator);
fprintf('(s/bz)noz = %12.4f\n', sbzstator);
fprintf('pitch noz = %12.4f m\n', sstator);
fprintf('chord noz = %12.4f m\n', cstator);
fprintf('bz noz = %12.4f m\n', bzstator);
fprintf('beta_s noz = %12.4f deg\n', betas);

fprintf('\n----- Rotor Geometry -----\n');
fprintf('(s/c)0 rot = %12.4f\n', sc0);
fprintf('(s/c)1 rot = %12.4f\n', sc1);
fprintf('xi = %12.4f\n', xi_rot);
fprintf('(s/c)rot = %12.4f\n', scopt);
fprintf('(s/bz)rot = %12.4f\n', sbzrot);
fprintf('pitch rot = %12.4f m\n', srot);
fprintf('chord rot = %12.4f m\n', crot);
fprintf('bz rot = %12.4f m\n', bzrot);
fprintf('beta_s rot = %12.4f deg\n', betasrot);
fprintf('---------------------------------------------------\n');

%% --------------------- LOSS CALCULATIONS ------------------------ %%
fprintf('\n------------ LOSS Calculations -----------------------\n');

% ---- STATOR ----
fprintf('------- STATOR --------------------\n');
fprintf('Yp_0       = %12.5f\n', yp0);
fprintf('Yp_1       = %12.5f\n', yp1);
fprintf('xi         = %12.5f\n', xi);
fprintf('Yp_noz     = %12.5f\n', yp_noz);
fprintf('KRe_noz    = %12.5f\n', KRe);
fprintf('Ys_noz     = %12.5f\n', ys_noz);
fprintf('Ycl_noz    = %12.5f\n', ycl_noz);
fprintf('Y_noz      = %12.5f\n', y_noz);
fprintf('New Kloss_N= %12.5f\n', Kloss_N);

% ---- ROTOR ----
fprintf('\n------- ROTOR --------------------\n');
fprintf('Yp_0       = %12.5f\n', yp0_rot);
fprintf('Yp_1       = %12.5f\n', yp1_rot);
fprintf('xi         = %12.5f\n', xi_rot);
fprintf('Yp_rot     = %12.5f\n', yp_rot);
fprintf('KRe_rot    = %12.5f\n', KRe_rotor);
fprintf('Ys_rot     = %12.5f\n', ys_rot);
fprintf('Ycl_rot    = %12.5f\n', ycl_rot);
fprintf('Y_rot      = %12.5f\n', y_rot);
fprintf('New Kloss_R= %12.5f\n', Kloss_R);

fprintf('\n T-T Efficiency = %12.4f\n', ett);
fprintf(' Power          = %12.4f hp\n', Power);
fprintf('===============================================================\n');