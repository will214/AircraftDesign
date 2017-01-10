function [TOW, LW, OWE, Weights] = weights(N_pax, MTOW, W_fuel, W_eng, M, fixed_weights)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File : weights.m
%
% Function : Based on aircraft geometry and evaluated weights, this
%   function evaluates each component's weight for computation of total
%   aircraft weight. This function has to be evaluated in iterations with
%   mission profile (mp.m) in order to converge on Structure and Fuel
%   weights. Structure weights is calculated with Kroo's equations,
%   Raymer's are used for systems weights, and Torenbeek for hydraulics and
%   avionics. Coefficients are applied to account for composite materials
%   and advanced technology factors.
%
% Inputs : (1) Number of passengers (to differentiate aircraft in family)
%          (2) Maximum Takeoff Weight early estimate
%          (3) Total Fuel Weight (mission + reserve + trapped fuel)
%          (4) Engines Weight
%          (5) Cruise Speed (Mach)
%          (6) Fixed Weights : matrix of weights calculated on biggest
%            aircraft in order to fix those weights on smallest aircraft
%            for commonality between both in aircraft family
%
% Outputs : (1) Maximum Takeoff Weight (recalculated)
%           (2) Maximum Landing Weight (recalculated, OWE + max payload)
%           (3) Operational Weight Empty
%           (4) Weight matrix for each comp used later for aircraft balance
%
% Author : Martin Lafrance
%
% Date : January 15th, 2015
% Last Edit : April 14th, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V_stall = 199.07; %stall speed, in knots (193.46 for 76pax, 199.07 for 100pax)

global output Parameters

Bulk_LS3 = output(431);
output_ = vertcat(output(1:429), output(434:441), output(446:end));

if N_pax == 77
    max_pax = 85;   % High Density Configuration
    W_ops = 2110;   % Operational weights - HARD CODED CONSTANT
elseif N_pax == 101
    max_pax = 120;  % High Density Configuration
    W_ops = 2380;   % Operational weights - HARD CODED CONSTANT
else
    error('Nombre de passagers invalide');
end

% number of Business Class passengers (if not high density conf)
nb_biz = 0;
if Parameters.HighDensity == 0
    nb_biz = sum(output_([745 764]).*output_([746 765]));
    max_pax = N_pax;
end
nb_econ = max_pax-nb_biz;


%% 
L_cabin = convlength(Bulk_LS3-output_(423), 'in', 'ft');

% Variables for RAYMER calculations
Swet_nacelle = pi*output_(160)*output_(161)*convlength(1,'in', 'ft')^2;%RAYMER %pi*D*L
N_f = 7;%RAYMER             %"Number of functions performed by controls (typically 4-7)"
N_m = 0;%RAYMER             %"Number of mechanical functions (typically 0-2)"
S_ele = 0.9*0.25*output_(335)*convlength(1,'m','ft')^2;%RAYMER %0.9(span)*0.25(chord)*S_h
S_rud = 0.5*0.32*output_(340)*convlength(1,'m','ft')^2;%RAYMER %0.5�9span)*0.32(chord)*S_v
S_aileron = 0.25*0.45* output_(346)*convlength(1,'m','ft')^2;%RAYMER %0.45(span)*0.2(chord)*S_ref
W_apu_ui = 240;%RAYMER      % Honeywell RE220 (crj900)
R_kva = 60;%RAYMER          % system electrical rating (kV*A) typically 40-60 for transport
N_gen = 2;%RAYMER           % number of generators (typ = nb_eng)
W_uav = 1400;%RAYMER        % lb - uninstalled avionics (800-1400 lb)
W_c = N_pax*30;%RAYMER      % lb - cargo weight
L_ec = convlength(sum(output_([26 27 30 34])), 'in', 'ft')-14;%RAYMER %(output_(159) - output_(833))*nb_eng;   %FOR OPEN ROTOR (OR_apex-cockpit_seats)*nb_eng

%% HARD CODED VARIABLES
rho0 = 0.0023769;   % slug/ft^3 - SL ISA air density
g = 32.174;         % ft/s^2 - gravitational constant
h_v = 1;            % 0 for conventional - 1 for T tail
Kflap = 1;          % 1.02 for Fowler flaps - else 1
Kspoil = 1.02;      % 1.02 with spoilers - 1 without
Kgear_w = 0.95;     % 1 if attached to wing - 0.95 if not
Keng = 1;           % 0.95 if attached to wing - 1 if not
Kng = 1.017;        % 1.017 for pylon-mounted nacelle, 1 otherwise
Kh = 0.12;          % 0.11 hydrolics for flaps, brakes& retracts, 0.12 w/ flt ctrls
Nlim = 2.5;         % limit load
Nult = Nlim*1.5;    % ultimate load
Kgear = 1;          % Kgear = 1 for gear on wing, 1.07 elsewhere
Klg = 1;            % 1.12 for gear in fuselage, 1.0 elsewhere
Kmp = 1;            % 1.126 for kneeling gear, 1 otherwise
Knp = 1;            % 1.15 for kneeling gear, 1 otherwise
Kfloor = 1;         % Kfloor = 1.1 for cargo floor, else 1
Kgr = 1;            % 1 for low wing, else 1.08
Kuht = 1.143;       % 1.143 for all-moving H-tail, else 1
nb_eng = 2;         % number of engines
nb_tanks = 3;       % number of fuel tanks (2x wing+center tank)
nb_pilot = 2;       % number of pilots
nb_att = 3;         % number of attendants
nb_1stclass = 0;    % number of 1st class pax
l_tc_rad = 4;       % ft - longueur tail cap + rad�me
P_delta = 8;        % psi - cabin pressure differential
nb_mw = 4;          % nombre de roues sur le MLG
nb_nw = 2;          % nombre de roues sur le NLG
nb_mss = 2;         % nombre de shock struts sur le MLG
Kdoor = 1.12;       % 1 no cargo door; 1.06 1 side cargo door; 1.12 2 side cargo doors OR aft clamshell door; 1.25 2 side cargo doors AND aft clamshell door
Ktr = 1.18;         % 1.18 w/ thrust reversers, 1 w/o
Kp = 1.4;           % 1.4 engine w/ propeller, 1 w/o
Kr = 1;             % 1.133 for reciprocating engine (piston engine), 1 otherwise
Ktp = 1;            % 0.793 if turboprop, 1 otherwise
R_bar_z = 0.44;     % Raymer table 16.1, jet transport, fuselage mounted engines
Isc = 3.5;          % KROO - 3.5 for fully powered ctrl surfaces, 2.5 part powered, 1.7 fully aerodynamic

%% HARD CODED CONSTANTS
V_mo = 330;         % kts - MR&O
W_1pilot = 240;     % lb (190 + 50)
W_att = 210;        % lb (170 + 40)
k1stclass = 80;     % lb/siege
kbiz = 58*0.9;          % lb/siege
kecon = 23*0.9;         % lb/siege
kbin = 10;          % lb/pax
koxygen = 1;        % lb/pax
kgalley = 500;      % lb/galley
klav = 370;         % lb/lav
kwardrobe = 150;    % lb/wr
ktrolley = 90;      % lb - 60lb half, 90lb full
kseat_att = 25;     % lb/att
kfurn = 0.5;        % lb/ft^3
fuel_density = 810; % kg/m^3 - (approx) Jet A p22 http://dtic.mil/dtic/tr/fulltext/u2/a132106.pdf

%% Inputs
payload = N_pax*225+175*(max_pax-N_pax);
MZFW = MTOW-W_fuel;
MLW = 0.85*MTOW;
OWE = MTOW-payload-W_fuel;
W_dg = sqrt(MTOW*MZFW);

S_ref = output_(346)*convlength(1,'m','ft')^2;
b = 2*convlength(output_(65), 'in', 'ft');
c = S_ref/b;
c_root = convlength(output_(59), 'in', 'ft');
AR = b^2/S_ref;
lambda = output_(62)/output_(59);     % c_tip/c_root
thick = mean((output_(66:68)+output_(67:69)/2).*[output_(63);(output_(64:65)-output_(63:64))]/b);
thick_root = output_(66);
LAMBDA = output_(327);   %C/2 sweep
L_fuse = convlength(sum(output_([26 27 30 34])), 'in', 'ft')-l_tc_rad;
L_a = L_fuse + 2*output_(107);               % electrical routing distance (see p. 462), L_fuse + 2* eng lat pos
w_fuse = convlength(output_(15), 'in', 'ft');
h_fuse = convlength(output_(16), 'in', 'ft');
d_fuse = convlength((2*h_fuse+w_fuse*12)/3, 'in', 'ft');
S_h = output_(335)*convlength(1,'m','ft')^2;
S_wet_h = output_(306)*convlength(1, 'm', 'ft')^2;
b_h = 2*convlength(output_(99), 'in', 'ft');
AR_h = b_h^2/S_h;
L_h = convlength(output_(336), 'in', 'ft');
thick_h = mean(output_(101:102));
lambda_h = output_(98)/output_(97);
LAMBDA_h = output_(322); %C/2 sweep
HtHv = output_(96);    %vertical position ratio for H tail (0 for conventional, 1 for T tail)
S_v = output_(340)*convlength(1,'m','ft')^2;
L_v = convlength(output_(338), 'in', 'ft');
thick_v = mean(output_(111:112));
AR_v = output_(341);
lambda_v = output_(106)/output_(105);
LAMBDA_v = output_(339); %C/4 sweep
thick_v_root = output_(111);
Swet_fuse = output_(308)*convlength(1,'m','ft')^2;
W_wf = convmass(output_(353)*fuel_density,'kg','lbm');
L_mlg = output_(854);    % in!!!! (divis� par 12 dans l'�quation)
L_nlg = output_(845);    % in!!!! (divis� par 12 dans l'�quation)
W_ec = 2.331*W_eng^.901*Kp*Ktr;
V_i = output_(354)*264.172;                                  %264.172 gal/m^3
V_t = convmass(W_fuel, 'lbm', 'kg')/fuel_density*264.172;   %264.172 gal/m^3
S_cs = S_ele+S_rud+S_aileron;
Iy = (OWE/g)*(R_bar_z*(b+L_fuse)/4)^2;  %Raymer eq. 16.54
Nseats = max_pax*nb_pilot+nb_att;
nb_galley = sum(output_([819 823 827]));
nb_lav = sum(output_([787 791 795 799 803]));
nb_wardrobe = output_(815);
N_Lt = convlength(output_(161), 'in', 'ft');%nacelle length - FOR OPEN ROTOR, ELSE CHANGE REFERENCE
N_w = convlength(output_(160), 'in', 'ft'); %nacelle width (open rotor diameter)
N_t = 3;                    %number of fuel tanks (2 wing + 1 center)
N_c = nb_pilot+nb_att;      % number of crew

% Largeur de la jonction du H_stab
if h_v == 1 %if T-tail
    F_w = output_(112)*convlength(output_(106), 'in', 'ft');  %width = Vtail thickness
else    % Jonction avec le fuselage
    %on enl�ve la proportion de la largeur du fuselage associ�e � la
    %position longitudinale du H-stab
    %F_w = Fuse_Width - (Hstab_apex-Debut_tailcone)/Tailcone_length
    %               * (1-Tailcone_aft_width_ratio)*Fuse_Width
    F_w = convlength(output_(15)-(output_(94)-(output_(30)+output_(26)+output_(27)))...
        /output_(34) * (1-output_(39))*output_(15), 'in', 'ft');
end

% Volume cabine
h_floor = output_(359);%+output_(360); (floor height + p-e floor thickness)
d_fusei = 12*d_fuse;
fun = @(x) sqrt((d_fusei/2)^2-x.^2)-h_floor;
A_cargo = integral(fun, -sqrt(d_fusei^2/4-h_floor^2), sqrt(d_fusei^2/4-h_floor^2));
A_cabin = pi*d_fusei^2/4-A_cargo;
V_cabin = A_cabin*L_cabin/1728; %1728 = 12^3: conversion in^3 -> ft^3

nb_p = max_pax + nb_pilot + nb_att;

nb_trolley = 3*nb_galley;

V_D = (V_mo+30); 
q_Df = V_D^2/295.37;

Kf = 1.1*Kgear*Kfloor;
Kws = 0.75*(1+2*lambda)/(1+lambda)*b*tand(LAMBDA)/L_fuse;

%% Kroo Weight calculations
W_wingK = 1.3*Kflap*Kspoil*Kgear_w*Keng*(4.22*S_ref...
    + 1.642e-6*Nult*b^3*sqrt(MTOW*MZFW)*(1+2*lambda)...
                             / (thick*cosd(LAMBDA)^2*S_ref*(1+lambda)));
                         
W_htailK = 5.25*S_wet_h + 8e-7*(Nult*b_h^3*MTOW*c*sqrt(S_wet_h)...
    / (thick_h*cosd(LAMBDA_h)^2*L_h*S_h^1.5));

W_vtailK = 2.62*S_v + 8e-7*(Nult*b_h^3*(8+0.44*MTOW/S_ref)...
    / (thick_v*cosd(LAMBDA_v)^2));

W_rudderK = 1.6*W_vtailK*S_rud/S_v;

Ip = 1.5e-3 * (P_delta*convlength(1,'ft', 'in')^2) * w_fuse;
Ib = 1.91e-4 * Nlim*(MZFW-W_wingK)*(L_fuse-c_root)/h_fuse^2;
if Ip > Ib
    Ifuse = Ip;
else
    Ifuse = (Ip^2+Ib^2)/2/Ib;
end
W_fuseK = (1.051 + 0.102*Ifuse)*Swet_fuse;

W_gearK = 0.04*MTOW;

W_flt_ctrlsK = Isc*(S_h + S_v);

% W_propK = 1.6*W_eng*nb_eng; %includes everything (pylon, nacelle, reversers, etc)
W_prop = nb_eng*0.7*W_eng^.736 + nb_eng*W_eng; %GASP pylon + engine

W_apuK = 7*Nseats;

W_instK = 800;   %100 for business jet, 800 for domestic transport, 1200 for long range / overwater

W_hydK = 0.65*S_ref;

W_electricalK = 13*Nseats;

W_electronicsK = 900;   %300 biz, 900 RJ, 1500 LR

W_furnK = (43.7 - 0.037*Nseats)*Nseats + 46*Nseats;

W_aircondK = 15*Nseats;

W_pilots = nb_pilot*W_1pilot;

W_atts =  + nb_att*W_att;

Weights_K = [0.85*W_wingK; 0.83*W_htailK; 0.83*W_vtailK; 0.83*W_rudderK; 0.90*W_fuseK; 0.95*W_gearK;...
    W_flt_ctrlsK; W_prop; W_apuK; W_instK; W_hydK; W_electricalK;...
    W_electronicsK; W_furnK; W_aircondK; W_ops; W_pilots; W_atts];

if N_pax == 77
    Weights_K([1:4, 6:9, 11]) = fixed_weights([1:4, 6:9, 11]);
end

OWE_K = sum(Weights_K);%� inclure dans la matrice Weights_K
MTOW_K = OWE_K + payload + W_fuel;

%% sortie KROO
TOW = MTOW_K;
LW = OWE_K+payload;

%% Raymer Transport Aircraft Weight calculations
W_wingR = 0.0051*(W_dg*Nult)^.557*S_ref^.649*AR^.5*thick_root^-.4...
    *(1+lambda)^.1*cosd(LAMBDA)^-1*S_aileron^.1;

W_htailR = 0.0379*Kuht*(1+w_fuse/b_h)^-.25*W_dg^.639*Nult^.1*S_h^.75...
    / L_h*(0.3*L_h)^.704/cosd(LAMBDA_h)*AR_h^.166*(1+S_ele)^.1;

W_vtailR = 0.0026*(1+HtHv)^.225*W_dg^.556*Nult^.536*L_v^-.5*S_v^.5...
    * L_v^.875/cosd(LAMBDA_v)*AR_v^.35*thick_v_root^-.5;

W_fuseR = 0.328*Kdoor*Klg*(W_dg*Nult)^.5*L_fuse^.25*Swet_fuse^.302...
    * (1+Kws)^.04*(L_fuse/d_fuse)^-.1;

W_MLGR = 0.0106*Kmp*MLW^.888*Nult^.25*L_mlg^.4*nb_mw^.321*nb_mss^-.5...
   * V_stall^.1;

W_NLGR = 0.032*Knp*MLW^.646*Nult^.2*L_nlg^.5*nb_nw^.45;

% W_nacelleR = 0.6724*Kng*N_Lt^.1*N_w^.294*Nult^.119*W_ec^.611...
%     * nb_eng^.984*Swet_nacelle^.224;    %includes air induction

% W_eng_ctrlsR = 5*nb_eng + 0.8*L_ec;

% W_starterR = 49.19*(nb_eng*W_eng/1000)^.541;

W_fuel_sysR = 2.405*V_t^.606/(1+V_i/V_t)*N_t^.5;
W_fuel_sys = 135;

W_flt_ctrlsR = 145.9*N_f^.554/(1+N_m/N_f)*S_cs^.2*(Iy*1e-6)^.06;

W_apu_iR = 2.2*W_apu_ui;    %W_apu_installed = 2.2*W_apu_uninstalled
W_apu_i = 300;
W_instrR = 4.509*Kr*Ktp*N_c^.541*nb_eng*(L_fuse+b)^.5;

W_hydR = 0.2673*N_f*(L_fuse+b)^.937;
W_hydT = 0.544*(45+1.318*(S_ref+1.44*S_h));

W_electrR = 7.291*R_kva^.782*L_a^.346*N_gen^.1;
W_electr = 3500;
W_avionicsR = 1.73*W_uav^.983;
W_avionics = 1600;
MWE = MTOW-W_fuel-225*N_pax-175*(max_pax-N_pax)-W_ops-W_pilots-W_atts;
W_avionicsT = 0.575*MWE^(5/9)*Parameters.Perf.d_mission^.25;

W_furnishingR = 0;%.0577*N_c^.1*W_c^.393*Swet_fuse^.75;    %does not include cargo & seats

W_ecsR = 62.36*(nb_pilot+nb_att+max_pax)^.25*(V_cabin/1000)^.604*W_instrR^.1;
W_ecs = 300;
W_anti_iceR = 0.002*W_dg;
W_de_ice = 100;
W_handling_gearR = 3e-4*W_dg;

%Unincluded weights

W_outfitting = nb_biz*kbiz + nb_econ*kecon...
    + nb_galley*kgalley + nb_lav*klav + nb_wardrobe*kwardrobe...
    + nb_trolley*ktrolley + nb_att*kseat_att;

% PAS DE NACELLE (OPEN ROTOR) � RAJOUTER POUR TURBOFAN %W_eng_ctrlsR;
% W_starterR;  aussi
% Weights_R = [0.85*W_wingR; 0.83*W_htailR; 0.83*W_vtailR; 0.9*W_fuseR; 0.95*W_MLGR; 0.95*W_NLGR;...
%     W_fuel_sysR; W_flt_ctrlsR;...
%     W_apu_iR; W_instrR; W_hydR; W_electrR; W_avionicsR; W_furnishingR;...
%     W_ecsR; W_anti_iceR; W_handling_gearR; W_prop; W_outfitting; W_ops; W_pilots; W_att];
% 
% if N_pax == 77
%     Weights_R([1:3, 5:9, 11]) = fixed_weightsR([1:3, 5:9, 11]);
% end

% OWE_R = sum(Weights_R);
% MTOW_R = OWE_R+payload+W_fuel;

% TOW = MTOW_R;
% LW = OWE_R+payload;

%% Mixed Weights (Kroo for struct, Raymer for systems, Torenbeek for hydraulics and avionics)
Weights = [0.85*W_wingK; 0.83*W_htailK; 0.83*W_vtailK; 0.83*W_rudderK;...
    0.88*W_fuseK; 0.95*W_gearK; W_flt_ctrlsR; W_prop; W_apu_i;...
    0*W_instrR; 0.33*W_hydT; W_electr; W_fuel_sys; W_avionics; W_furnishingR;...
    W_ecs; W_de_ice; W_outfitting; W_ops; W_pilots; W_atts];

%Fixed common weights for small A/C based on big A/C comp weights
if N_pax == 77
    Weights([1:4, 6:9, 11]) = fixed_weights([1:4, 6:9, 11]);
end

% if Parameters.Perf.d_mission == 500
    Weights(1:end-4) = fixed_weights([1:end-4]);
    Weights(end-2:end) = fixed_weights(end-2:end);
% end

OWE = sum(Weights);
MTOW = OWE+payload+W_fuel;
Parameters.MLW = OWE+payload+4500;
TOW = MTOW;
LW = OWE+payload;


%% Unused weight equations

% %% Raymer - General Aviation
% W_wingR2 = 0.036*S_ref^0.758*W_wf^.0035...
%     * (AR/cosd(LAMBDA)^2)^0.6*q_Df^0.006*lambda^0.04...
%     * (100*thick/cosd(LAMBDA))^-0.3*(Nult*(MZFW*MTOW)^.5)^0.49;
% 
% W_htailR2 = 0.016*(Nult*MTOW)^.414*q_Df^.168*S_h^.896...
%     * (100*thick/cosd(LAMBDA))^.3*(AR/cosd(LAMBDA_h)^2)^.043*lambda_h^-.2;
% 
% W_vtailR2 = 0.073*(1+0.2*HtHv)*(Nult*MTOW)^.376*q_Df^.122*S_v^.873...
%     * (100*thick/cosd(LAMBDA_v))^-.49*(AR/cosd(LAMBDA_v)^2)^.357...
%     * lambda_v^.039;
% 
% W_press = 11.9*(V_cabin*P_delta)^.271;  %weight penalty due to pressurization
% W_fuseR2 = 0.052*Swet_fuse^1.086*(Nult*MTOW)^.177*L_h^-.051...
%     * (L_fuse/d_fuse)^-.072*q_Df^.241 + W_press;
% 
% W_MLGR2 = 0.095*(Nult*MLW)^.768*(L_mlg/12)^.409;
% W_NLGR2 = 0.125*(Nult*MLW)^.566*(L_nlg/12)^.845;
% 
% W_eng_instR2 = 2.575*W_eng^.922*nb_eng;  %includes engine mount, prop...:/
% 
% W_fuel_sysR2 = 2.49*W_fuel^.726*(1+V_i/W_fuel)^.363...
%     * nb_tanks^.242*nb_eng^.157;
% 
% W_avionicsR2 = 2.117*W_uav^.933;
% 
% W_flt_ctrlsR2 = 0.053*L_fuse^1.536*b^.371*(Nult*MTOW*1e-4)^.8;
% 
% W_hydR2 = Kh*d_fuse^.8*M^.5;
% 
% W_electrR2 = 12.57*(W_fuel_sysR+W_avionicsR)^.51;
% 
% W_ecs_aiR2 = 0.265*MTOW^.52*nb_p^.68*W_avionicsR^.17*M^.08;
% 
% W_furnR2 = 0.0582*MTOW-65;
% 
% Weights_R2 = [W_wingR2; W_htailR2; W_vtailR2; W_fuseR2; W_MLGR2; W_NLGR2;...
%     W_eng_instR2; W_fuel_sysR2; W_avionicsR2; W_flt_ctrlsR2; W_hydR2;...
%     W_electrR2; W_ecs_aiR2; W_furnR2; W_ops];
% 
% OWE_R2 = sum(Weights_R2);
% MTOW_R2 = OWE_R2 + payload + W_fuel;

%% Torenbeek

% Torenbeek
% b_ref = []; %from SHEVELL - Fundamentals of Flight, Chap 18
% OMEGA_s = []; %from SHEVELL - Fundamentals of Flight, Chap 18
% W_wingT = 0.0013*Nult*(MZFW*MTOW)^.5*eta_cp*b/b_ref*AR/thick/cosd(LAMBDA)^2....
%     + OMEGA_s*S_ref;

% W_htailT = 1.1*S_h*(3.81*S_h^.2*V_D/(1000*sqrt(cosd(LAMBDA_h)))-0.287);

% W_fuse2 = 0.02625*Kf*sqrt(V_D*L_h/(d_fuse+d_fuse))*Swet_fuse^1.2;

% W_MLG2 = 1.8*Kgr*(33 + 0.04*MTOW^.75 + 0.021*MTOW);

% W_NLG2 = 1.8*(12 + 0.06*MTOW^.75);

% W_hyd = 0.544*(45+1.318*(S_ref+1.44*S_h));

%% OTHER
% 
% % Nacelle - semaine5 p49
% W_nacelle = 0.4*W_eng;

% % GASP
% W_fuel_sys = -2.2e-4*S_ref^2 + 1.059*S_ref + 26.335;
% 
% % APU - Roskam
% W_APU = 0.007*MTOW;

% % General Dynamics
% W_flt_ctrls = 56.08*(MTOW*q_D*1e-5)^.576*1.225;

% % General Dynamics
% W_electr = 1163*((W_fuel_sys+W_instr)/1000)^.506;

% % �quipements de finition
% W_green_furn = 0.0577*nb_pilot^.1*W_pilot^.393*Swet_fuse^.75;

% % Peinture - Roskam
% W_paint = 0.004*MTOW;

% W_outfitting = nb_1stclass*k1stclass + nb_biz*kbiz + nb_econ*kecon...
%     + nb_pax*(kbin+koxygen) + nb_galley*kgalley + nb_lav*klav...
%     + nb_wardrobe*kwardrobe + nb_trolley*ktrolley + nb_att*kseat_att...
%     + V_cabin*kfurn;
% 
