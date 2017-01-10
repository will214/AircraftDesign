clear all; clc; %close all;
global Parameters Aero COC Airplane
Parameters = struct();
Parameters.Weights = struct();
Parameters.MRO = struct();
Parameters.Conditions = struct();
Parameters.Aero = struct();
%% SIMULATION PARAMETERS
Parameters.Aircraft = 1; %1- 77 pax | 2- 101 pax
Parameters.Conditions.Cruise.Mach = 0.68; %% CHANGE EACH ITERATION
Parameters.Conditions.Cruise.Alt = 36000; %% CHANGE EACH ITERATION
Parameters.Perf.d_mission = 1800; %%% CHANGE EACH ITERATION
Parameters.excel_column = 'AG'; %%% CHANGE EACH ITERATION
Parameters.Perf.clbsched = [210 210 Parameters.Conditions.Cruise.Mach];
Parameters.HighDensity = 1; 
Parameters.Hardcoded_Wing = 0; %1 - Yes, 2- No
Parameters.Condtions.Cruise.delta_ISA = 0;
Parameters.Conditions.Flight=2; %1-Cruise 2-Landing 3-Takeoff
Parameters.Perf.clbschedvmd = [175 175 Parameters.Conditions.Cruise.Mach]; %do not mind, does not change anything
%% MRO PARAMETERS
Parameters.MRO.passengers=[77 101];
Parameters.MRO.mission = [1800 1200]; % No effect on code, only for payload range, post treatment;
Parameters.Perf.h_ptrans = 30680; % Transition altitude [ft]
%% PERFORMANCE
% Weight, steps
Parameters.Perf.MTOW = 80000;
Parameters.MLW = 72000;
Parameters.Perf.LW = 0.85*Parameters.Perf.MTOW;
Parameters.Perf.delta_hp = 500;
Parameters.Perf.Delta_T_cruise = 60*5;
% Engine settings
Parameters.Perf.T_SLS = 13000;
Parameters.Perf.T_BASELINE = 15000;
Parameters.Perf.SF = Parameters.Perf.T_SLS/Parameters.Perf.T_BASELINE;
Parameters.Weights.Engine = 5850*Parameters.Perf.SF^0.75; 
Parameters.Engine.NacelleDiameter = 47.6*sqrt(Parameters.Perf.SF); % [ft]
Parameters.Engine.FWDDiameter = 10.2*sqrt(Parameters.Perf.SF); % [ft]
Parameters.Engine.AFTDiameter = 9.2*sqrt(Parameters.Perf.SF); % [ft]
% Oruise settings
Parameters.Perf.ICA = Parameters.Conditions.Cruise.Alt; % Initial Cruise Altitude [ft]
Parameters.Perf.M_target = Parameters.Conditions.Cruise.Mach; % Target Cruise Mach

% Takeoff settings
Parameters.Perf.h_pito = 0;
Parameters.Perf.h_pfto = 1500;
Parameters.Perf.fuel_to = 150; % fixed [lb] 
Parameters.Perf.t_to = 2*60; % fixed [s]

% Landing settings
%Parameters.Perf.fuel_ldg = 150; % fixed [lb] 
Parameters.Perf.t_ldg = 3*60; % fixed [s]

% Taxi settings
Parameters.Perf.t_taxi_in = 10*60; %[s]
Parameters.Perf.fuel_taxi_in = 148; % [lb]
Parameters.Perf.t_taxi_out = 5*60; %[s] 
Parameters.Perf.fuel_taxi_out = 90; % [b]

% Climb/Descent settings
Parameters.Perf.h_pi = 1500; % Initial altitude [ft]
Parameters.Perf.h_pf = Parameters.Perf.ICA; % Final altitude [ft]
% Parameters.Perf.clbsched = [160 200 Parameters.Perf.M_target]; % Speed schedule in VCAS [kts], low speed schedule (module 5, p.15)
Parameters.Perf.disa = 0; % (module12, p.6)
Parameters.Perf.V_wind_ft = 0; % no crosswind

% Fuel Reserve settings
Parameters.Perf.d_mission_fr = 100; % [nm]
Parameters.Perf.ICA_fr = 15000; %[ft]
Parameters.Perf.h_pi_fr = 0; %[ft] 
Parameters.Perf.h_pf_fr = Parameters.Perf.ICA_fr; 
%% WEIGHTS PARAMETERS
% Parameters.Weights.nb_eng;
% Parameters.Weights.nb_pilot;
% Parameters.Weights.nb_att;
% Parameters.Weights.nb_pax;
% Parameters.Weights.nb_1stclass;
% Parameters.Weights.nb_biz;
% Parameters.Weights.nb_econ;
% Parameters.Weights.nb_galley;
% Parameters.Weights.nb_lav;
% Parameters.Weights.nb_wardrobe;
% Parameters.Weights.nb_trolley;
% Parameters.Weights.Vmo; %(V_mo+30) kts
% Parameters.Weights.Kflap;
% Parameters.Weights.Kspoil;
% Parameters.Weights.Kgear_w;
% Parameters.Weights.Keng;
% Parameters.Weights.Nult; %3.75;
% Parameters.Weights.Kgear;
% Parameters.Weights.Kfloor=0.7;
% Parameters.Weights.Kgr;
% Parameters.Weights.W_pilot; %240;  %lb
% Parameters.Weights.k1stclass;%80; %lb/siege
% Parameters.Weights.kbiz;%58;  %lb/siege
% Parameters.Weights.kecon;%23; %lb/siege
% Parameters.Weights.kbin;%10;  %lb/pax
% Parameters.Weights.koxygen;%1;    %lb/pax
% Parameters.Weights.kgalley;%500;  %lb/galley.Weights
% Parameters.Weights.klav;%%lb/lav
% Parameters.Weights.kwardrobe;%150;    %lb/wr
% Parameters.Weights.ktrolley;%90;  %lb - 60lb half, 90lb full
% Parameters.Weights.kseat_att;%25; %lb/att
% Parameters.Weights.kfurn;%0.5;    %lb/ft^3%% AERO PARAMTERS
%% AERO PARAMTERS
Parameters.Aero.Trans = struct();
Parameters.Aero.drag = struct();
Parameters.Aero.airfoil = struct();
Parameters.Aero.Wing = struct();
Parameters.Aero.Flap = struct();
Parameters.Aero.Slat = struct();
Parameters.Aero.H_Tail = struct();
% DRAG
Parameters.Aero.Trans.wing      = 0.1;      
Parameters.Aero.Trans.engine    = 0;      
Parameters.Aero.Trans.pylon     = 0;      
Parameters.Aero.Trans.V_Tail    = 0.1;      
Parameters.Aero.Trans.H_Tail    = 0.1;      
Parameters.Aero.Trans.fuse      = 0;        
Parameters.Aero.Trans.Belly     = 0;
Parameters.Aero.drag.Fa2=1.0; %0.5 for area ruling fuse-nacelle interference
Parameters.Aero.drag.CDn=0.11; %Fuselage-Nacelle interference
Parameters.Aero.drag.SHP=0.1;
Parameters.Aero.drag.windshield = 0.011; %(See Roskam page 1792)
Parameters.Aero.Nacelle_swet = 19.354; %(m^2)
Parameters.Aero.Pylon_swet = 4.646; %(m^2)
Parameters.Aero.elevateur_drag = 0.008;
% AIRFOIL
Parameters.Aero.Mach=0.7;
Parameters.Aero.airfoil.ao=-3.5;
Parameters.Aero.airfoil.a_clmax = 12.25;
Parameters.Aero.airfoil.clmax = 1.5;
Parameters.Aero.airfoil.cla=1.8*pi;
Parameters.Aero.airfoil.LER = 0.3;
Parameters.Aero.airfoil.cmo = -0.1078;
% WING
Parameters.Aero.Wing.twist=0;
Parameters.Aero.Wing.delta_ao_twist=0;
Parameters.Aero.Wing.df=5; %ft
Parameters.Aero.Wing.alpha=0;
Parameters.Aero.Wing.a_CLmax=14*pi/180;
Aero.Airplane.CLcruise = 0.5;
Aero.Drag.Wing.CDo = 0.005;
Parameters.Aero.Wing_CL_Ctmax = 0.9;
Parameters.Aero.Wing_CLmax_corr = 0.05;
% FLAP
Parameters.Aero.Flap.cf_c=[0.24 0.24];
Parameters.Aero.Flap.deflection=[20 42];
Parameters.Aero.Flap.type='single-slotted';
Parameters.Aero.Flap.span_start=6; %ft
Parameters.Aero.Flap.span_end=36;   %ft
Parameters.Aero.Flap.s_wf_s=0.6;
Parameters.Aero.Flap.bf_b=0.7;
Parameters.Aero.Flap.cs_c=[1.15 1.45];
Parameters.Aero.Flap.cp_c=[1.15 1.45];
Parameters.Aero.Flap.Swf=400;
% SLAT
Parameters.Aero.Slat.cf_c=[0.15 0.15];
Parameters.Aero.Slat.deflection=[20 25];
Parameters.Aero.Slat.cp_c=[1.1 1.4];
Parameters.Aero.Slat.bp_b=0.8;
Parameters.Aero.Slat.cs_c=[1.1 1.4];
Parameters.Aero.Slat.Swf=0.7;
Parameters.Aero.Slat.type = 'fowler';


% H-TAIL
Parameters.Aero.H_Tail.cl_a=1.5*pi;
Parameters.Aero.H_Tail.eoh=2;
Parameters.Aero.H_Tail.delta_CLh=0.5;
Parameters.Aero.H_Tail.flap_SF_Sref = 0.95;
Parameters.Aero.H_Tail.volet_cf_c = 0.25;
% CM
Parameters.Aero.Airplane.CG=22;
Parameters.Aero.trim_init = 22;
Parameters.Aero.H_Tail.eta_h=1;
Parameters.Aero.XLE_corr = 0;

%% STARTING PROGRAM
readinput();
% cd aero/Lift
% lift_coefficient(1);
% cd ../..


if Parameters.Aircraft == 2
    load ./weight/weights101 Weights W_fuel fuel_reserve MTOW OWE max_pax MZFW
else
    load ./weight/weights77 Weights W_fuel fuel_reserve MTOW OWE max_pax MZFW
end
% load ./weight/weights101
fixed_weights = Weights;
%% Main Code
TOW = Parameters.Perf.MTOW;
LW = Parameters.Perf.LW;
MTOW = Parameters.Perf.MTOW-1000;
% options = optimset('FunValCheck', 'on', 'MaxFunEvals', 7, 'PlotFcns', @optimplotfval, 'TolX', 1);%options for fminsearch for fuel reserve
options = optimset('FunValCheck', 'on', 'MaxFunEvals', 7, 'TolX', 15);
trapped_fuel = 1.02;    % 2% trapped fuel
OWE = 50000;
while abs(MTOW - TOW) > 100
    TOW = MTOW;
    cd perf
    W_fuel = mp(TOW, LW);
    [a,b,c] = fminsearch(@(x) mp_reserve(x,OWE), [500, 240, 50], options);
    fuel_reserve = a(1)*1.05;
    lrc_alt = a(2)*100;
    lrc_mach = a(3)/100;
    [a,b,c] = fminsearch(@(x) mp_vmd(x,OWE,fuel_reserve),[28], options);
    mach_vmd = a(1);
    fuel_vmd = COC.VMD_fuel*1.05;
    fuel_tot = (W_fuel + fuel_reserve + fuel_vmd)*trapped_fuel;
    cd ../weight/
    [MTOW, LW, OWE, Weights] = weights(Parameters.MRO.passengers(Parameters.Aircraft), TOW, fuel_tot, Parameters.Weights.Engine, Parameters.Conditions.Cruise.Mach, fixed_weights);
    fprintf('\nMTOW=%.f LW = %.f W_fuel = %.f\n',MTOW,LW,W_fuel)
    cd ../
end

fprintf('LRC_time = %.f LRC_SAR=%.3f LRC_FUEL = %.f LRC_MACH = %.2f LRC_ALT = %.f\n',COC.LRC_time*60, COC.LRC_SAR,COC.LRC_fuel,COC.LRC_mach,COC.LRC_alt)
fprintf('VMD_MACH = %.2f VMD_fuel = %.f VMD_drag = %.f VMD_d = %.f\n',mach_vmd,fuel_vmd,COC.VMD_drag,COC.VMD_d)

% % Post-traitement : Balance
% filename = strcat('output_dt_', num2str(Parameters.MRO.passengers(Parameters.Aircraft)), 'pax_V1.xlsx');
% cd weight
% [cg, tow, max_pax, MZFW] = pass_loading(Weights, W_fuel, filename);
% cd ..

% if Parameters.Aircraft == 1
%     save ./weight/weights77hd Weights W_fuel fuel_reserve MTOW OWE max_pax
% else
%     save ./weight/weights101hd Weights W_fuel fuel_reserve MTOW OWE max_pax LW
% end

%Post-traitement : Stability
% if Parameters.Aircraft==1;
%     COC.W_TOC=7.1784e+04;
%     COC.W_EOC=6.5051e+04;
% elseif Parameters.Aircraft==2;
%     COC.W_TOC=8.1811e+04;
%     COC.W_EOC=7.7415e+04;
% end

% cd perf
% xcg_weight = 0.27; %101 = 22 77=27
% xnp_weight = 0.44;
% s_vt_real = 152; %ft2, ESDU mismatch
% [cm_alpha cn_b cn_b_vt cn_b_w cn_b_f l_vt s_vt x_np x_aft static_margin] = stability(Parameters.Conditions.Cruise.Alt,Parameters.Conditions.Cruise.Mach, COC.W_TOC, COC.W_EOC,xcg_weight,xnp_weight,s_vt_real);
% xp2(Parameters.Conditions.Cruise.Alt,Parameters.Conditions.Cruise.Mach, COC.W_TOC, COC.W_EOC,xcg_weight,xnp_weight,s_vt_real);
% fprintf('cm_alpha : %.2f cn_b : %.2f cn_b_vt : %.2f cn_b_w : %.2f cn_b_f : %.2f l_vt : %.f ft s_vt : %.f ft2 x_np : %0.2f x_aft = %0.2f static margin 0.05-0.1 : %.3f\n',cm_alpha,cn_b,cn_b_vt,cn_b_w,cn_b_f, l_vt, s_vt, x_np, x_aft,static_margin)
% cd ..

% % Post-traitement : BFL-TOW
% TOW = COC.TOW;
% cd perf
% [W_MAT, BFL_MAT, BFL25_MAT]=bfltow(TOW);
% figure(6)
% plot(W_MAT,BFL_MAT,'m',W_MAT,BFL25_MAT,'r')
% grid on
% xlim([min(W_MAT) max(W_MAT)]);
% line([TOW TOW], ylim)
% legend('BFL', 'BFL 5000ft ISA+25','Design TOW','location','best')
% xlabel('TOW (lbs)')
% ylabel('BFL (ft)')
% title('BFL vs TOW chart (77 PAX)')
% cd ..


% %% Post-traitement : Payload-Range
% if Parameters.Aircraft == 2
%     W_fuel_total = 6629;
%     W_fuel_reserve = 703+1300;
% else
%     W_fuel_total = 8766;
%     W_fuel_reserve = 496+1409;
% end
% W_fuel_ldg = COC.W_fuel_ldg;
% W_fuel_taxi_in = COC.W_fuel_taxi_in;
% cd perf
% [Npax_MAT, Range_MAT] = payrangef(Parameters.MRO.passengers(Parameters.Aircraft), MZFW, OWE,MTOW,LW,W_fuel,W_fuel_ldg,W_fuel_taxi_in, max_pax,W_fuel_reserve);
% if Parameters.Aircraft == 2
%     Range_MAT(1) = 1200;
% else
%     Range_MAT(1) = 1800;
% end
% figure(9)
% plot([0, Range_MAT],[Parameters.MRO.passengers(Parameters.Aircraft)*225,Npax_MAT])
% grid on
% xlabel('Range (nm)')
% ylabel('Payload (lbs)')
% if Parameters.Aircraft == 2
%     title('Payload Range Chart (101 PAX)')
% else
%     title('Payload Range Chart (77 PAX)')
% end
% 
% cd ..

% % Post-traitement SC
% cd aero/Lift
% [~,CG,SH, ~]=lift_coefficient_sc('takeoff',-7,-5);
% [~,CG2,SH2, ~]=lift_coefficient_sc('takeoff',6,5);
% [~,~,~, NP]=lift_coefficient_sc('cruise',0,0);
% 
% [~, ind] = min(abs(SH+50));
% CG1 = CG(1:ind);
% SH1 = SH(1:ind);
% [~, ind] = min(abs(SH2+50));
% CG2 = CG2(ind:end);
% SH2 = SH2(ind:end);
% 
% cd ../..
% close all
% figure(10)
% s_ht = Airplane.ESDU.H_Tail.Area*convlength(1,'m','ft')^2;
% hold on
% plot(CG1/100,SH1,'r')
% plot(CG2/100,SH2,'g')
% plot(NP/100,0:300, 'm')
% grid on
% 
% title('Longitudinal X-Plot (101 PAX)')
% xlabel('X_C_G/c')
% ylabel('Horizontal Tail Area (ft2)')
% ylim([0 300])
% line(xlim, [s_ht s_ht])
% legend('Takeoff Nose Up','Takeoff Nose Down','Neutral Point','Current H-Tail','location','best')