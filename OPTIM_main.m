%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File : OPTIM_main.m
%
% Function : This program in set to be run, load parameters from CATALIST
%   aircraft geometry and MR&O specifications, used to iterate on design
%   characteritics, converge toward viable design and compare them
%   together in terms of operation costs, fuel consumption, travel time and
%   such variables.
%
% No Inputs, No Outputs
%
% Authors : Tony LAC
%           Martin LAFRANCE
%           Matthieu PARENTEAU
%
% Date : January 15th, 2015
% Last Edit : March 10th, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;
global Parameters

%% Parameters to change for each different iteration
Parameters.Aircraft = 2; %1- 77 pax | 2- 101 pax
Parameters.Conditions.Cruise.Mach = 0.74; %% Cruise Speed
Parameters.Conditions.Cruise.Alt = 34000; %% Cruise Altitude

Parameters.excel_column = 'AB'; % Excel File Column to Write In
Parameters.Condtions.Cruise.delta_ISA = 0; % Atmosphere Condition

Parameters.Perf.T_SLS = 13000;

% Auto set Mission Range according to requirements for each aircraft
if Parameters.Aircraft == 1
    Parameters.Perf.d_mission = 1800;
else
    Parameters.Perf.d_mission = 1200;
end

%% Optimization
OPTIM_param();
options1 = optimset('MaxIter', 2);%options for fminsearch for fuel reserve
mini = fminsearch(@(x) OPTIM_loop(x), [140, 200, 0.7, 30000], options1);

Parameters.Conditions.Cruise.Mach = mini(3);
Parameters.Perf.M_target = mini(3);
Parameters.Perf.clbsched = mini(1:3);
Parameters.Conditions.Cruise.Alt = mini(4);

%%
TOW = Parameters.Optim.MTOW;
LW = Parameters.Optim.LW;
OWE = Parameters.Optim.OWE;
trapped_fuel = 1.02; % 2% trapped fuel in wing
a = Parameters.Optim.a;
load ./weight/weights100
fixed_weights = Weights;

cd perf
W_fuel = mp(TOW, LW);
% [a,b,c] = fminsearch(@(x) mp_reserve(x,OWE), [500, 100, 50], options);
%why run fminsearch when it is already done in OPTIM-loop (TONY)
fuel_reserve = a(1);
fprintf('Fuel reserve : %.0f\nGo-Around Altitude : %.0f\nGo-Around Mach : %.2f\n', a(1), a(2)*100, a(3)/100);
fuel_tot = (W_fuel + fuel_reserve)*trapped_fuel;
cd ../weight/
[MTOW, LW, OWE, Weights] = weights(Parameters.MRO.passengers(Parameters.Aircraft), TOW, fuel_tot, Parameters.Weights.Engine, Parameters.Conditions.Cruise.Mach, fixed_weights);
fprintf('\nMTOW=%.f LW = %.f W_fuel = %.f\n',MTOW,LW,W_fuel)
cd ../

fprintf('\nClimb Schedule=%.f %.f %.2f Altitude = %.f\n',Parameters.Perf.clbsched(1),Parameters.Perf.clbsched(2),Parameters.Perf.clbsched(3),Parameters.Conditions.Cruise.Alt)

% %% Post Teatment
% 
% % Balance (PAX loading weight scheme to file and figure)
% if Parameters.Aircraft == 2 % save common weights from 100pax for 76pax
%     save weights100 Weights Weights_R
% end
% filename = strcat('output_dt_', num2str(Parameters.MRO.passengers(Parameters.Aircraft)), 'pax_V1.xlsx');
% cd weight
% [cg, tow] = pass_loading(Weights, W_fuel, filename);
% figure(2)
% plot(cg.c(4:end), tow.c(4:end));
% hold on
% plot(cg.d(4:end), tow.d(4:end));
% legend('77PAX M = 0.7 ICA = 29 000 ft')
% cd ..

% % BFL-TOW
% TOW = 82795;
% cd perf
% [W_MAT, BFL_MAT, BFL25_MAT]=bfltow(TOW);
% plot(W_MAT,BFL_MAT,'b',W_MAT,BFL25_MAT,'r')
% legend('BFL', 'BFL 5000ft ISA+25','best')
% xlabel('TOW (lbs)')
% ylabel('BFL (ft)')
% title('BFL vs TOW chart (101 PAX)')
% cd ..

% %% Post-traitement : Payload-Range
% TOW = 82795;
% LW = 74161;
% W_fuel_total = 7977*1.06*1.02;
% W_fuel_ldg = 24;
% W_fuel_taxi_in = 90;
% OWE = LW - Parameters.MRO.passengers(Parameters.Aircraft)*225;
% MZFW = OWE + (Parameters.MRO.passengers(Parameters.Aircraft)+20)*225;
% cd perf
% [Npax_MAT2, Range_MAT2] = payrange2(Parameters.MRO.passengers(Parameters.Aircraft), MZFW, OWE,TOW,LW,W_fuel_total,W_fuel_ldg,W_fuel_taxi_in);
% [Npax_MAT3, Range_MAT3] = payrange3(Parameters.MRO.passengers(Parameters.Aircraft), MZFW, OWE,TOW,LW,W_fuel_total,W_fuel_ldg,W_fuel_taxi_in);
% Npax_MAT1 = linspace((Parameters.MRO.passengers(Parameters.Aircraft)+20)*225,(Parameters.MRO.passengers(Parameters.Aircraft)+20)*225,25);
% Range_MAT1 = linspace(0,Range_MAT2(1),25);
% plot([Range_MAT1, Range_MAT2, Range_MAT3],[Npax_MAT1,Npax_MAT2,Npax_MAT3])
% xlabel('Range (nm)')
% ylabel('Payload (lbs)')
% title('Payload Range Chart (101 PAX)')
% cd ..