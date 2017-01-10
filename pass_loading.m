function [CG, TOW, max_pax, MZFW] = pass_loading(Weights, W_fuel, filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File : pass_loading.m
%
% Function : Based on aircraft geometry and component weights, this
%   function evaluates the CG position of the aircraft at all time (empty,
%   during loading scenarios, full) and writes the result in an MS Excel
%   file for figure generation.
%
% Inputs : (1) Weights matrix with 21 elements : Wing, HTail, VStab,
%               Rudder, Fuselage, Landing Gear, Flight Controls, 
%               Propulsion, APU, Instriments, Hydrolics, Electrical 
%               Systems, Fuel Systems, Avionics, Furnishing, ECS, De-Ice, 
%               Outfitting, Operational Weights, Pilots, Flight Attendants
%          (2) Fuel Weight (total mission fuel)
%          (3) Excel Filename for post-analysis
%
% Outputs : (1) CG position matrix for each load step in chord percentage
%           (2) Weight value for each load step in lb
%           (3) High density configuration number of passengers
%           (4) Maximum Zero Fuel Weight
%
% Author : Martin LAFRANCE
%
% Date : January 15th, 2015
% Last Edit : April 10th, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W_fuel = W_fuel+500; % fuel account for taxi out (MR&O 500lb diff)
OWE = sum(Weights);
global output COC
%%
% Cargo ratio fwd/aft
fwd_ratio = 0.6;
aft_ratio = 0.4;
% Cargo position (lower), to change
%% suppression of added output values
output_ = vertcat(output(1:429), output(434:441), output(446:end));

%% Fuel Tank geometry
% Wing chord at root, BL1, BL2, tip
Wchord0 = convlength(output_(59), 'in', 'ft');
Wchord1 = convlength(output_(60), 'in', 'ft');
Wchord2 = convlength(output_(61), 'in', 'ft');
Wchord3 = convlength(output_(62), 'in', 'ft');
% Half-span at BL1, BL2, tip
WBL1 = convlength(output_(63), 'in', 'ft');
WBL2 = convlength(output_(64), 'in', 'ft');
WBL3 = convlength(output_(65), 'in', 'ft');
% Leading edge longitudinal position at root, BL1, BL2, tip, mac
X_LE_0 = convlength(output_(54), 'in', 'ft');
X_LE_1 = WBL1*tand(output_(56)) + X_LE_0;
X_LE_2 = (WBL2-WBL1)*tand(output_(57)) + X_LE_1;
X_LE_3 = (WBL3-WBL2)*tand(output_(58)) + X_LE_2;
X_LE_MAC = convlength(output_(323)-0.25*output_(326), 'in', 'ft');
% Half tank span at center tank junction and whole
ctrtank_span = convlength(output_(92)*output_(15), 'in', 'ft')/2;
tank_span = output_(91)*WBL3;
% Longitudinal position of mid tank point at root, BL1, BL2, tip
X_midtank_0 = X_LE_0 + (output_(83)+output_(87))*Wchord0/2;
X_midtank_1 = X_LE_1 + (output_(84)+output_(88))*Wchord0/2;
X_midtank_2 = X_LE_2 + (output_(85)+output_(89))*Wchord0/2;
X_midtank_3 = X_LE_3 + (output_(86)+output_(90))*Wchord0/2;
% Tank chord at root, BL1, BL2, tip
Tchord0 = (output_(87)-output_(83))*Wchord0;
Tchord1 = (output_(88)-output_(84))*Wchord1;
Tchord2 = (output_(89)-output_(85))*Wchord2;
Tchord3 = (output_(90)-output_(86))*Wchord3;
% Tank thickness at root, BL1, BL2, tip
Tthick0 = output_(66)*Wchord0;
Tthick1 = output_(67)*Wchord1;
Tthick2 = output_(68)*Wchord2;
Tthick3 = output_(69)*Wchord3;
% % Tank section area at root, BL1, BL2
% Tsarea0 = 2.09;
% Tsarea1 = 0.411;
% Tsarea2 = 0.148;
% Volume of sections
Tvol1 = WBL1/6*(Tchord0*Tthick0+(Tchord0+Tchord1)*(Tthick0+Tthick1)+Tchord1*Tthick1);
Tvol2 = (WBL2-WBL1)/6*(Tchord1*Tthick1+(Tchord1+Tchord2)*(Tthick1+Tthick2)+Tchord2*Tthick2);
% Fuel density
rho_fuel = convmass(810,'kg','lbm')/convlength(1,'m','ft')^3; % kg/m^3->lb/ft^3 - (approx) Jet A p22 http://dtic.mil/dtic/tr/fulltext/u2/a132106.pdf
fill_factor = 0.98;      %expansion space in tank
W_fuel_wing_max = output_(353)*convlength(1,'m', 'ft')^3*fill_factor*rho_fuel;
tank_scaling_factor = (Tvol1+Tvol2)/output_(354)/convlength(1,'m', 'ft')^3;


%% Inputs
MAC = convlength(output_(326), 'in', 'ft');
H_MAC = convlength(output_(321), 'in', 'ft');
X_HLE = convlength(output_(318)-0.25*H_MAC, 'in', 'ft');
V_MAC = convlength(output_(313), 'in', 'ft');
X_VLE = convlength(output_(310)-0.25*V_MAC, 'in', 'ft');
X_NLG = convlength(output_(840), 'in', 'ft');
X_MLG = convlength(output_(847), 'in', 'ft');
X_eng = convlength(output_(144)+output_(146)/2, 'in', 'ft');
X_APU = convlength(output_(94), 'in', 'ft'); %SKETCHY
l_nose = convlength(output_(30), 'in', 'ft');
l_fuse = convlength(sum(output_([26 27 30 34])), 'in', 'ft');
mid_cabin = convlength(sum(output_([26 27])), 'in', 'ft')/2 + l_nose;


X_wing = X_LE_MAC+0.3*MAC;
X_htail = X_HLE+0.3*H_MAC;
X_vtail = X_VLE+0.3*V_MAC;
X_fuse = 0.4*l_fuse;
X_flt_ctrls = X_LE_MAC+0.6*MAC;
X_hyd = X_LE_MAC+0.65*MAC;
X_instr = l_nose;
X_electr = 0.4*l_fuse;
X_fuel_sys = X_LE_MAC + 0.5*MAC;
X_ecs = output(54)/12;%0.4*l_fuse;
X_outfitting = 0.45*l_fuse;
X_anti_ice = X_LE_MAC + 0.2*MAC;
X_wops = mid_cabin;
X_gearK = 0.22*X_NLG + 0.78*X_MLG;
X_pilot = convlength(output_(833), 'in', 'ft');
X_atts = convlength(output_(816)/3 + 2*(output_(792)+30)/3, 'in', 'ft'); %1 flt att en avant, deux en arri�re
X_instK = X_pilot;

xlswrite(filename, output_(326), 'Sheet1', 'C4:C4');
xlswrite(filename, convlength(X_LE_MAC, 'ft', 'in'), 'Sheet1', 'C5:C5');
tipover = (X_MLG-convlength(output(863)+output(856), 'in', 'ft')*tand(13) - X_LE_MAC)/MAC*100;
xlswrite(filename, [tipover; tipover], 'Sheet1', 'U53:U54');
xlswrite(filename, [COC.cgneutre; COC.cgneutre], 'Sheet1', 'U55:U56');


%%
Pos = [X_wing; X_htail; X_vtail; X_vtail; X_fuse; X_gearK; X_flt_ctrls;...
    X_eng; X_APU; X_instK; X_hyd; X_electr; X_fuel_sys; X_instr;...
    X_outfitting; X_ecs; X_anti_ice; X_outfitting; X_wops; X_pilot; X_atts];

%% Kroo
xle35htail = output(94) + 0.35*output(99)*sind(output(95));
c35htail = output(97) + 0.35*(output(98)-output(97));
X_htailK = convlength(xle35htail + 0.3*c35htail, 'in', 'ft');
h35vtail = output(107) - output(110) - 0.65*output(316);
xle35vtail = output(104) + output(110)*sind(108) + h35vtail*sind(109);
c35vtail = output(317)*(1-0.35*output(342));
X_vtailK = convlength(xle35htail + 0.3*c35vtail, 'in', 'ft');
X_fuseK = 0.45*l_fuse;
X_flt_ctrlsK = X_LE_MAC+0.4*MAC;
X_instK = 0.4*l_nose;
X_hydK = 0.75*X_wing + 0.25*X_htail;
X_electrK = 0.75*X_fuseK + 0.25*X_eng;
X_ecsK = l_nose;
X_anti_iceK = l_nose;

PosK = [X_wing; X_htailK; X_vtailK; X_vtailK; X_fuseK; X_gearK;...
    X_flt_ctrlsK; X_eng; X_APU; X_instK; X_hydK; X_electrK; X_fuel_sys;...
    X_instK; X_outfitting; X_ecsK; X_anti_iceK; X_outfitting; X_wops; X_pilot; X_atts];

%% Passenger distribution (do once)
% Pour chaque section (colonne), les �l�ments sont
% ligne 1) position longitudinale du premier si�ge (in)
% ligne 2) nombre de si�ges abreast de la sous-section
% ligne 3) nombre de rang�es de la sous-section
% ligne 4) seat pitch (in)
sections(:,1) = output_(744:747);    %biz left
sections(:,2) = output_(750:753);    %econ left front
sections(:,3) = output_(756:759);    %econ left back
sections(:,4) = output_(763:766);    %biz right
sections(:,5) = output_(769:772);    %econ right front
sections(:,6) = output_(775:778);    %econ right back

% N_pax (utile seulement pour calculer W_cargo)
N_pax = sum(sections(2,:).*sections(3,:));
W_cargo = 30*N_pax;

% Formattage des positions longitudinales et seat pitch en pieds
sections([1 4], :) = convlength(sections([1 4], :), 'in', 'ft');

% Initialisation des matrices avec la bonne taille
if output(768) && output(749)
    biz_abreast = sections(2,1)+sections(2,4);
    biz_rows = max(sum(sections(3,1)),sum(sections(3,4)));
    biz = NaN(biz_rows, biz_abreast);
else
    biz = [];
    sections(:,[1 4]) = zeros(4,2);
end

if output(756) && output(762) && output(775) && output(781)
    econ_abreast = max(sections(2,2:3)+sections(2,5:6));
    econ_rows = max(sum(sections(3,2:3)),sum(sections(3,5:6)));
    econ = NaN(econ_rows, econ_abreast);
elseif output(756) && output(775) && ~output(762) && ~output(781)
    econ_abreast = sections(2,2)+sections(2,5);
    econ_rows = max(sections(3,2), sections(3,5));
    econ = NaN(econ_rows, econ_abreast);
    sections(:,[3 6]) = zeros(4,2);
end

% biz loop
lat = 1;
for i=[1 4] %colonnes biz de la matrice sections
        lon = 1;
    for ii = 1:sections(3,i)	% pour toutes les rang�es de la sous-section
        for iii = 1:sections(2,i)   % pour tous les si�ge de la rang�e (abreast)
            biz(lon,lat) = sections(1,i) + (ii-1)*sections(4,i);% position longitudinale du premier si�ge de la rang�e + nombre de rang�es * seatpitch
            lat = lat + 1;  %on se tasse lat�ralement
            if i < 4 && lat > sections(2,i) && lon~=sections(3,i) % si on est dans la section de gauche mais avec une trop grande position abreast
                lat = 1;    % on r�initialise � la gauche
            elseif lat > sections(2,1)+sections(2,i)% m�me condition mais � droite
                lat = sections(2,1)+1;   % on r�initialise au premier si�ge de la section de droite
            end
        end
        lon = lon + 1;  %on recule d'un si�ge
    end
end

% econ loop (idem)
lat = 1;
lon = 1;
for i=[2 3 5 6]
    for ii = 1:sections(3,i)
        for iii = 1:sections(2,i)
            econ(lon,lat) = sections(1,i) + (ii-1)*sections(4,i);
            lat = lat + 1;
            if i < 4 && lat > sections(2,i)
                lat = 1;
            elseif lat - max(sections(2,2:3)) > sections(2,i)
                lat = lat-sections(2,i);
            end
        end
        lon = lon + 1;
    end
    
    %repositionnement du "premier si�ge" en arrivant � la fin de la
    %pr�sente section
    switch i    %apr�s la section i
        case 2  %fin de la premi�re section econ
            lon = sections(3,i)+1;  %on est derri�re la premi�re section
            lat = 1;    %on est toujours � gauche
        case 3  %fin de la deuxi�me section econ
            lon = 1;    %on retourne � l'avant
            lat = max(sections(2,2:3)) + 1; %� droite des sections de gauche
        case 5  %fin de la troisi�me section econ
            lon = sections(3,i)+1;  %derri�re la 3e section
            lat = max(sections(2,2:3)) + 1; %� droite des sections de gauche
    end
end

%une matrice avec tous les si�ges dedans (pour load case tous ensemble)
seats = NaN(size(biz, 1)+size(econ, 1), size(econ,2));
seats(1:size(biz,1),1:ceil(size(biz,2)/2)) = biz(:,1:ceil(size(biz,2)/2));
seats(1:size(biz,1), end-floor(size(biz,2)/2)+1:end) = biz(:,1:ceil(size(biz,2)/2)+1:end);
seats(size(biz,1)+1:end, :) = econ;
X_outfitting = mean(seats(:,1));
PosK(15) = X_outfitting;
Pos = PosK;
%% CARGO POSITIONING
X_fwd_cargo = biz(end,1);
X_aft_cargo = econ(ceil(length(econ)/2)+6,1);

% %% load case 1 - back to front, biz then econ
% 
% % temp matrices to be able to reuse originals for different load cases
% biz_ = biz;
% econ_ = econ;
% Weights_ = Weights;
% Pos_ = Pos;
% 
% %biz loop first
% for i = 1:ceil(size(biz_,2)/2)  %deux "colonnes" � la fois (ext�rieures)
%     for ii = size(biz_,1):-1:1  %back to front
%         if ~isnan(biz_(ii,1)) && ~isnan(biz_(ii,end))   %si les deux si�ges existent
%             Pos_ = [Pos_; biz_(ii, 1); biz_(ii, end)];  %on ajoute les deux
%             Weights_ = [Weights_; 195; 195];
%         elseif ~isnan(biz_(ii,1))       %si seulement le premier existe
%             Pos_ = [Pos_; biz_(ii, 1)]; %on n'ajoute que le premier
%             Weights_ = [Weights_; 195];
%         elseif ~isnan(biz_(ii,end))         %si seulement le dernier existe
%             Pos_ = [Pos_; biz_(ii, end)];   %on n'ajoute que le dernier
%             Weights_ = [Weights_; 195];
%         end
%     end
%     biz_ = biz_(:,2:end-1); %une fois les colonnes compl�tement remplies on les �limine de la matrice
% end
% 
% %econ loop second
% for i = 1:ceil(size(econ_,2)/2)
%     for ii = size(econ_,1):-1:1
%         if size(econ_,2) > 1
%             if ~isnan(econ_(ii,1)) && ~isnan(econ_(ii,end))   %si les deux si�ges existent
%                 Pos_ = [Pos_; econ_(ii, 1); econ_(ii, end)];  %on ajoute les deux
%                 Weights_ = [Weights_; 195; 195];
%             elseif ~isnan(econ_(ii,1))   %si seulement le premier existe
%                 Pos_ = [Pos_; econ_(ii, 1)]; %on n'ajoute que le premier
%                 Weights_ = [Weights_; 195];
%             elseif ~isnan(econ_(ii,end)) %si seulement le dernier existe
%                 Pos_ = [Pos_; econ_(ii, end)];   %on n'ajoute que le dernier
%                 Weights_ = [Weights_; 195];
%             end
%         else
%             if ~isnan(econ_(ii))
%                 Pos_ = [Pos_; econ_(ii)];  %on ajoute les deux
%                 Weights_ = [Weights_; 195];
%             end
%         end
%     end
%     econ_ = econ_(:,2:end-1);%une fois les colonnes compl�tement remplies on les �limine de la matrice
% end
% 
% % Cargo loading (rear then front)
% Pos_ = [Pos_; X_aft_cargo; X_fwd_cargo];
% Weights_ = [Weights_; aft_ratio*W_cargo; fwd_ratio*W_cargo];
% 
% % w1 = sum(Weights_);
% % Wing Tank Loading
% nstep = 100;     %nstep = nombre de points = nombre de tranches + 1
% Tspan = linspace(ctrtank_span, tank_span, nstep);    %vecteur des BL de la tank tranch�e
% thick_step = Tspan(2)-Tspan(1);     %�paisseur d'une tranche de la tank
% Tspan = Tspan(1:end-1)+thick_step/2;    %vecteur des BL des milieux de tranches
% W_nofuel = sum(Weights_);
% for i=Tspan
%     if i <= WBL1
%         ratio = i/WBL1;
%         Tchord = ratio*(Tchord1-Tchord0) + Tchord0;
%         Tthick = ratio*(Tthick1-Tthick0) + Tthick0;
%         X_cp = ratio*(X_midtank_1-X_midtank_0) + X_midtank_0;
%     elseif i <= WBL2
%         ratio = (i-WBL1)/(WBL2-WBL1);
%         Tchord = ratio*(Tchord2-Tchord1) + Tchord1;
%         Tthick = ratio*(Tthick2-Tthick1) + Tthick1;
%         X_cp = ratio*(X_midtank_2-X_midtank_1) + X_midtank_1;
%     else
%         ratio = (i-WBL2)/(WBL3-WBL2);
%         Tchord = ratio*(Tchord3-Tchord2) + Tchord2;
%         Tthick = ratio*(Tthick3-Tthick2) + Tthick2;
%         X_cp = ratio*(X_midtank_3-X_midtank_2) + X_midtank_2;
%     end
%     vol = thick_step*Tchord*Tthick/tank_scaling_factor;
%     
%     Pos_ = [Pos_; X_cp];
%     Weights_ = [Weights_; vol*rho_fuel*2];
%     
%     if sum(Weights_)-W_nofuel >= W_fuel
%         tank_full = true;
%         break
%     else
%         tank_full = false;
%         missing_fuel = W_fuel + W_nofuel - sum(Weights_);
%     end
% end
% % w2 = sum(Weights_);
% 
% % Center Tank
% if ~tank_full
%     X_ctr_tank = ctrtank_span/WBL1*(X_midtank_1-X_midtank_0)+X_midtank_0;
%     Pos_ = [Pos_; X_ctr_tank];
%     Weights_ = [Weights_; missing_fuel];
% end
% %final computation for excel spreadsheet
% TOW.a = [];
% CG.a = [];
% for i = length(Weights):length(Weights_)    
%     Weight_tot = sum(Weights_(1:i));
%     Weighed_pos = Weights_(1:i).*Pos_(1:i)/Weight_tot;
%     cg = sum(Weighed_pos);
%     TOW.a = [TOW.a; Weight_tot];
%     CG.a = [CG.a; (cg-X_LE_MAC)/MAC*100];
% end
% 
% fwd_lim = min(CG.a);
% aft_lim = max(CG.a);
% 
% xlswrite(filename, CG.a, 'Sheet1', strcat('H53:H200'));%', num2str(fin)));
% xlswrite(filename, TOW.a, 'Sheet1', strcat('J53:J200'));%', num2str(fin)));
% %% load case 2 - front to back, biz then econ
% biz_ = biz;
% econ_ = econ;
% Weights_ = Weights;
% Pos_ = Pos;
% 
% %biz loop first
% for i = 1:ceil(size(biz_,2)/2)  %deux "colonnes" � la fois (ext�rieures)
%     for ii = 1:size(biz_,1)  %front to back
%         if ~isnan(biz_(ii,1)) && ~isnan(biz_(ii,end))   %si les deux si�ges existent
%             Pos_ = [Pos_; biz_(ii, 1); biz_(ii, end)];  %on ajoute les deux
%             Weights_ = [Weights_; 195; 195];
%         elseif ~isnan(biz_(ii,1))       %si seulement le premier existe
%             Pos_ = [Pos_; biz_(ii, 1)]; %on n'ajoute que le premier
%             Weights_ = [Weights_; 195];
%         elseif ~isnan(biz_(ii,end))         %si seulement le dernier existe
%             Pos_ = [Pos_; biz_(ii, end)];   %on n'ajoute que le dernier
%             Weights_ = [Weights_; 195];
%         end
%     end
%     biz_ = biz_(:,2:end-1); %une fois les colonnes compl�tement remplies on les �limine de la matrice
% end
% 
% %econ loop second
% for i = 1:ceil(size(econ_,2)/2)
%     for ii = 1:size(econ_,1)
%         if size(econ_,2) > 1
%             if ~isnan(econ_(ii,1)) && ~isnan(econ_(ii,end))   %si les deux si�ges existent
%                 Pos_ = [Pos_; econ_(ii, 1); econ_(ii, end)];  %on ajoute les deux
%                 Weights_ = [Weights_; 195; 195];
%             elseif ~isnan(econ_(ii,1))   %si seulement le premier existe
%                 Pos_ = [Pos_; econ_(ii, 1)]; %on n'ajoute que le premier
%                 Weights_ = [Weights_; 195];
%             elseif ~isnan(econ_(ii,end)) %si seulement le dernier existe
%                 Pos_ = [Pos_; econ_(ii, end)];   %on n'ajoute que le dernier
%                 Weights_ = [Weights_; 195];
%             end
%         else
%             if ~isnan(econ_(ii))
%                 Pos_ = [Pos_; econ_(ii)];  %on ajoute les deux
%                 Weights_ = [Weights_; 195];
%             end
%         end
%     end
%     econ_ = econ_(:,2:end-1);%une fois les colonnes compl�tement remplies on les �limine de la matrice
% end
% 
% 
% %final computation for excel spreadsheet
% TOW.b = [];
% CG.b = [];
% for i = length(Weights)-3:length(Weights_)    
%     Weight_tot = sum(Weights_(1:i));
%     Weighed_pos = Weights_(1:i).*Pos_(1:i)/Weight_tot;
%     cg = sum(Weighed_pos);
%     TOW.b = [TOW.b; Weight_tot];
%     CG.b = [CG.b; (cg-X_LE_MAC)/MAC*100];
% end
% 
% fwd_lim = min(CG.b);
% aft_lim = max(CG.b);
% 
% xlswrite(filename, CG.b, 'Sheet1', strcat('O53:O200'));
% xlswrite(filename, TOW.b, 'Sheet1', strcat('Q53:Q200'));


%% load case 3 back to front all aircraft
% temp matrices to be able to reuse originals for different load cases
Weights_ = Weights;
Pos_ = Pos;
seats_ = seats;


%1 loop for all seats
for i = 1:ceil(size(seats_,2)/2)
    for ii = size(seats_,1):-1:1 %back to front
        if size(seats_,2) > 1
            if ~isnan(seats_(ii,1)) && ~isnan(seats_(ii,end))   %si les deux si�ges existent
                Pos_ = [Pos_; seats_(ii, 1); seats_(ii, end)];  %on ajoute les deux
                Weights_ = [Weights_; 195; 195];
            elseif ~isnan(seats_(ii,1))   %si seulement le premier existe
                Pos_ = [Pos_; seats_(ii, 1)]; %on n'ajoute que le premier
                Weights_ = [Weights_; 195];
            elseif ~isnan(seats_(ii,end)) %si seulement le dernier existe
                Pos_ = [Pos_; seats_(ii, end)];   %on n'ajoute que le dernier
                Weights_ = [Weights_; 195];
            end
        else
            if ~isnan(seats_(ii))
                Pos_ = [Pos_; seats_(ii)];  %on ajoute les deux
                Weights_ = [Weights_; 195];
            end
        end
    end
    seats_ = seats_(:,2:end-1);%une fois les colonnes compl�tement remplies on les �limine de la matrice
end

% Cargo loading (front then rear)
Pos_ = [Pos_; X_fwd_cargo; X_aft_cargo];
Weights_ = [Weights_; fwd_ratio*W_cargo; aft_ratio*W_cargo];

%final computation for excel spreadsheet
TOW.c = [];
CG.c = [];
for i = length(Weights)-3:length(Weights_)    
    Weight_tot = sum(Weights_(1:i));
    Weighed_pos = Weights_(1:i).*Pos_(1:i)/Weight_tot;
    cg = sum(Weighed_pos);
    TOW.c = [TOW.c; Weight_tot];
    CG.c = [CG.c; (cg-X_LE_MAC)/MAC*100];
end

fwd_lim = min(CG.c);
aft_lim = max(CG.c);

xlswrite(filename, CG.c, 'Sheet1', strcat('H53:H200'));
xlswrite(filename, TOW.c, 'Sheet1', strcat('J53:J200'));

%% load case 4 front to back all aircraft
% temp matrices to be able to reuse originals for different load cases
Weights_ = Weights;
Pos_ = Pos;

seats_ = seats;

%1 loop for all seats
for i = 1:ceil(size(seats_,2)/2)
    for ii = 1:size(seats_,1) %front to back
        if size(seats_,2) > 1
            if ~isnan(seats_(ii,1)) && ~isnan(seats_(ii,end))   %si les deux si�ges existent
                Pos_ = [Pos_; seats_(ii, 1); seats_(ii, end)];  %on ajoute les deux
                Weights_ = [Weights_; 195; 195];
            elseif ~isnan(seats_(ii,1))   %si seulement le premier existe
                Pos_ = [Pos_; seats_(ii, 1)]; %on n'ajoute que le premier
                Weights_ = [Weights_; 195];
            elseif ~isnan(seats_(ii,end)) %si seulement le dernier existe
                Pos_ = [Pos_; seats_(ii, end)];   %on n'ajoute que le dernier
                Weights_ = [Weights_; 195];
            end
        else
            if ~isnan(seats_(ii))
                Pos_ = [Pos_; seats_(ii)];  %on ajoute les deux
                Weights_ = [Weights_; 195];
            end
        end
    end
    seats_ = seats_(:,2:end-1);%une fois les colonnes compl�tement remplies on les �limine de la matrice
end

% Cargo loading (rear then front)
Pos_ = [Pos_; X_aft_cargo; X_fwd_cargo];
Weights_ = [Weights_; aft_ratio*W_cargo; fwd_ratio*W_cargo];

% Wing Tank Loading
nstep = 150;     %nstep = nombre de points = nombre de tranches + 1
Tspan = linspace(ctrtank_span, tank_span, nstep);    %vecteur des BL de la tank tranch�e
thick_step = Tspan(2)-Tspan(1);     %�paisseur d'une tranche de la tank
Tspan = Tspan(1:end-1)+thick_step/2;    %vecteur des BL des milieux de tranches
W_nofuel = sum(Weights_);
for i=Tspan
    if i <= WBL1
        ratio = i/WBL1;
        Tchord = ratio*(Tchord1-Tchord0) + Tchord0;
        Tthick = ratio*(Tthick1-Tthick0) + Tthick0;
        X_cp = ratio*(X_midtank_1-X_midtank_0) + X_midtank_0;
    elseif i <= WBL2
        ratio = (i-WBL1)/(WBL2-WBL1);
        Tchord = ratio*(Tchord2-Tchord1) + Tchord1;
        Tthick = ratio*(Tthick2-Tthick1) + Tthick1;
        X_cp = ratio*(X_midtank_2-X_midtank_1) + X_midtank_1;
    else
        ratio = (i-WBL2)/(WBL3-WBL2);
        Tchord = ratio*(Tchord3-Tchord2) + Tchord2;
        Tthick = ratio*(Tthick3-Tthick2) + Tthick2;
        X_cp = ratio*(X_midtank_3-X_midtank_2) + X_midtank_2;
    end
    vol = thick_step*Tchord*Tthick/tank_scaling_factor;
    
    Pos_ = [Pos_; X_cp];
    Weights_ = [Weights_; vol*rho_fuel*2];
    
    if sum(Weights_)-W_nofuel >= W_fuel
        tank_full = true;
        break
    else
        tank_full = false;
        missing_fuel = W_fuel + W_nofuel - sum(Weights_);
    end
end

% Center Tank
if ~tank_full
    X_ctr_tank = ctrtank_span/WBL1*(X_midtank_1-X_midtank_0)+X_midtank_0;
    Pos_ = [Pos_; X_ctr_tank];
    Weights_ = [Weights_; missing_fuel];
end

%final computation for excel spreadsheet
TOW.d = [];
CG.d = [];
for i = length(Weights)-3:length(Weights_)    
    Weight_tot = sum(Weights_(1:i));
    Weighed_pos = Weights_(1:i).*Pos_(1:i)/Weight_tot;
    cg = sum(Weighed_pos);
    TOW.d = [TOW.d; Weight_tot];
    CG.d = [CG.d; (cg-X_LE_MAC)/MAC*100];
end

fwd_lim = min(CG.d);
aft_lim = max(CG.d);

xlswrite(filename, CG.d, 'Sheet1', strcat('O53:O200'));
xlswrite(filename, TOW.d, 'Sheet1', strcat('Q53:Q200'));
MTOW = max(TOW.d)-500;
xlswrite(filename, MTOW, 'Sheet1', 'R2');
MRW = max(TOW.d);
xlswrite(filename, MRW, 'Sheet1', 'R5');

% High density load
nrow_max = floor((max(seats(end,:))-min(seats(1,:)))/29*12);
if N_pax == 77
    max_pax = 85;
else
    max_pax = 120;
end

% max_pax = nrow_max*size(seats,2);
max_payload = 225*N_pax + 175*(max_pax-N_pax);
MZFW = OWE+max_payload;
xlswrite(filename, MZFW, 'Sheet1', 'R6');
% MLW = 0.85*MTOW;
MLW = MZFW+4500;
xlswrite(filename, MLW, 'Sheet1', 'R7');
