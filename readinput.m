function [] = readinput()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File : readinput.m
%
% Function : This function is used to read Excel CATALIST file and store
%   geometrical data into structures used by other functions.
%
% No Inputs, No Outputs
%
% Authors : Tony LAC
%           Martin LAFRANCE
%           Matthieu PARENTEAU
%
% Date : January 15th, 2015
% Last Edit : March 24th, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Airplane Conditions Parameters output

%%
Conditions = struct();
Airplane   = struct();
xcelfile = 'Output_DT.xlsx';


%% Cruise Conditions
% Conditions.Mach = simulation(4);
% Conditions.Alt = simulation(5);
column = Parameters.excel_column;
%% Airplane DATA
ranges = [column '14:' column '863'];
DATA = xlsread(xcelfile,'Feuil1',ranges);
%% DATA STRUCTURE
Airplane.engine         = struct();
Airplane.Fuse           = struct();
Airplane.Nose           = struct();
Airplane.FSF            = struct();
Airplane.Tailcone       = struct();
Airplane.Tail           = struct();
Airplane.CokeBottle     = struct();
Airplane.Belly_fairing  = struct();
Airplane.Bullet         = struct();
Airplane.Wing           = struct();
Airplane.Winglet        = struct();
Airplane.Fuel           = struct();
Airplane.H_Tail         = struct();
Airplane.V_Tail         = struct();
Airplane.Volume         = struct();
Airplane.Swet           = struct();
Airplane.Delta          = struct();
Airplane.Nacelle        = struct();
Airplane.Open_rotor     = struct();
Airplane.Pylon          = struct();
Airplane.ESDU           = struct();
Airplane.ESDU.Wing      = struct();
Airplane.ESDU.H_Tail    = struct();
Airplane.ESDU.V_Tail    = struct();
%% ========================================================================
%% FSF
Airplane.FSF.Engine = DATA(285);    %(in)
Airplane.FSF.MAC = DATA(286);   %(in)

%% FUSE
Airplane.Fuse.MRP = DATA(1);                % in
Airplane.Fuse.Max_Width = DATA(2);          % in
Airplane.Fuse.Max_height = DATA(3);         % in
Airplane.Fuse.Cross_Sect_1st_a0 = DATA(4);  % ?
Airplane.Fuse.Cross_Sect_1st_a1 = DATA(5);  % ?
Airplane.Fuse.Cross_Sect_1st_b1 = DATA(6);  % ?
Airplane.Fuse.Cross_Sect_2nd_a0 = DATA(7);  % ?
Airplane.Fuse.Cross_Sect_2nd_a1 = DATA(8);  % ?
Airplane.Fuse.Cross_Sect_2nd_b1 = DATA(9);  % ?
Airplane.Fuse.Cross_Sect_3rd_a0 = DATA(10); % ?
Airplane.Fuse.Cross_Sect_3rd_a1 = DATA(11); % ?
Airplane.Fuse.Cross_Sect_3rd_b1 = DATA(12); % ?
Airplane.Fuse.First_section_length = DATA(13);
Airplane.Fuse.Second_section_length = DATA(14);
Airplane.Fuse.length = DATA(13)+DATA(14)+DATA(17)+DATA(21);            %(in)
%% NOSE
Airplane.Nose.Downsweep_angle = DATA(15);   % deg
Airplane.Nose.FS = DATA(16);                % in
Airplane.Nosecone.Length = DATA(17);        % (in)
Airplane.Nose.Upper_Mid_Point = DATA(18);   % (in)
Airplane.Nose.Fwd_Upper_Tangent = DATA(20); % (deg)
Airplane.Nose.Half_Lateral_Ratio = DATA(31);
%% TAILCONE
Airplane.Tailcone.Length = DATA(21);        % (in)
Airplane.Tailcone.to_Fuse_Angle = DATA(22); % (deg)
Airplane.Tailcone.Upsweep = DATA(23);       %(deg)
%% TAIL
Airplane.Tail.Upper_Mid_Point = DATA(24);   % (in)
Airplane.Tail.Lower_Mid_Point = DATA(25);   % (in)
Airplane.Tail.Aft_End_Width_Ratio = DATA(26); 
Airplane.Tail.Aft_End_Height_Ratio = DATA(27);
Airplane.Tail.Aft_Lateral_Angle = DATA(28); % (deg)
Airplane.Tail.Pylon_l_pylon = DATA(296);    %(in)
%%
Airplane.CokeBottle.Lateral_Ratio = DATA(29);
Airplane.CokeBottle.Position = DATA(30);

Airplane.Belly_fairing.Width_over_Fuse_W = DATA(33);
Airplane.Belly_fairing.max_thick_ratio = DATA(34);
Airplane.Belly_fairing.aft_ratio = DATA(35);
Airplane.Belly_fairing.fwd_ratio = DATA(36);
Airplane.Belly_fairing.Vertical_Position = DATA(37);    %(in)

Airplane.Bullet.fairing.Length_Over_V_Tail_Tip_Chord = DATA(38);
Airplane.Bullet.slendeness_ratio = DATA(39);             
Airplane.Bullet.Longitudinal_Offset = DATA(40);         %(in)
%% WING
Airplane.Wing.apex_position = DATA(41);                 %(in) 54 204
Airplane.Wing.Vertical_Position = DATA(42);             %(in)
Airplane.Wing.first_panel_leading_edge_sweep = DATA(43);    % (deg)
Airplane.Wing.second_panel_leading_edge_sweep = DATA(44);   %(deg)
Airplane.Wing.third_panel_leading_edge_sweep = DATA(45);    %(deg)
Airplane.Wing.Root_chord = DATA(46);    %(in)
Airplane.Wing.first_kink_chord = DATA(47);  %(in)
Airplane.Wing.second_kink_chord =DATA(48);  %(in)
Airplane.Wing.tip_chord = DATA(49); %(in)
Airplane.Wing.first_kink_BL = DATA(50); %(in)
Airplane.Wing.second_kink_BL = DATA(51);    %(in)
Airplane.Wing.half_span = DATA(52); %(in)
Airplane.Wing.Root_tc_ratio = DATA(53);
Airplane.Wing.first_kink_tc_ratio = DATA(54);
Airplane.Wing.second_kink_tc_ratio = DATA(55);
Airplane.Wing.tip_tc_ratio = DATA(56);
Airplane.Wing.Root_Incidence_Angle = DATA(57);  %(deg)
Airplane.Wing.first_Kink_Incidence_Angle = DATA(58);    %(deg)
Airplane.Wing.second_Kink_Incidence_Angle = DATA(59);   %(deg)
Airplane.Wing.Tip_Incidence_Angle = DATA(60);   %(deg)
Airplane.Wing.Root_dihedral_angle = DATA(61);   %(deg)
Airplane.Wing.second_panel_dihedral_angle = DATA(62);   %(deg)
Airplane.Wing.third_panel_dihedral_angle = DATA(63);    %(deg)
Airplane.Wing.Nacelle_Lateral_Position = DATA(281); %(in)
Airplane.Wing.Pylon_l_pylon = DATA(282);    %(in)
Airplane.Wing.Front_Spar_Root = DATA(70);
Airplane.Wing.Front_Spar_Kink1 = DATA(71);
Airplane.Wing.Front_Spar_Kink2 = DATA(72);
Airplane.Wing.Front_Spar_Tip = DATA(73);
Airplane.Wing.Rear_Spar_Root = DATA(74);
Airplane.Wing.Rear_Spar_Kink1 = DATA(75);
Airplane.Wing.Rear_Spar_Kink2 = DATA(76);
Airplane.Wing.Rear_Spar_Tip = DATA(77);
Airplane.Wing.TC_fuse_junction = DATA(235);
Airplane.Wing.Gross_Area = DATA(336);   %(m2)
Airplane.Wing.Area_Exposed = DATA(337); %(m2)
%% WINGLET
Airplane.Winglet.span = DATA(64);   %(in)
Airplane.Winglet.root_chord = DATA(65); %(in)
Airplane.Winglet.tip_chord = DATA(66);  %(in)
Airplane.Winglet.leading_edge_sweep = DATA(67); %(deg)
Airplane.Winglet.cant_angle = DATA(68); %(deg)
Airplane.Winglet.incidence_angle = DATA(69);    %(deg)
Airplane.Winglet.Planform_Area = DATA(234); %(m2)
%% FUEL
Airplane.Fuel.Tank_Span = DATA(78);
Airplane.Fuel.Center_Tank_span_over_max_fuse_width = DATA(79);
Airplane.Fuel.K_Factor = DATA(80);
Airplane.Fuel.Center_Tank_Volume = DATA(339);   %(m3)
Airplane.Fuel.Wing_Tank_Volume = DATA(340); %(m3)
Airplane.Fuel.Total_Volume = DATA(341); %(m3)
%% H-TAIL
Airplane.H_Tail.Apex = DATA(81);    %(in)
Airplane.H_Tail.Leading_edge_sweep = DATA(82);  %(deg)
Airplane.H_Tail.vertical_position_ratio = DATA(83);
Airplane.H_Tail.Root_Chord = DATA(84);  %(in)
Airplane.H_Tail.Tip_Chord = DATA(85);   %(in)
Airplane.H_Tail.Half_Span = DATA(86);   %(in)
Airplane.H_Tail.Dihedral_Angle = DATA(87);   %(deg)
Airplane.H_Tail.tc_ratio_at_root = DATA(88);   
Airplane.H_Tail.tc_ratio_at_tip = DATA(89);
%% V-TAIL
Airplane.V_Tail.V_tail_Symmetry = DATA(90);
Airplane.V_Tail.Apex = DATA(91);    %(in)
Airplane.V_Tail.Root_Chord = DATA(92);  %(in)
Airplane.V_Tail.Tip_Chord = DATA(93);   %(in)
Airplane.V_Tail.Span = DATA(94);    %(in)
Airplane.V_Tail.first_panel_Leading_edge_sweep = DATA(95);  %(deg)
Airplane.V_Tail.second_Panel_Leading_Edge_Sweep = DATA(96); %(deg)
Airplane.V_Tail.Kink_Chord_Vert_Pos = DATA(97); %(in)
Airplane.V_Tail.tc_ratio_root = DATA(98);
Airplane.V_Tail.tc_ratio_tip = DATA(99);
%%
% Airplane.Canard.Canard_apex_position = DATA(98);    %(in)
% Airplane.Canard.Canard_Vertical_Position_ratio = DATA(99);
% Airplane.Canard.Canard_Root_chord = DATA(100);  %(in)
% Airplane.Canard.Canard_1st_kink_chord = DATA(101);  %(in)
% Airplane.Canard.Canard_2nd_kink_chord = DATA(102);  %(in)
% Airplane.Canard.Canard_tip_chord = DATA(103);   %(in)
% Airplane.Canard.Canard_1st_kink_BL = DATA(104); %(in)
% Airplane.Canard.Canard_2nd_kink_BL = DATA(105); %(in)
% Airplane.Canard.Canard_Half_Span = DATA(106);   %(in)
% Airplane.Canard.Canard_1st_panel_leading_edge_sweep = DATA(107); %(deg)
% Airplane.Canard.Canard_2nd_panel_leading_edge_sweep = DATA(108);(deg)
% Airplane.Canard.Canard_3rd_panel_leading_edge_sweep (deg)
% Airplane.Canard.Canard_1st_panel_dihedral_angle (deg)
% Airplane.Canard.Canard_2nd_panel_dihedral_angle (deg)
% Airplane.Canard.Canard_3rd_Panel_dihedral_angle (deg)
% Airplane.Canard.Canard_Root_Incidence_Angle (deg)
% Airplane.Canard.Canard_1st_kink_Incidence_Angle (deg)
% Airplane.Canard.Canard_2nd_kink_Incidence_Angle (deg)
% Airplane.Canard.Canard_Tip_Incidence_Angle (deg)
% Airplane.Canard.Canard_Root_tc_ratio
% Airplane.Canard.Canard_1st_kink_tc_ratio
% Airplane.Canard.Canard_2nd_kink_tc_ratio
% Airplane.Canard.Canard_tip_tc_ratio
%%
Airplane.Delta.Fin_Apex_Position = DATA(123);   %(in) 
Airplane.Delta.Fin_Vertical_Position_ratio = DATA(124);
Airplane.Delta.Fin_Root_Chord = DATA(125);  %(in)
Airplane.Delta.Fin_Tip_chord = DATA(126);   %(in)
Airplane.Delta.Fin_Length = DATA(127);  %(in)
Airplane.Delta.Fin_Leading_Edge_Sweep = DATA(128);  %(deg)
Airplane.Delta.Fin_Dihedral_angle = DATA(129);  %(deg)
Airplane.Delta.Fin_tc_ratio = DATA(130);
%% NACELLE
Airplane.Nacelle.Apex = DATA(131);  %(in)
Airplane.Nacelle.Diameter = DATA(132);  %(in)
Airplane.Nacelle.Length = DATA(133);    %(in)
Airplane.Nacelle.Cowl_length = DATA(134);   %(in)
Airplane.Nacelle.Incidence_Angle = DATA(135);   %(deg)
Airplane.Nacelle.Tow_in_Angle = DATA(136);  %(deg)
Airplane.Nacelle.Underwing_Engine_Vertical_Position_Ratio = DATA(137);
Airplane.Nacelle.Underwing_Nacelle_Lateral_Position = DATA(138);    %(in)
Airplane.Nacelle.second_Apex = DATA(139);   %(in)
Airplane.Nacelle.second_Diameter = DATA(140);   %(in)
Airplane.Nacelle.second_Length = DATA(141); %(in)
Airplane.Nacelle.second_Cowl_length = DATA(142);    %(in)
Airplane.Nacelle.second_Incidence_Angle = DATA(143);    %(deg)
Airplane.Nacelle.second_Tow_in_Angle = DATA(144);   %(deg)
Airplane.Nacelle.Underwing_2nd_Nacelle_lateral_position = DATA(145);    %(in)
%%
Airplane.Open_rotor.Apex = DATA(146);   %(in)
Airplane.Open_rotor.Diameter = DATA(147);   %(in)
Airplane.Open_rotor.Length = DATA(148); %(in)
Airplane.Open_rotor.Length = DATA(149); %(in)
Airplane.Open_rotor.CS_Diameter = DATA(150);    %(in)
Airplane.Open_rotor.CS_Length = DATA(151);  %(in)
Airplane.Open_rotor.End_Diameter = DATA(152);   %(in)
Airplane.Open_rotor.Incidence_Angle = DATA(153);    %(deg)
Airplane.Open_rotor.Tow_in_Angle = DATA(154);   %(deg)
Airplane.Open_rotor.Exhaust_Diameter = DATA(155);   %(in)
Airplane.Open_rotor.Exhaust_Length = DATA(156); %(in)
Airplane.Open_rotor.Blade_Root_chord = DATA(157);   %(in)
Airplane.Open_rotor.Blade_Kink_Chord = DATA(158);   %(in)
Airplane.Open_rotor.Blade_Tip_Chord = DATA(159);    %(in)
Airplane.Open_rotor.Blade_Kink_BL = DATA(160);  %(in)
Airplane.Open_rotor.Blade_Length = DATA(161);   %(in)
Airplane.Open_rotor.Blade_Kink_Fwd_Offset = DATA(162);  %(in)
Airplane.Open_rotor.Blade_Tip_Fwd_Offset = DATA(163);   %(in)
Airplane.Open_rotor.Blade_2nd_Ratio = DATA(164);    
Airplane.Open_rotor.Blade_Pitch = DATA(165);    %(deg)
Airplane.Open_rotor.Blade_Offset_ratio = DATA(166);
Airplane.Open_rotor.Blade_2nd_Offset_ratio = DATA(167);
Airplane.Open_rotor.Number_of_blades = DATA(168);
Airplane.Open_rotor.Number_of_blades_2 = DATA(169);
%% PYLON
Airplane.Pylon.Apex = DATA(170);    %(in)
Airplane.Pylon.Root_Chord = DATA(171);  %(in)
Airplane.Pylon.Tip_Chord = DATA(172);   %(in)
Airplane.Pylon.Leading_edge_sweep = DATA(173);  %(deg)
Airplane.Pylon.Vertical_Position_ratio = DATA(174);
Airplane.Pylon.Incidence_Angle = DATA(175);  %(deg)
Airplane.Pylon.Dihedral_Angle = DATA(176);  %(deg)
Airplane.Pylon.Tail_Lateral_Position = DATA(177);   %(in)
Airplane.Pylon.Wing_Pylon_Length = DATA(178);   %(in)
Airplane.Pylon.Apex2 = DATA(179);   %(in)
Airplane.Pylon.Root_Chord2 = DATA(180); %(in)
Airplane.Pylon.Tip_Chord2 = DATA(181);  %(in)
Airplane.Pylon.Leading_edge_sweep2 = DATA(182); %(deg)
Airplane.Pylon.Vertical_Position_ratio2 = DATA(183);
Airplane.Pylon.Incidence_Angle2 = DATA(184);    %(deg)
Airplane.Pylon.Dihedral_Angle2 = DATA(185); %(deg)
Airplane.Pylon.Tail_Lateral_Position = DATA(186);   %(in)
%%
% Airplane.Open_rotor_Pylon.Pylon_3_Apex (in)
% Airplane.Open_rotor_Pylon.Pylon_3_Root_Chord (in)
% Airplane.Open_rotor_Pylon.Pylon_3_Tip_Chord (in)
% Airplane.Open_rotor_Pylon.Pylon_3_Leading_edge_sweep (deg)
Airplane.Open_rotor_Pylon.Pylon_3_Vertical_Position_ratio = DATA(191);
% Airplane.Open_rotor_Pylon.Pylon_3_Incidence_Angle (deg) 
% Airplane.Open_rotor_Pylon.Pylon_3_Dihedral_Angle (deg)
% Airplane.Open_rotor_Pylon.Tail_Pylon_3_Lateral_Position (in)
% Airplane.Turboprop_Apex_Position (in)
% Airplane.Turboprop.Turboprop_Vertical_Position (in)
% Airplane.Turboprop_Lateral_Position (in)
% `Airplane.Turboprop.Cross_Section_Lengths.2nd_Length` (in)
% `Airplane.Turboprop.Cross_Section_Lengths.3rd_Length` (in)
% `Airplane.Turboprop.Cross_Section_Lengths.4th_Length` (in)
% `Airplane.Turboprop.Cross_Section_V-Offset.2nd_V-Offset` (in)
% `Airplane.Turboprop.Cross_Section_V-Offset.3rd_V-Offset` (in)
% `Airplane.Turboprop.Cross_Section_V-Offset.4th_V-Offset` (in)
% `Airplane.Turboprop.Nose_Cone.Nose_Cone_V-Offset` (in)
% Airplane.Turboprop.Nose_Cone.Nose_Cone_Fwd_offset (in)
% Airplane.Turboprop.Nose_Cone.Nose_Cone_Radius (in)
% Airplane.Turboprop.Nose_Cone.Nose_Cone_Lateral_radius (in)
% Airplane.Turboprop.Nose_Cone.Nose_Cone_Length (in)
% Airplane.Turboprop.Nose_Cone.Propeller.Propeller_Fwd_offset (in)
% Airplane.Turboprop.Nose_Cone.Propeller.Propeller_Radius (in)
% `Airplane.Turboprop.Intake.Intake_V-Offset` (in)
% Airplane.Turboprop.Intake.Intake_Fwd_offset (in)
% Airplane.Turboprop.Intake.Intake_Upper_radius (in)
% Airplane.Turboprop.Intake.Intake_Lower_radius (in)
% Airplane.Turboprop.Intake.Intake_Upper_limit (in)
% Airplane.Turboprop.Intake.Intake_Width_Ratio
% Airplane.Turboprop.Intake.Intake_Height_Ratio
% `Airplane.Turboprop.1st_Cross_section.1st_CS_Width` (in)
% `Airplane.Turboprop.1st_Cross_section.1st_CS_Heigth` (in)
% `Airplane.Turboprop.1st_Cross_section.1st_CS_Upper_radius` (in)
% `Airplane.Turboprop.1st_Cross_section.1st_CS_Lower_radius` (in)
% `Airplane.Turboprop.2nd_Cross_section.2nd_CS_Width` (in)
% `Airplane.Turboprop.2nd_Cross_section.2nd_CS_Heigth` (in)
% `Airplane.Turboprop.2nd_Cross_section.2nd_CS_Upper_radius` (in)
% `Airplane.Turboprop.2nd_Cross_section.2nd_CS_Lower_radius` (in)
% `Airplane.Turboprop.3rd_Cross_section.3rd_CS_Width` (in)
% `Airplane.Turboprop.3rd_Cross_section.3rd_CS_Heigth` (in)
% `Airplane.Turboprop.3rd_Cross_section.3rd_CS_Upper_radius` (in)
% `Airplane.Turboprop.3rd_Cross_section.3rd_CS_Lower_radius` (in)
% `Airplane.Turboprop.4th_Cross_section.4th_CS_Width` (in)
% `Airplane.Turboprop.4th_Cross_section.4th_CS_Heigth` (in)
% `Airplane.Turboprop.4th_Cross_section.4th_CS_Upper_radius` (in)
% `Airplane.Turboprop.4th_Cross_section.4th_CS_Lower_radius` (in)
%% VOLUME
Airplane.Volume.Fins = DATA(236);   %(m3)
Airplane.Volume.Wing2 = DATA(237);  %(m3)
Airplane.Volume.Tail_Nacelle2 = DATA(238);  %(m3)
Airplane.Volume.Winglets = DATA(239);   %(m3)
Airplane.Volume.Nacelle_Pylon = DATA(240);  %(m3)
Airplane.Volume.Empennage = DATA(241);  %(m3)
Airplane.Volume.Fuse_Belly = DATA(242); %(m3)
Airplane.Volume.Wing = DATA(243);   %(m3)
Airplane.Volume.Total_Plane_Volume = DATA(283);    %(m3)
%% SWET
Airplane.Swet.Wing2 = DATA(244);    %(m2)
Airplane.Swet.Tail_Nacelle2 = DATA(245);    %(m2)
Airplane.Swet.Ventral_Fins = DATA(246); %(m2)
Airplane.Swet.Bullet = DATA(247);   %(m2)
Airplane.Swet.Belly = DATA(248);    %(m2)
Airplane.Swet.Winglets = DATA(249); %(m2)
Airplane.Swet.Fairing = DATA(289);  %(m2)
Airplane.Swet.Fuse_Nacelle = DATA(290); %(m2)
Airplane.Swet.Fuse_Pylon = DATA(291);   %(m2)
Airplane.Swet.Vtail = DATA(292);    %(m2)
Airplane.Swet.Htail = DATA(293);    %(m2)
Airplane.Swet.Wing = DATA(294); %(m2)
Airplane.Swet.Fuse = DATA(295); %(m2)
Airplane.Swet.Wing_Nacelle = DATA(285); %(m2)
Airplane.Swet.Wing_Pylon = DATA(286);   %(m2)
Airplane.Swet.Total_Wetted_Area = DATA(284); %(m2)
Airplane.Swet.Turboprop_engine = DATA(338); %(m2)
%%
% Airplane.WArea.31 (m2)
% Airplane.WArea.30 (m2)
% Airplane.WArea.29 (m2)
% Airplane.WArea.28 (m2)
% Airplane.WArea.27 (m2)
% Airplane.WArea.26 (m2)
% Airplane.WArea.25 (m2)
% Airplane.WArea.24 (m2)
% Airplane.WArea.23 (m2)
% Airplane.WArea.22 (m2)
% Airplane.WArea.21 (m2)
% Airplane.WArea.20 (m2)
% Airplane.WArea.19 (m2)
% Airplane.WArea.18 (m2)
% Airplane.WArea.17 (m2)
% Airplane.WArea.16 (m2)
% Airplane.WArea.15 (m2)
% Airplane.WArea.14 (m2)
% Airplane.WArea.13 (m2)
% Airplane.WArea.12 (m2)
% Airplane.WArea.11 (m2)
% Airplane.WArea.10 (m2)
% Airplane.WArea.9 (m2)
% Airplane.WArea.8 (m2)
% Airplane.WArea.7 (m2)
% Airplane.WArea.6 (m2)
% Airplane.WArea.5 (m2)
% Airplane.WArea.4 (m2)
% Airplane.WArea.3 (m2)
% Airplane.WArea.2 (m2)
% Airplane.WArea.1 (m2)
%%  ESDU
Airplane.ESDU.V_Tail.Z25_MAC_Nose_Pos = DATA(297);  %(in)
Airplane.ESDU.V_Tail.MAC_LE_X_pos = DATA(298);  %(in)
Airplane.ESDU.V_Tail.MAC_Lat_Pos = DATA(299);   %(in)
Airplane.ESDU.V_Tail.MAC = DATA(300);   %(in)
Airplane.ESDU.V_Tail.HC_Sweep = DATA(301);  %(deg)
Airplane.ESDU.V_Tail.LE_Sweep = DATA(302);  %(deg)
Airplane.ESDU.V_Tail.Span = DATA(303);  %(in)
Airplane.ESDU.V_Tail.Root_Chord = DATA(304);    %(in)
Airplane.ESDU.V_Tail.Arm = DATA(325);   %(in)
Airplane.ESDU.V_Tail.Z25_Chord_Sweep = DATA(326);   %(deg)
Airplane.ESDU.V_Tail.Area = DATA(327);  %(m2)
Airplane.ESDU.V_Tail.Aspect_Ratio = DATA(328);
Airplane.ESDU.V_Tail.Taper_Ratio = DATA(329);

Airplane.ESDU.H_Tail.Z25_MAC_Nose_Pos = DATA(305);  %(in)
Airplane.ESDU.H_Tail.MAC_LE_X_pos = DATA(306);  %(in)
Airplane.ESDU.H_Tail.MAC_Lat_Pos = DATA(307);   %(in)
Airplane.ESDU.H_Tail.MAC = DATA(308);   %(in)
Airplane.ESDU.H_Tail.Half_Chord_Sweep = DATA(309);  %(deg)
Airplane.ESDU.H_Tail.Z25_Chord_Sweep = DATA(319);   %(deg)
Airplane.ESDU.H_Tail.taper_ratio = DATA(320);
Airplane.ESDU.H_Tail.Aspect_Ratio = DATA(321);
Airplane.ESDU.H_Tail.Area = DATA(322);  % (m2) 
Airplane.ESDU.H_Tail.Arm  = DATA(323);  % (in)
Airplane.ESDU.H_Tail.Taper_Ratio = DATA(324); 

Airplane.ESDU.Wing.Z25_MAC_Nose_Pos = DATA(310);    %(in)
Airplane.ESDU.Wing.MAC_LE_X_pos = DATA(311);    %(in)
Airplane.ESDU.Wing.MAC_Lat_Pos = DATA(312); %(in)
Airplane.ESDU.Wing.MAC = DATA(313); %(in)
Airplane.ESDU.Wing.Half_Chord_Sweep = DATA(314);    %(deg)
Airplane.ESDU.Wing.LE_Sweep = DATA(315);    %(deg)
Airplane.ESDU.Wing.Fuse_Jct = DATA(316);    %(in)
Airplane.ESDU.Wing.Fuse_Chord = DATA(317);  %(in)
Airplane.ESDU.Wing.RC = DATA(318);  %(in)
Airplane.ESDU.Wing.Quarter_Chord_Sweep = DATA(330);  %(deg)
Airplane.ESDU.Wing.Taper_Ratio = DATA(331);
Airplane.ESDU.Wing.Aspect_Ratio = DATA(332);
Airplane.ESDU.Wing.Surface = DATA(333); %(m2)

Airplane.ESDU.Orig_Planform_Aspect_Ratio = DATA(334);
Airplane.ESDU.Orig_Planform_Chord_at_Fuse_jct = DATA(335);   %(in)
%% STRUCRURE AND SEATS
% Structure.Structure_fwd_longitudinal (in)
% Structure.Structure_end_longitudinal (in)
% Structure.Wall_insulation_thickness (in)
% Structure.Wall_structure_thickness (in)
% Structure.Floor_heigth (in)
% Structure.Floor_thickness (in)
% Structure.Cargo_Definition.Cargo_Floor_offset (in)
% Structure.Cargo_Doors.Cargo_Door_1.Cargo_1_Presence
% Structure.Cargo_Doors.Cargo_Door_1.Cargo_1_Longitudinal (in)
% Structure.Cargo_Doors.Cargo_Door_1.Cargo_1_Vertical_Offset (in)
% Structure.Cargo_Doors.Cargo_Door_1.Cargo_1_Heigth (in)
% Structure.Cargo_Doors.Cargo_Door_1.Cargo_1_Width (in)
% Structure.Cargo_Doors.Cargo_Door_1.Cargo_1_Radius (in)
% Structure.Cargo_Doors.Cargo_Door_2.Cargo_2_Presence
% Structure.Cargo_Doors.Cargo_Door_2.Cargo_2_Longitudinal (in)
% Structure.Cargo_Doors.Cargo_Door_2.Cargo_2_Vertical_Offset (in)
% Structure.Cargo_Doors.Cargo_Door_2.Cargo_2_Heigth (in)
% Structure.Cargo_Doors.Cargo_Door_2.Cargo_2_Width (in)
% Structure.Cargo_Doors.Cargo_Door_2.Cargo_2_Radius (in)
% Structure.Emergency_Exits.LS_1.Door_LS_1_Presence
% Structure.Emergency_Exits.LS_1.Door_LS_1_Longitudinal (in)
% Structure.Emergency_Exits.LS_1.Door_LS_1_Vertical_Offset (in)
% Structure.Emergency_Exits.LS_1.Door_LS_1_Heigth (in)
% Structure.Emergency_Exits.LS_1.Door_LS_1_Width (in)
% Structure.Emergency_Exits.LS_1.Door_LS_1_Radius (in)
% Structure.Emergency_Exits.LS_2.Door_LS_2_Presence
% Structure.Emergency_Exits.LS_2.Door_LS_2_Longitudinal (in)
% Structure.Emergency_Exits.LS_2.Door_LS_2_Vertical_Offset (in)
% Structure.Emergency_Exits.LS_2.Door_LS_2_Heigth (in)
% Structure.Emergency_Exits.LS_2.Door_LS_2_Width (in)
% Structure.Emergency_Exits.LS_2.Door_LS_2_Radius (in)
% Structure.Emergency_Exits.LS_3.Door_LS_3_Presence
% Structure.Emergency_Exits.LS_3.Door_LS_3_Longitudinal (in)
% Structure.Emergency_Exits.LS_3.Door_LS_3_Vertical_Offset (in)
% Structure.Emergency_Exits.LS_3.Door_LS_3_Heigth (in)
% Structure.Emergency_Exits.LS_3.Door_LS_3_Width (in)
% Structure.Emergency_Exits.LS_3.Door_LS_3_Radius (in)
% Structure.Emergency_Exits.LS_4.Door_LS_4_Presence
% Structure.Emergency_Exits.LS_4.Door_LS_4_Longitudinal (in)
% Structure.Emergency_Exits.LS_4.Door_LS_4_Vertical_Offset (in)
% Structure.Emergency_Exits.LS_4.Door_LS_4_Heigth (in)
% Structure.Emergency_Exits.LS_4.Door_LS_4_Width (in)
% Structure.Emergency_Exits.LS_4.Door_LS_4_Radius (in)
% Structure.Emergency_Exits.RS_1.Door_RS_1_Presence
% Structure.Emergency_Exits.RS_1.Door_RS_1_Longitudinal (in)
% Structure.Emergency_Exits.RS_1.Door_RS_1_Vertical_Offset (in)
% Structure.Emergency_Exits.RS_1.Door_RS_1_Heigth (in)
% Structure.Emergency_Exits.RS_1.Door_RS_1_Width (in)
% Structure.Emergency_Exits.RS_1.Door_RS_1_Radius (in)
% Structure.Emergency_Exits.RS_2.Door_RS_2_Presence
% Structure.Emergency_Exits.RS_2.Door_RS_2_Longitudinal (in)
% Structure.Emergency_Exits.RS_2.Door_RS_2_Vertical_Offset (in)
% Structure.Emergency_Exits.RS_2.Door_RS_2_Heigth (in)
% Structure.Emergency_Exits.RS_2.Door_RS_2_Width (in)
% Structure.Emergency_Exits.RS_2.Door_RS_2_Radius (in)
% Structure.Emergency_Exits.RS_3.Door_RS_3_Presence
% Structure.Emergency_Exits.RS_3.Door_RS_3_Longitudinal (in)
% Structure.Emergency_Exits.RS_3.Door_RS_3_Vertical_Offset (in)
% Structure.Emergency_Exits.RS_3.Door_RS_3_Heigth (in)
% Structure.Emergency_Exits.RS_3.Door_RS_3_Width (in)
% Structure.Emergency_Exits.RS_3.Door_RS_3_Radius (in)
% Structure.Emergency_Exits.RS_4.Door_RS_4_Presence
% Structure.Emergency_Exits.RS_4.Door_RS_4_Longitudinal (in)
% Structure.Emergency_Exits.RS_4.Door_RS_4_Vertical_Offset (in)
% Structure.Emergency_Exits.RS_4.Door_RS_4_Heigth (in)
% Structure.Emergency_Exits.RS_4.Door_RS_4_Width (in)
% Structure.Emergency_Exits.RS_4.Door_RS_4_Radius (in)
% Structure.Bulkheads.Bulk_LS_1.Bulk_LS_1_Presence
% Structure.Bulkheads.Bulk_LS_1.Bulk_LS_1_Longitudinal (in)
% Structure.Bulkheads.Bulk_LS_1.Bulk_LS_1_Lateral (in)
% Structure.Bulkheads.Bulk_LS_1.Bulk_LS_1_Thickness (in)
% Structure.Bulkheads.Bulk_LS_2.Bulk_LS_2_Presence
% Structure.Bulkheads.Bulk_LS_2.Bulk_LS_2_Longitudinal (in)
% Structure.Bulkheads.Bulk_LS_2.Bulk_LS_2_Lateral (in)
% Structure.Bulkheads.Bulk_LS_2.Bulk_LS_2_Thickness (in)
% Structure.Bulkheads.Bulk_RS_1.Bulk_RS_1_Presence
% Structure.Bulkheads.Bulk_RS_1.Bulk_RS_1_Longitudinal (in)
% Structure.Bulkheads.Bulk_RS_1.Bulk_RS_1_Lateral (in)
% Structure.Bulkheads.Bulk_RS_1.Bulk_RS_1_Thickness (in)
% Structure.Bulkheads.Bulk_RS_2.Bulk_RS_2_Presence
% Structure.Bulkheads.Bulk_RS_2.Bulk_RS_2_Longitudinal (in)
% Structure.Bulkheads.Bulk_RS_2.Bulk_RS_2_Lateral (in)
% Structure.Bulkheads.Bulk_RS_2.Bulk_RS_2_Thickness (in)
% Structure.Bulkheads.Bulk_Middle_1.Bulk_Middle_1_Presence
% Structure.Bulkheads.Bulk_Middle_1.Bulk_Middle_1_Longitudinal (in)
% Structure.Bulkheads.Bulk_Middle_1.Bulk_Middle_1_Thickness (in)
% Structure.Frames.Frame_1.Frame_1_Presence
% Structure.Frames.Frame_1.Frame_1_Offset (in)
% Structure.Frames.Frame_1.Frame_1_Thickness (in)
% Structure.Frames.Frame_2.Frame_2_Presence
% Structure.Frames.Frame_2.Frame_2_Offset (in)
% Structure.Frames.Frame_2.Frame_2_Thickness (in)
% Structure.Frames.Frame_3.Frame_3_Presence
% Structure.Frames.Frame_3.Frame_3_Offset (in)
% Structure.Frames.Frame_3.Frame_3_Thickness (in)
% Structure.Frames.Frame_4.Frame_4_Presence
% Structure.Frames.Frame_4.Frame_4_Offset (in)
% Structure.Frames.Frame_4.Frame_4_Thickness (in)
% Structure.Frames.Frame_5.Frame_5_Presence
% Structure.Frames.Frame_5.Frame_5_Offset (in)
% Structure.Frames.Frame_5.Frame_5_Thickness (in)
% Structure.Frames.Frame_6.Frame_6_Presence
% Structure.Frames.Frame_6.Frame_6_Offset (in)
% Structure.Frames.Frame_6.Frame_6_Thickness (in)
% Structure.Frames.Frame_7.Frame_7_Presence
% Structure.Frames.Frame_7.Frame_7_Offset (in)
% Structure.Frames.Frame_7.Frame_7_Thickness (in)
% Structure.Frames.Frame_8.Frame_8_Presence
% Structure.Frames.Frame_8.Frame_8_Offset (in)
% Structure.Frames.Frame_8.Frame_8_Thickness (in)
% Structure.Frames.Frame_9.Frame_9_Presence
% Structure.Frames.Frame_9.Frame_9_Offset (in)
% Structure.Frames.Frame_9.Frame_9_Thickness (in)
% Structure.Frames.Frame_10.Frame_10_Presence
% Structure.Frames.Frame_10.Frame_10_Offset (in)
% Structure.Frames.Frame_10.Frame_10_Thickness (in)
% Structure.Frames.Frame_11.Frame_11_Presence
% Structure.Frames.Frame_11.Frame_11_Offset (in)
% Structure.Frames.Frame_11.Frame_11_Thickness (in)
% Structure.Frames.Frame_12.Frame_12_Presence
% Structure.Frames.Frame_12.Frame_12_Offset (in)
% Structure.Frames.Frame_12.Frame_12_Thickness (in)
% Structure.Frames.Frame_13.Frame_13_Presence
% Structure.Frames.Frame_13.Frame_13_Offset (in)
% Structure.Frames.Frame_13.Frame_13_Thickness (in)
% Structure.Frames.Frame_14.Frame_14_Presence
% Structure.Frames.Frame_14.Frame_14_Offset (in)
% Structure.Frames.Frame_14.Frame_14_Thickness (in)
% Structure.Frames.Frame_15.Frame_15_Presence
% Structure.Frames.Frame_15.Frame_15_Offset (in)
% Structure.Frames.Frame_15.Frame_15_Thickness (in)
% Structure.Frames.Frame_16.Frame_16_Presence
% Structure.Frames.Frame_16.Frame_16_Offset (in)
% Structure.Frames.Frame_16.Frame_16_Thickness (in)
% Structure.Frames.Frame_17.Frame_17_Presence
% Structure.Frames.Frame_17.Frame_17_Offset (in)
% Structure.Frames.Frame_17.Frame_17_Thickness (in)
% Structure.Frames.Frame_18.Frame_18_Presence
% Structure.Frames.Frame_18.Frame_18_Offset (in)
% Structure.Frames.Frame_18.Frame_18_Thickness (in)
% Structure.Frames.Frame_19.Frame_19_Presence
% Structure.Frames.Frame_19.Frame_19_Offset (in)
% Structure.Frames.Frame_19.Frame_19_Thickness (in)
% Structure.Frames.Frame_20.Frame_20_Presence
% Structure.Frames.Frame_20.Frame_20_Offset (in)
% Structure.Frames.Frame_20.Frame_20_Thickness (in)
% Structure.Frames.Frame_21.Frame_21_Presence
% Structure.Frames.Frame_21.Frame_21_Offset (in)
% Structure.Frames.Frame_21.Frame_21_Thickness (in)
% Structure.Frames.Frame_22.Frame_22_Presence
% Structure.Frames.Frame_22.Frame_22_Offset (in)
% Structure.Frames.Frame_22.Frame_22_Thickness (in)
% Structure.Frames.Frame_23.Frame_23_Presence
% Structure.Frames.Frame_23.Frame_23_Offset (in)
% Structure.Frames.Frame_23.Frame_23_Thickness (in)
% Structure.Frames.Frame_24.Frame_24_Presence
% Structure.Frames.Frame_24.Frame_24_Offset (in)
% Structure.Frames.Frame_24.Frame_24_Thickness (in)
% Structure.Frames.Frame_25.Frame_25_Presence
% Structure.Frames.Frame_25.Frame_25_Offset (in)
% Structure.Frames.Frame_25.Frame_25_Thickness (in)
% Structure.Frames.Frame_26.Frame_26_Presence
% Structure.Frames.Frame_26.Frame_26_Offset (in)
% Structure.Frames.Frame_26.Frame_26_Thickness (in)
% Structure.Frames.Frame_27.Frame_27_Presence
% Structure.Frames.Frame_27.Frame_27_Offset (in)
% Structure.Frames.Frame_27.Frame_27_Thickness (in)
% Structure.Frames.Frame_28.Frame_28_Presence
% Structure.Frames.Frame_28.Frame_28_Offset (in)
% Structure.Frames.Frame_28.Frame_28_Thickness (in)
% Structure.Frames.Frame_29.Frame_29_Presence
% Structure.Frames.Frame_29.Frame_29_Offset (in)
% Structure.Frames.Frame_29.Frame_29_Thickness (in)
% Structure.Frames.Frame_30.Frame_30_Presence
% Structure.Frames.Frame_30.Frame_30_Offset (in)
% Structure.Frames.Frame_30.Frame_30_Thickness (in)
% Structure.Frames.Frame_31.Frame_31_Presence
% Structure.Frames.Frame_31.Frame_31_Offset (in)
% Structure.Frames.Frame_31.Frame_31_Thickness (in)
% Structure.Frames.Frame_32.Frame_32_Presence
% Structure.Frames.Frame_32.Frame_32_Offset (in)
% Structure.Frames.Frame_32.Frame_32_Thickness (in)
% Structure.Frames.Frame_33.Frame_33_Presence
% Structure.Frames.Frame_33.Frame_33_Offset (in)
% Structure.Frames.Frame_33.Frame_33_Thickness (in)
% Structure.Frames.Frame_34.Frame_34_Presence
% Structure.Frames.Frame_34.Frame_34_Offset (in)
% Structure.Frames.Frame_34.Frame_34_Thickness (in)
% Structure.Frames.Frame_35.Frame_35_Presence
% Structure.Frames.Frame_35.Frame_35_Offset (in)
% Structure.Frames.Frame_35.Frame_35_Thickness (in)
% Structure.Frames.Frame_36.Frame_36_Presence
% Structure.Frames.Frame_36.Frame_36_Offset (in)
% Structure.Frames.Frame_36.Frame_36_Thickness (in)
% Structure.Frames.Frame_37.Frame_37_Presence
% Structure.Frames.Frame_37.Frame_37_Offset (in)
% Structure.Frames.Frame_37.Frame_37_Thickness (in)
% Structure.Frames.Frame_38.Frame_38_Presence
% Structure.Frames.Frame_38.Frame_38_Offset (in)
% Structure.Frames.Frame_38.Frame_38_Thickness (in)
% Structure.Frames.Frame_39.Frame_39_Presence
% Structure.Frames.Frame_39.Frame_39_Offset (in)
% Structure.Frames.Frame_39.Frame_39_Thickness (in)
% Structure.Frames.Frame_40.Frame_40_Presence
% Structure.Frames.Frame_40.Frame_40_Offset (in)
% Structure.Frames.Frame_40.Frame_40_Thickness (in)
% Structure.Frames.Frame_41.Frame_41_Presence
% Structure.Frames.Frame_41.Frame_41_Offset (in)
% Structure.Frames.Frame_41.Frame_41_Thickness (in)
% Structure.Frames.Frame_42.Frame_42_Presence
% Structure.Frames.Frame_42.Frame_42_Offset (in)
% Structure.Frames.Frame_42.Frame_42_Thickness (in)
% Structure.Frames.Frame_43.Frame_43_Presence
% Structure.Frames.Frame_43.Frame_43_Offset (in)
% Structure.Frames.Frame_43.Frame_43_Thickness (in)
% Structure.Frames.Frame_44.Frame_44_Presence
% Structure.Frames.Frame_44.Frame_44_Offset (in)
% Structure.Frames.Frame_44.Frame_44_Thickness (in)
% Structure.Frames.Frame_45.Frame_45_Presence
% Structure.Frames.Frame_45.Frame_45_Offset (in)
% Structure.Frames.Frame_45.Frame_45_Thickness (in)
% Structure.Frames.Frame_46.Frame_46_Presence
% Structure.Frames.Frame_46.Frame_46_Offset (in)
% Structure.Frames.Frame_46.Frame_46_Thickness (in)
% Structure.Frames.Frame_47.Frame_47_Presence
% Structure.Frames.Frame_47.Frame_47_Offset (in)
% Structure.Frames.Frame_47.Frame_47_Thickness (in)
% Structure.Frames.Frame_48.Frame_48_Presence
% Structure.Frames.Frame_48.Frame_48_Offset (in)
% Structure.Frames.Frame_48.Frame_48_Thickness (in)
% Structure.Frames.Frame_49.Frame_49_Presence
% Structure.Frames.Frame_49.Frame_49_Offset (in)
% Structure.Frames.Frame_49.Frame_49_Thickness (in)
% Structure.Frames.Frame_50.Frame_50_Presence
% Structure.Frames.Frame_50.Frame_50_Offset (in)
% Structure.Frames.Frame_50.Frame_50_Thickness (in)
% Structure.Frames.Frame_51.Frame_51_Presence
% Structure.Frames.Frame_51.Frame_51_Offset (in)
% Structure.Frames.Frame_51.Frame_51_Thickness (in)
% Structure.Frames.Frame_52.Frame_52_Presence
% Structure.Frames.Frame_52.Frame_52_Offset (in)
% Structure.Frames.Frame_52.Frame_52_Thickness (in)
% Structure.Windows.Windows_Vertical_Offset (in)
% Structure.Windows.Windows_Heigth (in)
% Structure.Windows.Windows_Width (in)
% Structure.Windows.Windows_Radius (in)
% Structure.Windows.Window_1.Win_1_Presence
% Structure.Windows.Window_1.Win_1_Offset (in)
% Structure.Windows.Window_2.Win_2_Presence
% Structure.Windows.Window_2.Win_2_Offset (in)
% Structure.Windows.Window_3.Win_3_Presence
% Structure.Windows.Window_3.Win_3_Offset (in)
% Structure.Windows.Window_4.Win_4_Presence
% Structure.Windows.Window_4.Win_4_Offset (in)
% Structure.Windows.Window_5.Win_5_Presence
% Structure.Windows.Window_5.Win_5_Offset (in)
% Structure.Windows.Window_6.Win_6_Presence
% Structure.Windows.Window_6.Win_6_Offset (in)
% Structure.Windows.Window_7.Win_7_Presence
% Structure.Windows.Window_7.Win_7_Offset (in)
% Structure.Windows.Window_8.Win_8_Presence
% Structure.Windows.Window_8.Win_8_Offset (in)
% Structure.Windows.Window_9.Win_9_Presence
% Structure.Windows.Window_9.Win_9_Offset (in)
% Structure.Windows.Window_10.Win_10_Presence
% Structure.Windows.Window_10.Win_10_Offset (in)
% Structure.Windows.Window_11.Win_11_Presence
% Structure.Windows.Window_11.Win_11_Offset (in)
% Structure.Windows.Window_12.Win_12_Presence
% Structure.Windows.Window_12.Win_12_Offset (in)
% Structure.Windows.Window_13.Win_13_Presence
% Structure.Windows.Window_13.Win_13_Offset (in)
% Structure.Windows.Window_14.Win_14_Presence
% Structure.Windows.Window_14.Win_14_Offset (in)
% Structure.Windows.Window_15.Win_15_Presence
% Structure.Windows.Window_15.Win_15_Offset (in)
% Structure.Windows.Window_16.Win_16_Presence
% Structure.Windows.Window_16.Win_16_Offset (in)
% Structure.Windows.Window_17.Win_17_Presence
% Structure.Windows.Window_17.Win_17_Offset (in)
% Structure.Windows.Window_18.Win_18_Presence
% Structure.Windows.Window_18.Win_18_Offset (in)
% Structure.Windows.Window_19.Win_19_Presence
% Structure.Windows.Window_19.Win_19_Offset (in)
% Structure.Windows.Window_20.Win_20_Presence
% Structure.Windows.Window_20.Win_20_Offset (in)
% Structure.Windows.Window_21.Win_21_Presence
% Structure.Windows.Window_21.Win_21_Offset (in)
% Structure.Windows.Window_22.Win_22_Presence
% Structure.Windows.Window_22.Win_22_Offset (in)
% Structure.Windows.Window_23.Win_23_Presence
% Structure.Windows.Window_23.Win_23_Offset (in)
% Structure.Windows.Window_24.Win_24_Presence
% Structure.Windows.Window_24.Win_24_Offset (in)
% Structure.Windows.Window_25.Win_25_Presence
% Structure.Windows.Window_25.Win_25_Offset (in)
% Structure.Windows.Window_26.Win_26_Presence
% Structure.Windows.Window_26.Win_26_Offset (in)
% Structure.Windows.Window_27.Win_27_Presence
% Structure.Windows.Window_27.Win_27_Offset (in)
% Structure.Windows.Window_28.Win_28_Presence
% Structure.Windows.Window_28.Win_28_Offset (in)
% Structure.Windows.Window_29.Win_29_Presence
% Structure.Windows.Window_29.Win_29_Offset (in)
% Structure.Windows.Window_30.Win_30_Presence
% Structure.Windows.Window_30.Win_30_Offset (in)
% Structure.Windows.Window_31.Win_31_Presence
% Structure.Windows.Window_31.Win_31_Offset (in)
% Structure.Windows.Window_32.Win_32_Presence
% Structure.Windows.Window_32.Win_32_Offset (in)
% Structure.Windows.Window_33.Win_33_Presence
% Structure.Windows.Window_33.Win_33_Offset (in)
% Structure.Windows.Window_34.Win_34_Presence
% Structure.Windows.Window_34.Win_34_Offset (in)
% Structure.Windows.Window_35.Win_35_Presence
% Structure.Windows.Window_35.Win_35_Offset (in)
% Structure.Windows.Window_36.Win_36_Presence
% Structure.Windows.Window_36.Win_36_Offset (in)
% Structure.Windows.Window_37.Win_37_Presence
% Structure.Windows.Window_37.Win_37_Offset (in)
% Structure.Windows.Window_38.Win_38_Presence
% Structure.Windows.Window_38.Win_38_Offset (in)
% Structure.Windows.Window_39.Win_39_Presence
% Structure.Windows.Window_39.Win_39_Offset (in)
% Structure.Windows.Window_40.Win_40_Presence
% Structure.Windows.Window_40.Win_40_Offset (in)
% Structure.Windows.Window_41.Win_41_Presence
% Structure.Windows.Window_41.Win_41_Offset (in)
% Structure.Windows.Window_42.Win_42_Presence
% Structure.Windows.Window_42.Win_42_Offset (in)
% Structure.Windows.Window_43.Win_43_Presence
% Structure.Windows.Window_43.Win_43_Offset (in)
% Structure.Windows.Window_44.Win_44_Presence
% Structure.Windows.Window_44.Win_44_Offset (in)
% Structure.Windows.Window_45.Win_45_Presence
% Structure.Windows.Window_45.Win_45_Offset (in)
% Structure.Windows.Window_46.Win_46_Presence
% Structure.Windows.Window_46.Win_46_Offset (in)
% Structure.Windows.Window_47.Win_47_Presence
% Structure.Windows.Window_47.Win_47_Offset (in)
% Structure.Windows.Window_48.Win_48_Presence
% Structure.Windows.Window_48.Win_48_Offset (in)
% Structure.Windows.Window_49.Win_49_Presence
% Structure.Windows.Window_49.Win_49_Offset (in)
% Structure.Windows.Window_50.Win_50_Presence
% Structure.Windows.Window_50.Win_50_Offset (in)
% Structure.Wing.Wing_Spars_Activation
% Structure.Wing.Wing_Rib_Root_Activation
% Structure.Wing.Wing_Rib_1st_Activation
% Structure.Wing.Wing_Rib_2nd_Activation
% Structure.Wing.Wing_Rib_Root_Number
% Structure.Wing.Wing_Rib_1st_Number
% Structure.Wing.Wing_Rib_2nd_Number
% Structure.Wing.Wing_Rib_Root_Spacing (in)
% Structure.Wing.Wing_Rib_1st_Spacing (in)
% Structure.Wing.Wing_Rib_2nd_Spacing (in)
% `Structure.H-Tail.H-Tail_Structure_Activation`
% `Structure.H-Tail.H-Tail_Span_Ratio`
% `Structure.H-Tail.H-Tail_Root_Front_Spar`
% `Structure.H-Tail.H-Tail_Root_Rear_Spar`
% `Structure.H-Tail.H-Tail_Tip_Front_Spar`
% `Structure.H-Tail.H-Tail_Tip_Rear_Spar`
% `Structure.H-Tail.H-Tail_Ribs_Angle` (deg)
% `Structure.H-Tail.H-Tail_1st_Rib_Offset` (in)
% `Structure.H-Tail.H-Tail_Ribs_Number`
% `Structure.H-Tail.H-Tail_Ribs_Spacing` (in)
% `Structure.V-Tail.V-Tail_Structure_Activation`
% `Structure.V-Tail.V-Tail_Middle_Spar_Activation`
% `Structure.V-Tail.V-Tail_Root_Front_Spar`
% `Structure.V-Tail.V-Tail_Root_Mid_Spar`
% `Structure.V-Tail.V-Tail_Root_Rear_Spar`
% `Structure.V-Tail.V-Tail_Tip_Front_Spar`
% `Structure.V-Tail.V-Tail_Tip_Mid_Spar`
% `Structure.V-Tail.V-Tail_Tip_Rear_Spar`
% `Structure.V-Tail.V-Tail_Ribs_Angle` (deg)
% `Structure.V-Tail.V-Tail_1st_Rib_Offset` (in)
% `Structure.V-Tail.V-Tail_Ribs_Number`
% `Structure.V-Tail.V-Tail_Ribs_Spacing` (in)

% Airplane.Seats.Seat_width (in)
% Airplane.Seats.Seat_recline_angle (deg)
% Airplane.Seats.Seat_heigth (in)
% Airplane.Seats.Armrest_width (in)
% Airplane.Seats.Biz_Seat_width (in)
% Airplane.Seats.Biz_Seat_recline_angle (deg)
% Airplane.Seats.Biz_Seat_heigth (in)
% Airplane.Seats.Biz_Armrest_width (in)
% Airplane.Seats.Presence_LS_Seat_1
% Airplane.Seats.LS_Seat_1_Biz
% Airplane.Seats.LS_Seat_1_Lateral_Position (in)
% Airplane.Seats.LS_Seat_1_Longitudinal_Position (in)
% Airplane.Seats.LS_Seat_1_LS_Seat_1_Abreast
% Airplane.Seats.LS_Seat_1_Number_of_rows
% Airplane.Seats.LS_Seat_1_Seat_pitch (in)
% Airplane.Seats.Presence_LS_Seat_2
% Airplane.Seats.LS_Seat_2_Lateral_Position (in)
% Airplane.Seats.LS_Seat_2_Longitudinal_Position (in)
% Interior_Control.Seats_Configuration.LS_Seat_2.LS_Seat_2_Abreast
% Interior_Control.Seats_Configuration.LS_Seat_2.Number_of_rows
% Interior_Control.Seats_Configuration.LS_Seat_2.Seat_pitch (in)
% Interior_Control.Seats_Configuration.LS_Seat_3.Presence_LS_Seat_3
% Interior_Control.Seats_Configuration.LS_Seat_3.Lateral_Position (in)
% Interior_Control.Seats_Configuration.LS_Seat_3.Longitudinal_Position (in)
% Interior_Control.Seats_Configuration.LS_Seat_3.LS_Seat_3_Abreast
% Interior_Control.Seats_Configuration.LS_Seat_3.Number_of_rows
% Interior_Control.Seats_Configuration.LS_Seat_3.Seat_pitch (in)
% Interior_Control.Presence_RS_Seat_1
% Interior_Control.Seats_Configuration.RS_Seat_1.RS_Seat_1_Biz
% Interior_Control.Seats_Configuration.RS_Seat_1.Lateral_Position (in)
% Interior_Control.Seats_Configuration.RS_Seat_1.Longitudinal_Position (in)
% Interior_Control.Seats_Configuration.RS_Seat_1.RS_Seat_1_Abreast
% Interior_Control.Seats_Configuration.RS_Seat_1.Number_of_rows
% Interior_Control.Seats_Configuration.RS_Seat_1.Seat_pitch (in)
% Interior_Control.Presence_RS_Seat_2
% Interior_Control.Seats_Configuration.RS_Seat_2.Lateral_Position (in)
% Interior_Control.Seats_Configuration.RS_Seat_2.Longitudinal_Position (in)
% Interior_Control.Seats_Configuration.RS_Seat_2.RS_Seat_2_Abreast
% Interior_Control.Seats_Configuration.RS_Seat_2.Number_of_rows
% Interior_Control.Seats_Configuration.RS_Seat_2.Seat_pitch (in)
% Interior_Control.Seats_Configuration.RS_Seat_3.Presence_RS_Seat_3
% Interior_Control.Seats_Configuration.RS_Seat_3.Lateral_Position (in)
% Interior_Control.Seats_Configuration.RS_Seat_3.Longitudinal_Position (in)
% Interior_Control.Seats_Configuration.RS_Seat_3.RS_Seat_3_Abreast
% Interior_Control.Seats_Configuration.RS_Seat_3.Number_of_rows
% Interior_Control.Seats_Configuration.RS_Seat_3.Seat_pitch (in)
% `Interior_Control.Cargo.LD3_Container_1.Presence_LD3-Container_1`
% Interior_Control.Cargo.LD3_Container_1.Container_Longitudinal_Position (in)
% `Interior_Control.LD3-Container_1_Spacing` (in)
% `Interior_Control.LD3-Container_1_Number`
% `Interior_Control.Cargo.LD3_Container_2.Presence_LD3-Container_2`
% Interior_Control.Cargo.LD3_Container_2.Container_Longitudinal_Position_2 (in)
% `Interior_Control.Cargo.LD3_Container_2.LD3-Container_2_Spacing` (in)
% `Interior_Control.Cargo.LD3_Container_2.LD3-Container_2_Number`
% Interior_Control.Lavatory.Lavatory_LS_1.Lavatory_LS_1_Presence
% Interior_Control.Lavatory.Lavatory_LS_1.Lavatory_LS_1_Longitudinal (in)
% Interior_Control.Lavatory.Lavatory_LS_1.Lavatory_LS_1_Lateral (in)
% Interior_Control.Lavatory.Lavatory_LS_1.Lavatory_LS_1_Width (in)
% Interior_Control.Lavatory.Lavatory_LS_2.Lavatory_LS_2_Presence
% Interior_Control.Lavatory.Lavatory_LS_2.Lavatory_LS_2_Longitudinal (in)
% Interior_Control.Lavatory.Lavatory_LS_2.Lavatory_LS_2_Lateral (in)
% Interior_Control.Lavatory.Lavatory_LS_2.Lavatory_LS_2_Width (in)
% Interior_Control.Lavatory.Lavatory_LS_3.Lavatory_LS_3_Presence
% Interior_Control.Lavatory.Lavatory_LS_3.Lavatory_LS_3_Longitudinal (in)
% Interior_Control.Lavatory.Lavatory_LS_3.Lavatory_LS_3_Lateral (in)
% Interior_Control.Lavatory.Lavatory_LS_3.Lavatory_LS_3_Width (in)
% Interior_Control.Lavatory.Lavatory_RS_1.Lavatory_RS_1_Presence
% Interior_Control.Lavatory.Lavatory_RS_1.Lavatory_RS_1_Longitudinal (in)
% Interior_Control.Lavatory.Lavatory_RS_1.Lavatory_RS_1_Lateral (in)
% Interior_Control.Lavatory.Lavatory_RS_1.Lavatory_RS_1_Width (in)
% Interior_Control.Lavatory.Lavatory_RS_2.Lavatory_RS_2_Presence
% Interior_Control.Lavatory.Lavatory_RS_2.Lavatory_RS_2_Longitudinal (in)
% Interior_Control.Lavatory.Lavatory_RS_2.Lavatory_RS_2_Lateral (in)
% Interior_Control.Lavatory.Lavatory_RS_2.Lavatory_RS_2_Width (in)
% Interior_Control.Overhead_Bin.Overhead_X_Position (in)
% Interior_Control.Overhead_Bin.Overhead_Y_Position (in)
% Interior_Control.Overhead_Bin.Overhead_H_Angle (deg)
% Interior_Control.Overhead_Bin.Overhead_V_Angle (deg)
% Interior_Control.Overhead_Bin.Overhead_Radius (in)
% Interior_Control.Overhead_Bin.Overhead_Longitudinal (in)
% Interior_Control.Overhead_Bin.Overhead_Length (in)
% Interior_Control.Overhead_Bin.Overhead_Roof_Radius (in)
% Interior_Control.Galleys.Galley_LS_1.Galley_LS_1_Presence
% Interior_Control.Galleys.Galley_LS_1.Galley_LS_1_Longitudinal (in)
% Interior_Control.Galleys.Galley_LS_1.Galley_LS_1_Lateral (in)
% Interior_Control.Galleys.Galley_LS_1.Galley_LS_1_Width (in)
% Interior_Control.Galleys.Galley_LS_2.Galley_LS_2_Presence
% Interior_Control.Galleys.Galley_LS_2.Galley_LS_2_Longitudinal (in)
% Interior_Control.Galleys.Galley_LS_2.Galley_LS_2_Lateral (in)
% Interior_Control.Galleys.Galley_LS_2.Galley_LS_2_Width (in)
% Interior_Control.Galleys.Galley_RS_1.Galley_RS_1_Presence
% Interior_Control.Galleys.Galley_RS_1.Galley_RS_1_Longitudinal (in)
% Interior_Control.Galleys.Galley_RS_1.Galley_RS_1_Lateral (in)
% Interior_Control.Galleys.Galley_RS_1.Galley_RS_1_Width (in)
% Interior_Control.Galleys.Galley_RS_2.Galley_RS_2_Presence
% Interior_Control.Galleys.Galley_RS_2.Galley_RS_2_Longitudinal (in)
% Interior_Control.Galleys.Galley_RS_2.Galley_RS_2_Lateral (in)
% Interior_Control.Galleys.Galley_RS_2.Galley_RS_2_Width (in)
% Interior_Control.Cockpit_Seats.Presence_Cockpit_Seats
% Interior_Control.Cockpit_Seats.Lateral_Position (in)
% Interior_Control.Cockpit_Seats.Longitudinal_Position (in)
% Interior_Control.Cockpit_Seats.Seats_Lateral_Separation (in)
% Interior_Control.Cockpit_Seats.Seat_width (in)
% Interior_Control.Cockpit_Seats.Seat_recline_angle (deg)
% Interior_Control.Cockpit_Seats.Seat_heigth (in)
% Interior_Control.Cockpit_Seats.Armrest_width (in)
% LG_Control.Landing_Gear_Position
% LG_Control.Nose_Landing_Gear.NLG_Longitudinal_Position (in)
% LG_Control.Nose_Landing_Gear.NLG_Vertical_Position (in)
% LG_Control.Nose_Landing_Gear.NLG_Tire_Diameter (in)
% LG_Control.Nose_Landing_Gear.NLG_Rim_Diameter (in)
% LG_Control.Nose_Landing_Gear.NLG_Tire_Width (in)
% LG_Control.Nose_Landing_Gear.NLG_Support_Length (in)
% LG_Control.Nose_Landing_Gear.NLG_Heigth (in)
% LG_Control.Main_Landing_Gear.MLG_Longitudinal_Position (in)
% LG_Control.Main_Landing_Gear.MLG_Vertical_Position (in)
% LG_Control.Main_Landing_Gear.MLG_Lateral_Position (in)
% LG_Control.MLG_retracted_angle (deg)
% LG_Control.Main_Landing_Gear.MLG_Tire_Diameter (in)
% LG_Control.Main_Landing_Gear.MLG_Rim_Diameter (in)
% LG_Control.Main_Landing_Gear.MLG_Tire_Width (in)
% LG_Control.Main_Landing_Gear.MLG_Support_Length (in)
% LG_Control.Main_Landing_Gear.MLG_Heigth (in)

ranges = [column '1:' column '863'];
[num, ~, raw] = xlsread(xcelfile,'Feuil1',ranges);
output = [NaN(13,1); num];


for i = 1:size(raw,1)
    if strcmp(raw(i),'true ')
        output(i) = 1;
    elseif strcmp(raw(i),'false ')
        output(i) = 0;
    end
end


end