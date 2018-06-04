clc; clear all;
%% inputs

p_out = 250e3;    %W, output power
v_in = 400;   %V, line to line
v_in_ph = v_in/sqrt(3);   %V, phase
n = 742;    %rpm, rated speed

%% design parameters

p = 8;  %pole numbers
f = 50;  %Hz, frequency

mag_load = 0.5;   %T, average flux density in airgap
el_load = 49e3;   %A/m,  electrical loading
b_teeth_peak = 1.5;  %T, teeth flux density
b_yoke_peak = 1.5;  %T, yoke flux density

cur_dens = 6e6;   %A/m2,  stator current density

pf_init = 0.84;   %initial power factor assumption
eff_init = 0.95;  %initial efficiency assumption
kw_1 = 0.966;    %initial winding factor assumption
stack = 0.96;   % stacking factor
dia_strand = 1.45e-3;   %m, AWG 12 wire diameter
ff = 0.5; %fill factor
slot_mount = 2e-3;   %m, slot mount

q_s = 2;  %number of slots per pole per phase
m = 3;  %phase number


%% aspect ratio

asp_ratio = pi/p * (p/2)^(1/3);   % ratio of l_eff/dia_bore approximately


%% torque calculation

torque = p_out / (n*pi/30);  %Nm, rated torque

%% tangential stress calculation

mag_load_peak = mag_load * pi/2;  %T, magnetic loading peak
el_load_peak  = el_load * sqrt(2);  %A/m, electrical loading peak

stress = el_load_peak * mag_load_peak * pf_init / 2;



%% bore diameter and effective length

dia_bore = (  2*torque / (stress*pi*asp_ratio)  )^(1/3);   %m, bore diameter
l_eff = dia_bore * asp_ratio;  %m, effective axial length

%% air gap calculation

g_unrounded = 0.18e-3 + 0.006*(p_out)^(0.4)*1e-3; %m, air gap
g = round(g_unrounded,3);
g=3e-3;

%% stator design

Q_s = q_s * m * p;  % stator slot number
dia_in_stator = dia_bore+g;  %m, stator inner diameter
slot_pitch = pi * dia_in_stator / Q_s;  %m, stator slot pitch

%% teeth width calculation

% b_teeth_avg = b_teeth_peak * 2/pi;   %T, average teeth flux density
l_act = l_eff - 2*g;   %m, actual axial length
width_teeth = l_eff * slot_pitch * mag_load_peak / (stack * l_act * b_teeth_peak );
%% winding factor 

kp = 1; %pitch factor, single layer
slot_ang = pi * p / Q_s;   %radians, slot angle

kd = sin(q_s*slot_ang/2) / (q_s * sin(slot_ang/2));  %dist factor
kw_act = kp * kd;  %winding factor

%% number of series turns

area_pole_eff = pi * dia_bore * l_eff / p;   %m, pole area
flux_pp = mag_load * area_pole_eff;   %Wb, flux per pole
number_of_series_turns = v_in_ph / (4.44*f*flux_pp*kw_act);  %number of series turns

number_of_cond_per_slot = number_of_series_turns * 2 /(Q_s/m);   %conductors per slot

%% winding

s_in = p_out /eff_init /pf_init;  %VA, apparent power
area_strand = pi*(dia_strand/2)^2;   %m^2, strand cross sect area
cur_stator = s_in / 3 / v_in_ph;    %Arms, stator current
area_wire_bulk = cur_stator/cur_dens;  %m^2, bir bütün conductor kesit alaný
number_of_strands_per_cond = area_wire_bulk / area_strand;  %number of strands per conductor
number_of_strands_per_slot = number_of_strands_per_cond * number_of_cond_per_slot; 

%% slot sizing

slot_area = number_of_strands_per_slot * area_strand / ff; %m2, slot area
slot_width = slot_pitch - width_teeth;  %m, slot width
slot_depth = slot_area / slot_width;   %m, slot depth


%% stator yoke thickness

depth_yoke = flux_pp / (2*stack*l_act*b_yoke_peak);   %m,depth of yoke

dia_out_stator = dia_in_stator + 2*slot_mount + 2*(slot_depth+depth_yoke); %m
dia_out_rotor = dia_bore-g; %m

%% torque calculation

rotor_vol = pi*(dia_out_rotor/2)^2*l_eff; 
stress = el_load_peak * mag_load_peak * pf_init/2;
torque = 2*stress*rotor_vol;













































