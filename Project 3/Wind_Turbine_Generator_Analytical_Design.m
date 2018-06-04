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

temp = 75;  %deg, operation temperature
ro_cu_20 = 1.68e-8;      %ohm*m, resistivity of copper at 20 deg
temp_ref = 20;    %deg, temp at ro_cu is defined
temp_coef = 0.004041;   %temperature dependency of copper resistance


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

%% stator design

Q_s = q_s * m * p;  % stator slot number
dia_in_stator = dia_bore+g;  %m, stator inner diameter
slot_pitch = pi * dia_in_stator / Q_s;  %m, stator slot pitch

%% teeth width calculation

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

%% effective air gap

k = ( width_teeth / g ) / ( 5 + width_teeth / g ); 
k_cs = slot_pitch / ( slot_pitch - k*width_teeth );  % carter's coefficient
g_eff = k_cs * g;  %m, effective air gap ignoring rotor teeths


%% magneting inductance 

u0 = pi*4e-7;  %vacuum permeability
Lm_ph = 2*u0 * dia_in_stator * l_eff * ( kw_1 * number_of_series_turns )^2 / ( pi * (p/2)^2 * g_eff );  %H, magnetizing inductance per phase
Lm = Lm_ph * m/2;  %H, total magnetizing inductance

%% resistance

l_end_winding = ( dia_in_stator + slot_depth ) * pi * 1.1 / p; %m, end winding length
slots_per_phase = Q_s / m;  % slots per phase
l_total_per_turn = 2*( l_end_winding + l_act );  %m, total length per turn 

ro_cu = ro_cu_20 * (1+temp_coef*(temp-temp_ref));     %ohm*m, resistivity of copper at operating temperature

res_per_strand_per_turn = ro_cu * l_total_per_turn / area_strand;  %ohm, resistance per strand per turn

res_per_cond_per_turn = res_per_strand_per_turn / number_of_strands_per_cond;  %ohm

res_ph = number_of_series_turns * res_per_cond_per_turn; %ohm, phase resistance


%% winding factors for harmonics

coil_span = pi;  %rad, full pitch coils

for h = 1:1:15
    
    kp(h) = sin(h*coil_span/2);  %pitch factor is equal to one due to full pitch
    kd(h) = sin(h*q_s*slot_ang/2) / (q_s*sin(h*slot_ang/2));  %dist factor
    kw(h) = kp(h) * kd(h) * mod(h,2);    %winding factor
end


%% efficiency calculation

p_cu_loss_stator = cur_stator^2 * res_ph * 3;  %W, stator copper losses
p_loss_fw = 2e3; %W, friction and windage losses
p_loss_stray = 2.5e3;  %W, stray losses
p_loss_core = 1e3;  %W, iron core losses
p_loss_rotor = 3e3;  %W, rotor ohmic losses
p_loss = p_cu_loss_stator+p_loss_fw+p_loss_stray+p_loss_core+p_loss_rotor;  %W, total losses

eff = p_out / (p_out+p_loss)*100;  %efficiency


%% mass calculation

dens_cu = 8900;   %kg/m3, copper density
dens_core = 7650;     %kg/m3, core density
dia_in_rotor = 100e-3;  %m, rotor inner diameter

vol_stator = pi * l_act * stack * ( (dia_out_stator/2)^2 - (dia_in_stator/2+20e-3)^2 );  %m3, stator volume
vol_rotor = pi * l_act * stack * ( (dia_out_rotor/2)^2 - (dia_in_rotor/2)^2 );  %m3, stator volume
vol_cu = area_strand * l_total_per_turn * number_of_strands_per_cond * number_of_series_turns * m;  %m3, copper volume

mass_stator = dens_core * vol_stator;  %kg, stator mass
mass_rotor = dens_core * vol_rotor;   %kg, rotor mass
mass_cu = dens_cu * vol_cu;   %kg, copper mass

mass_total = mass_stator + mass_rotor + mass_cu;  %kg, total motor mass




























