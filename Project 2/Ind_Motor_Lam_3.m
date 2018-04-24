% EE 564 - Project 2
% Gökhan Çakal 2332120
% Lamimation 3 is chosen.

clc; clear all; 

%% motor parameters

p = 4; %number of poles
m = 3; %number of phases
s_stator = 36;  %number of stator slots
s_rotor = 34;  %number of rotor slots
area_s = 126e-6; %m2, stator slot area
f = 50;  %Hz, rated frequency
g = 0.35e-3;   %m, air gap clearence
b_avg = 0.35;   %T, magnetic loading, avg air gap B
ff = 0.7;  %fill factor
elec_loading = 24e3;    %A/m, electrical loading

dia_out_s = 200e-3;   %m,  stator outer diameter
dia_in_s = 135e-3;   %m,  stator inner diameter
dia_in_r = 44e-3;   %m,  stator inner diameter
dia_out_r = dia_in_s-2*g;  %m, rotor outer diameter
length = 135e-3; %m, depth of motor
v_in_ll = 380;  %Vrms, line to line input voltage
v_in_ph = v_in_ll/sqrt(3);   %Vrms, phase input voltage

%% winding factor calculation

h = 1:1:7;   %harmonics to be included

cs = 180;   %deg, coil span
kp = sin(h*cs*pi/180/2);  %pitch factor

alpha = 360/s_stator * p/2;   %deg, electrical slot angle
q = s_stator/m/p;   %number of slots per pole per phase
kd = (sin(q*alpha*pi/180*h/2)) ./ (q*sin(h*alpha*pi/180/2));   %dist factor

kw = kp.*kd; %winding factor

%% N calculation

dia_gap = (dia_in_s+dia_out_r) / 2;  %m, air gap diameter
area_pole = pi*dia_gap*length/p;   %m^2, pole area
flux_pp = b_avg*area_pole;  %Wb, peak flux per pole

e_ind = v_in_ph;  %V, induced voltage for phase
turns_phase = e_ind / (4.44*f*flux_pp*kw(1));  %number of turns per phase
turns_slot = turns_phase / (s_stator/m/2)  %number of conductors in a stator slot

%% wire size 

area_wire = area_s*ff/turns_slot;  %m2, wire cross section area
dia_wire = 2*sqrt(area_wire/pi)  %m, wire diameter

%% current calculation

cur = elec_loading*pi*dia_gap/turns_slot/s_stator;  %A, phase current
cur_density = cur/area_wire;   %A/mm2, current density


%% output power

pf = 0.82;   %power factor expectation
eff = 0.9;   %efficiency expectation
s_in = sqrt(3)*v_in_ll*cur;   %VA,  input apperant power
p_out = s_in*pf*eff  %W, output power

%% MMF waveform

ic = 1;  %A, phase A
ia = -0.5;   %A, phase B
ib = -0.5;  %A, phase C

winding_config = [ia ia ia -ic -ic -ic ib ib ib -ia -ia -ia ic ic ic -ib -ib -ib];  % winding configuration
mmf_slot = turns_slot*winding_config;  %mmf of slots

mmf(1) = mmf_slot(1);
for i=2:1:numel(winding_config)
mmf(i) = mmf(i-1)+mmf_slot(i);
end
mmf = mmf-sum(mmf)/numel(mmf);   %average cýkartýldý

elect_angle = 0:20:340;
% bar (elect_angle,mmf);

% xlabel('Electrical Angle (deg)')
% ylabel('MMF (A)')
% title ('MMF Waveform for i_c= 1 A   i_a=i_b= -0.5 A')



%% core losses


kh = 173.3;  %hysteresis
kc = 0.086; %eddy current
ke = 2.068; %excess core losses
f = 50;  %Hz, Frequency
bm = 1.8;  %T, peak flux density
arm_volume = 1.35e-3;  %m3, armature volume

p_core = arm_volume*(kh*f*bm^2 + kc*f^2*bm^2 + ke*(f*bm)^1.5)  %W/m3, core loss calculation


%% mutual inductance 

u0 = 4*pi*10^-7;  %permeability constant

% lm = m*u0*dia_gap*length*turns_phase^2*kw(1)/pi/g/p^2; %H, magneting inductance

lm = 4*m*u0*dia_gap*length*turns_phase^2*kw(1)^2/pi/g/p^2; %H, magneting inductance

xm = 2*pi*f*lm;  %ohm, magnetizing reactance










