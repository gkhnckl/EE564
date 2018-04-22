% EE 564 - Project 2
% Gökhan Çakal 2332120
% Lamimation 3 is chosen.

%% motor parameters

p = 4; %number of poles
m = 3; %number of phases
s_stator = 36;  %number of stator slots
s_rotor = 34;  %number of rotor slots
area_s = 114.3e-6; %m2, stator slot area
f = 50;  %Hz, rated frequency
g = 0.5e-3;   %m, air gap clearence
b_avg = 0.6;   %T, magnetic loading, avg air gap B

dia_out_s = 200e-3;   %m,  stator outer diameter
dia_in_s = 135e-3;   %m,  stator inner diameter
dia_in_r = 44e-3;   %m,  stator inner diameter
dia_out_r = dia_in_s-2*g;  %m, rotor outer diameter
length = 135e-3; %m, depth of motor
v_in_ll = 380;  %Vrms, line to line input voltage
v_in_ph = v_in_ll/sqrt(3);   %Vrms, phase input voltage


%% N calculation

dia_gap = (dia_in_s+dia_out_r) / 2;  %m, air gap diameter
area_pole = pi*dia_gap*length/p;   %m^2, pole area
flux_pp = b_avg*sqrt(2)*area_pole;  %Wb, peak flux per pole

e_ind = v_in_ph*0.9;  %V, induced voltage for phase
turns_phase = e_ind / (4.44*f*flux_pp);  %number of turns per phase
turns_slot = turns_phase / (s_stator/m/2);  %number of conductors in a stator slot

% kw ekle buraya























