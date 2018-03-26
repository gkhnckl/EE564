%% EE 564 Project 1 - Inductance and Transformer Modeling

% Gökhan Çakal 
% Öðrenci No: 2332120

clc; clear; close all;
%% Q1 - Inductor Design

% core selected = 0077739A7 u=90 Mag-Inc Kool Mu Toroid Core

or = 74.1e-3/2;  %m, outer radius
ir = 45.3e-3/2;  %m, inner radius
height = 35e-3;    %m, height
cross_sect = 497e-6;   %m^2, cross section area
path_eff = 184e-3;  %m, effective path
u0 = pi*4e-7;   %permeability constant
ur = 83.82;    % relative permeability for chosen operation point
n = 57;    % number of turns for chosen operation point


% BH Curve

a=4.182e-2;
b=2.99e-2;
c=7.826e-4;
d=6.542e-2;
e=7.669e-4;
x=1.569;

H = 0:1:3000;   %A*T/cm,  magnetic intensity
B = ((a+b*H+c*H.^2) ./ (1+d*H+e*H.^2)).^x;  %T, magnetic field density
plot(H,B);
title('BH Curve')
xlabel('H(A*T/cm)')
ylabel('B(T)')

% H vs Initial Relative Permeability

a1=0.01;
b1=2.698e-5;
c1=1.558;

initial_perm = 1./(a1+b1*(H).^c1);
figure;
plot(H,initial_perm);
title('Initial Permeability vs DC Bias')
xlabel('H(A*T/cm)')
ylabel('Initial Relative Permeability')


%% Part A, Homogeneous linear core

dia_mean = or+ir;  %m, mean diameter
path_mean = pi*dia_mean; %m, mean path
reluct_hl = path_mean/(ur*u0*cross_sect);  %1/H, reluctance for homogenous linear core
ind_hl = n^2/reluct_hl;   % H, inductance for hom. linear core



%% Part A, Non-homogeneous linear core

num_path = 4e4;
delta_r = (or-ir)/(num_path-1);     %m, radial length of each part
delta_cs = delta_r*height;   %m^2, cross section area of each part
radi_delta = linspace(ir,or,num_path);   %m, radius of each path
path_each = 2*pi*radi_delta;    %m, length of the each path
reluct_each = path_each/(ur*u0*delta_cs);  %1/H, reluctance of each path
reluct_nl = (sum(1./reluct_each))^-1;   %1/H, equivalent reluctance for non homg, linear core
ind_nl = n^2/reluct_nl;    %H, inductance for non-homoge. linear core

%% Part A, Homogeneous non-linear core

cur = 7.5;   %A, current
ur_hn = 73.97;  % relative permeability for homogenous non linear case
reluct_hn = path_mean/(ur_hn*u0*cross_sect);     %1/H, reluctance 
ind_hn = n^2/reluct_hn;    %H, inductance

%% Part A, Non-homogeneous non-linear core

ur_nn=ur_hn;   %relative permeability
num_path_nn = 4e4;
delta_r_nn = (or-ir)/(num_path_nn-1);     %m, radial length of each part
delta_cs_nn = delta_r_nn*height;   %m^2, cross section area of each part
radi_delta = linspace(ir,or,num_path_nn);   %m, radius of each path
path_each = 2*pi*radi_delta;    %m, length of the each path
reluct_each_nn = path_each/(ur_nn*u0*delta_cs_nn);  %1/H, reluctance of each path
reluct_nn = (sum(1./reluct_each_nn))^-1;   %1/H, equivalent reluctance for non homg, linear core
ind_nn = n^2/reluct_nn;    %H, inductance for non-homoge. linear core

%% Part A, Air gapped homogeneous linear core

g = 2e-3;  %m, air gap distance
reluct_gap = g/(u0*cross_sect);   %1/H,  reluctance of gap ignoring fringing
reluct_core = (path_mean-g)/(ur*u0*cross_sect);  %1/H, reluct of core
reluct_total = reluct_gap+reluct_core;   %1/H, total reluctance of system
ind_gapped = n^2/reluct_total;   %H, inductance of gapped linear core

%% Part A, Air gapped homogeneous linear core with fringing effect

gap_cross_sec = (or-ir+g)*(height+g);  %m^2 gap cross section area
reluct_gap_fr = g/(u0*gap_cross_sec);  %1/H,  reluctance of gap with fringing
reluct_total_fr = reluct_gap_fr+reluct_core;   %1/H, total reluctance of system
ind_gapped_fr = n^2/reluct_total_fr;   %H, inductance of gapped linear core








