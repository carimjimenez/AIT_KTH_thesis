%Grid-Forming VSC with droop control

clear all
clc
close all

T_s=1e-4;     % Simulation time step [s]  Ts highers --> makes simulation run faster, but it may increase the numerical error

Tend=35;      %Endtime of the simulation  [s]

T_loss=Tend;

T_load=12;    %Load-step time [s]

Tsc_on=20;    %Short-circuit time ON [s]

Tsc_off=20.25;  %Short-circuit time OFF [s]

T_en=T_load-0.5;%Enabling the DC source saturation after the initial synchronization

T_ms=0.001;


%% Basea values
S_b=100*(10^6);  % Base Power (VA)

V=230*(10^3);    %Base Voltage L-L (V)

f_b=50;          % Base frequency (Hz)

w_ref=2*pi*f_b;  %Base speed (rad/s)


% High voltage Zone
P_b=S_b;          %Base power (W)

Q_b=S_b;
V_b=V;            %Base voltage (V)

I_b=S_b/(sqrt(3)*V_b); %Base current (A)

Z_b=(V_b^2)/S_b;    %Base impedance (Ohm)

L_b=Z_b/w_ref;     %Base inductance (H) 

C_b=1/(w_ref*Z_b); %Base capacitance (F)

% Medium voltage Zone
V2_rms=13800;
I2_b=S_b/(sqrt(3)*V2_rms);
Z2_b=(V2_rms^2)/S_b;
L2_b=Z2_b/w_ref;
C2_b=1/(w_ref*Z2_b);

% Low voltage Zone
V1_rms=1000;
I1_b=S_b/(sqrt(3)*V1_rms);
Z1_b=(V1_rms^2)/S_b;
L1_b=Z1_b/w_ref;
C1_b=1/(w_ref*Z1_b);

%% Branch parameters 
R_lines=[0,0.017,0.039,0,0.0119,0.0085,0,0.032,0.01]*Z_b;  %Resistance Value of each Line (Ohm)

L_lines=[0.0576,0.092,0.17,0.0586,0.1008,0.072,0.0625,0.161,0.085]*L_b;  %Inductance Value of each Line (Ohm)

R_14=R_lines(1); L_14=L_lines(1);
R_45=R_lines(2); L_45=L_lines(2);
R_56=R_lines(3); L_56=L_lines(3);
R_36=R_lines(4); L_36=L_lines(4);
R_67=R_lines(5); L_67=L_lines(5);
R_78=R_lines(6); L_78=L_lines(6);
R_82=R_lines(7); L_82=L_lines(7);
R_89=R_lines(8); L_89=L_lines(8);
R_94=R_lines(9); L_94=L_lines(9);

%% Shunt parameters ([pu]*base)
C_4=0.1670*C_b; %Capacitance Value of each shunt (F)

C_6=0.2835*C_b;

C_8=0.2275*C_b;

%% VSC Parameters
VSC_Pn = 1* 1e6;            % Nom. VSC power (W)  

V1_rms=1000;      % Nom. phase-to-phase VSC voltage (V)

V_m=sqrt(2/3)*V1_rms;  % peak phase-to-neutral VSC voltage (V)

I_b_LV=S_b/(sqrt(3)*V1_rms);

Vdc_n=3*V_m;     % Nom. DC voltage (V)

n=200;             %number of parallel-inverters 

C_dc=0.008*(n);

%Filter Parameters

R_f=0.001/n;
L_f=(1/n)*200*10^-6;
C_f=n*300*10^-6;

R_f_pu=0.005;   % Low-pass filter resistance Rf [p.u.]
L_f_pu=0.15;    % Low-pass filter inductance Lf [p.u.]
C_f_pu=0.066;   % Low-pass filter capacitance Cf [p.u.]


R_f=R_f_pu*Z1_b;  % Low-pass filter resistance Rf [Ohm]
L_f=L_f_pu*L1_b;  % Low-pass filter inductance Lf [H]
C_f=C_f_pu*C1_b;  % Low-pass filter capacitance Cf [F]


R_dc=(Vdc_n/(0.05*(S_b)/Vdc_n));


%%  transformer parameters
%LV/MV
m=100;
Pn=210* 1e6;%Nominal power and frequency  [ Pn(VA)
V2_rms=13800;     %Medium voltage side (v)

R1_pu=1*0.00734/m;
L1_pu=1*0.0186/m;
R2_pu=R1_pu;
L2_pu=L1_pu;

Rm_pu=10*347.82/m; %Magnetization resistance  Rm (pu)
Lm_pu=10*156.68/m; %Magnetization inductance  Lm (pu)

%MV/HV
R1_mh=0.0027;  %resistance primary side [p.u.]
L1_mh=0.08 ;    %%inductance primary side [p.u.]

R2_mh=0.0027;  %resistance secondary side [p.u.]
L2_mh=0.08 ;    %%inductance secondary side [p.u.]

Rm_mh=500;  %Magnetization resistance  Rm (pu)
Lm_mh=500;  %Magnetization inductance  Lm (pu)
%% Virtual Impedance 
% R_v=1*10^-6;
% L_v=1*10^-6;

XR=3;  % X/R ratio
R_v=1*10^-5; %virtual resistance 
L_v=(XR*R_v)/(2*pi*f_b);  %virtual inductance

%or
% by using the following values, the behaviour will not be adequate
%R_vpu=0.1;  
%L_vpu=0.3;
%R_v=R_vpu*Z1_b;
%L_v=L_vpu*L1_b;
%% Control parameters

%DC source and governor-turbine time constants
tau_dc=0.05;tau_g=5;

% Strategy
w_f=7.5; %cut-off frequency  (Hz)
w_c=w_ref/4.95;

%defining SM governer gain----------------------
droop_percentage=5;

% Active power droop parameter 
%m_p=(2*pi*0.5)/(S_b);
%m_p=(50*(1.05)-  50*(0.95))/(2*S_b); 
m_p=(50*(1.025)-  50*(0.975))/(2*S_b);
%m_p=((10/100)*f_b)/(0.75*S_b);

% grid-forming converter control----------------
I_b_dc=S_b/Vdc_n;

i_loss_dc=Vdc_n/R_dc;

i_ul=1.15*(S_b/Vdc_n)+i_loss_dc;%dc source saturation limits

i_ll=-1.15*(S_b/Vdc_n)-i_loss_dc;%dc source saturation limits


% DC voltage control
eta_1= w_ref/Vdc_n;
k_dc=eta_1/(Vdc_n*m_p);
K_p=(1/Vdc_n);
K_r=1/R_dc;


% AC voltage control
%ki_v_ac=2*0.25;
ki_v_ac=10;    %Integral gain AC voltage control 
%kp_v_ac=0.001;
kp_v_ac=0.05;  %Proportional gain AC voltage control 

% Voltage loop
Kp_v =0.52;  %Proportional gain voltage control loop
Ki_v =(n)*1.161022; %Integral gain voltage control loop
Kff_v = 1;
Ti_v = Kp_v/Ki_v; 


% Current loop
Kp_i =0.738891;   %Proportional gain current control loop
Ki_i =(1/n)*1.19; %Integral gain current control loop
Kff_i = 1;
Ti_i = Kp_i / Ki_i;

%% Network loading and set-points
base=2.25; % base load
load_change=0.75;% load disturbance
ps=base/3; %set-ponit in [MW]
pl=S_b*ps; %loads in [W]
load_step=S_b*load_change; %disturbance in [W]
