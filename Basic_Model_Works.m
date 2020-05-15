%% Grid-Forming VSC with droop control

clear all
clc
close all


T_s=1e-4;     % Simulation time step [s]  Ts highers --> makes simulation run faster, but it may increase the numerical error

% Faults-Events times

T_conn=13;    %Load-step increase time [s]

%1-phase short-circuit
T1sc_on=13;     %1ph Short-circuit time ON [s]

T1sc_off=13.15; %15.25; %1ph Short-circuit time OFF [s]

%phase-phase short-circuit
Tsc_on=21;     %1ph Short-circuit time ON [s]

Tsc_off=21.15; %1ph Short-circuit time OFF [s]


%3-phase short-circuit
T3sc_on=30;     %3ph Short-circuit time ON [s]

T3sc_off=30.15; %3ph Short-circuit time OFF [s]

T_discon=32;    %Loss of Load time [s]

T_en=T_conn-0.5; %Enabling the DC source saturation after the initial synchronization

T_ms=0.001;

Tend=25;      %Endtime of the simulation  [s]

%% Basea values
S_b=100*(10^6);  % Base Power (VA)

V=230*(10^3);    %Base Voltage L-L (V)

f_b=50;          % Base frequency (Hz)

w_ref=2*pi*f_b;  %Base speed (rad/s)
w_b=1;

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

%Migrate and L2EP Values
R_f_pu=0.005;   % Low-pass filter resistance Rf [p.u.]
L_f_pu=0.15;    % Low-pass filter inductance Lf [p.u.]
C_f_pu=0.1; %0.066;  % Low-pass filter capacitance Cf [p.u.]


R_f=R_f_pu*Z1_b;  % Low-pass filter resistance Rf [Ohm]
L_f=L_f_pu*L1_b;  % Low-pass filter inductance Lf [H]
C_f=C_f_pu*C1_b;  % Low-pass filter capacitance Cf [F]


R_dc=(Vdc_n/(0.05*(S_b)/Vdc_n));


%%  transformer parameters
%LV/MV
m=100; %
Pn=210* 1e6;%Nominal power  Pn(VA)
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
%% Transient Virtual Impedance 
wc_HPF=2*2*pi();
XR=3;  % X/R ratio

% %Alborg VALUES
R_v_pu=0.2; %[pu] 
L_v_pu=0.3;  %virtual inductance

R_v=R_v_pu*Z1_b; %virtual resistance 
L_v=L_v_pu*L1_b;  %virtual inductance
%% Overcurrent detection and Virtual impedance Current Limiting

%%Current Limiting Virtual Impedance
I_max=1.2;
I_thresh=1;
I_threshold=I_thresh;% Output current limiting threshold to activate DZ0_vi (pu)
X_R=5;       %DX_vi/DR_vi
Rc_ZVI0 =R1_pu; %transformer resistor pu 
Lc_ZVI0 = L1_pu;%transformer inductor pu (Migrate values)
w0=1;

a=(I_max-I_thresh)^2*(1+(X_R)^2);
b=2*(I_max-I_thresh)*(Rc_ZVI0+w0*Lc_ZVI0*X_R);
c=(Rc_ZVI0)^2+(w0*Lc_ZVI0)^2-(1/I_max)^2;

kp_VI=(-b+sqrt(b^2-4*a*c))/(2*a);

% Virtual resistor
K_RV = 0.09;
W_RV = 16.66;


%% Control parameters

%DC source and governor-turbine time constants
tau_dc=0.05;tau_g=0.25;
%reduce it tau_g  0.5
%time constant models turbine behavior, as they are slow  SM will be more
%like a VSC

% Strategy
w_f=7.5; %cut-off frequency  (Hz)
w_c=w_ref/4.95;
T1=0.109;
T2=0.0182;

%defining SM governer gain----------------------
droop_percentage=1;

% Active power droop parameter 
%m_p=(2*pi*0.5)/(S_b);
%m_p=(50*(1.05)-  50*(0.95))/(n*VSC_Pn); 
m_p=(50*(1.025)-  50*(0.975))/(n*VSC_Pn);
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

% % Voltage loop
Kp_v =0.52;  %Proportional gain voltage control loop
Ki_v =(n)*1.16;%1022; %Integral gain voltage control loop
Kff_v = 1;
Ti_v = Kp_v/Ki_v; 


% Current loop
Kp_i =0.73;%8891;   %Proportional gain current control loop
Ki_i =(1/n)*1.19; %Integral gain current control loop
Kff_i = 1;
Ti_i = Kp_i / Ki_i;


%Saturation
idppu=0.9;
isatpu=1.4;
K_dP=2.3*S_b;

%% Network loading and set-points

%VSC(with Virtual Impedance) and Synchronous machine  IEEE 9-BUS SYSTEM
Ptot_Load=1.5; % base load (p.u.)

Qtot_Load=0.8; % base load (p.u.)

n_gen=2; %number of energy sources

n_load= 1;  %number of loads

ps=0.67*Ptot_Load/n_gen; %set-ponit for each Power Unit [p.u.]

pl=S_b*(Ptot_Load/n_load); %loads in [W]

ql=S_b*(Qtot_Load/n_load); %loads in [VAr]

