% Code for the reduced variant of the Maleckar et al. model under control
% conditions using ode15s

% 1.38 is the minimum multiplier to clamped [Na^{+}]_i required to trigger
% Ca^{2+} oscillations under control conditions (implemented in the code)

function MaleckarCaSubsystemControl

%%% y(1) = Ca_i;
%%% y(2) = Ca_d;

% y(3) = O_C;
% y(4) = O_TC;
% y(5) = O_TMgC;
% y(6) = O_TMgMg;
% y(7) = O;

% y(8) = Ca_c;

% y(9) = F1;
% y(10) = F2;

% y(11) = O_Calse;

%%% y(12) = Ca_up;
%%% y(13) = Ca_rel;

%%============ Parameters

% Physical constants

F = 96487;
R = 8314;
T = 306.15;

% Geometry

Vol_i = 0.005884;
Vol_d = 0.00011768;
Vol_rel = 0.0000441;
Vol_up = 0.0003969;
Vol_c = 0.000800224;

% I_NaCa constants

K_NaCa = 0.0374842;
gamma_Na = 0.45;
d_NaCa = 0.0003;

% I_di constants

tau_di = 0.01;

% [Ca^{2+}] buffering

Mg_i = 2.5;

% [Ca^{2+}] handling constants

k_rel_d = 0.003;
k_rel_i = 0.0003;
r_recov = 0.815;
alpha_rel = 200000;
I_up_max = 2800;
k_cyca = 0.0003;
k_xcs = 0.4;
k_srca = 0.5;
tau_tr = 0.01;

Ca_b = 1.8;
tau_Ca = 24.7;

%%============ Clamped ionic concentrations

Nai_clamp = 8.12952418640886*1.38;

Nac_clamp = 130.016896185515;

Vclamp = -72.7393769138908;

% Time discretization and simulation time

tspan = [0;5000];

%===== Initial conditions

Ca_i0 = 6.5e-5;
Ca_d0 = 7.1e-5;

O_C0 = 0.026766;
O_TC0 = 0.012922;
O_TMgC0 = 0.190369;
O_TMgMg0 = 0.714463;
O0 = 1.38222;

Ca_c0 = 1.815768;

F10 = 0.470055;
F20 = 0.002814;

O_Calse0 = 0.431547;
Ca_up0 = 0.649195;
Ca_rel0 = 0.632613;

y0 = [Ca_i0 Ca_d0 O_C0 O_TC0 O_TMgC0 O_TMgMg0 O0 Ca_c0 F10 F20 O_Calse0 Ca_up0 Ca_rel0];

% Load initial conditions (commented)

% load('MaleckarCasubsystem'); % load output of previous simulation saved as yfinal.mat
% y0 = yfinal;

% ODE15S
options = odeset('RelTol',1e-5,'MaxStep',1,'Stats','on'); 
[t,y] = ode15s(@(t,y) f(t,y,F,R,T,Vol_i,Vol_d,Vol_rel,Vol_up,Vol_c,K_NaCa,gamma_Na,d_NaCa,tau_di,Mg_i,k_rel_d,k_rel_i,r_recov,alpha_rel,I_up_max,k_cyca,k_xcs,k_srca,tau_tr,Ca_b,tau_Ca,Nai_clamp,Nac_clamp,Vclamp), tspan, y0, options);

% Save the last values for the variables in the simulation (commented)

% yfinal = y(end,:);
% output = yfinal;
% save 'MaleckarCasubsystem' 'yfinal'

% Plot

figure(1)
plot(t,y(:,1),'LineWidth',5)
ylabel('\fontsize{45} [Ca^{2+}]_i (mM)')
xlabel('\fontsize{45} Time(ms)')
title('\fontsize{45} [Ca^{2+}]_i vs time')
set(gca,'fontsize',45)

figure(2)
plot(t,y(:,2),'LineWidth',5)
ylabel('\fontsize{45} [Ca^{2+}]_d (mM)')
xlabel('\fontsize{45} Time(ms)')
title('\fontsize{45} [Ca^{2+}]_d vs time')
set(gca,'fontsize',45)

figure(12)
plot(t,y(:,12),'LineWidth',5)
ylabel('\fontsize{45} [Ca^{2+}]_u_p (mM)')
xlabel('\fontsize{45} Time(ms)')
title('\fontsize{45} [Ca^{2+}]_u_p vs time')
set(gca,'fontsize',45)

figure(4)
plot(t,y(:,13),'LineWidth',5)
ylabel('\fontsize{45} [Ca^{2+}]_r_e_l (mM)')
xlabel('\fontsize{45} Time(ms)')
title('\fontsize{45} [Ca^{2+}]_r_e_l vs time')
set(gca,'fontsize',45)

% Save the simulation time-course of the variables (commented)

% save MaleckarCasubsystemdata.mat y t;  % Save all vectors in y matrix. 

function ydot = f(t,y,F,R,T,Vol_i,Vol_d,Vol_rel,Vol_up,Vol_c,K_NaCa,gamma_Na,d_NaCa,tau_di,Mg_i,k_rel_d,k_rel_i,r_recov,alpha_rel,I_up_max,k_cyca,k_xcs,k_srca,tau_tr,Ca_b,tau_Ca,Nai_clamp,Nac_clamp,Vclamp)
ydot = zeros(size(y));

% I_NaCa
I_NaCa = (K_NaCa*((Nai_clamp^3)*y(8)*(exp((F*Vclamp*gamma_Na)/(R*T)))-(Nac_clamp^3)*y(1)*(exp(((gamma_Na-1)*Vclamp*F)/(R*T)))))/(1+d_NaCa*((Nac_clamp^3)*y(1)+(Nai_clamp^3)*y(8)));

% Intracellular ion concentration flux

I_di = ((y(2)-y(1))*2*Vol_d*F)/tau_di;

% Intracellular [Ca^{2+}] buffering

J_O_C = 200000*y(1)*(1-y(3))-(476*y(3));
J_O_TC = 78400*y(1)*(1-y(4))-(392*y(4));
J_O_TMgC = 200000*y(1)*(1-y(5)-y(6))-(6.6*y(5));
J_O_TMgMg = 2000*Mg_i*(1-y(5)-y(6))-(666*y(6));

J_O = 0.08*J_O_TC + 0.16*J_O_TMgC + 0.045*J_O_C;

% [Ca^{2+}] handling by the SR

r_Ca_d_term = y(2)/(y(2)+k_rel_d);
r_Ca_i_term = y(1)/(y(1)+k_rel_i);
r_Ca_d_factor = (r_Ca_d_term)^4;
r_Ca_i_factor = (r_Ca_i_term)^4;
r_act = 203.8*(r_Ca_i_factor+r_Ca_d_factor);
r_inact = 33.96+339.6*r_Ca_i_factor;

I_rel_f2 = y(10)/(y(10)+0.25);
I_rel_factor = (I_rel_f2)^2;
I_rel = alpha_rel*I_rel_factor*(y(13)-y(1));

I_up = (I_up_max*((y(1)/k_cyca)-(((k_xcs^2)*y(12))/k_srca)))/(((y(1)+k_cyca)/k_cyca)+((k_xcs*(y(12)+k_srca))/k_srca));
I_tr = ((y(12)-y(13))*2*Vol_rel*F)/tau_tr;

J_O_Calse = 480*y(13)*(1-y(11))-(400*y(11));

%%======= ODEs

ydot(1) = ((-(I_up-(I_di+I_rel+2*I_NaCa)))/(2*Vol_i*F))-(J_O);
ydot(2) = (-(I_di))/(2*Vol_d*F);

ydot(3) = J_O_C;
ydot(4) = J_O_TC;
ydot(5) = J_O_TMgC;
ydot(6) = J_O_TMgMg;
ydot(7) = J_O;

ydot(8) = ((Ca_b-y(8))/tau_Ca) + ((-(2*I_NaCa))/(2*Vol_c*F));

ydot(9) = r_recov*(1-y(9)-y(10))-(r_act*y(9));
ydot(10) = r_act*y(9)-(r_inact*y(10));

ydot(11) = J_O_Calse;

ydot(12) = ((I_up-I_tr)/(2*Vol_up*F));
ydot(13) = ((I_tr-I_rel)/(2*Vol_rel*F)) - (31*J_O_Calse);
