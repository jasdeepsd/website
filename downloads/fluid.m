close all;
clc;
clear all;
lpip = (0.01:0.024:0.20)'
% Parameters
g = 9.8; % [m/s^2]
h0 = 0.15; % initial height in [m]
H3 = 0.044; % height from the bottom of the bottle, where the tube in located in [m]
K1 = 164;
u = 1.003 * 10 ^ (-3); % viscosity in [N-s/m^2]
p = 998; % density in [kg/m^3] at room temp
K2 = 1; % for a re-entrant tube
e = 0.0000015; % roughness factor for plastic from table 6.1 in [m]
D = 0.109;
 
% Variables
L = 0.04; % va ries between 0.01 and 0.25 [m]
d = 0.004; % varies between 0.004 and 0.01 [m]; does not have to vary, but it is good to show the involvement of this variable with laminar vs turbulent
 
% Solve for h
[t h] = call_dstate();
figure(1);
hold on;
title('Height Vs time');
xlabel('Time (s)');
ylabel('Height [h] (m)');
grid on;
h = real(h);
plot(t,h,'LineWidth',3);
saveas(1,'h_vs_t','png');
hold off;
 
% ODE
A = ((D^2)/(p*d^4))*(K1*u*d + 64*u*L);
B = (D^4/((p^2)*(d^8)))*((K1*u*d + 64*u*L)^2);
C = 8*(D^4)*((K2 + 1)/(d^4))*g*(H3 - h);
G = (2*(K2+1)*(D^4/d^4));
 
% Ignnoring the Flow acceleration term in Bernoullis equation
A2 = ((D^2)/(p*d^4))*(K1*u*d + 64*u*L);
B2 = (D^4/((p^2)*(d^8)))*((K1*u*d + 64*u*L)^2);
C2 = 8*(D^4)*((K2)/(d^4))*g*(H3 - h);
G2 = (2*(K2)*(D^4/d^4));
 
% Ignnoring the Entry Friction term in Bernoullis equation
%K1 = 0;
%K2 = 0;
A3 = ((D^2)/(p*d^4))*(K1*u*d + 64*u*L);
B3 = (D^4/((p^2)*(d^8)))*((K1*u*d + 64*u*L)^2);
C3 = 8*(D^4)*((K2 + 1)/(d^4))*g*(H3 - h);
G3 = (2*(K2+1)*(D^4/d^4));
 
% Calculate dh/dt
dhdt = (A - sqrt(B - C))/G;
% Velocity
vt = -((D/d)^2)*(dhdt);
figure(2);
hold on;
grid on;
title('Velocity');
xlabel('Time (s)');
ylabel('Velocity [v(t)] (m/s)');
vt = real(vt);
plot(t,vt,'LineWidth',3);
saveas(2,'v_vs_t','png');
hold off;
% Reynolds number
Rey = (p*vt*d/u) ;
figure(3);
hold on;
grid on;
title('Reynolds');
xlabel('Time (s)');
ylabel('Reynolds Number [Re(t)]');
Rey = real(Rey);
plot(t,Rey,'LineWidth',3);
saveas(3,'Re_vs_t','png');
hold off;
 
% Energy analysis
PE = g*(h-H3);
E_flow_acc = (vt.^2)/2;
% Acceleration
ax = ((vt.^2)./2).*(p/L);
%ax = diff(vt);
%t1 = t(1:length(t)-1);
figure(4);
hold on;
grid on;
title('Flow Acceleration Factor');
xlabel('Time (s)');
ylabel('Flow Acceleration Factor [ax(t)] (N/m^3)');
ax = real(ax);
ax_abs = abs(ax);
%ax_abs = p*ax_abs;
plot(t,ax_abs,'LineWidth',3);
%plot(t1,ax_abs,'LineWidth',3);
saveas(4,'a_vs_t','png');
hold off;
% dp/dx
figure(5);
hold on;
title('Pressure Gradient Vs time');
xlabel('Time (s)');
ylabel('dp/dx (N/m^3)');
dPdX = (p*g*(h-H3))/L; 
grid on;
dPdX = real(dPdX);
plot(t,dPdX,'LineWidth',3);
saveas(5,'dPdX_vs_t','png');
hold off;
 
% Entry friction equivalent
figure(6);
hold on;
title('Entry Friction Factor Vs time');
xlabel('Time (s)');
ylabel('Entry Friction Factor (N/m^3)');
F_entry_1 = ((vt.^2)./2);
F_entry_2 = (K1/Rey)';
F_entry_2 = F_entry_2 + K2;
F_entry = F_entry_1 .* F_entry_2;
friction_factor = F_entry .*(p/L);
grid on;
friction_factor = real(friction_factor);
plot(t,friction_factor,'LineWidth',3);
saveas(6,'f_vs_t','png');
hold off;
