clear all
close all
clc

s = tf('s')
G = 40 / ((s+ 5.72) * (s - 1.72))
Ts = 0.02
Kc = 2
C_SS = Kc /s
G_ZOH = 1/ (1 + s*Ts/2)
L = C_SS * G
L11 = L * G_ZOH
Tp = 1.76
Sp = mag2db(1.55)
M_T_LF = mag2db(0.01 / 0.1)
%zeta = 0.46
figure
nichols(L11, 'b')
hold on 
T_grid(Tp)
T_grid(M_T_LF)
S_grid(Sp)

wc_des = 7  %6.3

w_norm = 3
wz = wc_des / w_norm
C_Z = (1  + s/wz)^2
L2 = L11 * C_Z
nichols(L2, 'r')

K = 10 ^(0/10)
L3 = L2 * K
nichols(L3, 'g')

%Closure pole is one pole
wp = 90 / 2.5
C_P = 1 / (1  + s/wp)
L4  = L3 * C_P
nichols(L4, 'm')
%simulation 
C0 = C_SS  * K * C_Z * C_P
Cd = c2d(C0 , Ts, 'Tustin')
%%
%simulation for controlling the max overshoot nad rise time
rho = 1
delta_a = 0
delta_y = 0 
delta_t = 0 
wt = 0 

out = sim("my1sim.slx")
stepinfo(out.y.data, out.y.time, 1, 0,  'RiseTimeLimits', [0 1], 'SettlingTimeTreshold' , 0.05)
%%
%simulation of the ramp input e_r_inf
rho = 2

delta_a = 0
delta_y = 0 
delta_t = 0 
wt = 0 

out = sim("my1sim.slx")

figure
plot(out.e.time, out.e.data, LineWidth=0.5)
%%
%simulation in presence of sinusoidal disturbance
rho = 0
delta_a = 0
delta_y = 0 
delta_t = 0.1
wt = 90

out = sim("my1sim.slx")

figure
plot(out.y.time, out.y.data, LineWidth=0.5)

return 
