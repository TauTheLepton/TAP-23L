%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                              %
%%    Technika Automatyzacji Procesów - projekt                                 %
%%    Zadanie 1 - linearyzacja                                                  %
%%    Autorzy: Jakub Ostrysz, Mateusz Palczuk, Tymon Kobylecki                  %
%%    Wielkości regulowane: h, T_out                                            %
%%    Wielkości sterujące:  F_h, F_cin                                          %
%%    Zakłócenia: T_d, F_d, T_h, T_c                                            %
%%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables;
close all;

fig_position = [10,10,800,800];

%% parametry zbiornika
global c alpha
c = 0.4;
alpha = 18;

%% punkt pracy (warunki początkowe)
global T_c T_h T_d F_c F_h F_d F_cin Tau_c Tau h T T_out
T_c = 20;
T_h = 72;
T_d = 30;
F_c = 31;
F_cin = F_c;
F_h = 25;
F_d = 15;
Tau_c = 100;
Tau = 120;
h = 15.56;
T = 40.42;
T_out = T;
% ZMIENNYCH GLOBALNYCH NIE ZMIENIAMY

%T_out(t) = T(t - Tau);
%F_c(t) = F_cin(t - Tau_c);

%% okres próbkowania
Tp = 10;
t0 = 0;
tk = 1000;

%% wyznaczenie modelu liniowego w postaci równań stanu
dhdt = (F_h + F_c + F_d - F(h))/(3*c*h^2);
dT_outdt = (F_h*T_h + F_c*T_c + F_d*T_d - F(h)*T_out)/(c*h^3);

dhdh = (3*alpha*sqrt(h) - 4*(F_h + F_c + F_d))/(6*c*h^3);
dhdT_out = 0;
dhdF_h = 1/(3*c*h^2);
dhdF_cin = 1/(3*c*h^2);

dT_outdh = (5*alpha*sqrt(h)*T_out - 6*(F_h*T_h + F_c*T_c + F_d*T_d))/(2*c*h^4);
dT_outdT_out = (-alpha*sqrt(h))/(c*h^3);
dT_outdF_h = T_h/(c*h^3);
dT_outdF_cin = T_c/(c*h^3);

global A B
A = [dhdh dhdT_out; dT_outdh dT_outdT_out];
B = [dhdF_h dhdF_cin; dT_outdF_h dT_outdF_cin];
C = eye(2);
D = zeros(2);

model_lin_ciag = ss(A,B,C,D);

%% rozwiązanie równań różniczkowych zmiana F_h
h_start = h;
T_out_start = T_out;
options = odeset('RelTol',1e-5);

u = [F_h+0; F_cin+0];
[t,x] = ode23(@(t,x) model_ciagly_ode(t,x,u), [t0 tk], [h_start; T_out_start],options);
[tL,xL] = ode23(@(t,x) model_liniowy_ode(t,x,u), [t0 tk], [h_start; T_out_start],options);

u = [F_h+10; F_cin+0];
[t1,x1] = ode23(@(t,x) model_ciagly_ode(t,x,u), [t0 tk], [h_start; T_out_start],options);
[tL1,xL1] = ode23(@(t,x) model_liniowy_ode(t,x,u), [t0 tk], [h_start; T_out_start],options);

u = [F_h+20; F_cin+0];
[t2,x2] = ode23(@(t,x) model_ciagly_ode(t,x,u), [t0 tk], [h_start; T_out_start],options);
[tL2,xL2] = ode23(@(t,x) model_liniowy_ode(t,x,u), [t0 tk], [h_start; T_out_start],options);

u = [F_h-10; F_cin-0];
[t3,x3] = ode23(@(t,x) model_ciagly_ode(t,x,u), [t0 tk], [h_start; T_out_start],options);
[tL3,xL3] = ode23(@(t,x) model_liniowy_ode(t,x,u), [t0 tk], [h_start; T_out_start],options);

u = [F_h-20; F_cin-0];
[t4,x4] = ode23(@(t,x) model_ciagly_ode(t,x,u), [t0 tk], [h_start; T_out_start],options);
[tL4,xL4] = ode23(@(t,x) model_liniowy_ode(t,x,u), [t0 tk], [h_start; T_out_start],options);

color1 = "#0072BD";
color2 = "#D95319";

figure('Position', fig_position);
subplot(211);
title("Skoki F_{H}");
hold on;
plot([0;t],[x(1,1);x(:,1)], 'Color', color1);
plot([0;tL],[xL(1,1);xL(:,1)], 'Color', color2);
plot([0;t1],[x1(1,1);x1(:,1)], 'Color', color1);
plot([0;tL1],[xL1(1,1);xL1(:,1)], 'Color', color2);
plot([0;t2],[x2(1,1);x2(:,1)], 'Color', color1);
plot([0;tL2],[xL2(1,1);xL2(:,1)], 'Color', color2);
plot([0;t3],[x3(1,1);x3(:,1)], 'Color', color1);
plot([0;tL3],[xL3(1,1);xL3(:,1)], 'Color', color2);
plot([0;t4],[x4(1,1);x4(:,1)], 'Color', color1);
plot([0;tL4],[xL4(1,1);xL4(:,1)], 'Color', color2);
legend('h', 'h zlinearyzowane');
xlabel('Czas');
ylabel('h');
xlim([0,tk]);

subplot(212);
title("Skoki F_{H}");
hold on;
plot([0;t+Tau],[x(1,2);x(:,2)], 'Color', color1);
plot([0;tL+Tau],[xL(1,2);xL(:,2)], 'Color', color2);
plot([0;t1+Tau],[x1(1,2);x1(:,2)], 'Color', color1);
plot([0;tL1+Tau],[xL1(1,2);xL1(:,2)], 'Color', color2);
plot([0;t2+Tau],[x2(1,2);x2(:,2)], 'Color', color1);
plot([0;tL2+Tau],[xL2(1,2);xL2(:,2)], 'Color', color2);
plot([0;t3+Tau],[x3(1,2);x3(:,2)], 'Color', color1);
plot([0;tL3+Tau],[xL3(1,2);xL3(:,2)], 'Color', color2);
plot([0;t4+Tau],[x4(1,2);x4(:,2)], 'Color', color1);
plot([0;tL4+Tau],[xL4(1,2);xL4(:,2)], 'Color', color2);
legend('T_{out}', 'T_{out} zlinearyzowane');
xlabel('Czas');
ylabel('T_{out}');
xlim([0,tk]);
exportgraphics(gcf, 'skoki_fh.pdf');

%% rozwiązanie równań różniczkowych zmiana F_cin
h_start = h;
T_out_start = T_out;
options = odeset('RelTol',1e-5);

u = [F_h+0; F_cin+0];
[t,x] = ode23(@(t,x) model_ciagly_ode(t,x,u), [t0 tk], [h_start; T_out_start],options);
[tL,xL] = ode23(@(t,x) model_liniowy_ode(t,x,u), [t0 tk], [h_start; T_out_start],options);

u = [F_h+0; F_cin+10];
[t1,x1] = ode23(@(t,x) model_ciagly_ode(t,x,u), [t0 tk], [h_start; T_out_start],options);
[tL1,xL1] = ode23(@(t,x) model_liniowy_ode(t,x,u), [t0 tk], [h_start; T_out_start],options);

u = [F_h+0; F_cin+20];
[t2,x2] = ode23(@(t,x) model_ciagly_ode(t,x,u), [t0 tk], [h_start; T_out_start],options);
[tL2,xL2] = ode23(@(t,x) model_liniowy_ode(t,x,u), [t0 tk], [h_start; T_out_start],options);

u = [F_h-0; F_cin-10];
[t3,x3] = ode23(@(t,x) model_ciagly_ode(t,x,u), [t0 tk], [h_start; T_out_start],options);
[tL3,xL3] = ode23(@(t,x) model_liniowy_ode(t,x,u), [t0 tk], [h_start; T_out_start],options);

u = [F_h-0; F_cin-20];
[t4,x4] = ode23(@(t,x) model_ciagly_ode(t,x,u), [t0 tk], [h_start; T_out_start],options);
[tL4,xL4] = ode23(@(t,x) model_liniowy_ode(t,x,u), [t0 tk], [h_start; T_out_start],options);

color1 = "#0072BD";
color2 = "#D95319";

figure('Position', fig_position);
subplot(211);
title("Skoki F_{Cin}");
hold on;
plot([0;t+Tau_c],[x(1,1);x(:,1)], 'Color', color1);
plot([0;tL+Tau_c],[xL(1,1);xL(:,1)], 'Color', color2);
plot([0;t1+Tau_c],[x1(1,1);x1(:,1)], 'Color', color1);
plot([0;tL1+Tau_c],[xL1(1,1);xL1(:,1)], 'Color', color2);
plot([0;t2+Tau_c],[x2(1,1);x2(:,1)], 'Color', color1);
plot([0;tL2+Tau_c],[xL2(1,1);xL2(:,1)], 'Color', color2);
plot([0;t3+Tau_c],[x3(1,1);x3(:,1)], 'Color', color1);
plot([0;tL3+Tau_c],[xL3(1,1);xL3(:,1)], 'Color', color2);
plot([0;t4+Tau_c],[x4(1,1);x4(:,1)], 'Color', color1);
plot([0;tL4+Tau_c],[xL4(1,1);xL4(:,1)], 'Color', color2);
legend('h', 'h zlinearyzowane');
xlabel('Czas');
ylabel('h');
xlim([0,tk]);

subplot(212);
title("Skoki F_{Cin}");
hold on;
plot([0;t+Tau+Tau_c],[x(1,2);x(:,2)], 'Color', color1);
plot([0;tL+Tau+Tau_c],[xL(1,2);xL(:,2)], 'Color', color2);
plot([0;t1+Tau+Tau_c],[x1(1,2);x1(:,2)], 'Color', color1);
plot([0;tL1+Tau+Tau_c],[xL1(1,2);xL1(:,2)], 'Color', color2);
plot([0;t2+Tau+Tau_c],[x2(1,2);x2(:,2)], 'Color', color1);
plot([0;tL2+Tau+Tau_c],[xL2(1,2);xL2(:,2)], 'Color', color2);
plot([0;t3+Tau+Tau_c],[x3(1,2);x3(:,2)], 'Color', color1);
plot([0;tL3+Tau+Tau_c],[xL3(1,2);xL3(:,2)], 'Color', color2);
plot([0;t4+Tau+Tau_c],[x4(1,2);x4(:,2)], 'Color', color1);
plot([0;tL4+Tau+Tau_c],[xL4(1,2);xL4(:,2)], 'Color', color2);
legend('T_{out}', 'T_{out} zlinearyzowane');
xlabel('Czas');
ylabel('T_{out}');
xlim([0,tk]);
exportgraphics(gcf, 'skoki_fcin.pdf');

%% model liniowy ciągły w postaci transmitancji (zawiera opóźnienia)
Gs = tf(model_lin_ciag);
Gs.InputDelay = [0 Tau_c];
Gs.OutputDelay = [0 Tau];
% Gs11 = Gs(1,1);     %h/F_h
% Gs12 = Gs(1,2);     %h/F_cin
% Gs21 = Gs(2,1);     %T_out/F_h
% Gs22 = Gs(2,2);     %T_out/F_cin

%% model liniowy dyskretny w postaci równań stanu (nie zawiera opóźnień)
model_lin_dysk = c2d(model_lin_ciag,Tp,'zoh');

%% model liniowy dyskretny w postaci transmitancji (zawiera opóźnienia)
Gz = c2d(Gs,Tp,'zoh'); % używamy transmitancji ciągłej, bo zawiera już opóźnienia

%% wykresy modelu liniowego i liniowego dyskretnego w postaci transmitancji
figure('Position', fig_position);
subplot(221)
hold on;
step(Gs(1,1)); step(Gz(1,1));

subplot(222)
hold on;
step(Gs(1,2)); step(Gz(1,2));

subplot(223)
hold on;
step(Gs(2,1)); step(Gz(2,1));

subplot(224)
hold on;
step(Gs(2,2)); step(Gz(2,2));

%% rozwiązanie równań różniczkowych dla dyskretnego
u = [F_h+10; F_cin+0];
[tL1, xL1] = ode23(@(t,x) model_liniowy_ode(t,x,u), [t0 tk], [h_start; T_out_start], options);
conf = RespConfig('Amplitude', [10, 0], 'Delay', 0);%'InputOffset', [F_h, F_cin], 'InitialState', [h_start, T_out_start]);
[xd1, td1] = step(model_lin_dysk, tk, conf);
xd1(:, 1, 1) = xd1(:, 1, 1) + h;
xd1(:, 2, 1) = xd1(:, 2, 1) + T_out;

u = [F_h+0; F_cin+10];
[tL2, xL2] = ode23(@(t,x) model_liniowy_ode(t,x,u), [t0 tk], [h_start; T_out_start], options);
conf = RespConfig('Amplitude', [0, 10], 'Delay', Tau_c);%'InputOffset', [F_h, F_cin], 'InitialState', [h_start, T_out_start]);
[xd2, td2] = step(model_lin_dysk, tk, conf);
xd2(:, 1, 2) = xd2(:, 1, 2) + h;
xd2(:, 2, 2) = xd2(:, 2, 2) + T_out;

%% wykresy modelu liniowego ciągłego i liniowego dyskretnego w postaci równań stanu
figure('Position', fig_position);
subplot(211);
title("Skoki F_{H}");
hold on;
plot([0; tL1], [xL1(1, 1); xL1(:, 1)], 'Color', color1);
stairs(td1, xd1(:, 1, 1), 'Color', color2);
legend('Model liniowy ciagly', 'Model liniowy dyskretny');
xlabel('Czas');
ylabel('h');
xlim([0,tk]);

subplot(212);
title("Skoki F_{H}");
hold on;
plot([0; tL1+Tau], [xL1(1, 2); xL1(:, 2)], 'Color', color1);
stairs([0; td1+Tau], [xd1(1, 2, 1); xd1(:, 2, 1)], 'Color', color2);
legend('Model liniowy ciagly', 'Model liniowy dyskretny');
xlabel('Czas');
ylabel('T_{out}');
xlim([0,tk]);
exportgraphics(gcf, 'lin_ciag_lin_dysk_fh.pdf');

figure('Position', fig_position);
subplot(211);
title("Skoki F_{Cin}");
hold on;
plot([0; tL2+Tau_c], [xL2(1, 1); xL2(:, 1)], 'Color', color1);
stairs(td2, xd2(:, 1, 2), 'Color', color2);
legend('Model liniowy ciagly', 'Model liniowy dyskretny');
xlabel('Czas');
ylabel('h');
xlim([0, tk]);

subplot(212);
title("Skoki F_{Cin}");
hold on;
plot([0; tL2+Tau+Tau_c], [xL2(1, 2); xL2(:, 2)], 'Color', color1);
stairs([0; td2+Tau], [xd2(1, 2, 2); xd2(:, 2, 2)], 'Color', color2);
legend('Model liniowy ciagly', 'Model liniowy dyskretny');
xlabel('Czas');
ylabel('T_{out}');
xlim([0, tk]);
exportgraphics(gcf, 'lin_ciag_lin_dysk_fcin.pdf');

%% sprawdzanie jakości dyskretyzacji (źle, wzięte są pojedyncze tory, przez co pozostałe nie są wgl liczone)
conf_no_delay = RespConfig('InputOffset',0,'Amplitude',10,'Delay',0);

figure('Position', fig_position);
subplot(211);
conf = RespConfig('InputOffset',0,'Amplitude',10,'Delay',0);
hold on;
step(model_lin_ciag(1,1), conf);
step(Gz(1,1), conf_no_delay);
xlabel('Czas');
ylabel('h');
title('Skok F_{H}');

subplot(212);
conf = RespConfig('InputOffset',0,'Amplitude',10,'Delay',Tau_c);
hold on;
step(model_lin_ciag(1,2), conf);
step(Gz(1,2), conf_no_delay);
xlabel('Czas');
ylabel('h');
title('Skok F_{Cin}');

figure('Position', fig_position);
subplot(211);
conf = RespConfig('InputOffset',0,'Amplitude',10,'Delay',Tau);
hold on;
step(model_lin_ciag(2,1), conf);
step(Gz(2,1), conf_no_delay);
xlabel('Czas');
ylabel('T_{out}');
title('Skok F_{H}');

subplot(212);
conf = RespConfig('InputOffset',0,'Amplitude',10,'Delay',Tau+Tau_c);
hold on;
step(model_lin_ciag(2,2), conf);
step(Gz(2,2), conf_no_delay);
xlabel('Czas');
ylabel('T_{out}');
title('Skok F_{Cin}');

%% porównanie transmitancji dyskretnej i równań stanu zlinearyzowanych dyskretnych
conf_no_delay = RespConfig('Amplitude', [10, 0], 'Delay', 0);
figure('Position', fig_position);
subplot(211);
conf = RespConfig('Amplitude', [10, 0], 'Delay', 0);
[xd, td] = step(model_lin_dysk, tk, conf);
[xg, tg] = step(Gz, tk, conf_no_delay);
xd(:,1,1) = xd(:,1,1) + h;
xg(:,1,1) = xg(:,1,1) + h;
hold on;
stairs(td, xd(:,1,1), 'Color', color1);
stairs(tg, xg(:,1,1), '--', 'Color', color2);
legend('Model liniowy dyskretny', 'Transmitancja dyskretna');
xlabel('Czas');
ylabel('h');
title('Skok F_{H}');

subplot(212);
conf = RespConfig('Amplitude', [10, 0], 'Delay', Tau_c+Tp*2);
[xd, td] = step(model_lin_dysk, tk, conf);
[xg, tg] = step(Gz, tk, conf_no_delay);
xd(:,2,1) = xd(:,2,1) + T_out;
xg(:,2,1) = xg(:,2,1) + T_out;
hold on;
stairs(td, xd(:,2,1), 'Color', color1);
stairs(tg, xg(:,2,1), '--', 'Color', color2);
legend('Model liniowy dyskretny', 'Transmitancja dyskretna');
xlabel('Czas');
ylabel('T_{out}');
title('Skok F_{H}');
exportgraphics(gcf, 'lin_dysk_trans_dysk_fh.pdf');

conf_no_delay = RespConfig('Amplitude', [0, 10], 'Delay', 0);
figure('Position', fig_position);
subplot(211);
conf = RespConfig('Amplitude', [0, 10], 'Delay', Tau-Tp*2);
[xd, td] = step(model_lin_dysk, tk, conf);
[xg, tg] = step(Gz, tk, conf_no_delay);
xd(:,1,2) = xd(:,1,2) + h;
xg(:,1,2) = xg(:,1,2) + h;
hold on;
stairs(td, xd(:,1,2), 'Color', color1);
stairs(tg, xg(:,1,2), '--', 'Color', color2);
legend('Model liniowy dyskretny', 'Transmitancja dyskretna');
xlabel('Czas');
ylabel('h');
title('Skok F_{Cin}');

subplot(212);
conf = RespConfig('Amplitude', [0, 10], 'Delay', Tau+Tau_c);
[xd, td] = step(model_lin_dysk, tk, conf);
[xg, tg] = step(Gz, tk, conf_no_delay);
xd(:,2,2) = xd(:,2,2) + T_out;
xg(:,2,2) = xg(:,2,2) + T_out;
hold on;
stairs(td, xd(:,2,2), 'Color', color1);
stairs(tg, xg(:,2,2), '--', 'Color', color2);
legend('Model liniowy dyskretny', 'Transmitancja dyskretna');
xlabel('Czas');
ylabel('T_{out}');
title('Skok F_{Cin}');
exportgraphics(gcf, 'lin_dysk_trans_dysk_fcin.pdf');

%% funkcja opisująca równania różniczkowe dla modelu ciągłego
function dxdt = model_ciagly_ode(t,x,u)
global T_c T_h T_d F_c F_h F_d c
dxdt(1,1) = (u(1) + u(2) + F_d - F(x(1)))/(3*c*x(1)^2);
dxdt(2,1) = (u(1)*T_h + u(2)*T_c + F_d*T_d - F(x(1))*x(2))/(c*x(1)^3);
end

%% funkcja opisująca równania różniczkowe dla modelu liniowego
function dxdtL = model_liniowy_ode(t,x,u)
global F_h F_cin h T_out A B
dxdtL = A*(x - [h; T_out]) + B*(u - [F_h; F_cin]);
end

%% funkcje zmiennych
function f = F(h)
global alpha;
f = alpha*sqrt(h);
end

function v = V(h)
global c;
v = c*h^3;
end
