%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                        %
%%    Technika Automatyzacji Procesów - projekt                           %
%%    Zadanie 2 - regulacja                                               %
%%    Autorzy: Jakub Ostrysz, Mateusz Palczuk, Tymon Kobylecki            %
%%    Wielkości regulowane: h, T_out                                      %
%%    Wielkości sterujące:  F_h, F_cin                                    %
%%    Zakłócenia: T_d, F_d, T_h, T_c                                      %
%%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%% Modele
model_lin_ciag = modelLinCont();
Gs = gs(model_lin_ciag);
model_lin_dysk = modelLinDisc(model_lin_ciag,Tp);
Gz = gz(Gs,Tp);

ks = tfLimit(Gs,2,2,0)
kz = tfLimit(Gz,2,2,1)

rga_ciag = RGA(ks)
rga_dysk = RGA(kz)

%% PID bez odsprzegania
% classPID(K, Ti, Kd, Td, Tp, Hlim, Llim, Dir, AutoMan, ManVal) 
pid_y2u1 = classPID(4, 10, 1, 100, 10, 100, -100, 1, 1, 0);
pid_y1u2 = classPID(4, 10, 1, 100, 10, 100, -100, 1, 1, 0);

t = 0:Tp:Tp*200;
t = t';
u21 = zeros(size(t));
u12 = zeros(size(t));
len = size(t);
len = len(1);

% [y,tout,x] = lsim(model_lin_dysk(2,1),u,t);

% reTune(obj, K, Ti, Kd, Td)
pid_y2u1.reTune(0.2, 100, 0, 0);
pid_y1u2.reTune(0.2, 100, 0, 0);

delay = 12;

stpt21 = ones(size(t))*17;
stpt21(1:10) = 0;
stpt12 = ones(size(t))*17;
stpt12(1:10) = 0;

lazy_start = 2;  % potrzebne do ustawienia wyjść procesu potrzebnych do PIDa

y21 = lsim(Gz(2,1), u21(1:lazy_start), t(1:lazy_start));
y12 = lsim(Gz(1,2), u12(1:lazy_start), t(1:lazy_start));

for i = lazy_start+1:len
%     u0 = pid_y2u1.calc(y(end), stpt)
%     u = [u; u0];
%     t = [t; t(end)+Tp];
%     [y_new,tout,x] = lsim(model_lin_dysk(2,1), u(end-1-delay:end-delay), t(end-1:end), x(end,:));
%     y = [y; y_new(end)];

    % jestesmy w czasie 'i-1'
    u21(i) = pid_y2u1.calc(y21(end), stpt21(i));
    u12(i) = pid_y1u2.calc(y12(end), stpt12(i));
    % jestesmy w czasie 'i'
    y21 = lsim(Gz(2,1), u21(1:i), t(1:i));
    y12 = lsim(Gz(1,2), u12(1:i), t(1:i));
end

figure;
subplot(211)
hold on
stairs(t,y21)
stairs(t,u21)
stairs(t,stpt21)
legend('y21', 'u21', 'stpt21')
subplot(212)
hold on
stairs(t,y12)
stairs(t,u12)
stairs(t,stpt12)
legend('y12', 'u12', 'stpt12')


% [y_new(end),~,~] = lsim(model_lin_dysk(1,1),u(end),t_new(end),x(end,:));
% [y_new(end),~,~] = lsim(sys                ,u,t_new(end),x(end,:),t_new(end));

% figure;
% hold on;
% stairs(t,u);
% stairs(t_new,y_new);

% for i = 2:len
%     u = pid_y2u1.calc(cv(i-1), stpt);
%     mv(i) = u;
% 
% %     cv(i) = o.calc(mv(i-19));  % TU ZAMIENIĆ NA COŚ MĄDREGO
%     % Zrobić równanie na wyjście od sterowania (może być jakiegoś bardziej
%     % przeszłego, wtedy będzie opóźnienie)
% 
% end

%% Odsprzeganie
Gz22 = Gz(2,2);
Gz22.InputDelay = 0;
Gz22.OutputDelay = 0;
Gz(2,2)
D21 = -Gz(2,1) / Gz22;
D21.InputDelay = 0;
D21.OutputDelay = 0;
D21

zera21 = zero(D21)
bieguny21 = pole(D21)

D12 = -Gz(1,2) / Gz(1,1)

zera12 = zero(D12)
bieguny12 = pole(D12)

[num,den] = tfdata(Gz(2,2))
[z,p,k] = tf2zpk(num{1},den{1})
