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

SAVE_FIG = false;

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

F_d = F_d + 5;
model_lin_ciag2 = modelLinCont();
Gs2 = gs(model_lin_ciag2);
model_lin_dysk2 = modelLinDisc(model_lin_ciag2,Tp);
Gz2 = gz(Gs2,Tp);

ks = tfLimit(Gs,2,2,0);
kz = tfLimit(Gz,2,2,1);

rga_ciag = RGA(ks);
rga_dysk = RGA(kz);

%% PID bez odsprzegania
% classPID(K, Ti, Kd, Td, Tp, Hlim, Llim, Dir, AutoMan, ManVal) 
pid_y2u1 = classPID(4, 10, 1, 100, 10, 100, -100, 1, 1, 0);
pid_y1u2 = classPID(4, 10, 1, 100, 10, 100, -100, 1, 1, 0);

t = 0:Tp:Tp*300;
t = t';
len = size(t);
len = len(1);
u = zeros(len,2);

% reTune(obj, K, Ti, Kd, Td)
pid_y2u1.reTune(0.2, 100, 0, 0);
pid_y1u2.reTune(0.2, 100, 0, 0);

% X = [0.5, 180, 0, 50 1.7833 161.7775 0 8.5290];
% pid_y2u1.reTune(X(1), X(2), X(3), X(4));
% pid_y1u2.reTune(X(5), X(6), X(7), X(8));

pid_y2u1.reTune(0.1, 120, 0, 50);
pid_y1u2.reTune(1.5, 400, 0, 100); % niezle

stpt = zeros(len,2);
stpt(20:end,1) = 10;
stpt(40:end,2) = 5;

F_d_t = zeros(len,1);
F_d_t(1:150) = 15;
F_d_t(151:end) = 20;


lazy_start = 2;  % potrzebne do ustawienia wyjść procesu potrzebnych do PIDa

y = lsim(Gz, u(1:lazy_start,:), t(1:lazy_start));

for i = lazy_start+1:len
    if i == 150
        Gz = Gz2;
        y2 = y;
    end
%     u0 = pid_y2u1.calc(y(end), stpt)
%     u = [u; u0];
%     t = [t; t(end)+Tp];
%     [y_new,tout,x] = lsim(model_lin_dysk(2,1), u(end-1-delay:end-delay), t(end-1:end), x(end,:));
%     y = [y; y_new(end)];

    % jestesmy w czasie 'i-1'
    u(i,1) = pid_y2u1.calc(y(end,2), stpt(i,2));
    u(i,2) = pid_y1u2.calc(y(end,1), stpt(i,1));
    % jestesmy w czasie 'i'
    y = lsim(Gz, u(1:i,:), t(1:i));
end

figure('Position', fig_position);
subplot(211)
xlabel('t [s]')
ylabel('F_H [cm^3/s]')
hold on
stairs(t,u(:,1)+F_h)
legend('Wartość sterowania','Location','southeast')

subplot(212)
hold on
xlabel('t [s]')
ylabel('T_{out} [\circC]')
stairs(t,y(:,2)+T_out)
% ytmp = y;
% ytmp(1:149,:) = y2;
% stairs(t,ytmp(:,2)+T_out)
stairs(t,stpt(:,2)+T_out)
stairs(t, F_d_t+25, '--')
legend('Wartość wyjściowa', 'Wartość zadana','Zakłócenie','Location','southeast')
if SAVE_FIG; saveas(gcf, 'y2u1_zak.png'); end

figure('Position', fig_position);
subplot(211)
xlabel('t [s]')
ylabel('F_{Cin} [cm^3/s]')
hold on
stairs(t,u(:,2)+F_cin)
legend('Wartość sterowania','Location','southeast')

subplot(212)
xlabel('t [s]')
ylabel('h [cm]')
hold on
stairs(t,y(:,1)+h)
stairs(t,stpt(:,1)+h)
stairs(t, F_d_t, '--')
legend('Wartość wyjściowa', 'Wartość zadana','Zakłócenie' ,'Location','southeast')
if SAVE_FIG; saveas(gcf, 'y1u2_zak.png'); end


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
% D12 tutaj to w notatkach D22
% D21 tutaj to w notatkach D11
Gz = gz(Gs,Tp);
Gz22 = Gz(2,2);
Gz22.InputDelay = 0;
Gz22.OutputDelay = 0;
Gz(2,2)
D21 = -Gz(2,1) / Gz22;
% trzeba wyzerować opóźnienia, bo z dzielenia wyżej by wyszło
% przyspieszenie, co jest nierealizowalne, więc trzeba odrzucić
D21.InputDelay = 0;
D21.OutputDelay = 0;
D21 = stabilizeTf(D21);
% pole(D21)
% zero(D21)

D12 = -Gz(1,2) / Gz(1,1);
D12 = stabilizeTf(D12);
% pole(D12)
% zero(D12)

% NOWE
Gz12 = Gz(1,2);
Gz12.InputDelay = 0;
Gz12.OutputDelay = 0;
D11 = -Gz(1,1) / Gz12;
D11 = stabilizeTf(D11);

Gz21 = Gz(2,1);
Gz22 = Gz(2,2);
Gz22.InputDelay = Gz22.InputDelay - Gz21.InputDelay;
Gz22.OutputDelay = Gz22.OutputDelay - Gz21.OutputDelay;
Gz21.InputDelay = 0;
Gz21.OutputDelay = 0;
D22 = -Gz22 / Gz21;
D22 = stabilizeTf(D22);

D12 = D22;
D21 = D11;

% classPID(K, Ti, Kd, Td, Tp, Hlim, Llim, Dir, AutoMan, ManVal) 
pid_y2u1 = classPID(4, 10, 1, 100, 10, 100, -100, 1, 1, 0);
pid_y1u2 = classPID(4, 10, 1, 100, 10, 100, -100, 1, 1, 0);

t = 0:Tp:Tp*300;
t = t';
len = size(t);
len = len(1);
u = zeros(len,2);

%% dostrajanie algorytmem genetycznym
% options = gaoptimset('StallGenLimit', 10, 'PopulationSize', 100);
% [X] = ga(@PID_func2,8,[],[],[],[],[0,160,0,40,1.7833,161.7775,0,8.5290],[1,200,0,60,1.7833,161.7775,0,8.5290],[],[],options);
% fprintf('\nK=%f; Ti=%f; Kd=%f; Td=%f;\n K=%f; Ti=%f; Kd=%f; Td=%f;\n', X)

% X = [0.5855 160.0047 0 40.5876 1.7833 161.7775 0 8.5290]; %better
% X = [0.5999  160.0072  0 52.9825 1.8791 152.4977 0 13.6712]; min AE
% X = [0.5861 160.0007 0 53.0921 1.7833 161.7775 0 8.5290]; % best y1u2

% X = [0.5, 180, 0, 50 1.7833 161.7775 0 8.5290]; best
% pid_y2u1.reTune(X(1), X(2), X(3), X(4));
% pid_y1u2.reTune(X(5), X(6), X(7), X(8));

%% dostrajanie ręczne 
% reTune(obj, K, Ti, Kd, Td)
pid_y2u1.reTune(0.5, 180, 0, 50);
pid_y1u2.reTune(2.5, 150, 0, 5); % niezle

stpt = zeros(len,2);
stpt(20:end,1) = 10;
stpt(40:end,2) = 5;

lazy_start = 2;  % potrzebne do ustawienia wyjść procesu potrzebnych do PIDa

y = lsim(Gz, u(1:lazy_start,:), t(1:lazy_start));

for i = lazy_start+1:len
    % jestesmy w czasie 'i-1'
    ur(i,1) = pid_y2u1.calc(y(end,2), stpt(i,2));
    ur(i,2) = pid_y1u2.calc(y(end,1), stpt(i,1));
    u1 = lsim(D12, ur(1:i,2), t(1:i)) + ur(1:i,1);
    u2 = lsim(D21, ur(1:i,1), t(1:i)) + ur(1:i,2);
    % jestesmy w czasie 'i'
    y = lsim(Gz, [u1, u2], t(1:i));
end

% figure('Position', fig_position);
% subplot(211)
% hold on
% stairs(t,y(:,2))
% stairs(t,u1)
% stairs(t,stpt(:,2))
% legend('y2', 'u1', 'stpt2')
% 
% subplot(212)
% hold on
% stairs(t,y(:,1))
% stairs(t,u2)
% stairs(t,stpt(:,1))
% legend('y1', 'u2', 'stpt1')

figure('Position', fig_position);
subplot(211)
xlabel('t [s]')
ylabel('F_H [cm^3/s]')
hold on
stairs(t,u1+F_h)
legend('Wartość sterowania','Location','southeast')

subplot(212)
hold on
xlabel('t [s]')
ylabel('T_{out} [\circC]')
stairs(t,y(:,2)+T_out)
stairs(t,stpt(:,2)+T_out)
legend('Wartość wyjściowa', 'Wartość zadana','Location','southeast')
if SAVE_FIG; saveas(gcf, 'y2u1_odsprz.png'); end

figure('Position', fig_position);
subplot(211)
xlabel('t [s]')
ylabel('F_{Cin} [cm^3/s]')
hold on
stairs(t,u2+F_cin)
legend('Wartość sterowania','Location','southeast')

subplot(212)
xlabel('t [s]')
ylabel('h [cm]')
hold on
stairs(t,y(:,1)+h)
stairs(t,stpt(:,1)+h)
legend('Wartość wyjściowa', 'Wartość zadana','Location','southeast')
if SAVE_FIG saveas(gcf, 'y1u2_odsprz.png'); end

%% Zmiana zakłócenia z odsprzęganiem 

% D12 tutaj to w notatkach D22
% D21 tutaj to w notatkach D11
Gz = gz(Gs,Tp);
Gz22 = Gz(2,2);
Gz22.InputDelay = 0;
Gz22.OutputDelay = 0;
Gz(2,2)
D21 = -Gz(2,1) / Gz22;
% trzeba wyzerować opóźnienia, bo z dzielenia wyżej by wyszło
% przyspieszenie, co jest nierealizowalne, więc trzeba odrzucić
D21.InputDelay = 0;
D21.OutputDelay = 0;
D21 = stabilizeTf(D21);
% pole(D21)
% zero(D21)

D12 = -Gz(1,2) / Gz(1,1);
D12 = stabilizeTf(D12);
% pole(D12)
% zero(D12)

% NOWE
Gz12 = Gz(1,2);
Gz12.InputDelay = 0;
Gz12.OutputDelay = 0;
D11 = -Gz(1,1) / Gz12;
D11 = stabilizeTf(D11);

Gz21 = Gz(2,1);
Gz22 = Gz(2,2);
Gz22.InputDelay = Gz22.InputDelay - Gz21.InputDelay;
Gz22.OutputDelay = Gz22.OutputDelay - Gz21.OutputDelay;
Gz21.InputDelay = 0;
Gz21.OutputDelay = 0;
D22 = -Gz22 / Gz21;
D22 = stabilizeTf(D22);

D12 = D22;
D21 = D11;

% classPID(K, Ti, Kd, Td, Tp, Hlim, Llim, Dir, AutoMan, ManVal) 
pid_y2u1 = classPID(4, 10, 1, 100, 10, 100, -100, 1, 1, 0);
pid_y1u2 = classPID(4, 10, 1, 100, 10, 100, -100, 1, 1, 0);

t = 0:Tp:Tp*300;
t = t';
len = size(t);
len = len(1);
u = zeros(len,2);

%% dostrajanie algorytmem genetycznym
% options = gaoptimset('StallGenLimit', 10, 'PopulationSize', 100);
% [X] = ga(@PID_func,8,[],[],[],[],[0,160,0,45,1,130,0,0],[1,200,0,55,4,170,0,15],[],[],options);
% fprintf('\nK=%f; Ti=%f; Kd=%f; Td=%f;\n K=%f; Ti=%f; Kd=%f; Td=%f;\n', X)

% X = [0.5855 160.0047 0 40.5876 1.7833 161.7775 0 8.5290]; %better
% X = [0.5999  160.0072  0 52.9825 1.8791 152.4977 0 13.6712]; min AE
% X = [0.5861 160.0007 0 53.0921 1.7833 161.7775 0 8.5290]; % best y1u2

% X = [0.5, 180, 0, 50 1.7833 161.7775 0 8.5290]; best
% pid_y2u1.reTune(X(1), X(2), X(3), X(4));
% pid_y1u2.reTune(X(5), X(6), X(7), X(8));

%% dostrajanie ręczne 
pid_y2u1.reTune(0.5, 180, 0, 50);
pid_y1u2.reTune(2.5, 150, 0, 5); % niezle

stpt = zeros(len,2);
stpt(20:end,1) = 10;
stpt(40:end,2) = 5;

lazy_start = 2;  % potrzebne do ustawienia wyjść procesu potrzebnych do PIDa

y = lsim(Gz, u(1:lazy_start,:), t(1:lazy_start));

for i = lazy_start+1:len
    if i == 150
        Gz = Gz2;
    end
    % jestesmy w czasie 'i-1'
    ur(i,1) = pid_y2u1.calc(y(end,2), stpt(i,2));
    ur(i,2) = pid_y1u2.calc(y(end,1), stpt(i,1));
    u1 = lsim(D12, ur(1:i,2), t(1:i)) + ur(1:i,1);
    u2 = lsim(D21, ur(1:i,1), t(1:i)) + ur(1:i,2);
    % jestesmy w czasie 'i'
    y = lsim(Gz, [u1, u2], t(1:i));
end

%% DMC analityczny



%% plots
% figure('Position', fig_position);
% subplot(211)
% hold on
% stairs(t,y(:,2))
% stairs(t,u1)
% stairs(t,stpt(:,2))
% legend('y2', 'u1', 'stpt2')
% 
% subplot(212)
% hold on
% stairs(t,y(:,1))
% stairs(t,u2)
% stairs(t,stpt(:,1))
% legend('y1', 'u2', 'stpt1')

figure('Position', fig_position);
subplot(211)
xlabel('t [s]')
ylabel('F_H [cm^3/s]')
hold on
stairs(t,u1+F_h)
legend('Wartość sterowania','Location','southeast')

subplot(212)
hold on
xlabel('t [s]')
ylabel('T_{out} [\circC]')
stairs(t,y(:,2)+T_out)
stairs(t,stpt(:,2)+T_out)
stairs(t, F_d_t+25, '--')
legend('Wartość wyjściowa', 'Wartość zadana','Zakłócenie' ,'Location','southeast')
if SAVE_FIG; saveas(gcf, 'y2u1_odsprz_zak.png'); end

figure('Position', fig_position);
subplot(211)
xlabel('t [s]')
ylabel('F_{Cin} [cm^3/s]')
hold on
stairs(t,u2+F_cin)
legend('Wartość sterowania','Location','southeast')

subplot(212)
xlabel('t [s]')
ylabel('h [cm]')
hold on
stairs(t,y(:,1)+h)
stairs(t,stpt(:,1)+h)
stairs(t, F_d_t, '--')
legend('Wartość wyjściowa', 'Wartość zadana','Zakłócenie' ,'Location','southeast')
if SAVE_FIG saveas(gcf, 'y1u2_odsprz_zak.png'); end
