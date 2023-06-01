clear all;
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

%% Odpowiedzi
Tprob = 10;

length = 300;

skokU = 1; %0 jesli badamy skoki Z
K = 1; %wartosc skoku

hvect = zeros(length, 101);
Tvect = zeros(length, 101);

conf = RespConfig('Amplitude', [10, 0], 'Delay', 0);
% model_lin_dysk = c2d(model_lin_ciag,Tp,'zoh');

if skokU == 1
    u = ones(length, 1)*K;
    for k=8:length
        [xd1, taxis] = step(model_lin_dysk, tk, conf);
        xd1(:, 1, 1) = xd1(:, 1, 1) + h;
        xd1(:, 2, 1) = xd1(:, 2, 1) + T_out;
        hvect(k,:) = transpose(xd1(:,1,1));
        Tvect(k,:) = transpose(xd1(:,2,1));
    end
else
    z = ones(length, 1)*K;
    for k=3:length
        % [h(k), T(k)]=step(model_lin_dysk, tk, conf);
    end
end
hvect
Tvect
%program napisany tak zeby od razu dawal tez S-y dla dmc


hold on
subplot(3,1,1);
plot(h)
subplot(3,1,2);
stairs(u)
subplot(3,1,3);
stairs(z)
hold off