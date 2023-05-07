function E_abs = PID_func(X)
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

    Gz = gz(Gs,Tp);
    Gz22 = Gz(2,2);
    Gz22.InputDelay = 0;
    Gz22.OutputDelay = 0;
    Gz(2,2);
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
    
    % % reTune(obj, K, Ti, Kd, Td)
    pid_y2u1.reTune(X(1), X(2), X(3), X(4));
    pid_y1u2.reTune(X(5), X(6), X(7), X(8));
    
    stpt = zeros(len,2);
    stpt(20:end,1) = 10;
    stpt(40:end,2) = 5;
    E_abs = 0;
    
    lazy_start = 2;  % potrzebne do ustawienia wyjść procesu potrzebnych do PIDa
    
    y = lsim(Gz, u(1:lazy_start,:), t(1:lazy_start));
    
        for i = lazy_start+1:len
            if i == 150
                Gz = Gz2;
            end
            % jestesmy w czasie 'i-1'
            ur(i,1) = pid_y2u1.calc(y(end,2), stpt(i,2));
            ur(i,2) = pid_y1u2.calc(y(end,1), stpt(i,1));
            E_abs = E_abs + abs(pid_y2u1.err(1)) + abs(pid_y2u1.err(2)) + abs(pid_y1u2.err(1)) + abs(pid_y1u2.err(2));
            u1 = lsim(D12, ur(1:i,2), t(1:i)) + ur(1:i,1);
            u2 = lsim(D21, ur(1:i,1), t(1:i)) + ur(1:i,2);
            % jestesmy w czasie 'i'
            y = lsim(Gz, [u1, u2], t(1:i));
        end
    
end