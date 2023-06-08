function E_abs = DMC_func2(X)
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

h = 15.56;
T = 40.42;
T_out = T;

%% Modele
model_lin_ciag = modelLinCont();
Gs = gs(model_lin_ciag);
model_lin_dysk = modelLinDisc(model_lin_ciag,Tp);
Gz = gz(Gs,Tp);

conf = RespConfig('Amplitude', [1, 0], 'Delay', 0);%'InputOffset', [F_h, F_cin], 'InitialState', [h_start, T_out_start]);
[xd1, ~] = step(model_lin_dysk, tk, conf);
xd1(:, 1, 1) = xd1(:, 1, 1) + h;
xd1(:, 2, 1) = xd1(:, 2, 1) + T_out;

conf = RespConfig('Amplitude', [0, 1], 'Delay', Tau_c);%'InputOffset', [F_h, F_cin], 'InitialState', [h_start, T_out_start]);
[xd2, ~] = step(model_lin_dysk, tk, conf);
xd2(:, 1, 2) = xd2(:, 1, 2) + h;
xd2(:, 2, 2) = xd2(:, 2, 2) + T_out;

%% set parameters
use_plot = true; % for optimalization with algorithm should be 'false'
% dmc params
N = X(1);
Nu = X(2);
D = X(4);
lambda = X(3);
psi = 1;
umax = 100;
umin = 0;
% goal values
Tp = 10;
t = 1:Tp:Tp*300;
t = t';
len = size(t, 1);
err = zeros(len, 2);
yzad = zeros(len, 2);
yzad(len/15:end,1) = 10;
yzad(len/2:end,2) = 5;

%% read s and d
shFh = xd1(1:D,1,1); % launch zadanie1.m for xd1 and xd2
sToutFh = xd1(1:D,2,1);
shFcin = xd2(1:D,1,2);
sToutFcin = xd2(1:D,2,2);
S = zeros(2,2,D);
S(1,1,:) = shFh;
S(1,2,:) = shFcin;
S(2,1,:) = sToutFh;
S(2,2,:) = sToutFcin;
%% calculate offline
% generate mp
mp = zeros(2*N, 2*(D-1));
for i=1:N
    for j=1:D-1
        idx = i+j;
        if idx > D
            idx = D;
        end
        mp(2*i-1:2*i, 2*j-1:2*j) = S(:,:,idx) - S(:,:,j);
    end
end
% generate m
m = zeros(2*N, 2*Nu);
for i=1:N
    for j=1:Nu
        if i >= j
            idx = i-j+1;
            if idx > D
                idx = D;
            end
            m(2*i-1:2*i, 2*j-1:2*j) = S(:,:,idx);
        end
    end
end
% calculate K
K = (m'*(psi*eye(N*2))*m+(lambda*eye(Nu*2)))\m'*(psi*eye(N*2));
%% main loop, online calculations
dup = zeros(2*(D-1), 2);
u = zeros(len, 2);
y = zeros(len, 2);    
for k=2:len
    %% process simulation
    yvect = lsim(Gz, u(1:k,:), t(1:k));
    y(k,:) = yvect(end,:);
    %% calculate error
    err(k,:) = yzad(k,:) - y(k,:);
    %% calculate control values
    % calculate dup
    for i=1:D-1
        if k-i-1 >= 1
            dup(i, :) = u(k-i, :) - u(k-i-1, :);
        end
    end
    dup;
    % calculate y0
    y0 = y(k,:) + mp*dup;
    du = K(1:2,:)*(yzad(k,:) - y0);
    % change u increment to u value
    u(k,1) = u(k-1,1) + du(2,2);
    u(k,2) = u(k-1,2) + du(2,1); 
    if u(k,1) > umax
        u(k,1) = umax;
    end
    if u(k,1) < umin
        u(k,1) = umin;
    end
    if u(k,2) > umax
        u(k,2) = umax;
    end
    if u(k,2) < umin
        u(k,2) = umin;
    end
end
%% calculate error
error = [0; 0];
error(1) = err(:,1)' * err(:,1);
error(2) = err(:,2)'*err(:,2);
E_abs = error(1) + error(2);
%% plot
end