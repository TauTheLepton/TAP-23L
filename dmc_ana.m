close all;
h = 15.56;
T = 40.42;
T_out = T;
%% set parameters
use_plot = true; % for optimalization with algorithm should be 'false'
% dmc params
N = 17;
Nu = 15;
D = 100;
lambda = 70;
psi = 1;
umax = 50;
umin = 0;
% goal values
Tp = 10;
t = 1:Tp:Tp*1000;
t = t';
len = size(t, 1);
err = zeros(len, 2);
yzad = zeros(len, 2);
yzad(len/4:end,1) = 35;
yzad(len/2:end,2) = 20;

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
%% plot
if use_plot
    subplot(4, 1, 1)
    hold on
    plot(yzad(:,1)+h)
    plot(y(:,1)+h)
    hold off
    legend('hzad', 'h')

    subplot(4, 1, 2)
    plot(u(:,1))
    legend('Fh')

    subplot(4, 1, 3)
    hold on
    plot(yzad(:,2)+T_out)
    plot(y(:,2)+T_out)
    hold off
    legend('Tzad', 'T')

    subplot(4, 1, 4)
    plot(u(:,2))
    legend('Fc')
end