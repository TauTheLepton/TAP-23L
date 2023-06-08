close all;

h = 15.56;
T = 40.42;
T_out = T;

%% Set parameters
use_plot = true; % For optimization with the algorithm, set to 'false'
% DMC params
N = 3;
Nu = 2;
D = 100;
lambda = 6;
psi = 1;
umax = 150;
umin = -50;
ymin = 0;
ymax = 100;
ymin2 = [ymin, ymin];
ymax2 = [ymax, ymax];
umax2 = [umax, umax];
umin2 = [umin, umin];

% Goal values
Tp = 10;
t = (1:Tp:Tp*1000)';
len = size(t, 1);
err = zeros(len, 2);
dy = zeros(len, 2);
yzad = zeros(len, 2);
yzad(len/4:end, 1) = 35;
yzad(len/2:end, 2) = 20;

%% Read s and d
shFh = xd1(1:D, 1, 1); % Launch zadanie1.m for xd1 and xd2
sToutFh = xd1(1:D, 2, 1);
shFcin = xd2(1:D, 1, 2);
sToutFcin = xd2(1:D, 2, 2);
S = zeros(2, 2, D);
S(1, 1, :) = shFh;
S(1, 2, :) = shFcin;
S(2, 1, :) = sToutFh;
S(2, 2, :) = sToutFcin;

%% Calculate offline
% Generate mp
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
% Generate m
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

%% Main loop, online calculations
dup = zeros(2*(D-1), 2);
u = zeros(len, 2);
y = zeros(len, 2);
H = 2 * (m' * psi * m + lambda * eye(Nu*2));
J = tril(ones(2*Nu));
A = [-J; J; -m; m];

for k = 2:len
    %% Process simulation
    yvect = lsim(Gz, u(1:k, :), t(1:k));
    y(k, :) = yvect(end, :);

    %% Calculate error
    err(k, :) = yzad(k, :) - y(k, :);

    %% Calculate dup
    for i=1:D-1
        if k-i-1 >= 1
            dup(i, :) = u(k-i, :) - u(k-i-1, :);
        end
    end
    dup;

    %% Calculate y0
    y0 = y(k,:) + mp*dup;
    dy = yzad(k, :) - y0;

    %% Calculate control values
    if k > D

        % Calculate control signal by minimizing the objective function
        F1 = -2 * m' * psi * dy(:,1);
        F2 = -2 * m' * psi * dy(:,2);
        b1 = [];
        b2 = [];
        for i = 1:2*Nu
            b1 = [b1 ; -umin2 + u(k-1,:)];
        end
        for i = 1:2*Nu
            b2 = [b2 ; umax2 - u(k-1,:)];
        end
        b3 = -ymin2 + y0;
        b4 = ymax2 - y0;

        b = [b1; b2; b3; b4];
        delta_u1 = quadprog(H, F1, A, b(:,1), [], [], umin * ones(2*Nu, 1), umax * ones(2*Nu, 1));
        delta_u2 = quadprog(H, F2, A, b(:,2), [], [], umin * ones(2*Nu, 1), umax * ones(2*Nu, 1));
        du = [delta_u2(2), delta_u1(2)];

        % Update control signal
        u(k, :) = u(k-1, :) + du;
    else
        % Use previous control signal if not enough history
        u(k, :) = u(k-1, :);
    end

    % Saturate control signals
    u(k, 1) = max(min(u(k, 1), umax), umin);
    u(k, 2) = max(min(u(k, 2), umax), umin);
end

%% Calculate error
error = [0; 0];
error(1) = err(:, 1)' * err(:, 1);
error(2) = err(:, 2)' * err(:, 2);

%% Plot
if use_plot
    subplot(4, 1, 1)
    hold on
    plot(yzad(:, 1) + h)
    plot(y(:, 1) + h)
    hold off
    legend('hzad', 'h')

    subplot(4, 1, 2)
    plot(u(:, 1))
    legend('Fh')

    subplot(4, 1, 3)
    hold on
    plot(yzad(:, 2) + T_out)
    plot(y(:, 2) + T_out)
    hold off
    legend('Tzad', 'T')

    subplot(4, 1, 4)
    plot(u(:, 2))
    legend('Fc')
end
