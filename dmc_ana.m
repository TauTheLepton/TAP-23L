% function error = dmc_ana(params)
close all;
    h = 15.56;
    T = 40.42;
    T_out = T;
    params=[10, 6, 50];
    %% set parameters
    use_plot = true; % for optimalization with algorithm should be 'false'
    use_z = 0; %0-brak zaklocen, 1- zaklocenie skokowe, 2-zaklocenie sinusoidalne,3-zaklocenie skokowe z szumem
    dB = 40; %wielkosc szumu: 10 - duza, 20 - srednia, 40 - mala
    calc_z = false; %czy probuje kontrowac zaklocenie (mzp*dzp)
    % dmc params
    N = params(1);
    Nu = params(2);
    lambda = params(3);
    N = 5;
    Nu = 3;
    lambda = 7;
    psi = 0.5;
    % goal values
    Tp = 10;
    t = 1:Tp:Tp*1000;
    t = t';
    len = size(t, 1);
    err = zeros(len, 2);
    yzad = zeros(len, 2);
    yzad(len/4:end,1) = 2;
    yzad(len/4:end,2) = 2;
    
    %% read s and d
    % temp = load('S_u.mat');
    D = 7;
    shFh = xd1(1:D,1,1); %TODO: pobierać z pliku
    sToutFh = xd1(1:D,2,1);
    shFcin = xd2(1:D,1,2);
    sToutFcin = xd2(1:D,2,2);
    S = zeros(2,2,D);
    S(1,1,:) = shFh;
    S(1,2,:) = shFcin;
    S(2,1,:) = sToutFh;
    S(2,2,:) = sToutFcin;
    % S = [[shFh shFcin]; [sToutFh sToutFcin]]
    % temp = load('S_z.mat');
    % sz = temp.S_z;
    % dz = 150;
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
    % generate mzp
    % mzp = zeros(n, dz);
    % for i=1:dz
    %     if i==1
    %        mzp(1:n,1)=sz(1:n);
    %     else
    %        mzp(1:n,i)=sz(i:n+i-1)-sz(i-1);
    %     end
    % end
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
    K = (m'*(psi*eye(N*2))*m+(lambda*eye(Nu*2)))\m';
    %% main loop, online calculations
    dup = zeros(2*(D-1), 2);
    % dzp = zeros(dz,1);
    u = zeros(len, 2);
    y = zeros(len, 2);
    % z = zeros(len, 1);
    % if use_z == 1
    %     z(350:len,1)=1;
    % end
    % if use_z == 2
    %     z(350:len, 1) = sin(linspace(0,1,len-350+1)*20)/3; % sinusoidal disturbance
    % end
    % if use_z == 3
    %     z(350:len,1)=1;
    %     z=awgn(z, dB);
    % end
    
    for k=2:len
        %% process simulation
        yvect = lsim(Gz, u(1:k,:), t(1:k));
        last = size(yvect,1);
        y(k,:) = yvect(last,:);
        %% calculate error
        err(k,:) = yzad(k,:) - y(k,:);
        %% calculate control values
        % calculate dup
        for i=1:D-1
            if k-i-1 >= 1
                dup(i, :) = u(k-i, :) - u(k-i-1, :);
            end
        end
        % calculate y0
        y0 = y(k,:) + mp*dup;
        du = K*(yzad(k,:) - y0);
        % change u increment to u value
        u(k,:) = u(k-1,:) + du(1,:);
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
        legend('yzad', 'y')
    
        subplot(4, 1, 2)
        plot(u(:,1))
        legend('u')

        subplot(4, 1, 3)
        hold on
        plot(yzad(:,2)+T_out)
        plot(y(:,2)+T_out)
        hold off
        legend('yzad', 'y')

        subplot(4, 1, 4)
        plot(u(:,2))
        legend('u')
    end
% end