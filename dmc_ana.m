function error = dmc_ana(params)
    % params=[70, 6, 50];
    %% set parameters
    use_plot = true; % for optimalization with algorithm should be 'false'
    use_z = 0; %0-brak zaklocen, 1- zaklocenie skokowe, 2-zaklocenie sinusoidalne,3-zaklocenie skokowe z szumem
    dB = 40; %wielkosc szumu: 10 - duza, 20 - srednia, 40 - mala
    calc_z = false; %czy probuje kontrowac zaklocenie (mzp*dzp)
    % dmc params
    n = params(1);
    nu = params(2);
    lambda = params(3);
    % goal values
    len = 600;
    yzad = ones(len, 1);
    yzad(1:len/4) = 0;
    
    %% read s and d
    temp = load('S_u.mat');
    s = temp.S_u;
    d = 170;
    temp = load('S_z.mat');
    sz = temp.S_z;
    dz = 150;
    %% calculate offline
    % generate mp
    mp = zeros(n, d-1);
    for i=1:n
        for j=1:d-1
            first = i+j;
            if first > d
                first = d;
            end
            mp(i, j) = s(first) - s(j);
        end
    end
    % generate mzp
    mzp = zeros(n, dz);
    for i=1:dz
        if i==1
           mzp(1:n,1)=sz(1:n);
        else
           mzp(1:n,i)=sz(i:n+i-1)-sz(i-1);
        end
    end
    % generate m
    m = zeros(n, nu);
    for i=1:n
        for j=1:nu
            if i >= j
                idx = i-j+1;
                if idx > d
                    idx = d;
                end
                m(i, j) = s(idx);
            end
        end
    end
    % calculate K
    K = (m'*m+(lambda*eye(nu)))\m';
    %% main loop, online calculations
    dup = zeros(d-1, 1);
    dzp = zeros(dz,1);
    u = zeros(len, 1);
    y = zeros(len, 1);
    z = zeros(len, 1);
    if use_z == 1
        z(350:len,1)=1;
    end
    if use_z == 2
        z(350:len, 1) = sin(linspace(0,1,len-350+1)*20)/3; % sinusoidal disturbance
    end
    if use_z == 3
        z(350:len,1)=1;
        z=awgn(z, dB);
    end
    
    err = zeros(len, 1);
    for k=13:len
        %% process simulation
        y(k)=symulacja_obiektu8y_p2(u(k-6),u(k-7),z(k-1),z(k-2),y(k-1),y(k-2));
        %% calculate error
        err(k) = yzad(k) - y(k);
        %% calculate control values
        % calculate dup
        for i=1:d-1
            if k-i-1 >= 1
                dup(i) = u(k-i) - u(k-i-1);
            end
        end
        % calculate y0
        y0 = (y(k)*ones(n, 1)) + mp*dup;
        if calc_z
            y0=y0+mzp*dzp;
        end
        % calculate du
        du = K*(yzad(k)*ones(n, 1) - y0);
        % change u increment to u value
        u(k) = u(k-1) + du(1);
        %calculate dzp
        for i = dz:-1:2
            dzp(i)=dzp(i-1);
        end
        dzp(1)=z(k)-z(k-1);
    end
    %% calculate error
    error = err' * err;
    %% plot
    if use_plot
        subplot(3, 1, 1)
        hold on
        plot(yzad)
        plot(y)
        hold off
        legend('yzad', 'y')
    
        subplot(3, 1, 2)
        plot(u)
        legend('u')
        subplot(3, 1, 3)
        plot(z)
        legend('z')
    end
end