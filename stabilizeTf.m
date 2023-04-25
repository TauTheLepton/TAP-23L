function G = stabilizeTf(G)
%STABILIZETF Summary of this function goes here
%   Detailed explanation goes here

idel = G.InputDelay;
odel = G.OutputDelay;
iodel = G.IODelay;
Tp = G.Ts;

[num,den] = tfdata(G);
[z,p,k] = tf2zpk(num{1},den{1});
up = p(abs(p) >= 1);
z = [z; up];

minzp = min([z;p]) / 10;

while len(z) > len(p)
    p = [p; minzp];
end

G = zpk(z, p, k, Tp);
G = tf(G);

G.InputDelay = idel;
G.OutputDelay = odel;
G.IODelay = iodel;

end

function l = len(x)
s = size(x);
l = max(s);
end
