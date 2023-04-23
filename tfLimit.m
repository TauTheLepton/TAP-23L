function lim = tfLimit(G, u_num, y_num, a)
%TFLIMIT Summary of this function goes here
%   Detailed explanation goes here
lim = zeros(y_num, u_num);
for u = 1:u_num
    for y = 1:y_num
        [num, den] = tfdata(G(y,u),'v');
        syms s;
        G_syms = poly2sym(num,s)/poly2sym(den,s);
        lim(y,u) = limit(G_syms, s, a);
    end
end
end
