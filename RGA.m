function lambda = RGA(K)
%RGA Summary of this function goes here
%   Detailed explanation goes here
K_dash = inv(K)';
lambda = K .* K_dash;
end

