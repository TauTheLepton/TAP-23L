function Gs = gs(model_lin_ciag)
%GS Summary of this function goes here
%   Detailed explanation goes here
global Tau_c Tau
Gs = tf(model_lin_ciag);
Gs.InputDelay = [0 Tau_c];
Gs.OutputDelay = [0 Tau];
% Gs11 = Gs(1,1);     %h/F_h
% Gs12 = Gs(1,2);     %h/F_cin
% Gs21 = Gs(2,1);     %T_out/F_h
% Gs22 = Gs(2,2);     %T_out/F_cin
end

