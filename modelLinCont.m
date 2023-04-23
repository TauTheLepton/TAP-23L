function model_lin_ciag = modelLinCont()
%MODELLINCONT Summary of this function goes here
%   Detailed explanation goes here

%% wyznaczenie modelu liniowego w postaci równań stanu
global F_h F_c F_d h c T_out T_h T_c T_d alpha
dhdt = (F_h + F_c + F_d - F(h))/(3*c*h^2);
dT_outdt = (F_h*T_h + F_c*T_c + F_d*T_d - F(h)*T_out)/(c*h^3);

dhdh = (3*alpha*sqrt(h) - 4*(F_h + F_c + F_d))/(6*c*h^3);
dhdT_out = 0;
dhdF_h = 1/(3*c*h^2);
dhdF_cin = 1/(3*c*h^2);

dT_outdh = (5*alpha*sqrt(h)*T_out - 6*(F_h*T_h + F_c*T_c + F_d*T_d))/(2*c*h^4);
dT_outdT_out = (-alpha*sqrt(h))/(c*h^3);
dT_outdF_h = T_h/(c*h^3);
dT_outdF_cin = T_c/(c*h^3);

global A B
A = [dhdh dhdT_out; dT_outdh dT_outdT_out];
B = [dhdF_h dhdF_cin; dT_outdF_h dT_outdF_cin];
C = eye(2);
D = zeros(2);

model_lin_ciag = ss(A,B,C,D);
end
