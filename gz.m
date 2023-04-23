function Gz = gz(Gs,Tp)
%GZ Summary of this function goes here
%   Detailed explanation goes here
Gz = c2d(Gs,Tp,'zoh'); % używamy transmitancji ciągłej, bo zawiera już opóźnienia
end
