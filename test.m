%Ten przyk³ad pokazuje jak u¿yæ klasê regulatora PID

clear all;

%classPID(K, Ti, Kd, Td, Tp, Hlim, Llim, Dir, AutoMan, ManVal) 
%DIR: 1 - direct (SP-PV), 0 - indirect (PV-SP)
% AutoMan = 0 regulator na wyjœciu podaje wartoœæ sterowania rêcznego ManVal
% AutoMan = 1 regulator na wyjœciu podaje wartoœæ wyliczon¹ z prawa regulacji
%Przyklad: deklaracja regulatora D
p=classPID(4, 10, 1, 100, 1, 100, -100, 1, 1, 0)
o=classLEADLAG(0.38, 0, 90, 1, 100, 0)


%PID
%reTune(obj, K, Ti, Kd, Td)
% funkcja umo¿liwia zmianê nastaw regulatora
%p.reTune(1, 30, 0.5, 30)

stpt = 1;


cv=zeros(1,300);

for i=20:1:300
  %metoda wyliczaj¹ca prawo PID w oparciu o PV i STPT
  u =  p.calc(cv(i-1),stpt);
  mv(i)=u;
  cv(i) = o.calc(mv(i-19));
  
end
figure
plot(cv);
figure
plot(mv)



%test LEADLAG

o2=classLEADLAG(0.38, 0, 90, 1, 100, 0)
for i=1:1:300
  cv(i) = o2.calc(1);
  end
figure
plot(cv);



o3=classLEADLAG(0.38, 80, 90, 1, 100, 0)
for i=1:1:300
  cv(i) = o3.calc(1);
  end
figure
plot(cv);
