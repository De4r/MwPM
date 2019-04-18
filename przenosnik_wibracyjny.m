function dX=przenosnik_wibracyjny(t, X, t1)
% ------------------------ parametry --------------------------------
m1 = 3.86; %kg masy nie wywazone
e = 0.02; %cm
mk= 55.2323; %kg mkorpusu+2*mstoj+2*mwalow+2*mrotorow
Jc = 1.3336; %kgm2 moment bezwladnosci korpusu
Jc1 = 4.8403*10^-05; %kgm2 moment bezwldnosci masy nie wywazonej wzgledem jej sr ciez.
Js1 = 0.0013/2 ; %kgm2 moment bezwaldnosci stojanu
Jw1 = 0.0013/2; %kgm2 moment bezwladnosci wirnika i rotora
R=0.01; % wspol restytucji
K=1.5*10^8; % cos tam szwarca
bx = 6.882; % wspol tlum na X
by = 4.2571; % wspol tlum na Y
kx = 1.2714*10^4; % wpol sprez na X
ky = 1.0144*10^4; % wspol sprez na Y
l = 0.5; %m
H = 0.15; %m
a = 0.225; %m
beta = pi/6;
% --------------------nadawa------------------
mn1 = 0.72; %kg % dla pola 0.25m^2 i wys. 2 mm dla piasku 1440 kq/m^3
mn2 = 0.72; %kg
u = 0.4; % wspol tracia rynna 1 warstaw
u2 = 0.7; % wspol tarcia miedzy nadawami
g = 9.81;
% ----------------------  zmienne  -------------------------------------
walfa = X(1); vx = X(2); vy = X(3); wfi1 = X(4); wfi2 = X(5); %predkosci
alfa = X(6); x = X(7); y = X(8); fi1 = X(9); fi2 = X(10); % przemieszcenia
vx1 = X(11); vy1 = X(12); vx2 = X(13); vy2 = X(14); % predkosci
x1 = X(15); y1 = X(16); x2 = X(17); y2 = X(18); % przemieszczenia
if t<t1
    %  -----------------------------  macierz mas  --------------------------
    diagonala = [Jc+2*m1*a^2, 2*m1+mk, 2*m1+mk, Jc1+Js1+Jw1+m1*e^2, Jc1+Js1+Jw1+m1*e^2, 1, 1, 1, 1, 1, mn1, mn1, mn2, mn2, 1, 1, 1, 1];
    M = diag(diagonala);
    M(1,4) = -1*a*m1*e*cos(beta-alfa-fi1); M(1,5) = a*m1*e*cos(beta-alfa+fi2);
    M(2,4) = m1*e*sin(fi1); M(2,5) = -m1*e*sin(fi2);
    M(3,4) = m1*e*cos(fi1); M(3,5) = m1*e*cos(fi2);
    M(4,2) = m1*e*sin(fi1); M(4,3) = m1*e*cos(fi1);
    M(5,2) = -m1*e*sin(fi2); M(5,3) = m1*e*cos(fi2);
    % ------------------------dane silnikowe----------------
    n_nom=1400;     
    n_synchr=1500; 
    p=2.5;   % przeciazalnosc      
    Prob_nom=0.16;  %kW
    nrob_nom=1400;  

    % --------wyznaczenie sta³ych dla silnika asynchronicznego---------
    w_znam=(pi*n_nom)/30; %omega znamionowa
    w_synchr=(pi*n_synchr)/30;  %omega synchroniczna
    Mrob_el=Prob_nom*1000/(w_znam);  %moment silnika roboczego
    Mrob_ut=p*Mrob_el;            %moment utyku silnika roboczego
    srob_nom=(n_synchr-nrob_nom)/n_synchr;  %poœliz nominalny silnika roboczego
    srob_ut=srob_nom*(p+sqrt(p*p-1));       %poœlizg utyku nominalny silnika roboczego

    % ------------------obliczenie momentow elektrycznych ---------------------
    Mel1=((2*Mrob_ut*srob_ut*((w_synchr-wfi1)/w_synchr))/(srob_ut*srob_ut+((w_synchr-wfi1)/w_synchr)*((w_synchr-wfi1)/w_synchr)));
    Mel2=((2*Mrob_ut*srob_ut*((w_synchr-wfi2)/w_synchr))/(srob_ut*srob_ut+((w_synchr-wfi2)/w_synchr)*((w_synchr-wfi2)/w_synchr)));

    % ------------wyznacznie sil nadawy ------------------
%     F12 = (y1-y2)*K*(1-(1-R^2)/2*(1-sign(y1-y2)*sign(vy1-vy2)));
%     Fr1 = (y-y1)*K*(1-(1-R^2)/2*(1-sign(y-y1)*sign(vy-vy1)));
%     T12 = u*F12*sign(vx1-vx);
%     Tr1 = u*Fr1*sign(vx2-vx1);
    F12=0; Fr1=0; T12=0; Tr1=0;
    %  ------------------------- macierz wyrazów wolnych ----------------------
    X(1) = Mel1-Mel2-2*H*bx*vx-2*H*kx*x-2*H^2*kx*alfa-2*H^2*bx*walfa-2*ky*l^2*alfa-2*by*l^2*walfa+m1*e*a*(sin(beta-alfa-fi1)*wfi1^2+sin(beta-alfa+fi2)*wfi2^2); % po walfa
    X(2) = -Tr1-m1*e*cos(fi1)*wfi1^2+m1*e*cos(fi2)*wfi2^2-2*bx*vx-2*kx*x-2*H*kx*alfa-2*H*bx*walfa; % po vx
    X(3) = -Fr1-g*(2*m1+mk)+m1*e*sin(fi2)*wfi2^2+m1*e*sin(fi1)*wfi1^2-2*by*vy-2*ky*y; % po vy
    X(4) = Mel1-m1*e*g*cos(fi1)-a*e*m1*cos(beta-alfa-fi1)*walfa-a*e*m1*sin(beta-alfa-fi1)*(alfa*walfa-alfa*wfi1+walfa*wfi1); % po wfi1
    X(5) = Mel2-m1*e*g*cos(fi2)-a*e*m1*cos(beta-alfa+fi2)*walfa-a*e*m1*sin(beta-alfa+fi2)*(alfa*walfa-alfa*wfi2+walfa*wfi2); % po wfi2
    X(6) = walfa; % po alfa
    X(7) = vx; % po x
    X(8) = vy; % po y
    X(9) = wfi1; % po fi1
    X(10) = wfi2; % po fi2
%     X(11) = Tr1-T12; % po vx1
%     X(12) = Fr1-F12-mn1*g; % po vy1
%     X(13) = T12; % po vx2
%     X(14) = F12-mn2*g; % po vy2
%     X(15) = vx1; % po x1
%     X(16) = vy1; % po y1
%     X(17) = vx2; % po x2
%     X(18) = vy2; % po y2
    X(11) = 0; % po vx1
    X(12) = 0; % po vy1
    X(13) = 0; % po vx2
    X(14) = 0; % po vy2
    X(15) = 0; % po x1
    X(16) = 0; % po y1
    X(17) = 0; % po x2
    X(18) = 0; % po y2
    dX=M\X;  
else
    %  -----------------------------  macierz mas  --------------------------
    diagonala = [Jc+2*m1*a^2, 2*m1+mk, 2*m1+mk, Jc1+Js1+Jw1+m1*e^2, Jc1+Js1+Jw1+m1*e^2, 1, 1, 1, 1, 1, mn1, mn1, mn2, mn2, 1, 1, 1, 1];
    M = diag(diagonala);
    M(1,4) = -1*a*m1*e*cos(beta-alfa-fi1); M(1,5) = a*m1*e*cos(beta-alfa+fi2);
    M(2,4) = m1*e*sin(fi1); M(2,5) = -m1*e*sin(fi2);
    M(3,4) = m1*e*cos(fi1); M(3,5) = m1*e*cos(fi2);
    M(4,2) = m1*e*sin(fi1); M(4,3) = m1*e*cos(fi1);
    M(5,2) = -m1*e*sin(fi2); M(5,3) = m1*e*cos(fi2);
    % ------------------------dane silnikowe----------------
    n_nom=1400;     
    n_synchr=1500; 
    p=2.5;   % przeciazalnosc      
    Prob_nom=0.16;  %kW
    nrob_nom=1400;  

    % --------wyznaczenie sta³ych dla silnika asynchronicznego---------
    w_znam=(pi*n_nom)/30; %omega znamionowa
    w_synchr=(pi*n_synchr)/30;  %omega synchroniczna
    Mrob_el=Prob_nom*1000/(w_znam);  %moment silnika roboczego
    Mrob_ut=p*Mrob_el;            %moment utyku silnika roboczego
    srob_nom=(n_synchr-nrob_nom)/n_synchr;  %poœliz nominalny silnika roboczego
    srob_ut=srob_nom*(p+sqrt(p*p-1));       %poœlizg utyku nominalny silnika roboczego

    % ------------------obliczenie momentow elektrycznych ---------------------
    Mel1=((2*Mrob_ut*srob_ut*((w_synchr-wfi1)/w_synchr))/(srob_ut*srob_ut+((w_synchr-wfi1)/w_synchr)*((w_synchr-wfi1)/w_synchr)));
    Mel2=((2*Mrob_ut*srob_ut*((w_synchr-wfi2)/w_synchr))/(srob_ut*srob_ut+((w_synchr-wfi2)/w_synchr)*((w_synchr-wfi2)/w_synchr)));

    % ------------wyznacznie sil nadawy ------------------
    F12 = (y1-y2)*K*(1-((1-R^2)/2)*(1-sign(y1-y2)*sign(vy1-vy2)));
    Fr1 = (y-y1)*K*(1-((1-R^2)/2)*(1-sign(y-y1)*sign(vy-vy1)));
    T12 = -u2*F12*sign(vx2-vx1);
    Tr1 = -u*Fr1*sign(vx1-vx);
    %  ------------------------- macierz wyrazów wolnych ----------------------
    X(1) = Mel1-Mel2-2*H*bx*vx-2*H*kx*x-2*H^2*kx*alfa-2*H^2*bx*walfa-2*ky*l^2*alfa-2*by*l^2*walfa+m1*e*a*(sin(beta-alfa-fi1)*wfi1^2+sin(beta-alfa+fi2)*wfi2^2); % po walfa
    X(2) = -Tr1-m1*e*cos(fi1)*wfi1^2+m1*e*cos(fi2)*wfi2^2-2*bx*vx-2*kx*x-2*H*kx*alfa-2*H*bx*walfa; % po vx
    X(3) = -Fr1-g*(2*m1+mk)+m1*e*sin(fi2)*wfi2^2+m1*e*sin(fi1)*wfi1^2-2*by*vy-2*ky*y; % po vy
    X(4) = Mel1-m1*e*g*cos(fi1)-a*e*m1*cos(beta-alfa-fi1)*walfa-a*e*m1*sin(beta-alfa-fi1)*(alfa*walfa-alfa*wfi1+walfa*wfi1); % po wfi1
    X(5) = Mel2-m1*e*g*cos(fi2)-a*e*m1*cos(beta-alfa+fi2)*walfa-a*e*m1*sin(beta-alfa+fi2)*(alfa*walfa-alfa*wfi2+walfa*wfi2); % po wfi2
    X(6) = walfa; % po alfa
    X(7) = vx; % po x
    X(8) = vy; % po y
    X(9) = wfi1; % po fi1
    X(10) = wfi2; % po fi2
    X(11) = Tr1-T12; % po vx1
    X(12) = Fr1-F12-mn1*g; % po vy1
    X(13) = T12; % po vx2
    X(14) = F12-mn2*g; % po vy2
    X(15) = vx1; % po x1
    X(16) = vy1; % po y1
    X(17) = vx2; % po x2
    X(18) = vy2; % po y2

    dX=M\X;  
end
