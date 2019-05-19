clc
close all
nazwa = 'wyniki1525';
wyniki = load([nazwa, '.mat']);
wyniki = cell2mat(struct2cell(wyniki));

t = wyniki(:,1);
Y = wyniki(:,2:end);

% Parametry
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
g=9.81;

%zmienne
walfa=Y(:,1);
vx=Y(:,2);
vy=Y(:,3);
alfa=Y(:,6);
x=Y(:,7);
y=Y(:,8);

%podpora A - , polowa masy, tlumik, sprezyna na y i na alfie
yA = g*(2*m1+mk)/2 + by*vy + ky*y - by*l*walfa - ky*l*alfa;
xA = bx*vx + kx*x + H*kx*alfa + H*bx*walfa;

%podpora B - , polowa masy, tlumik, sprezyna na y i na alfie
yB = g*(2*m1+mk)/2 + by*vy + ky*y + by*l*walfa + ky*l*alfa;
xB = -bx*vx - kx*x - H*kx*alfa - H*bx*walfa;

figure(1)
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
plot(t,xA,'Linewidth',1.5)
title({'Reakcja w podporze $A$ na kierunku $X$'},'Interpreter','latex')
xlabel({'Czas [s]'},'Interpreter','latex')
ylabel({'Sila [N]'},'Interpreter','latex')
set(gca,'fontsize',14)
grid on

subplot(2,1,2)
plot(t,xB,'Linewidth',1.5)
title({'Reakcja w podporze $B$ na kierunku $X$'},'Interpreter','latex')
xlabel({'Czas [s]'},'Interpreter','latex')
ylabel({'Sila [N]'},'Interpreter','latex')
grid on
set(gca,'fontsize',14)
print(gcf,'ReakcjeNaX','-dpng','-r300'); 

figure(2)
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
plot(t,yA,'Linewidth',1.5)
title({'Reakcja w podporze $A$ na kierunku $Y$'},'Interpreter','latex')
xlabel({'Czas [s]'},'Interpreter','latex')
ylabel({'Sila [N]'},'Interpreter','latex')
set(gca,'fontsize',14)
grid on

subplot(2,1,2)
plot(t,yB,'Linewidth',1.5)
title({'Reakcja w podporze $B$ na kierunku $Y$'},'Interpreter','latex')
xlabel({'Czas [s]'},'Interpreter','latex')
ylabel({'Sila [N]'},'Interpreter','latex')
grid on
set(gca,'fontsize',14)
print(gcf,'ReakcjeNaY','-dpng','-r300'); 