clc
clear all
close all
% parpool

t0 = 0;
tk = 25;
t1 = 15;
% krok = 0.1;
% t = t0:krok:tk-krok;
t = [t0, tk];

y0 = zeros(1,18);
tic
[t,y] = ode113(@(t,x) przenosnik_wibracyjny(t,x,t1),t,y0);
toc
xopis = ['$Czas\ [s]$'];
yopis = {['$Predkosc\ katowa\ \frac{rad}{s}$'], ['$Predkosc\ \frac{m}{s}$'], ['$Przem.\ katowe\ [rad]$'], ['$Przem.\ [m]$']};
tytuly = {['$Predkosc\ katowa\ dla\ \alpha$'], ['$Przem.\ katowe\ dla\ \alpha$'], ['$Predkosc\ dla\ x$'], ['$Przemieszczenie\ dla\ x$']...
    ['$Predkosc\ dla\ y$'], ['$Przemieszczenie\ dla\ y$'], ['$Predkosc\ katowa\ dla\ \phi_1$'], ['$Przem.\ katowe\ dla\ \phi_1$']...
    ['$Predkosc\ katowa\ dla\ \phi_2$'], ['$Przem.\ katowe\ dla\ \phi_2$'], ['$Predkosc\ dla\ x_1$']...
    ['$Predkosc\ dla\ y_1$'], ['$Predkosc\ dla\ x_2$'], ['$Predkosc\ dla\ y_2$'],  ['$Przemieszcenie\ dla\ x_1$']...
    ['$Przemieszczenie\ dla\ y_1$'], ['$Przemieszczenie\ dla\ x_2$'], ['$Przemieszczenie\ dla\ y_2$']};
yopis = string(yopis);
tytuly = string(tytuly);
Y = [t, y];
nazwa = ['wyniki', num2str(t1), num2str(tk),'.mat'];
save(nazwa, "Y");
for i=1:5
    figure(i)
    subplot(211)
    plot(t,y(:,i))
    xlabel({xopis},'Interpreter','latex')
    grid on
    if i==1 || i==4 || i==5
        ylabel({yopis(1)},'Interpreter','latex')
    else
        ylabel({yopis(2)},'Interpreter','latex')
    end
    title({tytuly(2*i-1)},'Interpreter','latex')
    
    subplot(212)
    plot(t,y(:,i+5))
    xlabel({xopis},'Interpreter','latex')
    grid on
    if i==1 || i==4 || i==5
        ylabel({yopis(3)},'Interpreter','latex')
    else
        ylabel({yopis(4)},'Interpreter','latex')
    end
    title({tytuly(2*i)},'Interpreter','latex')
end

for i=1:4
    figure(5+i)
    subplot(211)
    plot(t,y(:,10+i))
    xlabel({xopis},'Interpreter','latex')
    grid on
    ylabel({yopis(2)},'Interpreter','latex')
    title({tytuly(10+i)},'Interpreter','latex')
    
    subplot(212)
    plot(t,y(:,i+14))
    xlabel({xopis},'Interpreter','latex')
    grid on
    ylabel({yopis(4)},'Interpreter','latex')
    title({tytuly(14+i)},'Interpreter','latex')
end
