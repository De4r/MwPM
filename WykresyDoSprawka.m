clc
close all
nazwa = 'wyniki1525';
wyniki = load([nazwa, '.mat']);
wyniki = cell2mat(struct2cell(wyniki));

t = wyniki(:,1);
y = wyniki(:,2:end);
xopis = ['$Czas\ [s]$'];
yopis = {['$Predkosc\ katowa\ \frac{rad}{s}$'], ['$Predkosc\ \frac{m}{s}$'], ['$Przem.\ katowe\ [rad]$'], ['$Przem.\ [m]$']};
tytuly = {['$Predkosc\ katowa\ dla\ \alpha$'], ['$Przem.\ katowe\ dla\ \alpha$'], ['$Predkosc\ dla\ x$'], ['$Przemieszczenie\ dla\ x$']...
    ['$Predkosc\ dla\ y$'], ['$Przemieszczenie\ dla\ y$'], ['$Predkosc\ katowa\ dla\ \phi_1$'], ['$Przem.\ katowe\ dla\ \phi_1$']...
    ['$Predkosc\ katowa\ dla\ \phi_2$'], ['$Przem.\ katowe\ dla\ \phi_2$'], ['$Predkosc\ dla\ x_1$']...
    ['$Predkosc\ dla\ y_1$'], ['$Predkosc\ dla\ x_2$'], ['$Predkosc\ dla\ y_2$'],  ['$Przemieszcenie\ dla\ x_1\ i\ x_2$']...
    ['$Przemieszczenie\ dla\ y_1\ i\ y_2$'], ['$Przemieszczenie\ dla\ x_2$'], ['$Przemieszczenie\ dla\ y_2$']};
legenda = {['$x_1$'], ['$y_1$'], ['$x_2$'], ['$y_2$']};
yopis = string(yopis);
tytuly = string(tytuly);
for i=1:5
    figure(i)
    figure('units','normalized','outerposition',[0 0 1 1])
    plot(t,y(:,i+5),'Linewidth',1.5)
    xlabel({xopis},'Interpreter','latex')
    grid on
    if i+5==6 || i+5==9 || i+5==10
        ylabel({yopis(3)},'Interpreter','latex')
    else
        ylabel({yopis(4)},'Interpreter','latex')
    end
    title({tytuly(2*i)},'Interpreter','latex')
    set(gca,'fontsize',14)
    
    nazwazapisu = [nazwa, num2str(i),'.png'];
    print(gcf,nazwazapisu,'-dpng','-r300'); 
end

for i=1:2
    figure(5+i)
    figure('units','normalized','outerposition',[0 0 1 1])
    plot(t,y(:,14+i),'Linewidth',1.5)
    hold on
    plot(t,y(:,16+i),'Linewidth',1)
    xlabel({xopis},'Interpreter','latex')
    grid on
    ylabel({yopis(4)},'Interpreter','latex')
    title({tytuly(14+i)},'Interpreter','latex')
    legend([string(legenda(i)),string(legenda(i+2))],'Interpreter','latex')
    set(gca,'fontsize',14)
    nazwazapisu = [nazwa, num2str(i),num2str(i),'.png'];
    print(gcf,nazwazapisu,'-dpng','-r300'); 
end

figure(10)
figure('units','normalized','outerposition',[0 0 1 1])
plot(y(:,7),y(:,8),'Linewidth',1)
title({'$Trajektoria\ ruchu\ srodka\ ciezkosci$'}, 'Interpreter','latex')
xlabel({'$x\ [m]$'}, 'Interpreter','latex')
ylabel({'$y\ [m]$'}, 'Interpreter','latex')
grid on
axis equal
set(gca,'fontsize',14)
nazwazapisu = [nazwa,'trajek.png'];
print(gcf,nazwazapisu,'-dpng','-r300'); 

figure(11)
figure('units','normalized','outerposition',[0 0 1 1])
plot(t,y(:,9)-y(:,10),'Linewidth',1.5)
title({'$Roznica\ katow\ fazowych$'}, 'Interpreter','latex')
xlabel({xopis}, 'Interpreter','latex')
ylabel({yopis(3)}, 'Interpreter','latex')
grid on
set(gca,'fontsize',14)
nazwazapisu = [nazwa,'katy1.png'];
print(gcf,nazwazapisu,'-dpng','-r300'); 

figure(12)
figure('units','normalized','outerposition',[0 0 1 1])
plot(t,y(:,4)-y(:,5),'Linewidth',1.5)
title({'$Roznica\ katow\ fazowych$'}, 'Interpreter','latex')
xlabel({xopis}, 'Interpreter','latex')
ylabel({yopis(1)}, 'Interpreter','latex')
grid on
set(gca,'fontsize',14)
nazwazapisu = [nazwa,'katy2.png'];
print(gcf,nazwazapisu,'-dpng','-r300'); 