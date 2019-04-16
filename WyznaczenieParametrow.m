clc
clear all


A = readtable('drgania.csv');
A = table2array(A);
% 1 kolumna - czas, 2 -kolumna drgania w pionie sama rynna
% 3 koluman drgania w pionie z masa dodatkowa md=7.4kg
% 4 kolumna sama rynna w poziomie
% Kolumny maja wymiar napiecia [V]
md= 7.4; % kg
S = (0.1^-1)*9.81; % [mV* m/s^2] to jest 9,81 na 0.1 V
%wyznaczenie przyspieszenia V/S
A(:,2:end) = A(:,2:end)/0.1*9.81;
A1 = movmean(A,30);
figure(1)

plot(A(:,1), A(:,4), 'b')
hold on
plot(A1(:,1), A1(:,4),'r')
hold on
zero = zeros(length(A(:,1)),1);
plot(A1(:,1), zero)

y11 = 2.695;
y21 = 2.127;
T1 = 0.35;
n1 = 10;
w1= 2*pi/T1
beta1 = -log(y21/y11)/(n1*T1)

y12 = 2.734;
y22 = 1.883;
T2 = 0.37;
n2 = 15;
w2 = 2*pi/T2
beta2 = -log(y22/y12)/(n2*T2)

y13 = 0.759;
y23 = 0.4344;
T3 = 0.85;
n3 = 6;
w3 = 2*pi/T3
beta3 = -log(y23/y13)/(n3*T3)

z_ky = md*((beta1.^2+w1^2)*(beta2.^2+w2.^2))/((beta1.^2+w1.^2)-(beta2.^2+w2.^2)) % zastepczy zky =2 ky
m = z_ky /(beta1.^2+w1^2) %masa
z_by = 2*m*beta1
z_kx = (w3.^3+beta3.^2)*m % wspol na x
z_bx = 2*m*beta3

kx = z_kx / 2
ky = z_ky / 2
bx = z_bx / 2
by = z_by / 2
