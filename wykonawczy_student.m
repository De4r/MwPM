clear all
clc
t = 0:0.0001:10;
X0 = zeros(4,1);
% X0(7,1)=-pi/6;
% X0(8,1)=-(pi-pi/6);
[t X] = ode45(@model_student, t, X0);
 
 
figure(1);
subplot(4,1,1); plot(t,X(:,1)); title('Prêdkoœæ masy m'), 
xlabel('Czas [s]'), ylabel('vy [m/s]');
grid on
subplot(4,1,2); plot(t,X(:,2)); title('Predksc obrotu w'),
xlabel('Czas [s]'), ylabel('vwfi [m]');
grid on
subplot(4,1,3); plot(t,X(:,3)); title('Przemieszczenie masy m'), 
xlabel('Czas [s]'), ylabel('vy [m/s]');
grid on
subplot(4,1,4); plot(t,X(:,4)); title('kat obrotu w'),
xlabel('Czas [s]'), ylabel('w [m]');
grid on
% figure(2);
% subplot(2,1,1); plot(t,X(:,2)); title('Prêdkoœæ k¹towa masy m1'),
% xlabel('Czas [s]'), ylabel('wfi [rad/s]');
% subplot(2,1,2); plot(t,X(:,4)); title('Przemieszczenie katowe masy m1'),
% xlabel('Czas [s]'), ylabel('fi [rad]');
 
 




