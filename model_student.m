function dX=model_student(t,X)

% dane



k=100000; 
b=600;
g=9.81;
    
m=1000;
m1=50;
e=0.06; 

J1=0.5;


% %% dane do obliczen momentu oporu na wale silników
% zred=0.0018;
% d=0.04;

%% dane silnikowe
n_nom=980;     
n_synchr=1000; 
p=2.5;   % przeciazalnosc      
Prob_nom=4;  %kW
nrob_nom=950;  
 
 
%% wyznaczenie sta³ych dla silnika asynchronicznego
w_synchr=(3.14*n_synchr)/30;  %omega synchroniczna
Mrob_el=9550*Prob_nom/n_nom;  %moment silnika roboczego
Mrob_ut=p*Mrob_el;            %moment utyku silnika roboczego
srob_nom=(n_synchr-nrob_nom)/n_synchr;  %poœliz nominalny silnika roboczego
srob_ut=srob_nom*(p+sqrt(p*p-1));       %poœlizg utyku nominalny silnika roboczego



%% ----------------------  zmienne  -------------------------------------
vy   = X(1,1);  %dy/dt
y    = X(3,1);  %y
wfi = X(2,1);
fi  = X(4,1);


%  -----------------------------  macierz mas  --------------------------
M=zeros(4,4);

M(1, 1) = m+m1;
M(1, 2) = m1*e*cos(fi);
M(2, 1) = m1*e*cos(fi);
M(2, 2) = m1*e.^2+J1;
M(3,3) = 1;
M(4,4) = 1;




% M01=-zred*d*m1*e1*wfi1*wfi1;

Mel=((2*Mrob_ut*srob_ut*((w_synchr-wfi)/w_synchr))/(srob_ut*srob_ut+((w_synchr-wfi)/w_synchr)*((w_synchr-wfi)/w_synchr)));

%% wy³¹czanie silnika 



%  ------------------------- macierz wyrazów wolnych ----------------------
Q(1,1)=m1*e*sin(fi)*wfi.^2-k*y-b*vy-(m1+m)*g ;  
Q(2,1)=Mel-m1*e*g*cos(fi);
Q(3,1)=vy;
Q(4,1)=wfi;



dX=inv(M)*Q;  
