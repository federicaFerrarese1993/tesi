% Problema pag. 504 (esempio 2.0.1 + esempio 2.6.2) 
clc
clear all
close all

% funzione esempio
f=@(t,y) [y(2);3/2*y(1)^2]; 

% estremi
a=0;
b1=1;

% valori agli estremi  
alfa=4;
beta=1; 

% trovo lo zero s=-8
s0=-10;
df1=@(y)3*y(1);
df2=@(y) 1;
[s,iter] = shootingNewton(f,df1,df2,s0,a,b1,alfa,beta);

% soluzione analitica 
yesatta=@(x) 4./(1+x).^2; 


% TABLEAU
% ordine 4 
A=diag([1/2,1/2,1],-1);
c=[0,1/2,1/2,1];
b=[1/6,1/3,1/3,1/6];

% ordine 2 
% A=diag(2/3,-1);
% c=[0,2/3];
% b=[1/4,3/4];


% CALCOLO LA SOLUZIONE DI RIFERIMENTO

m=2.^10; % numero di nodi 

% dato iniziale
w0=[alfa;s]; 

% soluzione problema con metodo di R-K in tableau
u=rk(f,a,b1,A,b,c,m,w0);

% soluzione di riferimento
urif=u(1,end);

% CALCOLO LE SOLUZIONI CON MENO PASSI
mrange = 2.^(4:8); % numero nodi (variabile)
cc = 0; % contatore 

for m = mrange + 1
   cc=cc+1;
   % trovo le soluzioni numeriche con meno passi
   U=rk(f,a,b1,A,b,c,m,w0);
   x0 = linspace(a,b1,m);
   sol=yesatta(x0);
%    % confronto con soluzione di riferimento
%    err(cc) = norm(U(1,end)-urif,inf); 
   % confronto con soluzione analitica 
   err(cc)=norm(U(1,:)-sol,inf);
end

% GRAFICO SOLUZIONI
figure
 plot(x0,U(1,:),'r+')
 hold on
plot(x0,yesatta(x0),'b')
legend('soluzione approssimata','soluzione analitica')


% GRAFICO ERRORE   
ord1 = err(1)*(mrange/mrange(1)).^(-1);
ord2 = err(1)*(mrange/mrange(1)).^(-2);
ord3 = err(1)*(mrange/mrange(1)).^(-3);
ord4 = err(1)*(mrange/mrange(1)).^(-4);

figure
loglog(mrange,err,'*r','linewidth',3)
hold on
loglog(mrange,ord1,'k','linewidth',2)
hold on
loglog(mrange,ord2,'linewidth',2)
loglog(mrange,ord3,'linewidth',2)
loglog(mrange,ord4,'linewidth',2)
xlabel('m')
ylabel('errore')
legend('Ordine','h','h^2','h^3','h^4')























