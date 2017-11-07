% Esempio 2.4+shooting modificato 

clc
clear all
close all

% funzione problema 
f=@(t,y) [y(2); 1/8*(32+2*t^3-y(1)*y(2))]; 

% estremi
a=1; 
b1=3;

% valori agli estremi 
alfa=17; 
beta=43/3; 
 
% approssimazioni iniziali per s 
s0=-1;
s1=-20; 

% derivate parziali di f 
 df1=@(y) -1/8*y(2);
 df2=@(y) -1/8*y(1); 
 
% METODO DI SHOOTING-NEWTON 
[sN,iterS] = shootingNewton(f,df1,df2,s0,a,b1,alfa,beta)
% METODO DI SHOOTING-BISEZIONE 
 [sB,iterB] = shootingbisezione(f,s0,s1,a,b1,alfa,beta)
 % METODO DI SHOOTING-SECANTI 
[sS,iterS] = shootingSecanti(f,s0,s1,a,b1,alfa,beta)
% METODO DI SHOOTING MODIFICATO 


% ESEMPIO METODO DI NEWTON 

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
w0=[alfa;sN]; 

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
   % confronto con soluzione di riferimento
   err(cc) = norm(U(1,end)-urif,inf); 
  
end

% GRAFICO ERRORE    
ord1 = err(1)*(mrange/mrange(1)).^(-1);
ord2 = err(1)*(mrange/mrange(1)).^(-2);
ord3 = err(1)*(mrange/mrange(1)).^(-3);
ord4 = err(1)*(mrange/mrange(1)).^(-4);

figure
loglog(mrange,err,'*r',mrange,ord1,mrange,ord2,mrange,ord3,mrange,ord4)
title('Ordine metodo')
legend('Ordine','h','h^2','h^3','h^4')


