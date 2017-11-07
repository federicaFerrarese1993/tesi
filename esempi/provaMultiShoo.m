% Esempio 2.4+esempio 3.3.1

% prova shooting semplice, shooting multiplo, shooting modificato,
% shooting con fsolve 

% confronto tra shooting con metodo di Newton e shooting con fsolve
% confronto tra metodo di shooting  e shooting modificato (mms) 
clc
clear all
close all

m=7; % numero di nodi (se m=1 shooting semplice) 
n=10;
% funzione problema 
f=@(t,y) [y(2); 1/8*(32+2*t^3-y(1)*y(2))]; 
 
% soluzione esatta 
yesatta=@(x)x.^2+16./x; 

alfa=17; 
beta=43/3; 

a=1; 
b=3; 
df1=@(y) -1/8*y(2);
df2=@(y) -1/8*y(1); 

s0=ones(1,m)*(beta-alfa)/(b-a);

% multi shooting con fsolve  (per ogni m) 
tic
s=solveMulti(f,a,b,alfa,beta,m)
toc
% multi shooting modificato (per ogni m) 
[sM,iterM]=mssm(f,a,b,alfa,beta,df1,df2,n);


% caso m>1 (multi shooting)
if m>1
    % multi shooting con Newton  (con m<=8) 
    tic
    [sN,iterN]=multiNew(f,a,b,alfa,beta,df1,df2,s0,m)
    toc
    % intervallo temporale per soluzione analitica 
    x0=linspace(a,b,(m-1)*m);

     x2=linspace(a,b,m);

    % soluzione del metodo usando i valori ottenuti con il metodo di shooting e
    % il metodo di Newton 

    j=1;
for k=1:m-1
         % metodo con Newton 
         
         y0=[sN(j);sN(j+1)];
         x1=linspace(x2(k),x2(k+1),m);
         [x,y]=ode45(f,x1,y0);

         hold on
         plot(x,y(:,1),'bo','linewidth',3)

         % metodo con fsolve 
         y0=[s(j);s(j+1)];
         x1=linspace(x2(k),x2(k+1),m);
         [x,y]=ode45(f,x1,y0);

         hold on
         plot(x,y(:,1),'r*','linewidth',2)
         hold on
     plot(x0,yesatta(x0),'g','linewidth',2)
    legend('Metodo di Newton','Algoritmo di trust-region','Soluzione analitica')
    xlabel('x')
 ylabel('y(x)')    
    j=j+2;
end
%      % metodo di shooting modificato 
%      y0=[alfa;sM];
%      [x,y]=ode45(f,x0,y0);
%      hold on
%      plot(x,y(:,1),'k+')
end


 % caso m=1 (shooting semplice) 
 if m==1
     
   
     x0=linspace(a,b,20);
     s0=1;
     % metodo di shooting semplice
     n=10;
     [sN,iter] = shootingNewton(f,df1,df2,s0,a,b,alfa,beta)
     
     % metodo con Newton
     y0=[alfa;sN];
     [x,y]=ode45(f,x0,y0);
     hold on 
     plot(x,y(:,1),'bo')
     % metodo con fsolve
      y0=[s(1);s(2)];
     [x,y]=ode45(f,x0,y0);
     hold on
     plot(x,y(:,1),'r*')
     % metodo di shooting modificato 
    y0=[alfa;sM];
    [x,y]=ode45(f,x0,y0);
    hold on
    plot(x,y(:,1),'k+')
 end
hold on
 plot(x0,yesatta(x0),'g')
 legend('Metodo di Newton','Algoritmo di trust-region','Soluzione analitica')