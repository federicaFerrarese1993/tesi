% BRACHISTOCRONA-->shooting multiplo
clc
clear all
close all

% funzione problema
f1=@(t,y) [y(2);-(1+(y(2))^2)/(2*y(1))];

% condizioni ai bordi
xa=0;
ya=1;

xb=1;
yb=1.1;

% trovo i valori di a,b,r,C1
f=@(y) [y(1)*(y(3)-sin(y(3)))/2+y(2)-xa;...
                y(1)*(1-cos(y(3)))/2-ya;...
                y(1)*(y(4)-sin(y(4)))/2+y(2)-xb;...
                y(1)*(1-cos(y(4)))/2-yb];
            
              
y0=[1;-1;2.5;3.6];

y=fsolve(f,y0);

r=y(1);
% intervallo temporale 
a=y(3);
b=y(4);  
% estremi y 
alfa=ya; 
beta=yb; 

% soluzione analitica in forma parametrica 
xesatta=@(t) r*(t-sin(t))/2+y(2); 
yesatta=@(t) r*(1-cos(t))/2;


% confronto tra soluzione analitica e soluzione calcolata con il metodo
% di shooting 

mrange=20:25;
%  mrange=10:2:20;
% solo con fsolve solo con Eulero esplicito 

% ordine 2 
A=diag(2/3,-1);
c=[0,2/3];
b1=[1/4,3/4];

counter=0;

for m=mrange+1
     counter=counter+1;
     s0=ones(1,m);
    s=solveMulti(f1,xa,xb,alfa,beta,m);

     x0=linspace(a,b,m*(m-1));

     x2=linspace(xa,xb,m);
    x3=linspace(a,b,m);
    % soluzione del metodo usando i valori ottenuti con il metodo di shooting e
    % il metodo di Newton 

    j=1;
    g=1;
    g1=m;
 for k=1:m-1

    y0=[s(j);s(j+1)];
    x1=linspace(x2(k),x2(k+1),m);
    x11=linspace(x3(k),x3(k+1),m);
    h=(x2(k+1)-x2(k))/(m-1);
    y1=zeros(2,m);
    y1(:,1)=y0;
    for n=1:m-1
        y1(:,n+1)=y1(:,n)+h*f1(1,y1(:,n));
    end

    Y1(g:g1)=y1(1,:);
    X(g:g1)=x1;
    g=g1+1; 
    g1=g+m-1;
    j=j+2;
 end
 err(counter)=norm(yesatta(x11)-y1(1,:),inf);
 %  err(counter)=norm(yesatta(x0)-Y1,inf);
end

% confronto soluzione analitica e numerica 
figure 
plot(X,Y1,'ro','linewidth',2)
  hold on
    plot(xesatta(x0),yesatta(x0),'b','linewidth',3)
axis equal
legend('soluzione numerica','soluzione analitica')

% grafico errore 

ord1 = err(1)*(mrange/mrange(1)).^(-1);
ord2 = err(1)*(mrange/mrange(1)).^(-2);
ord3 = err(1)*(mrange/mrange(1)).^(-3);
ord4 = err(1)*(mrange/mrange(1)).^(-4);

figure
loglog(mrange,err,'*r',mrange,ord1,mrange,ord2,mrange,ord3,mrange,ord4)
title('Ordine metodo')
legend('Ordine','h','h^2','h^3','h^4')

