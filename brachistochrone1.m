% brachistochrone problem
clc
clear all
close all


f1=@(t,y) [y(2);-(1+(y(2))^2)/(2*y(1))];

xa=0;
ya=1;

xb=1;
yb=1;

f=@(y) [y(1)*(y(3)-sin(y(3)))/2+y(2)-xa;...
                y(1)*(1-cos(y(3)))/2-ya;...
                y(1)*(y(4)-sin(y(4)))/2+y(2)-xb;...
                y(1)*(1-cos(y(4)))/2-yb];
            
% 
% jf=@(y) [(y(3)-sin(y(3)))/2,1,y(1)*(1-cos(y(3)))/2,0;...
%                   (1-cos(y(3)))/2,0,y(1)/2*sin(y(3)),0;...
%                   (y(4)-sin(y(4)))/2,1,0,y(1)*(1-cos(y(4)))/2;...
%                   (1-cos(y(4)))/2,0,0,y(1)*sin(y(4))/2];                  
  y0=[1;-1;2.5;3.6];

% delta=-jf(y0)\f(y0);
% tol=1e-8;
% y=y0;
% 
% 
% while norm(delta)>tol
%     y=y+delta;
%     delta=-jf(y)\f(y);
% end
% y=y+delta;
y=fsolve(f,y0)
r=y(1);
% intervallo temporale 

a=y(3);
b=y(4);  

% estremi y 
alfa=ya;%alfa=1; 
beta=yb; %beta=2*r+1; 
% options=bvpset('RelTol',1e-8);
% sol0=bvpinit(linspace(xa,xb),[1,0]);
% bc=@(za,zb)[za(1)-alfa;zb(1)-beta];
% sol=bvp4c(f1,bc,sol0,options);


% intervallo temporale per il metodo di bisezione in cui f cambia segno 


s0=1;
s1=5;

% % metodo di bisezione 
% [sb,iter] = shootingbisezione(f1,s0,s1,a,b,alfa,beta)
% % metodo secanti
% [ss,iter] = shootingSecanti(f1,s0,s1,a,b,alfa,beta)

% metodo di Newton
 
df1=@(y) (1+y(2)^2)/(2*y(1)^2);
df2=@(y) y(2)/y(1);
[sN,iter] = shootingNewton(f1,df1,df2,s0,xa,xb,alfa,beta)


% soluzione analitica in forma parametrica 
xesatta=@(t) r*(t-sin(t))/2+y(2); 
yesatta=@(t) r*(1-cos(t))/2;


% confronto tra soluzione analitica e soluzione calcolata con il metodo
% di shooting 

mrange=10:10:100;
iter=0;


for m=mrange
    iter=iter+1;
    t1=linspace(xa,xb,m)';
    t2=linspace(a,b,m)'; % 0<= t <=2*pi
    
    y0=[alfa;sN];
     
    Y1(:,1)=y0;
    
    k=(xb-xa)/(m-1); 
    
    t=xa;
   for n=1:m-1
       Y1(:,n+1)=Y1(:,n)+k*f1(t,Y1(:,n));
       t=t+k;
   end
    [x,Y]=ode45(f1,t1,y0);
    err(iter)=norm(yesatta(t2)-Y(:,1), inf); 
end

% grafico soluzione analitica piu soluzione approssimata 
figure
plot(t1,Y(:,1),'r*')
hold on
plot(xesatta(t2),yesatta(t2),'k')
hold on
plot(t1,Y1(1,:),'b')
% hold on
% plot(sol.x,sol.y(1,:),'+')
% hold off

axis square
legend('soluzione approssimata ode45','soluzione analitica','soluzione approssimata Eulero')


% grafico errore 
    
ord1 = err(1)*(mrange/mrange(1)).^(-1);
ord2 = err(1)*(mrange/mrange(1)).^(-2);
ord3 = err(1)*(mrange/mrange(1)).^(-3);
ord4 = err(1)*(mrange/mrange(1)).^(-4);

figure
loglog(mrange,err,'*r',mrange,ord1,mrange,ord2,mrange,ord3,mrange,ord4)
title('Ordine metodo')
legend('Ordine','h','h^2','h^3','h^4')

