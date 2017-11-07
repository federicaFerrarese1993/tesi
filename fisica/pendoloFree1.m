% PENDOLO CON FRONTIERA LIBERA
% DESCRIZIONE: trovo il tempo e la velocità iniziale del seguente
% problema ai limiti a frontiera libera mediante i metodi di shooting a
% frontiera libera. 
clc
clear all 
close all

l=1;
g=9.8;

% aprossimazione lineare valida quando l'angolo iniziale e' piccolo
f=@(t,y) [y(2); -g/l*y(1)];

% estremi in cui F(s) cambia di segno 
% (intervallo in cui F ha uno solo zero) 
s0=pi/2*sqrt(l/g); % istante finale approssimazione 1 
s1=pi*sqrt(l/g); % istante finale approssimazione 2

a=2/3*pi; % istante iniziale

alfa= pi; % valore della soluzione all'istante finale (angolo finale)
beta=0; % valore della velocita' all'istante finale 
gamma=pi/100; % angolo iniziale piccolo


% METODO DI BISEZIONE 
tic
[sB,vB,iterB]= shootingBisFree(f,alfa,beta,gamma,s0,s1,a);
toc

% METODO DI NEWTON
df1=@(y) -g/l;
df2=@(y) 1;
% dato iniziale 2*s0 vicino alla radice

% METODO DELLE SECANTI
tic
 [sN,vN,iterN]=shootingNewFree(f,df1,df2,alfa,beta,gamma,2*s0,a);
toc
tic 
 [sS,vS,iterS]= shootingSecFree(f,s0,s1,alfa,beta,gamma,a);
toc

% SOLUZIONE NUMERICA 

% tableau
% ordine 2 
A=diag(2/3,-1);
c=[0,2/3];
b1=[1/4,3/4];

mrange=20:10:50;
counter=0; 

for m=mrange+1
    counter=counter+1;
    t=linspace(a,sS,m); 
    
    % SOLUZIONE ANALITICA 
    omega=sqrt(g/l);

    f1=@(y) [y(1)*cos(omega*a+y(2))-gamma;
        -y(1)*omega*sin(omega*a+y(2))-vS];
    y0=[1;-1];
    y=fsolve(f1,y0);

    sol=@(t) y(1)*cos(omega*t+y(2));
    
    % SOLUZIONE NUMERICA 
    y0=[gamma;vS];
    Y=rk(f,a,sS,A,b1,c,m,y0);
    err(counter)=norm(sol(t)-Y(1,:),inf);

end

% GRAFICI SOLUZIONE ANALITICA E NUMERICA 
plot(t,Y(1,:),'ro','linewidth',2)
hold on 
plot(t,sol(t),'b','linewidth',3)
legend('soluzione numerica', 'soluzione analitica')
xlabel('x')
ylabel('y(x)')

% GRAFICO ERRORE  
ord1 = err(1)*(mrange/mrange(1)).^(-1);
ord2 = err(1)*(mrange/mrange(1)).^(-2);
ord3 = err(1)*(mrange/mrange(1)).^(-3);
ord4 = err(1)*(mrange/mrange(1)).^(-4);

figure
loglog(mrange,err,'*r',mrange,ord1,mrange,ord2,mrange,ord3,mrange,ord4)
title('Ordine metodo')
legend('Ordine','h','h^2','h^3','h^4')
xlabel('m')
ylabel('errore')

