%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INICIALIZACION

% Limpiamos memoria y graficas
clear all
close all

% Liquido de trabajo y gravedad
rho = 1e3;
mu  = 1e-3;
g   = 9.81;

% Correlacion de Colebrook
colebrook = @(Re,rr,la) 1/sqrt(la)+2*log10(rr/3.71+2.51/Re/sqrt(la));

% Embalses
hE1 = 200;
hE2 = 200;
hE3 = 500;

%Valvulas

K5 = 10;
K8 = 15;
K13 = 10;

% Elementos
L1 = 600;
L2 = 400;
L3 = 900;
L4 = 400;
L5 = 400;
L6 = 350;
L7 = 600;
L8 = 400;
L9 = 300;
L10 = 500;
L11 = 300;
L12 = 350;
L13 = 600;
L14 = 300;

D1 = 0.5;
D2 = 0.5;
D3 = 0.5;
D4 = 0.5;
D5 = 0.4;
D6 = 0.5;
D7 = 0.5;
D8 = 0.4;
D9 = 0.3;
D10 = 0.3;
D11 = 0.3;
D12 = 0.3;
D13 = 0.3;
D14 = 0.3;





v_L  = [L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11 L12 L13 L14];
v_D  = [D1 D2 D3 D4 D5 D6 D7 D8 D9 D10 D11 D12 D13 D14];
rugCem = 1e-3; 
v_rr = rugCem./v_D;



% Matriz de conectividad

E = [1 -1 0 0 0 0 0 0 0; 
     1 0 -1 0 0 0 0 0 0; 
     1 0 0 -1 0 0 0 0 0;
     0 1 0 -1 0 0 0 0 0;
     0 1 0 0 -1 0 0 0 0;
     0 0 1 0 0 -1 0 0 0;
     0 0 0 1 0 -1 0 0 0;
     0 0 0 1 -1 0 0 0 0;
     0 0 1 0 0 0 -1 0 0;
     0 0 0 0 0 1 -1 0 0;
     0 0 0 0 0 1 0 -1 0;
     0 0 0 0 0 0 1 -1 0;
     0 0 0 0 0 0 0 1 -1;
     0 0 0 1 0 0 0 0 -1];
 
 % Vector auxiliar consumo=0, altura=1
 
 y_ch = [0 1 1 0 1 0 0 0 0];
 
% Factores de friccion. Comenzamos suponiendo rugosidad dominante
v_la = 0.25./(log10(v_rr/3.71)).^2;


%Creando vector de alturas piezometricas
Patm=101325;

H=[0 200 200 0 500 0 0 0 0];
C=[-0.2 0 0 0 0 0 0 -0.5 -0.3];



%Creando vector de perdidas secundarias

K_valv=[0 0 0 0 10 0 0 15 0 0 0 0 10 0];


%MAIN

%inicializacion de Q
Q=[1 1 1 1 1 1 1 1 1 1 1 1 1 1];
Kper=K_perdidas (v_la, v_D, v_L, K_valv, Q);
M=MatrizM(E, Kper);
A=MatrizA(y_ch,M);
b=vectorB(y_ch,M,H,C);
v_sol=inv(A)*transpose(b);
H=redefH(v_sol,y_ch,H);
C=redefC(v_sol,y_ch,C);
Q=recalcQ(Kper,E,H);
Msol(:,1)=v_sol;
Mq(:,1)=Q;

% datos de diagnostico:
maxiter= 1000;
tol=1e-9;
inc=1;
res=1;
iter=0;
coef=0.5; %Hemos ajustado el coeficiente para optimizar el numero de iteraciones

while(inc>tol) && (res>tol) && (iter<maxiter)
    iter=iter+1;
    v_sol_old=v_sol;
    Re=defRe(Q,v_D);
    v_la=landaColebrook(Re,v_rr,v_la);
    Kper=K_perdidas (v_la, v_D, v_L, K_valv, Q);
    M=MatrizM(E,Kper);
    A=MatrizA(y_ch,M);
    b=vectorB(y_ch,M,H,C);
    v_sol=inv(A)*transpose(b);
    H=redefH(v_sol,y_ch,H);
    C=redefC(v_sol,y_ch,C);
    
    QANT = Q;
    
    Q=recalcQ(Kper,E,H);
    
    Q=Q*coef + QANT*(1-coef);
    
    inc = norm((v_sol-v_sol_old)./v_sol,inf);
    v_inc(iter) = inc;
    
    %Almacenando valores
    
    Msol(:,iter+1)=v_sol;
    Mq(:,iter+1)=Q;
    
    %Calculo del residuo

    v_res=E'*Q-C';
    res  = norm(v_res,inf);
    resPp(iter)=res;
end

figure(1)
 
subplot(7,2,1)
plot(1:iter+1,Mq(1,:))
box on,grid on
xlabel('iteracion')
ylabel('Q1')

subplot(7,2,2)
plot(1:iter+1,Mq(2,:))
box on,grid on
xlabel('iteracion')
ylabel('Q2')

subplot(7,2,3)
plot(1:iter+1,Mq(3,:))
box on,grid on
xlabel('iteracion')
ylabel('Q3')

subplot(7,2,4)
plot(1:iter+1,Mq(4,:))
box on,grid on
xlabel('iteracion')
ylabel('Q4')

subplot(7,2,5)
plot(1:iter+1,Mq(5,:))
box on,grid on
xlabel('iteracion')
ylabel('Q5')

subplot(7,2,6)
plot(1:iter+1,Mq(6,:))
box on,grid on
xlabel('iteracion')
ylabel('Q6')

subplot(7,2,7)
plot(1:iter+1,Mq(7,:))
box on,grid on
xlabel('iteracion')
ylabel('Q7')

subplot(7,2,8)
plot(1:iter+1,Mq(8,:))
box on,grid on
xlabel('iteracion')
ylabel('Q8')

subplot(7,2,9)
plot(1:iter+1,Mq(9,:))
box on,grid on
xlabel('iteracion')
ylabel('Q9')

subplot(7,2,10)
plot(1:iter+1,Mq(10,:))
box on,grid on
xlabel('iteracion')
ylabel('Q10')

subplot(7,2,11)
plot(1:iter+1,Mq(11,:))
box on,grid on
xlabel('iteracion')
ylabel('Q11')

subplot(7,2,12)
plot(1:iter+1,Mq(12,:))
box on,grid on
xlabel('iteracion')
ylabel('Q12')

subplot(7,2,13)
plot(1:iter+1,Mq(13,:))
box on,grid on
xlabel('iteracion')
ylabel('Q13')

subplot(7,2,14)
plot(1:iter+1,Mq(14,:))
box on,grid on
xlabel('iteracion')
ylabel('Q14')

figure(2)

subplot(7,2,1)
plot(1:iter+1,Msol(1,:))
box on,grid on
xlabel('iteracion')
ylabel('H1')

subplot(7,2,2)
plot(1:iter+1,Msol(2,:))
box on,grid on
xlabel('iteracion')
ylabel('C2')

subplot(7,2,3)
plot(1:iter+1,Msol(3,:))
box on,grid on
xlabel('iteracion')
ylabel('C3')

subplot(7,2,4)
plot(1:iter+1,Msol(4,:))
box on,grid on
xlabel('iteracion')
ylabel('H4')

subplot(7,2,5)
plot(1:iter+1,Msol(5,:))
box on,grid on
xlabel('iteracion')
ylabel('C5')

subplot(7,2,6)
plot(1:iter+1,Msol(6,:))
box on,grid on
xlabel('iteracion')
ylabel('H6')

subplot(7,2,7)
plot(1:iter+1,Msol(7,:))
box on,grid on
xlabel('iteracion')
ylabel('H7')

subplot(7,2,8)
plot(1:iter+1,Msol(8,:))
box on,grid on
xlabel('iteracion')
ylabel('H8')

subplot(7,2,9)
plot(1:iter+1,Msol(9,:))
box on,grid on
xlabel('iteracion')
ylabel('H9')

figure(3)

subplot(2,2,1)
semilogy(1:iter,resPp)
box on,grid on
xlabel('iteracion')
ylabel('res')

subplot(2,2,2)
semilogy(1:iter,v_inc)
box on,grid on
xlabel('iteracion')
ylabel('inc')


%Funcion para crear matriz K de perdidas totales

function f=K_perdidas (landa, D, L, K_valv, Q)
    for i=1:14
        g=9.81;
        Kper(i)=((2*g*((pi/4)*(D(i))^2)^2)^-1)*(landa(i)*L(i)/D(i) + K_valv(i))*abs(Q(i));
        
    end
    
    f=diag(Kper);

end

        
%Creando funcion delta
function f=deltaFunc(n,m)
    if n==m
        f=1;
    else
        f=0;
    end
end


%Creando funcion vector nuevo landa
function landa=landaColebrook(Re,rr,la)
 for i=1:14
    landa(i)=(-2*log10((rr(i)/3.71)+(2.51/(Re(i)*sqrt(la(i))))))^(-2);
    %landa(i)=0.25/(log10((rr(i)/3.71)+(2.51/(Re(i)*sqrt(la(i))))))^2; 
 end
end

%Creando matriz M
function f=MatrizM(E,K)
 f=transpose(E)*inv(K)*E;
end 

%Creando Matriz A

function A=MatrizA(Y,M)
    for i=1:9
        if Y(i)==1
            for j=1:9
                A(j,i)=-deltaFunc(j,i);
            end
        else
            for j=1:9
                A(j,i)=M(j,i);
            end
        end
    end
    
end


%Creando vector B

function b=vectorB(Y,M,H,C)
    for i=1:9
        b(i)=0;
        for j=1:9
            b(i)=b(i)-Y(j)*M(i,j)*H(j);
        end
        
        if Y(i)==0
        b(i)=b(i)+C(i);
        end
    end
end


%Creando funcion recalcular Q
function Q=recalcQ(K,E,H)
    Q=inv(K)*E*transpose(H);
end


%Creando funcion para crear Reynolds
function Re=defRe(Q,D)
    rho = 1e3;
    mu  = 1e-3;
    for i=1:14
        Re(i)=(4*rho*abs(Q(i)))/(pi*D(i)*mu);
    end 
   
end

%Creando funcion redefinir H
function f=redefH(sol,Y,H)
    for i=1:9
        if Y(i)==0
            f(i)=sol(i);
        else
            f(i)=H(i);
        end
    end 
end


%Creando funcion redefinir C
function f=redefC(sol,Y,C)
    for i=1:9
        if Y(i)==1
            f(i)=sol(i);
        else
            f(i)=C(i);
        end
    end 
end
