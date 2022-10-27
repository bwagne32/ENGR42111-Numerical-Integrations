%% D _____________________________________________________________________
% Consts
E=100*10^6; % Pa
A=.01; % m^2
P=1*10^6; % Pa
size_t=10; % default array length
L=1; % m
a=[.25 .5 .75]*L; % m

% Preallocated vectors
x=linspace(0,L,size_t);
ux=zeros(length(a),length(x));

syms u(n); %Numeric differential equation
for i=1:3
    d2u = diff(u,n,2) == (-P/(E*A))*dirac(n-a(i)*L);
    cond1 = u(0) == 0;
    cond2 = u(L) == 0;
    conds = [cond1 cond2];
    uSol(n) = dsolve(d2u,conds);
    uSol = simplify(uSol);
    ux(i,:)=uSol(x);
end

figure(1) %Graphing u(x)
plot(x,ux)
title('Displacement vs Distance')
xlabel('Distance (m)')
ylabel('Displacement (m)')
legend('a=.25','a=.5','a=.75')

figure(2) %Graphing eps(x)
eps=100*[diff(ux(1,:));diff(ux(2,:));diff(ux(3,:))];
plot(linspace(0,L,size_t-1),eps)
title('Strain vs Distance')
xlabel('Distance (m)')
ylabel('Strain')
legend('a=.25','a=.5','a=.75')

figure(3) %Graphing sig(x)
sig=eps*E/10^6;
plot(linspace(0,L,size_t-1),sig)
title('Stress vs Distance')
xlabel('Distance (m)')
ylabel('Stress (GPa)')
legend('a=.25','a=.5','a=.75')

clear;
%% E _____________________________________________________________________
% Consts
E=100*10^6; % Pa
P=1*10^6; % Pa
size_t=75; % default array length
L=1; % m
a=[.25 .5 .75]*L; % m
b=[.05 .06 .07]; % m^2
c=[-.08 -.15 -.24]; % m calculated in part B
n=[1 2 3];

% Preallocated vectors
x=linspace(0,L,size_t);
ux=zeros(9,length(x));
A=zeros(length(a),length(x));

for i=1:3
    A(i,:)=b(i)+c(i)*(L-x).^n(i);
end

syms u(n); %Numeric differential equation
for i=1:3 % a=.25
    d2u = diff(u,n,2) == (-P/(E*A(i))).*dirac(n-a(1)*L);
    cond1 = u(0) == 0;
    cond2 = u(L) == 0;
    conds = [cond1 cond2];
    uSol(n) = dsolve(d2u,conds);
    uSol = simplify(uSol);
    ux(i,:)=uSol(x);
end
for i=1:3 % a=.5
    d2u = diff(u,n,2) == (-P/(E*A(i))).*dirac(n-a(2)*L);
    cond1 = u(0) == 0;
    cond2 = u(L) == 0;
    conds = [cond1 cond2];
    uSol(n) = dsolve(d2u,conds);
    uSol = simplify(uSol);
    ux(i+3,:)=uSol(x);
end
for i=1:3 % a=.75
    d2u = diff(u,n,2) == (-P/(E*A(i))).*dirac(n-a(3)*L);
    cond1 = u(0) == 0;
    cond2 = u(L) == 0;
    conds = [cond1 cond2];
    uSol(n) = dsolve(d2u,conds);
    uSol = simplify(uSol);
    ux(i+6,:)=uSol(x);
end

% Graphs
figure(1) %Graphing u(x) for n=1,2,3 and a(1)
plot(x,ux(1:3,:))
title('Displacement vs Distance for a=.25L')
xlabel('Distance (m)')
ylabel('Displacement (m)')
legend('n=1','a=2','a=3')
figure(2) %Graphing u(x) for n=1,2,3 and a(2)
plot(x,ux(4:6,:))
title('Displacement vs Distance for a=.5L')
xlabel('Distance (m)')
ylabel('Displacement (m)')
legend('n=1','a=2','a=3')
figure(3) %Graphing u(x) for n=1,2,3 and a(3)
plot(x,ux(7:9,:))
title('Displacement vs Distance for a=.75L')
xlabel('Distance (m)')
ylabel('Displacement (m)')
legend('n=1','a=2','a=3')

for i=1:9 % calculating eps strain values from ux
    eps(i,:)=diff(ux(i,:));
end

figure(4) %Graphing eps(x) for n=1,2,3 and a(1)
plot(linspace(0,L,size_t-1),eps(1:3,:))
title('Strain vs Distance for a=.25L')
xlabel('Distance (m)')
ylabel('Strain')
legend('n=1','n=2','n=3')
figure(5) %Graphing eps(x) for n=1,2,3 and a(2)
plot(linspace(0,L,size_t-1),eps(4:6,:))
title('Strain vs Distance for a=.5L')
xlabel('Distance (m)')
ylabel('Strain')
legend('n=1','n=2','n=3')
figure(6) %Graphing eps(x) for n=1,2,3 and a(3)
plot(linspace(0,L,size_t-1),eps(7:9,:))
title('Strain vs Distance for a=.75L')
xlabel('Distance (m)')
ylabel('Strain')
legend('n=1','n=2','n=3')

sig=eps*(E/10^6); % sigma stress calculations from eps

figure(7) %Graphing sig(x) for n=1,2,3 and a(1)
plot(linspace(0,L,size_t-1),sig(1:3,:))
title('Stress vs Distance for a=.25L')
xlabel('Distance (m)')
ylabel('Stress (GPa)')
legend('n=1','n=2','n=3')
figure(8) %Graphing sig(x) for n=1,2,3 and a(2)
plot(linspace(0,L,size_t-1),sig(4:6,:))
title('Stress vs Distance for a=.5L')
xlabel('Distance (m)')
ylabel('Stress (GPa)')
legend('n=1','n=2','n=3')
figure(9) %Graphing sig(x) for n=1,2,3 and a(3)
plot(linspace(0,L,size_t-1),sig(7:9,:))
title('Stress vs Distance for a=.75L')
xlabel('Distance (m)')
ylabel('Stress (GPa)')
legend('n=1','n=2','n=3')





