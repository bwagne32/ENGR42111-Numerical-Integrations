%% D _____________________________________________________________________
% Consts
E=100*10^6; % Pa
A=.01; % m^2
P=1*10^6; % Pa
size_t=100000; % default array length
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





