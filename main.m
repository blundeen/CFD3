%Brandon Lundeen
clear
clc
close all

%inputs
N=31;
CFL=.5;
s=1400;


empty3=zeros(3,N);



gamma=1.4;

x=linspace(0,3,N);
dx=x(2)-x(1);

A=1+2.2*(x-1.5).^2;
dA=4.4*(x-1.5)./(1+2.2*(x-1.5).^2);

%initial cinditions

for c=1:length(x)
    
   if x(c)<.5
     rho(c)=1;
     T(c)=1;
   
   elseif x(c)<1.5
       rho(c)=1-3.66*(x(c)-.5);
       T(c)=1-0.167*(x(c)-.5);
       
   elseif x(c)<2.1
       rho(c)=.634-.0702*(x(c)-1.5);
       T(c)=.833-.4908*(x(c)-1.5);
       
       
   elseif x(c)<3.1
       rho(c)=.5896+.10228*(x(c)-2.1);
       T(c)=.93968+.0622*(x(c)-2.1);
       
   end
    
end

V=.59./(rho.*A);

all=struct('rho',rho,'T',T,'V',V,'U',empty3, 'F', empty3,'J', empty3);
dt=getdt(T,V,dx,CFL);

all.U=encodeU(rho,A,V,T, dA);
[all.F, all.J]=encode(all.U, dA);
data(1)=all;

%McCormack solver
