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
%natural log derivative
dA=4.4*(x-1.5)./(1+2.2*(x-1.5).^2);

%initial cinditions

for c=1:length(x)
    
   if x(c)<.5
     rho(c)=1;
     T(c)=1;
   
   elseif x(c)<1.5
       rho(c)=1-.366*(x(c)-.5);
       T(c)=1-0.167*(x(c)-.5);
       
       
       
   elseif x(c)<3.5
       rho(c)=.634-.3879*(x(c)-1.5);
       T(c)=.833-.3507*(x(c)-1.5);
       
   end
    
end

V=.59./(rho.*A);

all=struct('rho',rho,'T',T,'V',V,'M',V./sqrt(T),'U',empty3, 'F', empty3,'J', empty3);
dt=getdt(T,V,dx,CFL);

all.U=encodeU(rho,A,V,T, dA);


for z=1:100 %will be total iterations
[all.F, all.J]=encode(all.U, dA);
data(z)=all;

%McCormack solver

for c=1:N-1
   Up(:,c)=-(all.F(:,c+1)-all.F(:,c))/dx-all.J(:,c); 
   Ubar(:,c)=all.U(:,c)+Up(:,c)*dt;
   
end

[rhobar,Tbar,Vbar]=decodeU(Ubar, A(1:end-1));
Fbar=encode(Ubar, dA(1:end-1), A(1:end-1));

for c=2:N-1
   Upbar(:,c)=-(Fbar(:,c)-Fbar(:,c-1))/dx-all.J(:,c);    
   
end

Upavg=.5*(Upbar+Ubar);

Upavg(:,end+1)=[0,0,0];

all.U=all.U+Upavg*dt;

%BC 1
[trash,trash,V2]=decodeU(all.U(:,2),A(2));
[trash,trash,V3]=decodeU(all.U(:,3),A(3));

V1=2*V2-V3;
all.U(:,1)=encodeU(1,A(1),V1,1,dA);

%BC end
[rhoN2,TN2,VN2]=decodeU(all.U(:,end-2),A(end-2));
[rhoN1,TN1,VN1]=decodeU(all.U(:,end-1),A(end-1));

VN=2*VN1-VN2;
rhoN=2*rhoN1-rhoN2;
TN=2*TN1-TN2;
all.U(:,end)=encodeU(rhoN,A(end),VN,TN,dA);


%update data
[all.rho,all.T,all.V]=decodeU(all.U,dA);
plot(x,all.U(1,:))
titl=sprintf('Mass flow rate (iteration %i)',z);
title(titl);
pause(.1)

end

[rho,T,V]=decodeU(all.U,T);



