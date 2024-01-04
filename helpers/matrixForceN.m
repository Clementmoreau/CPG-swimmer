function F=matrixForceN(t,z,zdot,params)

% unpack the parameters.
N = params.N;
gamma=params.gamma;
L=params.L;

% z=(x,y,th,alpha_1, alpha_{N-1}
% zdot= \dot{z}
% This function fills the matrix F_i as defined in the text
% then calculate the normal component F_N at each segment
% Eq.(3.7) of Moreau et al. Interface (2018)

% Caution: var rod length might not be supported! (not sure about this
% comment -- CM 10/2023)

F=zeros(N-1,1);

% Recover the full coordinates.
z3=zeros(3*N,1);
z3(1)=z(1);
z3(N+1)=z(2);
z3(2*N+1)=z(3);
for i=2:N
    z3(2*N+i)=z3(2*N+i-1)+z(i+2);
    z3(i)=z3(i-1)+L(i-1)*cos(z3(2*N+i-1))/sum(L);
    z3(N+i)=z3(N+i-1)+L(i-1)*sin(z3(2*N+i-1))/sum(L);
end

% Recover the full velocities.
C1=zeros(N,N);C2=zeros(N,N);C3=zeros(N,N);
for i=1:N
    for j=1:N
        C3(i,j)=i>=j;
    end
end
for i=2:N
    C1(i,i-1)=-L(i-1)*sin(z3(2*N+i-1));
    C2(i,i-1)=L(i-1)*cos(z3(2*N+i-1));
end
for i=N:-1:3
    for j=i-2:-1:1
        C1(i,j)=C1(i,j+1)-L(j)*sin(z3(2*N+j));
        C2(i,j)=C2(i,j+1)+L(j)*cos(z3(2*N+j));
    end
end
C=[ones(N,1),zeros(N,1),C1/sum(L);zeros(N,1),ones(N,1),C2/sum(L);zeros(N,2),C3];

Xdot3N=C*zdot;

% Compute the normal force on each segment using RFT. 

% Update Jan 2024 : added /sum(L) in Xdoti calculation.
% Not useful anymore for signal function of CPG model, use matrixTorque instead. 

for i=1:N-1
    Lamb=zeros(3,2);
    Xdoti=[Xdot3N(i);Xdot3N(i+N);Xdot3N(i+2*N)*L(i)/sum(L)];
    cs=cos(z3(2*N+i+1));
    sn=sin(z3(2*N+i+1));
    Lamb(1,1)=-gamma*cs^2-sn^2;
    Lamb(1,2)=-(gamma-1)*cs*sn;
    Lamb(2,1)=Lamb(1,2);
    Lamb(2,2)=-cs^2-gamma*sn^2;
    Lamb(3,1)=0.5*sn;
    Lamb(3,2)=-0.5*cs;
    Fi=Lamb'*Xdoti; %L(i)*Lamb'*Xdoti;
    F(i)=+Fi(1)*sn-Fi(2)*cs;
end
    
end