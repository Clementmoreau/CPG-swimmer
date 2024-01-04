function T=matrixTorque(t,z,zdot,params)

% unpack the parameters.
N = params.N;
gamma=params.gamma;
L=params.L;

% z=(x,y,th,alpha_1, alpha_{N-1}
% zdot= \dot{z}
% This function returns local torque load at the i-th junction 
% for singnal function
% gamma=C_perp/C_para, the anisotropic drag ratio

% refer to Eq.(3.9) of Moreau, J R Soc Interface (2018)
% note that G11 and G12 were wrongly swapped in Eq.(3.9)
% note that \eta and \xi were wrongly swapped in Eq.(3.9)

% Caution: var rod length might not be supported!

F=zeros(N,1);

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


% Compute the torque on i-th rod, around j-th node
% j=1 for the front end of the first rod
% i=j-1 and i=j
for j=2:N % j starts at j=2
    Ftmp=zeros(N,1);
    for i=j-1:j
     G=zeros(3,3);
     Rel=zeros(1,3);
     Xdoti=[Xdot3N(i);Xdot3N(i+N);Xdot3N(i+2*N)*L(i)/sum(L)]; % L(i)->L(i)/N need change!!
     cs=cos(z3(2*N+i));
     sn=sin(z3(2*N+i));
    
     G(1,1)=0.5*sn;
     G(1,2)=-0.5*cs;
     G(1,3)=-1/3;
     G(2,1)= (1-gamma)*sn*cs;
     G(2,2)= -cs^2-gamma*sn^2;
     G(2,3)=-0.5*cs;
     G(3,1)=  gamma*cs^2+sn^2;
     G(3,2)= (-1+gamma)*sn*cs;
     G(3,3)=-0.5*sn;

    
     Rel(1)=L(i)/sum(L);
     Rel(2)=z3(i)-z3(j);
     Rel(3)=z3(N+i)-z3(N+j);

     Ftmp(i)= Rel*G*Xdoti*L(i)/sum(L); %multiplied the 'L(i)/N' factor
    end

    F(j)=sum(Ftmp);
end

% return local torque load: (N-1)-dim vector
    T=F(2:N);
    
end
    