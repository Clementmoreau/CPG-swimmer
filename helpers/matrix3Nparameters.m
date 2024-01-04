function M=matrix3Nparameters(t,z,params)
% This function fills the matrix defined as A in the text 
% (see Appendix VII-C, equation (20) and following)
% and returns Sp^4 * A

% unpack the required parameters.
N = params.N;
gamma=params.gamma;
Sp=params.Sp;
L=params.L;

x=z(1:N);
y=z(N+1:2*N);
th=z(2*N+1:3*N);
F=zeros(2,3*N);
T=zeros(N,3*N);

for i=1:N
    u=cos(th(i));
    v=sin(th(i));
    F(1,i)=-L(i)*(gamma*u^2+v^2);
    F(1,N+i)=-L(i)*(gamma-1)*u*v;
    F(2,i)=-L(i)*(gamma-1)*u*v;
    F(2,N+i)=-L(i)*(u^2+gamma*v^2);
    F(1,2*N+i)=L(i)*L(i)*1/2*v/sum(L);
    F(2,2*N+i)=-L(i)*L(i)*1/2*u/sum(L);
end

F=Sp^4*F/sum(L);

for i=1:N
    for j=i:N
        u=cos(th(j));
        v=sin(th(j));
        A=(x(j)-x(i));
        B=(y(j)-y(i));
        T(i,j)=L(j)*L(j)*1/2*v/sum(L)+...
            L(j)*A*(-gamma+1)*v*u+...
            L(j)*B*(gamma*u*u+v*v);
        T(i,N+j)=-L(j)*L(j)*1/2*u/sum(L)+...
            L(j)*B*(gamma-1)*v*u-...
            L(j)*A*(u*u+gamma*v*v);
        T(i,2*N+j)=-L(j)*L(j)*L(j)*1/3/sum(L)/sum(L)-...
            L(j)*L(j)*1/2*A*u/sum(L)...
            -L(j)*L(j)*1/2*B*v/sum(L);
    end
end

T=Sp^4*T/sum(L);

M=[F;T];
end