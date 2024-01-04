function [X,Y,TH]=coordinates_filament(z,params)
% This function computes the '3N coordinates' -- X_3N in the text
% from the 'N+2 coordinates' -- X in the text 

% --- input : N+2 coordinates, number of links N
% --- output : X, Y coordinates of the N links
% TH orientation of each link

N = params.N;
L = params.L;

X=zeros(N+1,1);
Y=zeros(N+1,1);
TH=zeros(N,1);

X(1)=z(1);
Y(1)=z(2);
TH(1)=z(3);

for i=2:N
    X(i)=X(i-1)+L(i-1)*cos(TH(i-1))/sum(L);
    Y(i)=Y(i-1)+L(i-1)*sin(TH(i-1))/sum(L);
    TH(i)=TH(i-1)+z(i+2);
end

X(N+1)=X(N)+L(N)*cos(TH(N))/sum(L);
Y(N+1)=Y(N)+L(N)*sin(TH(N))/sum(L);

end