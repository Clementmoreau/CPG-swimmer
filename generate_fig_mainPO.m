function [tps, traj]=generate_fig_mainPO(params, inits, opts)


% Motion of a swimmer driven by CPG

addpath(genpath('./helpers'))


% Unpack the physical parameters in the params structure.
N= params.N;
gamma=params.gamma;
Sp=params.Sp;
kd=params.kd;
L=params.L;
omega=params.omega;
tau=params.tau;
coup=params.coup;
coupR=params.coupR;
sigma=params.sigma;
psi=params.psi;

% Unpack the options from the opts structure.

t_switch=opts.t_switch;
switch_width=opts.switch_width;
T=opts.T;
tpnum=opts.tpnum;

% load initial data.
init=inits;

%% Set the parameters.

%sigma_amp=sigma;


% Simulation options
tps=linspace(0,T,tpnum); % time step vector
option = odeset('RelTol',1e-6,'AbsTol',1e-6);

% Input Optiions.
Isym=1; % 0: asym, 1:sym for the proprioceptice control


%Initial condition. (x,y,theta,alpha_1...alpha_N-1)

pp0= [z0(4:N+2);squeeze(p0(2:N-1))-p0(1)]; % shape states (i=1,2,...,N-1)+ CPG states (i=2,3,...,N-1)

% Simulation options
T=10; % final time (enough larger than 1)
if Isym==1
    dZ=@(t,z) cpg_sym(t,z,params);
else
    dZ=@(t,z) cpg(t,z,params);   
end
options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Events',@poincare_event_func);

dp=1e-4;

iter=1;
err=1;
errMax=1e-4;

while err>errMax && err <1e4 && iter<=40


% Poincare map: F:X ->X with X=pp0
[tp,Zp,Te,Ze,ie]=ode15s(dZ,[0 T],init, options);

ppmap=[Ze(4:N+2)';Ze(N+4:2*N+1)'-Ze(N+3)]; %F(X)
Tppe=Te;

% compute Jacobian matrix of Poincare map
Jabobi=zeros(length(ppmap));
for i=1:N-1 % derivative with shape angles (1,3,...N-1)
    
    ppJ=zeros(2*N+1,1);
    ppJ(i+3)=dp;
    initJ=init+ppJ;
    [tp,Zp,Te,Ze,ie]=ode15s(dZ,[0 T],initJ, options);
    ppmapJ=[Ze(4:N+2)';Ze(N+4:2*N+1)'-Ze(N+3)];
    Jacobi(:,i)=(ppmapJ-ppmap)/dp;
    Jacobi(i,i)=Jacobi(i,i)-1;
end
for i=2:N-1 % derivative with CPG phase angles (2,3,...N-1)
    
    ppJ=zeros(2*N+1,1);
    ppJ(N+2+i)=dp;
    initJ=init+ppJ;
    [tp,Zp,Te,Ze,ie]=ode15s(dZ,[0 T],initJ, options);
    ppmapJ=[Ze(4:N+2)';Ze(N+4:2*N+1)'-Ze(N+3)];
    Jacobi(:,i+N-2)=(ppmapJ-ppmap)/dp;
    Jacobi(i+N-2,i+N-2)=Jacobi(i+N-2,i+N-2)-1;
end

pp0next=pp0-Jacobi\(ppmap-pp0); % Xnext=X-J^{-1}(F(X)-X)

err=norm(pp0next-pp0);
[iter err]


pp0=pp0next;
init=[0;0;0;pp0(1:N-1);0;pp0(N:2*N-3)];
iter=iter+1;


end % end of while


if err<errMax
    disp('Newton-Raphson Scheme converged');
    Tppe
    [eigV,eigD]=eig(Jacobi);
    Lyap=diag(eigD);

    % plot eigen values
    figsize = 400;
    fig4=figure(4);clf;
    set(gcf, 'Position',  [figsize, 640, figsize, figsize])
    hold on

    scatter(real(Lyap),imag(Lyap),50,...
              'MarkerEdgeColor',[0  .4 .7],...
              'MarkerFaceColor',[0  .7  1],...
              'LineWidth',1.5)
    hold off
    xlabel('Re[$\lambda$]','FontSize',14,'Interpreter','latex');
    ylabel('Im[$\lambda$]','FontSize',14,'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex') 
    axis equal
    grid on    
    set(gca,'FontSize',20)
    
    %save('init_guess','init');
    
else
    disp('Newton-Raphson Scheme did not converge');
    %Itraj=0; % 1: to visualise trajecotry
    IMovie=0; % 1: to generate a movie
end



%% Solve the equation for periodic orbit
tpnum=100; % number of time steps
tps=linspace(0,Tppe,tpnum); % time step vector
opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
[tps,traj]=ode15s(dZ,tps,init,opts);

end
