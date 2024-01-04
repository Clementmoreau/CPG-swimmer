function Bp=cpg(t,zp,params)

% unpack the required parameters.
N = params.N;
kd=params.kd;
L=params.L;
omega=params.omega;
tau=params.tau;
coup=params.coup;
coupR=params.coupR;
sigma=params.sigma;
psi=params.psi;

z = zp(1:N+2); % geometrical state
s = z(4:N+2); % shape state
phi = zp(N+3:end); % CPG state

B1 = zeros(N+2,1); 

% the elasticity
Khat = -kd*sum(L)*eye(N-1);

% add the activity comping from the CPG oscillators
K = Khat*s+tau*cos(phi);

% call to the function that computes Sp^4*A*Q
M = matrixNparam(t,z,params);

% solving the linear system
B1(4:end) = -K ;
B = M\B1;

% The CPG part.
phiA = [diff(phi);0];
phiB = [0;-diff(phi)];

% Get the normal drag force.
Fn = matrixForceN(t,z,B,params);

% Get the curvature.

% follow Herault for the phi dynamics.
%PhiDot = omega*ones(N-1,1) + coup*sin(phiA) + coupR*sin(phiB) + sigma(t)*cos(phi+psi(t)).*(Fn(1:end-1)+Fn(2:end))/2; 
%PhiDot = omega*ones(N-1,1) + coup*sin(phiA) + coupR*sin(phiB) + sigma(t)*cos(phi+psi(t)).*Fn; 
PhiDot = omega*ones(N-1,1) + coup*sin(phiA) + coupR*sin(phiB) + sigma(t)*cos(phi+psi(t)).*Fn; 

% follow Herault + phase veloc modulation
OmegaMod = 0.5; % Choose: 0 to 1
%PhiDot=Omega*ones(N-1,1)-Omega*OmegaMod*(sin(phi).^2)+Coup*sin(phiA)+CoupR*sin(phiB)+sigma*cos(phi).*Fn; 

% another version?
%PhiDot=Omega*ones(N-1,1)+Coup*(sin(phiA)+sin(phiB))+sigma*sin(phi).*Fn;

% build the final output
Bp = [B;PhiDot];

end
