function Bp=cpg_sym(t,zp,params)

% unpack the required parameters.
N = params.N;
kd=params.kd;
L=params.L;

omega=params.omega;
tau=params.tau;
coup=params.coup;
coupR=params.coupR;

sigma_1 = params.sigma_1;
sigma_2 = params.sigma_2; 

alpha_omega = params.alpha_omega;
K_omega = params.K_omega;


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

% Get the normal drag torque.
Tn = matrixTorque(t,z,B,params)*N^2;

% follow Herault for the phi dynamics.
% (old) PhiDot = omega*ones(N-1,1) + coup*sin(phiA) + coupR*sin(phiB) + sigma(t)*cos(phi+psi(t)).*(Fn(1:end-1)+Fn(2:end))/2; 

% Intrisic frequency
phi_int = omega*ones(N-1,1);

% CPG coupling
phi_cpg = coup*sin(phiA) + coupR*sin(phiB);

% Sensory feedback for locomotion
phi_sens = sigma_1(t)*cos(phi+pi).*Tn;
% NB: can also add phi(t) if using main_CPG and psi description of
% proprioception

% Feedback for omega-turn
alpha_omega_var = alpha_omega * [-1;-1; 1;1;1;1;1; -1;-1];
psi_omega = 2*pi - acos(kd * alpha_omega_var / tau) - acos(omega / K_omega);
phi_omega = sigma_2(t) .* K_omega .* cos(phi - psi_omega);

PhiDot = phi_int + phi_cpg + phi_sens + phi_omega;

% build the final output
Bp = [B;PhiDot];

end
