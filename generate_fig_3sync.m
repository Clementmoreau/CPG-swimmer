function out=generate_fig_3sync()

clear all

% Cstar data 
% Filament parameters.
%N = 3; % number of links
%gamma = 1/2; % ratio between the hydrodynamic drag coefficients 
%Sp = 4; % Sperm number (typical range 1-10)
%kd = 1; % Bending stiffness (default = 1)
%L = ones(1,N); % Lengths of the segments (default = 1)

% CPG model parameters.
%omega = 2*pi; % T=1 unit
%tau = --; % Strength of the activity 
%coup = --; % Oscillators coupling strength 
%coupR = coup; % Differential coupling
%sigma_amp = 4; % Proprioceptive strength minus=phase accel


Cstar=[5.8, 0;...
       6.0, 0.037;...
       7.0, 0.097;...
       8.0, 0.142;...
       9.0, 0.181;...
      10.0, 0.215;...
      11.0, 0.245;...
      12.0, 0.273;...
      13.0, 0.299;...
      14.0, 0.323;...
      15.0, 0.350;...
       ];

figsize = 400;
fig1=figure(1);clf;
set(gcf, 'Position',  [1, 640, figsize*0.6, figsize*0.3]) 

plot(Cstar(:,1),Cstar(:,2),'-o','LineWidth',2);
xlabel('$\tau$','FontSize',14,'Interpreter','latex');
ylabel('$C$','FontSize',14,'Interpreter','latex');
xlim([0,15])
ylim([0,0.4])
set(gca,'TickLabelInterpreter','latex');  
%title('angle','Interpreter','latex');
%axis tight
grid on    
box on
set(gca,'FontSize',20)

dpi = '-r400';
saveas(fig1,'fig_3sync','epsc')

end
