function out=generate_fig_sims_draw2(tps, traj, params)


% Motion of a swimmer driven by CPG

addpath(genpath('./helpers'))

% Output options.
Itraj=1; % 1: to visualise trajecotry
IMovie=0; % 1: to generate a movie

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

%% post-processing

Xc=[];Yc=[];
X1=[];Y1=[];
for i = 1:length(tps)
     [X,Y,TH]=coordinates_filament(traj(i,:),params);
     Xc=[Xc sum(X)/(N+1)];Yc=[Yc sum(Y)/(N+1)];
     X1=[X1 X(1)];
     Y1=[Y1 Y(1)];
end

figsize = 400;

%-------- draw trajectories

if Itraj==1
%{
fig7=figure(7);clf;
set(gcf, 'Position',  [1, 1, 3*figsize, 3*figsize])
for i=1:4:length(tps)%i=length(tps)-1000:10:length(tps)
    plot(Xc(1:i),Yc(1:i),'k','LineWidth',2)
    hold on
    plot(X1(1:i),Y1(1:i),'k')
    [X,Y,TH]=coordinates_filament(traj(i,:),params);
    plot(X,Y,'b','LineWidth',6)
    hold off 
    box on
    axis equal
    set(gca,'FontSize',20)
    set(gca,'TickLabelInterpreter','latex') 
    title(['$T$ =',num2str(tps(i))],'Interpreter','latex')
    drawnow
end
%}

fig8=figure(8);clf;
%set(gcf, 'Position',  [figsize, 640, figsize, figsize])
set(gcf, 'Position',  [1, 640, figsize*1.2, figsize*0.5])  
    hold on
    timeseq=1:12:length(tps);%floor(0.8*length(tps)):6:length(tps);
for i=timeseq
    [X,Y,TH]=coordinates_filament(traj(i,:),params);
    %plot(traj(1:i,1)+8*(i-1)/length(tps),traj(1:i,2),'LineWidth',1, 'Color','black')
    plot(X+12*(i-1)/length(tps),Y,'LineWidth',2) 
end

    colormap parula
    cm=colormap;
    ax=axis;
    icod=floor(linspace(1,256,length(timeseq)));
    colororder(cm(icod,:))

    for i=timeseq
    [X,Y,TH]=coordinates_filament(traj(i,:),params);
    plot(X(1)+12*(i-1)/length(tps),Y(1),'o','MarkerSize',3,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor',[.2 .2 .2])
    end
    %set(ax,'ColorOrder', cm(icod,:));
    axis equal
    hold off
    xlabel('time','FontSize',14,'Interpreter','latex');
    ylabel('$y$','FontSize',14,'Interpreter','latex');
    set(gca,'FontSize',20)
    set(gca,'TickLabelInterpreter','latex')
    xlim([-0.1 12.1])
    xticks([0:4:12])
    ylim([-1.1 1.1])
    %yticks([-0.1:0.1:0.2])
    grid on
    box on
    
    cb=colorbar;
    set(cb,'Ticks',[])
    %set(cb,'Position',[0.92 0.52 0.03 0.4])
    %set(cb, 'Location','northoutside')
    set(cb,'FontSize',20)
    set(cb,'LineWidth',1)
    set(cb,'Ticks',[0.5],...
       'TickLabels',{'time'},...
       'TickLabelInterpreter','latex')


%{
fig9=figure(9);clf;
set(gcf, 'Position',  [2*figsize, 640, figsize, figsize])
    phi1=traj(:,N+3);
    phi2=traj(:,N+4);
    phiD=mod(diff(traj(:,N+3:end),1,2)+pi,2*pi)-pi;
    phiav=mean(phiD,2);
    plot(tps,phiD(:,1)/pi,tps, phiav/pi,'LineWidth',3)  
    hold on
    for i=1:N-2
        plot(tps,phiD(:,i)/pi,'LineWidth',1) 
    end
    hold off
    xlabel('$t$','FontSize',14,'Interpreter','latex');
    ylabel('$\phi_i/\pi$','FontSize',14,'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex')    
    axis tight
    grid on 
    box on
    set(gca,'FontSize',20)
    hold off 


    phi1=traj(:,N+3);
    phi2=traj(:,N+4);
    the1=traj(:,4);
    the2=traj(:,5);    
fig10=figure(10);clf;
set(gcf, 'Position',  [2*figsize, 640, figsize, figsize])
    hold on
    plot(phi1/pi,phi2/pi,'LineWidth',2)  
    hold off
    xlabel('$\phi_1/\pi$','FontSize',14,'Interpreter','latex');
    ylabel('$\phi_2/\pi$','FontSize',14,'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex')    
    axis tight
    grid on  
    box on
    set(gca,'FontSize',20)
    hold off 

fig11=figure(11);clf;
set(gcf, 'Position',  [2*figsize, 640, figsize, figsize])
    hold on
    scatter3(the1/pi,the2/pi, mod(phi2/pi,2),20,tps,'filled')  
    hold off
    xlabel('$\alpha_1/\pi$','FontSize',14,'Interpreter','latex');
    ylabel('$\alpha_2/\pi$','FontSize',14,'Interpreter','latex');
    zlabel('$\phi_2/\pi$','FontSize',14,'Interpreter','latex');    
    set(gca,'TickLabelInterpreter','latex')    
    axis equal
    cb=colorbar;
    colormap jet;
    set(cb, "TickLabelInterpreter",'latex');
    grid on    
    box on
    set(gca,'FontSize',20)
    hold off  
%}

%fig6=figure(6);clf;
%set(gcf, 'Position',  [figsize, 640, figsize*2, figsize])
    Larc=cumsum(L)/sum(L);
    Tarc=tps;
    [Tmesh,Lmesh]=meshgrid(Tarc,Larc(1:end-1));
    Adat=traj(:,4:N+2);
    Phdat=mod(traj(:,N+3:end),2*pi);

    Larc1=linspace(0,1,N);
    [Tmesh,Lmesh]=meshgrid(Tarc,Larc1);
    Adat=cat(2,traj(:,4:N+2),traj(:,N+2));
    Phdat=cat(2,mod(traj(:,N+3:end),2*pi),mod(traj(:,end),2*pi));

fig5=figure(5);clf;
set(gcf, 'Position',  [figsize, 640, figsize*1.2, figsize*0.5])    
    sA=surf(Tmesh,Lmesh,Adat');
    sA.EdgeColor = 'none';
    %image(PhData','CDataMapping','scaled');
    cbA=colorbar;
    colormap jet
    xlabel('time','FontSize',14,'Interpreter','latex');
    ylabel('length','FontSize',14,'Interpreter','latex');
    ylim([0,1])
    title('angle','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');    
    set(cbA,'TickLabelInterpreter','latex'); 
    axis tight
    view(0,90)
    grid on    
    box on
    set(gca,'FontSize',20)

fig6=figure(6);clf;
set(gcf, 'Position',  [figsize, 300, figsize*1.2, figsize*0.5])    
    sPh=surf(Tmesh,Lmesh,Phdat');
    sPh.EdgeColor = 'none';
    cbPh=colorbar;
    colormap hsv
    xlabel('time','FontSize',14,'Interpreter','latex');
    ylabel('length','FontSize',14,'Interpreter','latex');
    title('phase','Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')    
    set(cbPh,'TickLabelInterpreter','latex'); 
    set(cbPh,'Limits',[0 2*pi]); 
    set(cbPh,'Ticks',[0 pi 2*pi]);
    set(cbPh,'TickLabels',{'0', '$\pi$', '$2\pi$'});       
    axis tight
    view(0,90)
    grid on 
    box on
    set(gca,'FontSize',20)

     

dpi = '-r400';
saveas(fig5,'fig_sims_swimH_angle','epsc')
saveas(fig6,'fig_sims_swimH_phase','epsc')
%saveas(fig7,'fig7','epsc') 
saveas(fig8,'fig_sims_swimH_shape','epsc') 
%saveas(fig9,'fig9','epsc') 


%saveas(fig10,'fig10','epsc') 
%saveas(fig11,'fig11','epsc') 

end



end
