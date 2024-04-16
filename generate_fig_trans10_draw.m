function out=generate_fig_trans10_draw(tps, traj, params)


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


fig8=figure(8);clf;
%set(gcf, 'Position',  [figsize, 640, figsize, figsize])
set(gcf, 'Position',  [1, 640, figsize*1.5, figsize*0.5])  
    hold on
    timeseq=1:70:length(tps);%floor(0.8*length(tps)):6:length(tps);
for i=timeseq
    [X,Y,TH]=coordinates_filament(traj(i,:),params);
    plot(X,Y,'LineWidth',2) 
end 

    colormap parula
    cm=colormap;
    ax=axis;
    icod=floor(linspace(1,256,length(timeseq)));
    colororder(cm(icod,:))

    plot(traj(:,1),traj(:,2),'-','LineWidth',1, 'Color',[.5 .5 .5])
    plot(traj(timeseq,1),traj(timeseq,2),'o','MarkerSize',3,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor',[.2 .2 .2])

    %set(ax,'ColorOrder', cm(icod,:));
    axis equal
    hold off
    xlabel('$x$','FontSize',14,'Interpreter','latex');
    ylabel('$y$','FontSize',14,'Interpreter','latex');
    set(gca,'FontSize',20)
    set(gca,'TickLabelInterpreter','latex')
    xlim([-10 1])
    xticks([-10:2:0])
    ylim([-1 1])
    yticks([-1:1:1])
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
set(gcf, 'Position',  [figsize, 640, figsize*0.7, figsize*0.4])    
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
set(gcf, 'Position',  [figsize, 300, figsize*0.7, figsize*0.4])    
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
saveas(fig5,'fig_trans10_crawl_angle','epsc')
saveas(fig6,'fig_trans10_crawl_phase','epsc')
%saveas(fig7,'fig7','epsc') 
saveas(fig8,'fig_trans10_crawl_shape','epsc') 
%saveas(fig9,'fig9','epsc') 


%saveas(fig10,'fig10','epsc') 
%saveas(fig11,'fig11','epsc') 

end



end
