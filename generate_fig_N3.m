

load('generate_fig_N3/trajSPO.mat')
load('generate_fig_N3/trajUPO.mat')

N=2;

figsize = 400;

fig11=figure(11);clf;
phiU1=trajUPO(:,N+3);
    phiU2=trajUPO(:,N+4);
    theU1=trajUPO(:,4);
    theU2=trajUPO(:,5); 
phiS1=trajSPO(:,N+3);
    phiS2=trajSPO(:,N+4);
    theS1=trajSPO(:,4);
    theS2=trajSPO(:,5);  


set(gcf, 'Position',  [2*figsize, 120, figsize*0.8, figsize*0.8])
    hold on
    plot(theS1/pi,theS2/pi,'-c','LineWidth',3)  
    plot(theU1/pi,theU2/pi,'-r','LineWidth',3)  
    hold off
    xlabel('$\alpha_1/\pi$','FontSize',14,'Interpreter','latex');
    ylabel('$\alpha_2/\pi$','FontSize',14,'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex')  
    xlim([-0.6,0.6])
    ylim([-0.6,0.6])
    xticks(-0.5:0.5:0.5)
    yticks(-0.5:0.5:0.5)
    %axis tight
    grid on 
    box on
    set(gca,'FontSize',20)
    hold off
saveas(fig11,'fig_N3','epsc')