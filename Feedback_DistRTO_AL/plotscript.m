function plotscript(data)

if nargin<1
    disp('Plotting preloaded data.')
    load('Data/data.mat')
end

%% Public variables

figure(1)
set(gcf,'position',[60,275,1000,667])
subplot(221)
set(gca,'FontSize',14) 
hold all
plot(data.public.qGLMax,'k:','linewidth',2)
plot(data.public.w_gl1_opt+data.public.w_gl2_opt,'k','linewidth',2)
ylabel('Total Gas Lift [kg/s]','Interpreter','Latex')
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'latex';

subplot(222)
set(gca,'FontSize',14)
hold all
plot(data.public.w_gl1_opt,'b','linewidth',2)
plot(data.public.w_gl2_opt,'r','linewidth',2)
legend('Company 1','Company 2','Interpreter','Latex')
ylabel('Gas Lift [kg/s]','Interpreter','Latex')
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'latex';

subplot(223)
set(gca,'FontSize',14)
hold all
% plot(data.public.QgMax,'k:','linewidth',2)
plot(data.public.lambda_wgl,'k','linewidth',2)
ylabel('Shadow price $\lambda$','Interpreter','Latex')
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'latex';


subplot(224)
set(gca,'FontSize',14)
hold all
plot(data.public.residual,'k','linewidth',2)
ylabel('Primal Residual','Interpreter','Latex')
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'latex';

%% private variables company 1

figure(2)
set(gcf,'position',[60,275,1000,667])
subplot(321)
set(gca,'FontSize',14) 
hold all
plot(zeros(size(data.data1.AL(1,:))),'k:','linewidth',2)
plot(data.data1.AL(1,:),'b','linewidth',2)
title('well 1-A','Interpreter','Latex')
ylabel('CV: $\nabla \mathcal{L}$ [kg/s]','Interpreter','Latex')
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'latex';

subplot(322)
set(gca,'FontSize',14) 
hold all
plot(data.data1.u(1,:),'b','linewidth',2)
title('well 1-A','Interpreter','Latex')
ylabel('Gas lift rate [kg/s]','Interpreter','Latex')
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'latex';

subplot(323)
set(gca,'FontSize',14) 
hold all
plot(zeros(size(data.data1.AL(2,:))),'k:','linewidth',2)
plot(data.data1.AL(2,:),'b','linewidth',2)
title('well 1-B','Interpreter','Latex')
ylabel('CV: $\nabla \mathcal{L}$ [kg/s]','Interpreter','Latex')
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'latex';

subplot(324)
set(gca,'FontSize',14) 
hold all
plot(data.data1.u(2,:),'b','linewidth',2)
title('well 1-B','Interpreter','Latex')
ylabel('Gas lift rate [kg/s]','Interpreter','Latex')
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'latex';

subplot(325)
set(gca,'FontSize',14) 
hold all
plot(zeros(size(data.data1.AL(3,:))),'k:','linewidth',2)
plot(data.data1.AL(3,:),'b','linewidth',2)
title('well 1-C','Interpreter','Latex')
ylabel('CV: $\nabla \mathcal{L}$ [kg/s]','Interpreter','Latex')
xlabel('Time [h]','Interpreter','Latex')
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'latex';

subplot(326)
set(gca,'FontSize',14) 
hold all
plot(data.data1.u(3,:),'b','linewidth',2)
title('well 1-C','Interpreter','Latex')
ylabel('Gas lift rate [kg/s]','Interpreter','Latex')
xlabel('Time [h]','Interpreter','Latex')
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'latex';

%% Private variables company 2

figure(3)
set(gcf,'position',[60,275,1000,667])
subplot(321)
set(gca,'FontSize',14) 
hold all
plot(zeros(size(data.data2.AL(1,:))),'k:','linewidth',2)
plot(data.data2.AL(1,:),'r','linewidth',2)
title('well 2-A','Interpreter','Latex')
ylabel('CV: $\nabla \mathcal{L}$ [kg/s]','Interpreter','Latex')
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'latex';

subplot(322)
set(gca,'FontSize',14) 
hold all
plot(data.data2.u(1,:),'r','linewidth',2)
title('well 2-A','Interpreter','Latex')
ylabel('Gas lift rate [kg/s]','Interpreter','Latex')
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'latex';

subplot(323)
set(gca,'FontSize',14) 
hold all
plot(zeros(size(data.data2.AL(2,:))),'k:','linewidth',2)
plot(data.data2.AL(2,:),'r','linewidth',2)
title('well 2-B','Interpreter','Latex')
ylabel('CV: $\nabla \mathcal{L}$ [kg/s]','Interpreter','Latex')
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'latex';

subplot(324)
set(gca,'FontSize',14) 
hold all
plot(data.data2.u(2,:),'r','linewidth',2)
title('well 2-B','Interpreter','Latex')
ylabel('Gas lift rate [kg/s]','Interpreter','Latex')
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'latex';

subplot(325)
set(gca,'FontSize',14) 
hold all
plot(zeros(size(data.data2.AL(3,:))),'k:','linewidth',2)
plot(data.data2.AL(3,:),'r','linewidth',2)
title('well 2-C','Interpreter','Latex')
ylabel('CV: $\nabla \mathcal{L}$ [kg/s]','Interpreter','Latex')
xlabel('Time [h]','Interpreter','Latex')
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'latex';

subplot(326)
set(gca,'FontSize',14) 
hold all
plot(data.data2.u(3,:),'r','linewidth',2)
title('well 2-C','Interpreter','Latex')
ylabel('Gas lift rate [kg/s]','Interpreter','Latex')
xlabel('Time [h]','Interpreter','Latex')
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'latex';




%%
figure(4)
subplot(321)
hold all
plot(data.public.qGLMax,'k:','linewidth',2)
plot(data.public.w_gl1_opt+data.public.w_gl2_opt,'k','linewidth',2)
ylabel('Total Gas Lift [kg/s]','Interpreter','Latex')
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'latex';

subplot(322)
hold all
plot(data.public.w_gl1_opt,'linewidth',2)
plot(data.public.w_gl2_opt,'linewidth',2)
legend('Company 1','Company 2','Interpreter','Latex')
ylabel('Gas Lift [kg/s]','Interpreter','Latex')
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'latex';

subplot(323)
hold all
% plot(data.public.QgMax,'k:','linewidth',2)
plot(data.public.w_pg1_opt + data.public.w_pg2_opt,'k','linewidth',2)
ylabel('Total produced gas [kg/s]','Interpreter','Latex')
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'latex';
ylim([15,30])

subplot(324)
hold all
plot(data.public.w_pg1_opt,'linewidth',2)
plot(data.public.w_pg2_opt,'linewidth',2)
legend('Company 1','Company 2','Interpreter','Latex')
ylabel('Produced gas [kg/s]','Interpreter','Latex')
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'latex';

subplot(325)
hold all
plot(sum(data.data1.y(7:9,:)) + sum(data.data2.y(7:9,:)),'k','linewidth',2)
ylabel('Total produced oil [kg/s]','Interpreter','Latex')
xlabel('Time [h]','Interpreter','Latex')
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'latex';
ylim([160,170])

subplot(326)
hold all
plot(sum(data.data1.y(7:9,:)),'linewidth',2)
plot(sum(data.data2.y(7:9,:)),'linewidth',2)
legend('Company 1','Company 2','Interpreter','Latex')
ylabel('Produced oil [kg/s]','Interpreter','Latex')
xlabel('Time [h]','Interpreter','Latex')
box on
grid on
axs = gca;
axs.TickLabelInterpreter = 'latex';

