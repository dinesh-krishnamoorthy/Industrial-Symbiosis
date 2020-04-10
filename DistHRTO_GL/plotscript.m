load('Data/data.mat')

w_pg_tot = sum(data.data1.y(10,:)') + sum(data.data2.y(10:12,:)');
w_gl_tot = sum(data.data1.u) + sum(data.data2.u);
w_po_tot = sum(data.data1.y(7:9,:)') + sum(data.data2.y(7:9,:)');

w_gl_opt = sum(data.data1.u_opt) + sum(data.data2.u_opt);
w_pg_opt = sum(data.data1.w_pg_sp) + sum(data.data2.w_pg_sp);

n = 1:nIter;
t = n*TimeRTO/3600;

figure(1)
clf
hold all
plot(data.data1.d(1:end-1,:)','linewidth',1.5)
plot(data.data2.d(1:end-1,:)','linewidth',1.5)
title('Disturbance in GOR','Interpreter','Latex')
legend('well 1-A','well 1-B','well 1-C',...
    'well 2-A','well 2-B','well 2-C','Interpreter','Latex')
xlabel('Time [h]','Interpreter','Latex')
box on
grid on

figure(2)
clf
subplot(321)
hold all
plot(8.5.*ones(size(w_gl_tot)),'k:','linewidth',1.5)
plot(w_gl_tot,'k','linewidth',1.5)
ylabel('Total Gas Lift [kg/s]','Interpreter','Latex')
box on
grid on

subplot(322)
hold all
plot(sum(data.data1.u),'linewidth',1.5)
plot(sum(data.data2.u),'linewidth',1.5)
legend('Company 1','Company 2','Interpreter','Latex')
ylabel('Gas Lift [kg/s]','Interpreter','Latex')
box on
grid on

subplot(323)
hold all
plot(27.5.*ones(size(w_gl_tot)),'k:','linewidth',1.5)
plot(sum(data.data1.w_pg) + sum(data.data2.w_pg),'k','linewidth',1.5)
ylabel('Total produced gas [kg/s]','Interpreter','Latex')

subplot(324)
hold all
plot(sum(data.data1.w_pg),'linewidth',1.5)
plot(sum(data.data2.w_pg),'linewidth',1.5)
legend('Company 1','Company 2','Interpreter','Latex')
ylabel('Produced gas [kg/s]','Interpreter','Latex')
box on
grid on

subplot(325)
hold all
plot(sum(data.data1.w_po) + sum(data.data2.w_po),'k','linewidth',1.5)
ylabel('Total produced oil [kg/s]','Interpreter','Latex')
xlabel('Time [h]','Interpreter','Latex')
box on
grid on

subplot(326)
hold all
plot(sum(data.data1.w_po),'linewidth',1.5)
plot(sum(data.data2.w_po),'linewidth',1.5)
legend('Company 1','Company 2','Interpreter','Latex')
ylabel('Produced oil [kg/s]','Interpreter','Latex')
xlabel('Time [h]','Interpreter','Latex')
box on
grid on

figure(3)
clf
subplot(321)
hold all
stairs(data.data1.w_pg_sp(1,:),'k','linewidth',1.5)
plot(data.data1.w_pg(1,:),'b','linewidth',1.5)
title('well 1-A','Interpreter','Latex')
ylabel('Produced gas [kg/s]','Interpreter','Latex')
box on
grid on

subplot(322)
hold all
plot(data.data1.u(1,:),'b','linewidth',1.5)
title('well 1-A','Interpreter','Latex')
ylabel('Gas lift rate [kg/s]','Interpreter','Latex')
box on
grid on

subplot(323)
hold all
stairs(data.data1.w_pg_sp(2,:),'k','linewidth',1.5)
plot(data.data1.w_pg(2,:),'b','linewidth',1.5)
title('well 1-B','Interpreter','Latex')
ylabel('Produced gas [kg/s]','Interpreter','Latex')
box on
grid on

subplot(324)
hold all
plot(data.data1.u(2,:),'b','linewidth',1.5)
title('well 1-B','Interpreter','Latex')
ylabel('Gas lift rate [kg/s]','Interpreter','Latex')
box on
grid on

subplot(325)
hold all
stairs(data.data1.w_pg_sp(3,:),'k','linewidth',1.5)
plot(data.data1.w_pg(3,:),'b','linewidth',1.5)
title('well 1-C','Interpreter','Latex')
ylabel('Produced gas [kg/s]','Interpreter','Latex')
xlabel('Time [h]','Interpreter','Latex')
box on
grid on

subplot(326)
hold all
plot(data.data1.u(3,:),'b','linewidth',1.5)
title('well 1-C','Interpreter','Latex')
ylabel('Gas lift rate [kg/s]','Interpreter','Latex')
xlabel('Time [h]','Interpreter','Latex')
box on
grid on


figure(4)
clf
subplot(321)
hold all
stairs(data.data2.w_pg_sp(1,:),'k','linewidth',1.5)
plot(data.data2.w_pg(1,:),'b','linewidth',1.5)
title('well 2-A','Interpreter','Latex')
ylabel('Produced gas [kg/s]','Interpreter','Latex')
box on
grid on

subplot(322)
hold all
plot(data.data2.u(1,:),'b','linewidth',1.5)
title('well 2-A','Interpreter','Latex')
ylabel('Gas lift rate [kg/s]','Interpreter','Latex')
box on
grid on

subplot(323)
hold all
stairs(data.data2.w_pg_sp(2,:),'k','linewidth',1.5)
plot(data.data2.w_pg(2,:),'b','linewidth',1.5)
title('well 2-B','Interpreter','Latex')
ylabel('Produced gas [kg/s]','Interpreter','Latex')
box on
grid on

subplot(324)
hold all
plot(data.data2.u(2,:),'b','linewidth',1.5)
title('well 2-B','Interpreter','Latex')
ylabel('Gas lift rate [kg/s]','Interpreter','Latex')
box on
grid on

subplot(325)
hold all
stairs(data.data2.w_pg_sp(3,:),'k','linewidth',1.5)
plot(data.data2.w_pg(3,:),'b','linewidth',1.5)
title('well 2-C','Interpreter','Latex')
ylabel('Produced gas [kg/s]','Interpreter','Latex')
xlabel('Time [h]','Interpreter','Latex')
box on
grid on

subplot(326)
hold all
plot(data.data2.u(3,:),'b','linewidth',1.5)
title('well 2-C','Interpreter','Latex')
ylabel('Gas lift rate [kg/s]','Interpreter','Latex')
xlabel('Time [h]','Interpreter','Latex')
box on
grid on

