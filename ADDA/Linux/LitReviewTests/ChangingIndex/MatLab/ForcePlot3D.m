fileDirectory5=['..' filesep '2_0i' filesep 'RadForce-Y.csv'];
File5=csvread(fileDirectory5,1,0);
fileDirectory6=['..' filesep '2_0.5i' filesep 'RadForce-Y.csv'];
File6=csvread(fileDirectory6,1,0);
fileDirectory7=['..' filesep '2_1i' filesep 'RadForce-Y.csv'];
File7=csvread(fileDirectory7,1,0);
figure
hold on;
quiver3(File5(:,1),File5(:,2),File5(:,3),File5(:,5),File5(:,6),File5(:,7),'LineWidth',2);
quiver3(File6(:,1),File6(:,2),File6(:,3),File6(:,5),File6(:,6),File6(:,7),'LineWidth',2);
quiver3(File7(:,1),File7(:,2),File7(:,3),File7(:,5),File7(:,6),File7(:,7),'LineWidth',2);
set(gca,'TickLabelInterpreter','latex')
xlabel('x \(\left(\mu m\right)\)','Interpreter','latex')
ylabel('y \(\left(\mu m\right)\)','Interpreter','latex')
zlabel('z \(\left(\mu m\right)\)','Interpreter','latex')
mylegend=legend({'$m=2 \left(10^{-13} N\right)$','$m=2+0.5i \left(10^{-13} N\right)$','$m=2+i \left(10^{-13} N\right)$'});
set(mylegend,'Interpreter','latex')
hold off;
view(45,45)
clear