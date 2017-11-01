fileDirectory5=['..' filesep 'Lambda_5' filesep 'RadForce-Y.csv'];
File5=csvread(fileDirectory5,1,0);
fileDirectory6=['..' filesep 'Lambda_6' filesep 'RadForce-Y.csv'];
File6=csvread(fileDirectory6,1,0);
fileDirectory7=['..' filesep 'Lambda_7' filesep 'RadForce-Y.csv'];
File7=csvread(fileDirectory7,1,0);
figure
hold on;
quiver3(File5(:,1),File5(:,2),File6(:,3),File5(:,5),File5(:,6),File6(:,7),'LineWidth',2);
quiver3(File6(:,1),File6(:,2),File6(:,3),File6(:,5),File6(:,6),File6(:,7),'LineWidth',2);
quiver3(File7(:,1),File7(:,2),File6(:,3),File7(:,5),File7(:,6),File6(:,7),'LineWidth',2);
set(gca,'TickLabelInterpreter','latex')
xlabel('x \(\left(\mu m\right)\)','Interpreter','latex')
ylabel('y \(\left(\mu m\right)\)','Interpreter','latex')
zlabel('z \(\left(\mu m\right)\)','Interpreter','latex')
mylegend=legend({'$\lambda = 5\mu m \left(\mu N\right)$','$\lambda = 6\mu m \left(\mu N\right)$','$\lambda = 7\mu m \left(\mu N\right)$'});
set(mylegend,'Interpreter','latex')
hold off;
view(45,45)
clear