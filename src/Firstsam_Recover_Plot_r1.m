%------------------------------------------------------------------------------------%
% Description: Recover info and plot, case M = N*log(N), rad = 1, First sampling. 
% Programer: Juan Manuel Cardenas
% Date: July 16 - 2019   / Last update: July 16 - 2019
%------------------------------------------------------------------------------------%
%--- Set up ---%
clear all;clc;close all
Trials = 10;
d = [2 3 5 10]; %15 %
Nvalues = zeros(length(d),Trials);
Mvalues = zeros(length(d),Trials);


%------------------------------------------------------------------------------------%
%----------------------         d = 2  to d = 10              -----------------------%                         
%------------------------------------------------------------------------------------%

load('Firstsam_d2to10_NlogN_r1.mat')

Err_d2 = zeros(Trials,length(N_VALUES));
Err_d3 = zeros(Trials,length(N_VALUES));
Err_d5 = zeros(Trials,length(N_VALUES));
Err_d10 = zeros(Trials,length(N_VALUES));

for i = 1:length(N_VALUES)
    
    %--- save Errors ---%
    
    Err_d2(:,i) = ERROR_BOX(:,i,1);
    Err_d3(:,i) = ERROR_BOX(:,i,2);  
    Err_d5(:,i) = ERROR_BOX(:,i,3);  
    Err_d10(:,i) = ERROR_BOX(:,i,4);      
end

%--- save N values ---%

Nvalues = N_VALUES;
Mvalues = M_VALUES;

%------------------------------------------------------------------------------------%
%----------------------                d = 15                 -----------------------%                           
%------------------------------------------------------------------------------------%


load('Firstsam_d15_NlogN_r1.mat')

%--- save results ---%

Err_d15 = zeros(Trials,length(N_values));

for i = 1:length(N_values)
    Err_d15(:,i) = Error_box(:,i);                    
end

%--- save N values ---%

Nvalues15 = N_values;
Mvalues15 = M_values;


%------------------------------------------------------------------------------------%
%----------------------              Median Plot              -----------------------%                        
%------------------------------------------------------------------------------------%

%--- Compute Median  ---%

d_total = [2 3 5 10 15];
Error_median = [];

for i = 1:length(Nvalues(:,1))
    Error_median(i,1) = median(Err_d2(:,i));
    Error_median(i,2) = median(Err_d3(:,i));
    Error_median(i,3) = median(Err_d5(:,i)); 
    Error_median(i,4) = median(Err_d10(:,i));
end


for j = 1:length(Nvalues15)
    Error_median15(j) = median(Err_d15(:,j));
end

%--- Plot ---%

d_values = [2 3 5 10];
ms = 13;
lw = 1.5;

default_color = [0    0.4470    0.7410 ;
    0.8500    0.3250    0.0980 ;
    0.9290    0.6940    0.1250 ;
    0.4940    0.1840    0.5560 ;
    0.4660    0.6740    0.1880 ;
    0.3010    0.7450    0.9330 ;
    0.6350    0.0780    0.1840 ] ;

markers= {'-*','-o','-s','-^','-v','-+','-x','--*','--o','--s','--^','--v'};

fig = figure(1);
for l = 1:length(d_values)
    loglog(Mvalues(:,l),Error_median(:,l),markers{l},'markersize',ms,'MarkerFaceColor',default_color(l,:),'MarkerEdgeColor',default_color(l,:),'LineWidth',lw);
    hold on
end
%Add d=15

loglog(Nvalues15,Error_median15,markers{5},'markersize',ms,'MarkerFaceColor',default_color(5,:),'MarkerEdgeColor',default_color(5,:),'LineWidth',lw);
hold off 

h=legend('$d = 2$','$d = 3$','$d = 5$','$d = 10$','$d = 15$','location','southwest');
set(h,'Interpreter','latex');
X = xlabel('M values');
Y = ylabel('Error');
T = title('First strategy');

set(h,'Interpreter','latex');
set(X,'Interpreter','latex');
set(Y,'Interpreter','latex');
set(T,'Interpreter','latex');

ax = gca;
ax.YMinorTick = 'off';
ax.YMinorGrid = 'off';
ax.FontSize = 15;
ax.LineWidth = 2;
ax.YGrid = 'on';
ax.XGrid = 'on';
ax.XMinorTick = 'on';
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'off';
ax.YMinorTick = 'off';

ax.XTick = [1e0 1e1 1e2 1e3 1e4 1e5];
xlim([1e0 1e4]);

ax.YTick = [1e-14 1e-12 1e-10 1e-8 1e-6 1e-4 1e-2  1e0 1e2 1e4];
ylim([1e-14,1e2]);
