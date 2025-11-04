clear;
%% Initial conditions
parameters_spin;

%% creating output directory to save plots or data
name=strcat('/data/Omega',num2str(Omega),'N',num2str(N),'Ntr',num2str(Ntrc),'dt',num2str(dt),'gamma',num2str(gamma_d),'is_simplified',num2str(is_simplified),'is_noise',num2str(is_noise));
foldername=strcat('decay/',name);
[status, msg, msgID] = mkdir(foldername);


%Sampling of initial conditions
% Here, you can do this within the loop directly, which saves memery,
% or sampling all initial conditions in once
randomNumbersx = (2*randi([0, 1], [N*Ntrc,1])-1);
randomNumbersy = (2*randi([0, 1], [N*Ntrc,1])-1);
randomNumbersz(1:N*Ntrc,1)=-1;


% initialization of memory for output arrays
% Here, I save all trajectories and then do sums afterwards; 
% but the same result can be obtained by sumation of SSm ets directly
% within the loop over Ntrc (depends only on how much memory you have)

SSm(1:Mev,1:N,1:Ntrc)=0;
SSz(1:Mev,1:N,1:Ntrc)=0;
SSm2(1:Mev,1:N,1:Ntrc)=0;
SSz2(1:Mev,1:N,1:Ntrc)=0;
SSpm(1:Mev,1:N,1:Ntrc)=0;
%%

J=vanderWaalscouplings(J0,Alpha,N); % Coupling constnats
h=sum(J,2)/4; % effective magnetic field

for i=1:Ntrc

    % initial conditions: all spins pointed down
    sx0=randomNumbersx((i-1)*N+1:i*N,1);
    sy0=randomNumbersy((i-1)*N+1:i*N,1);
    sz0=randomNumbersz((i-1)*N+1:i*N,1);

% Rotation of the initial state if any

sxp=(sx0*cos(phi)-sy0*sin(phi));
syp=(sx0*sin(phi)+sy0*cos(phi));

sx=sz0*sin(theta)+sxp*cos(theta);
sy=syp;
sz=sz0*cos(theta)-sxp*sin(theta);
sm=(sx-1i*sy)/2;


% solution for a single trajectory
[t,Sm,Sz]=solver_spins_Rydbergs_real_loop(sm',sz',J,h);

%saving data for each trajectory

     SSm(1:Mev,1:N,i)=Sm';
     SSz(1:Mev,1:N,i)=Sz';
     SSm2(1:Mev,1:N,i)=(Sm').^2;
     SSpm(1:Mev,1:N,i)=abs(Sm').^2;
     SSz2(1:Mev,1:N,i)=(Sz').^2;
% i
end

%% Evaluation of averages
% average over different trajectories
sm_av=mean( SSm,3);
sz_av=mean( SSz,3);
sm2_av=mean(SSm2,3);
sz2_av=mean(SSz2,3);
spm_av=mean(SSpm,3);

%% Large spins
clear Sm Sz SSz SSm SSm2 SSz2 SSpm
% collective spin 
Sm=sum(sm_av,2);
Sz=sum(sz_av,2);
SPM=sum(spm_av,2);
SZ2=sum(sz2_av,2);
SMM=sum(sm2_av,2);

Sx=2*real(Sm);
Sy=2*imag(Sm);

sxsx=(real(SMM)+real(SPM))/2;
sysy=(real(SPM)-real(SMM))/2;

clear SSM SPM SMM spm_av sz2_av sm_av sm2_av sz_av

%%
%% plotting results


fig=figure; 
set(0,'DefaultTextFontSize',16,'DefaultTextFontName','Times New Roman'); 
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontName','Times New Roman');
hold on;
% plot(t,real(mean(Sm,1)),t,-imag(mean(Sm,1)),t,real(mean(Sz,1)),'LineWidth',2);
% % plot(t,sqrt(abs(mean(Sm,1)).^2+abs(mean(Sz,1)).^2),'LineWidth',2,'Color','k');
% plot(t,sqrt(mean(abs((Sm)).^2+abs((Sz)).^2,1)),'LineWidth',2,'Color','k');

plot(t,Sx,t,Sy,t,real((Sz)),'LineWidth',2);

% plot(t,sqrt(abs(mean(Sm,1)).^2+abs(mean(Sz,1)).^2),'LineWidth',2,'Color','k');
% plot(t,sqrt((abs((Sm)).^2+abs((Sz)).^2)),'LineWidth',2,'Color','k');

plot(t,sqrt((abs((Sx)).^2+abs((Sy)).^2+abs((Sz)).^2)),'LineWidth',2,'Color','k');



legend('$s^x$','$s^y$','$s^z$','interpreter','latex');
xlabel('$t$','Interpreter','latex');
ylabel('$\langle s\rangle$','Interpreter','latex');
box on;
saveas(fig, fullfile(foldername, 'individual spins'), 'epsc');
saveas(fig, fullfile(foldername, 'individual spins'), 'jpg');

%%
% fig=figure; 
% set(0,'DefaultTextFontSize',16,'DefaultTextFontName','Times New Roman'); 
% set(0,'DefaultAxesFontSize',16,'DefaultAxesFontName','Times New Roman');
% hold on;
% Nex=((Sz)+N/2)/N;%% number of excitations 
% plot(t(2:end),abs(diff(Nex)')./diff(t),'LineWidth',2);
% xlabel('$t$','Interpreter','latex');
% ylabel('$\langle ds^z\rangle/dt$','Interpreter','latex');
% box on;
% 
% saveas(fig, fullfile(foldername, 'intencity'), 'epsc');
% saveas(fig, fullfile(foldername, 'intensity'), 'jpg');




%% saving data if needed
% close all;
%   filename = strcat(foldername,'/data.mat');
%   save(filename)

