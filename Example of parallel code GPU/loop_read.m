function [] = loop_read(J_l,Gamma_l,V_l,mask_l,J_p,Gamma_p,V_p,mask_p)
parameters_spin;

if DO_PARALLEL
    J_l     = gather(J_l);
    Gamma_l = gather(Gamma_l);
    V_l     = gather(V_l);
    mask_l  = gather(mask_l);

    J_p     = gather(J_p);
    Gamma_p = gather(Gamma_p);
    V_p     = gather(V_p);
    mask_p  = gather(mask_p);
end

%reading results and evaluation averages over all trajectories
gamma_d=Gamma_l(1,1);
gamma_u=Gamma_p(1,1);
% Gamma_l

%output folder to save data and plots
name=strcat('/data/delay',num2str(tau0),'_Omega',num2str(Omega),'omega0',num2str(omega0),'eps',num2str(Eps),'dim',num2str(dim),'d',num2str(d),'kappa',num2str(kappa),'J',num2str(J0),'alpha',num2str(Alpha),'N',num2str(N),'Ntr',num2str(Ntrc),'dt',num2str(dt),'gamma_d',num2str(gamma_d),'gamma_u',num2str(gamma_u));
foldername=strcat('collective/',name);
[status, msg, msgID] = mkdir(foldername);
% 

% some matrices I need to precalculate to evaluate some correlation
% functions below
G     = Gamma_l;                 % N×N
Gt    = G.';                     % N×N
B     = G .* G.';                % N×N  (static; move outside the trace loop)


% reserve memory 
SSm(1:Mev,1:N)=0;
SSz(1:Mev,1:N)=0;
SSM2(1:Mev,1:N)=0;
SSz2(1:Mev,1:N)=0;
SSPM(1:Mev,1:N)=0;
obs1(1:Mev,1)=0;
obs2(1:Mev,1)=0;
obs3(1:Mev,1)=0;
obs4(1:Mev,1)=0;
obs5(1:Mev,1)=0;
rate(1:Mev,1)=0;
g2_prep(1:Mev,1)=0;
SSpmtau0(1:Mev,1)=0;
deltaSp(1:Mev,1)=0;
SSpmtau02(1:Mev,1)=0;
deltaSp2(1:Mev,1)=0;
rt(1:Mev,1)=0;

%read out loop
for i=1:Ntrc

name_file = strcat(foldername_date, strcat('/data_', num2str(i), '.mat'));
load(name_file);

Sm=S_m;
Sz=S_z;
Smp=S_mp';
Szp=S_zp;
sm_0=sigma_m_0';

sm=Sm'; sz=Sz';


  
SSm=SSm+(sm);
SSz=SSz+(sz);
SSM2=SSM2+(sm.^2);
SSz2=SSz2+(sz).^2;
SSPM=SSPM+abs(sm).^2;
SM0=(sum(sm(1,1:N),2));

Spt=conj(sum(sm,2));

SSpmtau0(1:Mev,1) = SSpmtau0(1:Mev,1) + Spt .* SM0; %main contribution to two-time corr fucntion
eltaSp(1:Mev,1)=deltaSp+sum(conj(Smp(:,1:N)-sm(:,1:N)),2);% quantum correction to two-time corr fucntion (responce part)

%some checks that can be deleted

SSpmtau02(1:Mev,1) =SSpmtau02+ conj(sm(:,1)) .* (sm_0(1)); 
deltaSp2(1:Mev,1)=deltaSp+conj(Smp(:,1)-sm(:,1));

%SSpmtau0(1:Mev,1) =SSpmtau0+ conj(sm(:,1)) .* (sm_0); 



    % obs1 = real((sm * Gamma_l * Sm));  % Compute the quadratic form for all rows
obs1(1:Mev,1)=0;
obs2(1:Mev,1)=0;
obs3(1:Mev,1)=0;
obs4(1:Mev,1)=0;
obs5(1:Mev,1)=0;

smc   = conj(sm);                % Mev×N


    SG    = sm * Gt;                 % Mev×N where row r: sum_j sm(r,j)*G(i,j)
obs1  = sum(SG .* smc, 2);       % Mev×1

% obs5: z*B*z'
ZB    = sz * B.';                % Mev×N
obs5  = sum(ZB .* sz, 2);        % Mev×1

% obs3: dot(z, (G*conj(s)') .* (G*s'))
A     = smc * Gt;                % Mev×N, A(r,i) = sum_k conj(sm(r,k))*G(i,k)
Bv    = sm  * Gt;                % Mev×N, Bv(r,i) = sum_j sm(r,j)*G(i,j)
obs3  = sum((A .* Bv) .* sz, 2); % Mev×1
    % 
    % rate=rate+obs1+0.5*Gamma_l(1,1)*sum(sz(1:Mev,1:N),2);
    % 
    % rt=rt+abs(sm(:, 1)).^2+0.5*sz(:,1);
    % 
    % obs4=obs1.^2;
    % 
    % obs2=Gamma_l(1,1)*sum(sz(1:Mev,1:N),2).*obs1;
    % 
    % 
    % 
    % g2_prep=g2_prep+obs4+obs2+obs3-1./3*Gamma_l(1,1).^2*sum(sz,2)+0.25*obs5+0.25*Gamma_l(1,1).^2*(sum(sz,2)).^2-4./3*Gamma_l(1,1)*obs1;

 sum_sz = sum(sz, 2);                    % Mev×1, reuse repeatedly
rate   = rate + obs1 + 0.5*G(1,1)*sum_sz;
rt     = rt   + abs(sm(:,1)).^2 + 0.5*sz(:,1);
obs4   = obs1.^2;
obs2   = G(1,1) * sum_sz .* obs1;

g2_prep = g2_prep + obs4 + obs2 + obs3 ...
        - (1/3)*G(1,1)^2*sum_sz + 0.25*obs5 ...
        + 0.25*G(1,1)^2*(sum_sz).^2 - (4/3)*G(1,1)*obs1;


if mod(i,200) == 0   % true every 100th step
        fprintf('i = %d\n', i);
end


clear S_m S_z S_mp S_zp sigma_m_0

end

sm_av=SSm/Ntrc;
sz_av=SSz/Ntrc;
sm2_av=SSM2/Ntrc;
sz2_av=SSz2/Ntrc;
spm_av=SSPM/Ntrc;
%  nex_av=mean(Ex_num,3);
%% Large spins
clear Sm Sz SSm SSz SSm2 SSz2 SSpm SPM SZ2
%  Nex=sum(nex_av,2);
%  figure; plot(t,Nex);
Sm=sum(sm_av,2);
Sz=sum(sz_av,2);
SPM=sum(spm_av,2);
SZ2=sum(sz2_av,2);
SMM=sum(sm2_av,2);

g2_av=g2_prep/Ntrc;
Rate=rate/Ntrc;
Rt=rt/Ntrc;
sxsx=(real(SMM)+real(SPM))/2;
sysy=(real(SPM)-real(SMM))/2;
Stau=SSpmtau0/Ntrc;%+0.5*1i*mean(deltaSp,2)/Eps/N;
dStau=0.5*1i*(deltaSp/Ntrc)/Eps;

Stau2=SSpmtau02/Ntrc;%+0.5*1i*mean(deltaSp,2)/Eps/N;
dStau2=0.5*1i*(deltaSp2/Ntrc)/Eps;

figure; plot(t,abs(Stau),t,(dStau));
title('local')
box on;
omega_max=2*pi/(t(2)-t(1));
omega_grid=linspace(-omega_max/2,omega_max/2,Mev);

dSomega=2*real(fftshift(fft(dStau)));
Somega=2*real(fftshift(fft(Stau)));

dSomega2=2*real(fftshift(fft(dStau2)));
Somega2=2*real(fftshift(fft(Stau2)));



fig=figure;
set(0,'DefaultTextFontSize',16,'DefaultTextFontName','Times New Roman'); 
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontName','Times New Roman');
hold on;
plot(omega_grid,abs(Somega),'LineWidth',2);

plot(omega_grid,abs(dSomega),'LineWidth',2);
% xlim([-20,20]);
xlabel('$\omega$','interpreter','latex')
ylabel('$S(\omega)_{corr}$','interpreter','latex')
box on;
saveas(fig, fullfile(foldername, 'Somega'), 'epsc');
saveas(fig, fullfile(foldername, 'Somega'), 'jpg');


%%
fig=figure;
set(0,'DefaultTextFontSize',16,'DefaultTextFontName','Times New Roman'); 
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontName','Times New Roman');
hold on;
plot(omega_grid,abs(Somega+dSomega),'LineWidth',2);
% xlim([-20,20]);
xlabel('$\omega$','interpreter','latex')
ylabel('$S(\omega)$','interpreter','latex')
box on;
saveas(fig, fullfile(foldername, 'Somega_corr_full'), 'epsc');
saveas(fig, fullfile(foldername, 'Somega_corr_full'), 'jpg');


%%

fig=figure; 
plot(t,g2_av./Rate.^2,'linewidth',2);
% ylim([0,3])
xlabel('$t$','interpreter','latex')
ylabel('$g^{(2)}(0)$','interpreter','latex')
box on;
saveas(fig, fullfile(foldername, 'g2'), 'epsc');
 saveas(fig, fullfile(foldername, 'g2'), 'jpg');

 fig=figure; 
plot(t,g2_av/N^2,'linewidth',2);
% ylim([0,3])
xlabel('$t$','interpreter','latex')
ylabel('$G^{(2)}(0)/N^2$','interpreter','latex')
box on;
saveas(fig, fullfile(foldername, 'GG2'), 'epsc');
 saveas(fig, fullfile(foldername, 'GG2'), 'jpg');


fig=figure; 
hold on;
plot(t,Rate/N,'linewidth',2);%ylim([0,10])
plot(t,Rt/N,'linewidth',2);%ylim([0,10])

xlabel('$t$','interpreter','latex')
ylabel('$\langle S^+ S^-\rangle/N$','interpreter','latex')
box on;
saveas(fig, fullfile(foldername, 'g1'), 'epsc');
saveas(fig, fullfile(foldername, 'g1'), 'jpg');


%% Evaluation of averages

%%
% t=linspace(t_in,t_fin,M);
%% 

fig=figure; 
set(0,'DefaultTextFontSize',16,'DefaultTextFontName','Times New Roman'); 
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontName','Times New Roman');
hold on;
plot(t,2*real((Sm)),t,2*imag((Sm)),t,real((Sz)),'LineWidth',2);

% plot(t,sqrt(abs(mean(Sm,1)).^2+abs(mean(Sz,1)).^2),'LineWidth',2,'Color','k');
% plot(t,sqrt((abs((Sm)).^2+abs((Sz)).^2)),'LineWidth',2,'Color','k');
Sx=2*real(Sm);
Sy=2*imag(Sm);
plot(t,sqrt((abs((Sx)).^2+abs((Sy)).^2+abs((Sz)).^2)),'LineWidth',2,'Color','k');


legend('$s^x$','$s^y$','$s^z$','interpreter','latex');
xlabel('$t$','Interpreter','latex');
ylabel('$\langle s\rangle$','Interpreter','latex');
box on;
saveas(fig, fullfile(foldername, 'individual spins '), 'epsc');
saveas(fig, fullfile(foldername, 'individual spins'), 'jpg');
% figure;  plot(t,mean(rate,2))




%%

% sum_val = 1+(sum(sum(Gamma_l .* Gamma_l.'))-2*N*gamma_d^2)/N^2/gamma_d^2




%%

% sum_val = 1+(sum(sum(Gamma_l .* Gamma_l.'))-2*N*gamma_d^2)/N^2/gamma_d^2

close all;


name_file=strcat(foldername, '/data.mat');
save(name_file,"t","omega_grid",'dSomega','Somega','g2_av','Rate','N','Stau','dStau','Sz','Sx','Sy','SPM',"SZ2",'SMM')


end