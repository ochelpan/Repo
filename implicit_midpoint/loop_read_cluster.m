function [] = loop_read(J_l,Gamma_l,V_l,mask_l,J_p,Gamma_p,V_p,mask_p,Alpha)
%reading results and evaluation averages over all trajectories
parameters_spin;
gamma_d=Gamma_l(1,1);
gamma_u=Gamma_p(1,1);
% Gamma_l

%output folder to save data and plots
name=strcat('/delay',num2str(tau0),'_Omega',num2str(Omega),'omega0',num2str(omega0),'eps',num2str(Eps),'dim',num2str(dim),'d',num2str(d),'kappa',num2str(kappa),'J',num2str(J0),'alpha',num2str(Alpha),'N',num2str(N),'Ntr',num2str(Ntrc),'dt',num2str(dt),'gamma_d',num2str(gamma_d),'gamma_u',num2str(gamma_u));
foldername=strcat('far field/',name);
[status, msg, msgID] = mkdir(foldername);
% 

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
% name=strcat('/data/delay',num2str(tau0),'_Omega',num2str(Omega),'omega0',num2str(omega0),'eps',num2str(Eps),'dim',num2str(dim),'d',num2str(d),'kappa',num2str(kappa),'J',num2str(J0),'alpha',num2str(Alpha),'N',num2str(N),'Ntr',num2str(Ntrc),'dt',num2str(dt),'gamma_d',num2str(gamma_d),'gamma_u',num2str(gamma_u));
% foldername=strcat('collective/',name);
% [status, msg, msgID] = mkdir(foldername);
% 

% some matrices I need to precalculate to evaluate some correlation
% functions below
G     = Gamma_l;                 % N×N
Gt    = G.';                     % N×N
B     = G .* G.';                % N×N  (static; move outside the trace loop)


% reserve memory 
SSm(1:Mev,1:N)=0;
SSz_k(1:Mev,1:N)=0;
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


  
% SSm=SSm+(sm);
% SSz=SSz+(sz);

x=d:d:N*d; %% 
vec=exp(-1i*2*pi*x);
    
SSm=SSm+(sm);
SSz=SSz+(sz);
sz_fft_col = fftshift(fft(ifftshift(sz,2), [], 2), 2);
SSz_k=SSz_k+sz_fft_col;



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

mom_range=linspace(-N/2,N/2-1,N);
mom_range=2*pi/N*mom_range;
sz_k=real(SSz_k)/N/Ntrc;


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

omega_max=2*pi/(t(2)-t(1));
omega_grid=linspace(-omega_max/2,omega_max/2,Mev);

dSomega=2*real(fftshift(fft(dStau)));
Somega=2*real(fftshift(fft(Stau)));

dSomega2=2*real(fftshift(fft(dStau2)));
Somega2=2*real(fftshift(fft(Stau2)));





%% Evaluation of averages

%%
% t=linspace(t_in,t_fin,M);
%% 


% plot(t,sqrt(abs(mean(Sm,1)).^2+abs(mean(Sz,1)).^2),'LineWidth',2,'Color','k');
% plot(t,sqrt((abs((Sm)).^2+abs((Sz)).^2)),'LineWidth',2,'Color','k');
Sx=2*real(Sm);
Sy=2*imag(Sm);
% plot(t,sqrt((abs((Sx)).^2+abs((Sy)).^2+abs((Sz)).^2)),'LineWidth',2,'Color','k');


name_file=strcat(foldername, '/data.mat');
% save(name_file,"t","omega_grid",'dSomega','Somega','g2_av','Rate','N','Stau','dStau','Sz','Sx','Sy','SPM',"SZ2",'SMM')
save(name_file,"t","sz_k","mom_range","omega_grid",'dSomega','Somega','g2_av','Rate','N','Stau','dStau','Sz','Sx','Sy','SPM',"SZ2",'SMM')


end