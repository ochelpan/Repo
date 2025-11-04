function []= solver_spins_corr_emission_real_loop_par(J_l,Gamma_l,V_l,mask_l,J_p,Gamma_p,V_p,mask_p,Eps,indx_par)
% Here we evalute dynamics for a single trajectory
parameters_spin;
%% trying spherical mid_point method
%output folder to save intermediate results
% foldername_date=strcat('/Users/oksanachelpanova/Matlab_code/parallel_code/');
[status, msg, msgID] = mkdir(foldername_date); % create folder to save data if it does not exist

% gamma_d=Gamma;
% gamma_u=W;

%%
%sample initial conditions
 sx0= (2*randi([0, 1], [N,1])-1);
 sy0=(2*randi([0, 1], [N,1])-1);
 sz0(1:N,1)=1;

 %rotation if needed
sxp=(sx0*cos(phi)-sy0*sin(phi));
syp=(sx0*sin(phi)+sy0*cos(phi));

sx=sz0*sin(theta)+sxp*cos(theta);
sy=syp;
sz=sz0*cos(theta)-sxp*sin(theta);
clear sm
sm=(sx-1i*sy)/2;

% in this code, I am interested in the steady state properties of the
% results. This function does evolve system for time tau0 to prepare it 
% in steady state. 
[sm_0,sz_0] =  solver_spins_corr_emission_real_loop_prep(sm',sz',J_l,Gamma_l,V_l,mask_l,J_p,Gamma_p,V_p,mask_p);
%initial conditions
sigma_m_0=sm_0;
sigma_z_0=sz_0;
%%
   clear sm sz sy sx sxp syp sxp sy0 sz0 sx0
% [sm_0,sz_0] =  solver_spins_corr_emission_real_loop_prep(J_l,Gamma_l,V_l,mask_l,J_p,Gamma_p,V_p,mask_p);



    
    n=1;
    mu=adjust_mu-1;
    L=length(sigma_m_0);
%reserve meory for time evolution
S_m(1:L,1:Mev)=0;
S_z(1:L,1:Mev)=0;
S_mp(1:L,1:Mev)=0;
S_zp(1:L,1:Mev)=0;
t(1:Mev)=0;

%% preparation
rng('shuffle','v5normal');

%initial conditions in steady state 
sx(1:L,1)=2*real(sigma_m_0);
sy(1:L,1)=2*imag(sigma_m_0);
sz(1:L,1)=sigma_z_0;

% a bit rotated initial conditions to evaluate dynamical responce.
% it is used to calculate two-time correlation function
sx_p=(sx*cos(Eps)-sy*sin(Eps));
sy_p=(sx*sin(Eps)+sy*cos(Eps));
sz_p=sz;

for m=2:M
    mu=mu+1;
% noise that comes from the correlated loss     
noise_gamma_d_diag1=is_noise*mask_l.*(randn(L,1))*sqrt(dt); % noise in diagonal basis
noise_gamma_d_diag2=is_noise*mask_l.*(randn(L,1))*sqrt(dt); % noise in diagonal basis
Gamma_l=is_noise*Gamma_l;
noiseGammaPrev_real_d=(noise_gamma_d_diag1.'*V_l).';
noiseGammaPrev_imag_d=(noise_gamma_d_diag2.'*V_l).';

% noise that comes from correlated pump
noise_gamma_p_diag1=is_noise*mask_p.*(randn(L,1))*sqrt(dt); % noise in diagonal basis
noise_gamma_p_diag2=is_noise*mask_p.*(randn(L,1))*sqrt(dt); % noise in diagonal basis
Gamma_p=is_noise*Gamma_p;
noiseGammaPrev_real_p=(noise_gamma_p_diag1.'*V_p).';
noiseGammaPrev_imag_p=(noise_gamma_p_diag2.'*V_p).';
 %% previous step

 %dissipative dynamics 
    A11=-2*omega0*sy+sz.*(J_l*sy)+0.5*sz.*(Gamma_l*sx)+sz.*(J_p*sy)-0.5*sz.*(Gamma_p*sx);
    A12=2*omega0*sx-2*Omega*sz-sz.*(J_l*sx)+0.5*sz.*(Gamma_l*sy)-sz.*(J_p*sx)-0.5*sz.*(Gamma_p*sy);
    A13=2*Omega*sy+sy.*(J_l*sx)-sx.*(J_l*sy)-0.5*real(sx.*(Gamma_l*sx)+sy.*(Gamma_l*sy)) +sy.*(J_p*sx)-sx.*(J_p*sy)+0.5*real(sx.*(Gamma_p*sx)+sy.*(Gamma_p*sy));
   % dissipative dynamics for shifted initial conditions
    A11_p=-2*omega0*sy_p+sz_p.*(J_l*sy_p)+0.5*sz_p.*(Gamma_l*sx_p)+sz_p.*(J_p*sy_p)-0.5*sz_p.*(Gamma_p*sx_p);
    A12_p=2*omega0*sx_p-2*Omega*sz_p-sz_p.*(J_l*sx_p)+0.5*sz_p.*(Gamma_l*sy_p)-sz_p.*(J_p*sx_p)-0.5*sz_p.*(Gamma_p*sy_p);
    A13_p=2*Omega*sy_p+sy_p.*(J_l*sx_p)-sx_p.*(J_l*sy_p)-0.5*real(sx_p.*(Gamma_l*sx_p)+sy_p.*(Gamma_l*sy_p)) +sy_p.*(J_p*sx_p)-sx_p.*(J_p*sy_p)+0.5*real(sx_p.*(Gamma_p*sx_p)+sy_p.*(Gamma_p*sy_p));
    
    %% new step predictor
    sx_tilde=sx+A11*dt+noiseGammaPrev_imag_d.*sz-noiseGammaPrev_imag_p.*sz;
    sy_tilde=sy+A12*dt+noiseGammaPrev_real_d.*sz-noiseGammaPrev_real_p.*sz;
    sz_tilde=sz+A13*dt-(noiseGammaPrev_real_d.*sy+noiseGammaPrev_imag_d.*sx)+(noiseGammaPrev_real_p.*sy+noiseGammaPrev_imag_p.*sx); 
    % predictor for shifted initial conditions
    sx_tilde_p=sx_p+A11_p*dt+noiseGammaPrev_imag_d.*sz_p-noiseGammaPrev_imag_p.*sz_p;
    sy_tilde_p=sy_p+A12_p*dt+noiseGammaPrev_real_d.*sz_p-noiseGammaPrev_real_p.*sz_p;
    sz_tilde_p=sz_p+A13_p*dt-(noiseGammaPrev_real_d.*sy_p+noiseGammaPrev_imag_d.*sx_p)+(noiseGammaPrev_real_p.*sy_p+noiseGammaPrev_imag_p.*sx_p); 

    for i=1:n_iter
        
% dissipative dynamics for estimated state at time m
    A11_tilde=-2*omega0*sy_tilde+sz_tilde.*(J_l*sy_tilde)+0.5*sz_tilde.*(Gamma_l*sx_tilde)+sz_tilde.*(J_p*sy_tilde)-0.5*sz_tilde.*(Gamma_p*sx_tilde);
    A12_tilde=2*omega0*sx_tilde-2*Omega*sz_tilde-sz_tilde.*(J_l*sx_tilde)+0.5*sz_tilde.*(Gamma_l*sy_tilde)-sz_tilde.*(J_p*sx_tilde)-0.5*sz_tilde.*(Gamma_p*sy_tilde);
    A13_tilde=2*Omega*sy_tilde+sy_tilde.*(J_l*sx_tilde)-sx_tilde.*(J_l*sy_tilde)-0.5*real(sx_tilde.*(Gamma_l*sx_tilde)+sy_tilde.*(Gamma_l*sy_tilde)) +sy_tilde.*(J_p*sx_tilde)-sx_tilde.*(J_p*sy_tilde)+0.5*real(sx_tilde.*(Gamma_p*sx_tilde)+sy_tilde.*(Gamma_p*sy_tilde));
 % the same for shifted initial conditions
    A11_tilde_p=-2*omega0*sy_tilde_p+sz_tilde_p.*(J_l*sy_tilde_p)+0.5*sz_tilde_p.*(Gamma_l*sx_tilde_p)+sz_tilde_p.*(J_p*sy_tilde_p)-0.5*sz_tilde_p.*(Gamma_p*sx_tilde_p);
    A12_tilde_p=2*omega0*sx_tilde_p-2*Omega*sz_tilde_p-sz_tilde_p.*(J_l*sx_tilde_p)+0.5*sz_tilde_p.*(Gamma_l*sy_tilde_p)-sz_tilde_p.*(J_p*sx_tilde_p)-0.5*sz_tilde_p.*(Gamma_p*sy_tilde_p);
    A13_tilde_p=2*Omega*sy_tilde_p+sy_tilde_p.*(J_l*sx_tilde_p)-sx_tilde_p.*(J_l*sy_tilde_p)-0.5*real(sx_tilde_p.*(Gamma_l*sx_tilde_p)+sy_tilde_p.*(Gamma_l*sy_tilde_p)) +sy_tilde_p.*(J_p*sx_tilde_p)-sx_tilde_p.*(J_p*sy_tilde_p)+0.5*real(sx_tilde_p.*(Gamma_p*sx_tilde_p)+sy_tilde_p.*(Gamma_p*sy_tilde_p));
 
    
    
    %% new step predictor
    sx_new=sx+0.5*(A11+A11_tilde)*dt+0.5*noiseGammaPrev_imag_d.*(sz+sz_tilde)-0.5*noiseGammaPrev_imag_p.*(sz+sz_tilde);
    sy_new=sy+0.5*(A12+A12_tilde)*dt+0.5*noiseGammaPrev_real_d.*(sz+sz_tilde)-0.5*noiseGammaPrev_real_p.*(sz+sz_tilde);
    sz_new=sz+0.5*(A13+A13_tilde)*dt-0.5*(noiseGammaPrev_real_d.*(sy+sy_tilde)+noiseGammaPrev_imag_d.*(sx+sx_tilde))+ 0.5*(noiseGammaPrev_real_p.*(sy+sy_tilde)+noiseGammaPrev_imag_p.*(sx+sx_tilde)); 

    sx_new_p=sx_p+0.5*(A11_p+A11_tilde_p)*dt+0.5*noiseGammaPrev_imag_d.*(sz_p+sz_tilde_p)-0.5*noiseGammaPrev_imag_p.*(sz_p+sz_tilde_p);
    sy_new_p=sy_p+0.5*(A12_p+A12_tilde_p)*dt+0.5*noiseGammaPrev_real_d.*(sz_p+sz_tilde_p)-0.5*noiseGammaPrev_real_p.*(sz_p+sz_tilde_p);
    sz_new_p=sz_p+0.5*(A13_p+A13_tilde_p)*dt-0.5*(noiseGammaPrev_real_d.*(sy_p+sy_tilde_p)+noiseGammaPrev_imag_d.*(sx_p+sx_tilde_p))+ 0.5*(noiseGammaPrev_real_p.*(sy_p+sy_tilde_p)+noiseGammaPrev_imag_p.*(sx_p+sx_tilde_p)); 

    if max(abs(sz_new-sz_tilde))<epsilon_max
        if max(abs(sy_new-sy_tilde))<epsilon_max
            if max(abs(sx_new-sx_tilde))<epsilon_max
%                 i
                break;
            end
        end
    end

    

    sx_tilde=sx_new;
    sy_tilde=sy_new;
    sz_tilde=sz_new;

    sx_tilde_p=sx_new_p;
    sy_tilde_p=sy_new_p;
    sz_tilde_p=sz_new_p;

    end


    sx=sx_new;
    sy=sy_new;
    sz=sz_new;


    sx_p=sx_new_p;
    sy_p=sy_new_p;
    sz_p=sz_new_p;

    if mu==adjust_mu
       S_m(1:L,n)= 0.5*(sx-1i*sy);
       S_z(1:L,n)= sz;

       S_mp(1:L,n)= 0.5*(sx_p-1i*sy_p);
       S_zp(1:L,n)= sz_p;

       t(n)=t_in+(m-1)*dt;

        mu=0;
        n=n+1;
    end

end

% saving results for this trajectory.
close all;
name_file=strcat(foldername_date, strcat('/data_',num2str(indx_par),'.mat'));
save(name_file,'t','S_m', 'S_z', 'S_mp','S_zp','sigma_m_0')

end