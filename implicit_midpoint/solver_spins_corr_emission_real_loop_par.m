function []= solver_spins_corr_emission_real_loop_par(J_l,Gamma_l,V_l,mask_l,J_p,Gamma_p,V_p,mask_p,Eps,indx_par,Alpha)
% Here we evalute dynamics for a single trajectory
parameters_spin;
%% trying spherical mid_point method
%output folder to save intermediate results
% foldername_date=strcat('/Users/oksanachelpanova/Matlab_code/parallel_code/');

% gamma_d=Gamma;
% gamma_u=W;

%%
%sample initial conditions
if DO_PARALLEL
 sx0= (2*randi([0, 1], [N,ncopies],'gpuArray')-1);
 sy0=(2*randi([0, 1], [N,ncopies],'gpuArray')-1);
 sz0=ones(N,ncopies,'gpuArray');
 s_m=zeros(N,ncopies,Mev,'gpuArray');
 s_z=zeros(N,ncopies,Mev,'gpuArray');
s_mp=zeros(N,ncopies,Mev,'gpuArray');
s_zp=zeros(N,ncopies,Mev,'gpuArray');
t=zeros(1,Mev,'gpuArray');

else
    % L=N;
    sx0= (2*randi([0, 1], [N,ncopies])-1);
 sy0=(2*randi([0, 1], [N,ncopies])-1);
 sz0=ones(N,ncopies,1);
s_m=zeros(N,ncopies,Mev);
s_z=zeros(N,ncopies,Mev);
s_mp=zeros(N,ncopies,Mev);
s_zp=zeros(N,ncopies,Mev);
t=zeros(1,Mev);

end

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
[sm_0,sz_0] =  solver_spins_corr_emission_real_loop_prep(sm,sz,J_l,Gamma_l,V_l,mask_l,J_p,Gamma_p,V_p,mask_p);
%initial conditions
sigma_m_0=sm_0;
sigma_z_0=sz_0;
%%
   clear sm sz sy sx sxp syp sxp sy0 sz0 sx0
% [sm_0,sz_0] =  solver_spins_corr_emission_real_loop_prep(J_l,Gamma_l,V_l,mask_l,J_p,Gamma_p,V_p,mask_p);



    
    n=1;
    mu=adjust_mu-1;
    % L=length(sigma_m_0);
%reserve meory for time evolution


%% preparation
if DO_PARALLEL
    gpurng('shuffle');
else
    rng('shuffle','v5normal');
end
%initial conditions in steady state 
sx(1:N,1:ncopies,1)=2*real(sigma_m_0);
sy(1:N,1:ncopies,1)=2*imag(sigma_m_0);
sz(1:N,1:ncopies,1)=sigma_z_0;

% a bit rotated initial conditions to evaluate dynamical responce.
% it is used to calculate two-time correlation function
sx_p=(sx*cos(Eps)-sy*sin(Eps));
sy_p=(sx*sin(Eps)+sy*cos(Eps));
sz_p=sz;

mask_l=is_noise*mask_l;
mask_p=is_noise*mask_p;


% 
% 
% 
% noise_gamma_d_diag1=mask_l.*(randn(N,ncopies,Knoise,'like',sx))*sqrt(dt); % noise in diagonal basis
% noise_gamma_d_diag2=mask_l.*(randn(N,ncopies,Knoise,'like',sx))*sqrt(dt); % noise in diagonal basis
% % Gamma_l=is_noise*Gamma_l;
% NoiseGammaPrev_real_d=(noise_gamma_d_diag1.'*V_l).';
% NoiseGammaPrev_imag_d=(noise_gamma_d_diag2.'*V_l).';
% 
% % noise that comes from correlated pump
% noise_gamma_p_diag1=mask_p.*(randn(N,ncopies,Knoise,'like',sx))*sqrt(dt); % noise in diagonal basis
% noise_gamma_p_diag2=mask_p.*(randn(N,ncopies,Knoise,'like',sx))*sqrt(dt); % noise in diagonal basis
% % Gamma_p=is_noise*Gamma_p;
% NoiseGammaPrev_real_p=(noise_gamma_p_diag1.'*V_p).';
% NoiseGammaPrev_imag_p=(noise_gamma_p_diag2.'*V_p).';


    ng_d1 = randn(N, ncopies, Knoise, 'like', sx) .* mask_l * sqrt(dt);
    ng_d2 = randn(N, ncopies, Knoise, 'like', sx) .* mask_l * sqrt(dt);

    ng_p1 = randn(N, ncopies, Knoise, 'like', sx) .* mask_p * sqrt(dt);
    ng_p2 = randn(N, ncopies, Knoise, 'like', sx) .* mask_p * sqrt(dt);

    % Transform noise to real space:
    % Each slice: (N×ncopies)' * (N×N) = (ncopies×N),
    % then transpose back to N×ncopies
    noise_real_d = pagemtimes(permute(ng_d1,[2 1 3]), V_l); 
    noise_real_d = permute(noise_real_d,[2 1 3]);
    noise_imag_d = pagemtimes(permute(ng_d2,[2 1 3]), V_l); 
    noise_imag_d = permute(noise_imag_d,[2 1 3]);

    noise_real_p = pagemtimes(permute(ng_p1,[2 1 3]), V_p); 
    noise_real_p = permute(noise_real_p,[2 1 3]);
    noise_imag_p = pagemtimes(permute(ng_p2,[2 1 3]), V_p); 
    noise_imag_p = permute(noise_imag_p,[2 1 3]);

knoise=1;
   Aa11=@(x,y,z)-2*omega0*y+z.*(J_l*y)+0.5*z.*(Gamma_l*x)+z.*(J_p*y)-0.5*z.*(Gamma_p*x);
    Aa12=@(x,y,z)2*omega0*x-2*Omega*z-z.*(J_l*x)+0.5*z.*(Gamma_l*y)-z.*(J_p*x)-0.5*z.*(Gamma_p*y);
    Aa13=@(x,y,z)2*Omega*y+y.*(J_l*x)-x.*(J_l*y)-0.5*real(x.*(Gamma_l*x)+y.*(Gamma_l*y)) +y.*(J_p*x)-x.*(J_p*y)+0.5*real(x.*(Gamma_p*x)+y.*(Gamma_p*y));


for m=2:M
    mu=mu+1;

    noiseGammaPrev_real_d(1:N,1:ncopies)=noise_real_d(:,:,knoise);
    noiseGammaPrev_imag_d(1:N,1:ncopies)=noise_imag_d(:,:,knoise);
    noiseGammaPrev_real_p(1:N,1:ncopies)=noise_real_p(:,:,knoise);
    noiseGammaPrev_imag_p(1:N,1:ncopies)=noise_imag_p(:,:,knoise);

    knoise=knoise+1;

    if knoise==Knoise+1
 ng_d1 = randn(N, ncopies, Knoise, 'like', sx) .* mask_l * sqrt(dt);
    ng_d2 = randn(N, ncopies, Knoise, 'like', sx) .* mask_l * sqrt(dt);

    ng_p1 = randn(N, ncopies, Knoise, 'like', sx) .* mask_p * sqrt(dt);
    ng_p2 = randn(N, ncopies, Knoise, 'like', sx) .* mask_p * sqrt(dt);

    % Transform noise to real space:
    % Each slice: (N×ncopies)' * (N×N) = (ncopies×N),
    % then transpose back to N×ncopies
    noise_real_d = pagemtimes(permute(ng_d1,[2 1 3]), V_l); 
    noise_real_d = permute(noise_real_d,[2 1 3]);
    noise_imag_d = pagemtimes(permute(ng_d2,[2 1 3]), V_l); 
    noise_imag_d = permute(noise_imag_d,[2 1 3]);

    noise_real_p = pagemtimes(permute(ng_p1,[2 1 3]), V_p); 
    noise_real_p = permute(noise_real_p,[2 1 3]);
    noise_imag_p = pagemtimes(permute(ng_p2,[2 1 3]), V_p); 
    noise_imag_p = permute(noise_imag_p,[2 1 3]);

knoise=1;
    end

% noise that comes from the correlated loss  
% if DO_PARALLEL
% noise_gamma_d_diag1=mask_l.*(randn(N,ncopies,'gpuArray'))*sqrt(dt); % noise in diagonal basis
% noise_gamma_d_diag2=mask_l.*(randn(N,ncopies,'gpuArray'))*sqrt(dt); % noise in diagonal basis
% % Gamma_l=is_noise*Gamma_l;
% noiseGammaPrev_real_d=(noise_gamma_d_diag1.'*V_l).';
% noiseGammaPrev_imag_d=(noise_gamma_d_diag2.'*V_l).';
% 
% % noise that comes from correlated pump
% noise_gamma_p_diag1=mask_p.*(randn(N,ncopies,'gpuArray'))*sqrt(dt); % noise in diagonal basis
% noise_gamma_p_diag2=mask_p.*(randn(N,ncopies,'gpuArray'))*sqrt(dt); % noise in diagonal basis
% % Gamma_p=is_noise*Gamma_p;
% noiseGammaPrev_real_p=(noise_gamma_p_diag1.'*V_p).';
% noiseGammaPrev_imag_p=(noise_gamma_p_diag2.'*V_p).';
% else
% 
% noise_gamma_d_diag1=mask_l.*(randn(N,ncopies))*sqrt(dt); % noise in diagonal basis
% noise_gamma_d_diag2=mask_l.*(randn(N,ncopies))*sqrt(dt); % noise in diagonal basis
% % Gamma_l=is_noise*Gamma_l;
% noiseGammaPrev_real_d=(noise_gamma_d_diag1.'*V_l).';
% noiseGammaPrev_imag_d=(noise_gamma_d_diag2.'*V_l).';
% 
% % noise that comes from correlated pump
% noise_gamma_p_diag1=mask_p.*(randn(N,ncopies))*sqrt(dt); % noise in diagonal basis
% noise_gamma_p_diag2=mask_p.*(randn(N,ncopies))*sqrt(dt); % noise in diagonal basis
% % Gamma_p=is_noise*Gamma_p;
% noiseGammaPrev_real_p=(noise_gamma_p_diag1.'*V_p).';
% noiseGammaPrev_imag_p=(noise_gamma_p_diag2.'*V_p).';
% end
 %% previous step

 %dissipative dynamics 
   %  A11=-2*omega0*sy+sz.*(J_l*sy)+0.5*sz.*(Gamma_l*sx)+sz.*(J_p*sy)-0.5*sz.*(Gamma_p*sx);
   %  A12=2*omega0*sx-2*Omega*sz-sz.*(J_l*sx)+0.5*sz.*(Gamma_l*sy)-sz.*(J_p*sx)-0.5*sz.*(Gamma_p*sy);
   %  A13=2*Omega*sy+sy.*(J_l*sx)-sx.*(J_l*sy)-0.5*real(sx.*(Gamma_l*sx)+sy.*(Gamma_l*sy)) +sy.*(J_p*sx)-sx.*(J_p*sy)+0.5*real(sx.*(Gamma_p*sx)+sy.*(Gamma_p*sy));
   % % dissipative dynamics for shifted initial conditions
   %  A11_p=-2*omega0*sy_p+sz_p.*(J_l*sy_p)+0.5*sz_p.*(Gamma_l*sx_p)+sz_p.*(J_p*sy_p)-0.5*sz_p.*(Gamma_p*sx_p);
   %  A12_p=2*omega0*sx_p-2*Omega*sz_p-sz_p.*(J_l*sx_p)+0.5*sz_p.*(Gamma_l*sy_p)-sz_p.*(J_p*sx_p)-0.5*sz_p.*(Gamma_p*sy_p);
   %  A13_p=2*Omega*sy_p+sy_p.*(J_l*sx_p)-sx_p.*(J_l*sy_p)-0.5*real(sx_p.*(Gamma_l*sx_p)+sy_p.*(Gamma_l*sy_p)) +sy_p.*(J_p*sx_p)-sx_p.*(J_p*sy_p)+0.5*real(sx_p.*(Gamma_p*sx_p)+sy_p.*(Gamma_p*sy_p));
   % 

    A11=Aa11(sx,sy,sz);
    A12=Aa12(sx,sy,sz);
    A13=Aa13(sx,sy,sz);
    % dissipative dynamics for shifted initial conditions
    % A11_p=-2*omega0*sy_p+sz_p.*(J_l*sy_p)+0.5*sz_p.*(Gamma_l*sx_p)+sz_p.*(J_p*sy_p)-0.5*sz_p.*(Gamma_p*sx_p);
    % A12_p=2*omega0*sx_p-2*Omega*sz_p-sz_p.*(J_l*sx_p)+0.5*sz_p.*(Gamma_l*sy_p)-sz_p.*(J_p*sx_p)-0.5*sz_p.*(Gamma_p*sy_p);
    % A13_p=2*Omega*sy_p+sy_p.*(J_l*sx_p)-sx_p.*(J_l*sy_p)-0.5*real(sx_p.*(Gamma_l*sx_p)+sy_p.*(Gamma_l*sy_p)) +sy_p.*(J_p*sx_p)-sx_p.*(J_p*sy_p)+0.5*real(sx_p.*(Gamma_p*sx_p)+sy_p.*(Gamma_p*sy_p));
    A11_p=Aa11(sx_p,sy_p,sz_p);
    A12_p=Aa12(sx_p,sy_p,sz_p);
    A13_p=Aa13(sx_p,sy_p,sz_p);
    %% new step predictor
    sx_tilde=sx+A11*dt+noiseGammaPrev_imag_d.*sz-noiseGammaPrev_imag_p.*sz;
    sy_tilde=sy+A12*dt+noiseGammaPrev_real_d.*sz-noiseGammaPrev_real_p.*sz;
    sz_tilde=sz+A13*dt-(noiseGammaPrev_real_d.*sy+noiseGammaPrev_imag_d.*sx)+(noiseGammaPrev_real_p.*sy+noiseGammaPrev_imag_p.*sx); 
    % predictor for shifted initial conditions
    sx_tilde_p=sx_p+A11_p*dt+noiseGammaPrev_imag_d.*sz_p-noiseGammaPrev_imag_p.*sz_p;
    sy_tilde_p=sy_p+A12_p*dt+noiseGammaPrev_real_d.*sz_p-noiseGammaPrev_real_p.*sz_p;
    sz_tilde_p=sz_p+A13_p*dt-(noiseGammaPrev_real_d.*sy_p+noiseGammaPrev_imag_d.*sx_p)+(noiseGammaPrev_real_p.*sy_p+noiseGammaPrev_imag_p.*sx_p); 

    for i=1:n_iter
        
        sx_mid=0.5*(sx+sx_tilde);
        sy_mid=0.5*(sy+sy_tilde);
        sz_mid=0.5*(sz+sz_tilde);

        sx_mid_p=0.5*(sx_p+sx_tilde_p);
        sy_mid_p=0.5*(sy_p+sy_tilde_p);
        sz_mid_p=0.5*(sz_p+sz_tilde_p);
% dissipative dynamics for estimated state at time m
    A11_tilde=Aa11(sx_mid,sy_mid,sz_mid);
    A12_tilde=Aa12(sx_mid,sy_mid,sz_mid);
    A13_tilde=Aa13(sx_mid,sy_mid,sz_mid);

    A11_tilde_p=Aa11(sx_mid_p,sy_mid_p,sz_mid_p);
    A12_tilde_p=Aa12(sx_mid_p,sy_mid_p,sz_mid_p);
    A13_tilde_p=Aa13(sx_mid_p,sy_mid_p,sz_mid_p);
 %    A12_tilde=2*omega0*sx_tilde-2*Omega*sz_tilde-sz_tilde.*(J_l*sx_tilde)+0.5*sz_tilde.*(Gamma_l*sy_tilde)-sz_tilde.*(J_p*sx_tilde)-0.5*sz_tilde.*(Gamma_p*sy_tilde);
 %    A13_tilde=2*Omega*sy_tilde+sy_tilde.*(J_l*sx_tilde)-sx_tilde.*(J_l*sy_tilde)-0.5*real(sx_tilde.*(Gamma_l*sx_tilde)+sy_tilde.*(Gamma_l*sy_tilde)) +sy_tilde.*(J_p*sx_tilde)-sx_tilde.*(J_p*sy_tilde)+0.5*real(sx_tilde.*(Gamma_p*sx_tilde)+sy_tilde.*(Gamma_p*sy_tilde));
 % % the same for shifted initial conditions
    % A11_tilde_p=-2*omega0*sy_tilde_p+sz_tilde_p.*(J_l*sy_tilde_p)+0.5*sz_tilde_p.*(Gamma_l*sx_tilde_p)+sz_tilde_p.*(J_p*sy_tilde_p)-0.5*sz_tilde_p.*(Gamma_p*sx_tilde_p);
    % A12_tilde_p=2*omega0*sx_tilde_p-2*Omega*sz_tilde_p-sz_tilde_p.*(J_l*sx_tilde_p)+0.5*sz_tilde_p.*(Gamma_l*sy_tilde_p)-sz_tilde_p.*(J_p*sx_tilde_p)-0.5*sz_tilde_p.*(Gamma_p*sy_tilde_p);
    % A13_tilde_p=2*Omega*sy_tilde_p+sy_tilde_p.*(J_l*sx_tilde_p)-sx_tilde_p.*(J_l*sy_tilde_p)-0.5*real(sx_tilde_p.*(Gamma_l*sx_tilde_p)+sy_tilde_p.*(Gamma_l*sy_tilde_p)) +sy_tilde_p.*(J_p*sx_tilde_p)-sx_tilde_p.*(J_p*sy_tilde_p)+0.5*real(sx_tilde_p.*(Gamma_p*sx_tilde_p)+sy_tilde_p.*(Gamma_p*sy_tilde_p));
 
    
    
    %% new step predictor
    sx_new=sx+A11_tilde*dt+0.5*noiseGammaPrev_imag_d.*(sz+sz_tilde)-0.5*noiseGammaPrev_imag_p.*(sz+sz_tilde);
    sy_new=sy+A12_tilde*dt+0.5*noiseGammaPrev_real_d.*(sz+sz_tilde)-0.5*noiseGammaPrev_real_p.*(sz+sz_tilde);
    sz_new=sz+A13_tilde*dt-0.5*(noiseGammaPrev_real_d.*(sy+sy_tilde)+noiseGammaPrev_imag_d.*(sx+sx_tilde))+ 0.5*(noiseGammaPrev_real_p.*(sy+sy_tilde)+noiseGammaPrev_imag_p.*(sx+sx_tilde)); 

    sx_new_p=sx_p+A11_tilde_p*dt+0.5*noiseGammaPrev_imag_d.*(sz_p+sz_tilde_p)-0.5*noiseGammaPrev_imag_p.*(sz_p+sz_tilde_p);
    sy_new_p=sy_p++A12_tilde_p*dt+0.5*noiseGammaPrev_real_d.*(sz_p+sz_tilde_p)-0.5*noiseGammaPrev_real_p.*(sz_p+sz_tilde_p);
    sz_new_p=sz_p+A13_tilde_p*dt-0.5*(noiseGammaPrev_real_d.*(sy_p+sy_tilde_p)+noiseGammaPrev_imag_d.*(sx_p+sx_tilde_p))+ 0.5*(noiseGammaPrev_real_p.*(sy_p+sy_tilde_p)+noiseGammaPrev_imag_p.*(sx_p+sx_tilde_p)); 

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
       s_m(1:N,1:ncopies,n)= 0.5*(sx-1i*sy);
       s_z(1:N,1:ncopies,n)= sz;

       s_mp(1:N,1:ncopies,n)= 0.5*(sx_p-1i*sy_p);
       s_zp(1:N,1:ncopies,n)= sz_p;

       t(n)=t_in+(m-1)*dt;

        mu=0;
        n=n+1;
    end

end

% saving results for this trajectory.
close all;

if DO_PARALLEL
t=gather(t);
end

for j=1:ncopies
idx=Nrep*(j-1)+indx_par;


if DO_PARALLEL
    S_m  = gather(s_m(:,j,:));
    S_z  = gather(s_z(:,j,:));
    S_mp = gather(s_mp(:,j,:));
    S_zp = gather(s_zp(:,j,:));

    S_m  = reshape(S_m,  [N,Mev]);
    S_z  = reshape(S_z,  [N,Mev]);
    S_mp = reshape(S_mp, [N,Mev]);
    S_zp = reshape(S_zp, [N,Mev]);

    name_file=strcat(foldername_date, strcat('/data_',num2str(idx),'.mat'));
save(name_file,'t','S_m', 'S_z', 'S_mp','S_zp','sigma_m_0')

else

S_m=reshape(s_m(:,j,:),[N,Mev]);
S_z=reshape(s_z(:,j,:),[N,Mev]);
S_mp=reshape(s_mp(:,j,:),[N,Mev]);
S_zp=reshape(s_zp(:,j,:),[N,Mev]);
name_file=strcat(foldername_date, strcat('/data_',num2str(idx),'.mat'));
save(name_file,'t','S_m', 'S_z', 'S_mp','S_zp','sigma_m_0')
end
end


end