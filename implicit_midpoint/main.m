clear;

%% Initial conditions
parameters_spin;
disp('start')

if DO_PARALLEL==0
    % parpool(4); % Create a parallel pool with 3 workers
p = gcp('nocreate'); % Get current parallel pool without creating a new one
if isempty(p)
    % If no pool exists, create a new one with 4 workers
    parpool(4);
else
    % If a pool exists, delete it and create a new one
    delete(p);
    parpool(4);
end
end

% wmax=1.1*N*gamma_d;;
% mmin=0;
  pump_rates=0.;%[0:0.1*gamma_d:1.1*gamma_d*N];%[1,2,5,10];%[0.1:0.1:1,2:1:5,10:5:25,50:25:100];
% pump_rates=1;
% NNp=[1,5,10,25,50];
aalpha=0;%:0.02:1.5;


 for ind2=1:length(aalpha)

for ind=1:length(pump_rates)

mu=1;
tic;


Gamma=gamma_d;

W=pump_rates(ind);
Alpha=aalpha(ind2);
gamma_u=W;


%% preparation
if DO_PARALLEL 

Gamma_p=prepare_W(Alpha,W);%can be correlated/local/collective spin pump

Gamma_p=gpuArray(Gamma_p);

J_p=zeros(N,'gpuArray');
[V_p,D_p]=eig(Gamma_p);
D_p=diag(D_p);
V_p=inv(V_p);
coeff_p=D_p;
mask_p=sqrt(abs(coeff_p));

Gamma_l=ones(N,'gpuArray')*gamma_d;%can be correlated/local/collective spin loss
% Gamma_l=gamma_d*eye(N,'gpuArray');

J_l=zeros(N,'gpuArray');
[V_l,D_l]=eig(Gamma_l);
D_l=diag(D_l);
V_l=inv(V_l);
coeff_l=D_l;
mask_l=sqrt(abs(coeff_l));

else




Gamma_p=prepare_W(Alpha,W);

J_p=zeros(N);
[V_p,D_p]=eig(Gamma_p);
D_p=diag(D_p);
V_p=inv(V_p);
coeff_p=D_p;
mask_p=sqrt(abs(coeff_p));
Gamma_l=gamma_d*ones(N);%can be correlated/local/collective spin loss
% Gamma_l=gamma_d*eye(N);
J_l=zeros(N);
[V_l,D_l]=eig(Gamma_l);
D_l=diag(D_l);
V_l=inv(V_l);
coeff_l=D_l;
mask_l=sqrt(abs(coeff_l));
end


%% parallel computation
if DO_PARALLEL
    for i = 1:Nrep 
        solver_spins_corr_emission_real_loop_par(J_l,Gamma_l,V_l,mask_l,J_p,Gamma_p,V_p,mask_p,Eps,i,Alpha)
    end
else

    parfor i = 1:Nrep
    solver_spins_corr_emission_real_loop_par(J_l,Gamma_l,V_l,mask_l,J_p,Gamma_p,V_p,mask_p,Eps,i,Alpha)
    end
end
toc;

tic;
% loop_read_cluster(J_l,Gamma_l,V_l,mask_l,J_p,Gamma_p,V_p,mask_p,Alpha);
loop_read(J_l,Gamma_l,V_l,mask_l,J_p,Gamma_p,V_p,mask_p,Alpha);

toc;
end
 end








