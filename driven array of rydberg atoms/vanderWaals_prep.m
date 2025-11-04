function [h,Lqcl,Lqq,Mqcl,Mqq] = vanderWaals_prep()

% parameters_spins;
parameters_spin;
J=vanderWaalscouplings(J0,Alpha,N);
h=sum(J,2)/4;

 A=zeros(4*N,4*N);
 Rp=zeros(4*N,4*N); %Keldysh rotation for Hubbard-Stratonovich fields



for i=1:N

    for j=1:N
        if i==j
            A((i-1)*(4)+1,(i-1)*(4)+1)=1i*kappa;
            A((i-1)*(4)+2,(i-1)*(4)+1)=2i*kappa;
            A((i-1)*(4)+2,(i-1)*(4)+2)=1i*kappa;
            A((i-1)*(4)+3,(i-1)*(4)+3)=1i*gamma_d;
            A((i-1)*(4)+4,(i-1)*(4)+3)=2i*gamma_d;
            A((i-1)*(4)+4,i*4)=1i*gamma_d;
        end
        if i~=j
            A((i-1)*4+1,(j-1)*4+1)=J(i,j)/8;
            A((i-1)*4+2,(j-1)*4+2)=-J(i,j)/8;
%             A((i-1)*4+1,(j-1)*4+1)=J(i,j)/4;
%             A((i-1)*4+2,(j-1)*4+2)=-J(i,j)/4;
        end

        
    end

    Rp((i-1)*4+1,(i-1)*4+1)=1;
    Rp((i-1)*4+1,(i-1)*4+2)=-1;
    Rp((i-1)*4+2,(i-1)*4+1)=1;
    Rp((i-1)*4+2,(i-1)*4+2)=1;
    Rp((i-1)*4+3,(i-1)*4+3)=1;
    Rp((i-1)*4+3,(i-1)*4+4)=-1;
    Rp((i-1)*4+4,(i-1)*4+3)=1;
    Rp((i-1)*4+4,(i-1)*4+4)=1;



    % A(i,j) 

end
% R=R./sqrt(2);
Rp=Rp./sqrt(2);

AA=inv(Rp)*inv(A)*Rp; %% must be Block-Diagonal
% Lcc=AA(3:4:end,3:4:end)

Lqq(1:N,1:N)=AA(2:4:end,2:4:end);
% Lclq(1:N,1:N)=AA(2:4:end,1:4:end);
Lqcl(1:N,1:N)=AA(1:4:end,2:4:end);

Mqq(1:N,1:N)=AA(4:4:end,4:4:end);
% Mclq(1:N,1:N)=AA(4:4:end,3:4:end);
Mqcl(1:N,1:N)=AA(3:4:end,4:4:end);

%% clcl components must be zero!




end