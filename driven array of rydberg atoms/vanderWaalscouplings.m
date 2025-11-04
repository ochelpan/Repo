function [Jmatr] = vanderWaalscouplings(J0,Alpha,N)
% power-law interaction between nearest atoms assuming periodic BC
Jmatr=zeros(N,N);

 for i=1:N
     for j=i+1:N
         r=min(abs(i-j),abs(N-j+i));
         Jmatr(i,j)=J0/(r^Alpha);
         Jmatr(j,i)=Jmatr(i,j);
     end
 end
end