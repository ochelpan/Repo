function [Gamma_u] = prepare_W(Alpha,w)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    Alpha
    w
end
parameters_spin;
Gamma_u=zeros(N);
x=1:N;
for i=1:N
    for j=1:N
Gamma_u(i,j)=1/(1+abs(x(i)-x(j)))^Alpha;

    end
end

Gamma_u=w*Gamma_u;
end