function [data, A] = three_source_model(T, type)
%if type = 0, one non-interacting source else source 1 -> source 2 ->
%source 3
M = 3; %number of sources;
P=2; %order of the model
T0=1000; %length of ignored start 
Rmax=50;
%Generate stable AR matrix
while Rmax>3
lambdamax=10;

while lambdamax > 1 || lambdamax < 0.9
  A=[];
  for k=1:P
    aloc = zeros(M);
    aloc([1 5 9]) = -0.9; %diagonal elements
    if type==0
        aloc([2]) = randn(1, 1);
    else
        aloc([2 6]) = randn(2, 1);
    end
    A=[A,aloc];
  end
  E=eye(M*P);AA=[A;E(1:end-M,:)];lambda=eig(AA);lambdamax=max(abs(lambda));
end


%Generate MVAR data
[data, A] = gen_mvar_data(A, T, T0, P, M);
r= [norm(data(1,:)) norm(data(2,:)) norm(data(3,:))];
if r(1) > r(2)
R(1,1) = r(1)/r(2);
 else
   R(1,1) = r(2)/r(1);
end

if r(2) > r(3)
R(1,2) = r(2)/r(3);
 else
   R(1,2) = r(3)/r(2);
end


if r(1) > r(3)
R(1,3) = r(1)/r(3);
 else
   R(1,3) = r(3)/r(1);
end

Rmax = max(R);
end



 

 









