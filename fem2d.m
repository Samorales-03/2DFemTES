clc ; load('project.mat') ;
pts =  TRI.Points ;
te2p = TRI.ConnectivityList ;
[Nt, y] = size(te2p);[Np, y] = size(pts) ;
%% Assemble A matrix
col=ones([9*Nt 1]);row=ones([9*Nt 1]);val=zeros([9*Nt 1]);
for i=1:Nt
coeffi=inv([pts(te2p(i,:),1),pts(te2p(i,:),2),ones([3 1])]);
areai=abs(det([pts(te2p(i,:),1),pts(te2p(i,:),2),ones([3 1])]))/2;
col(9*(i-1)+1:9*i)=[te2p(i,1) te2p(i,1) te2p(i,1) te2p(i,2) te2p(i,2) te2p(i,2) te2p(i,3) te2p(i,3) te2p(i,3)];
row(9*(i-1)+1:9*i)=[te2p(i,1) te2p(i,2) te2p(i,3) te2p(i,1) te2p(i,2) te2p(i,3) te2p(i,1) te2p(i,2) te2p(i,3)];
val(9*(i-1)+1:9*i)=-conductivity(i)*areai*(coeffi(1:2,:)'*coeffi(1:2,:));
end
A=sparse(col,row,val,Np,Np); 
%% Assemble RHS
b = zeros([Np 1]) ; b(P_anode) = 1 ; b(P_cathode) = -1 ;
%% SET GROUND
A = A(1:end-1 , 1:end-1); b = b(1:end-1);
sol =  minres(A,b,10^(-15),2000); 
sol(end+1) = 0;