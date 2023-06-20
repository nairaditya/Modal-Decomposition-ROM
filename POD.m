function [D,FourCoef,uPOD,vPOD] = POD(nt,nModes,n,m,dx,dy,u,v)
%% Input: nt = number of snapshots 
%         nModes = number of modes
%         n = x-number of points
%         m = y-number of points
%         dx = spacing in x-dir
%         dy = spacing in y-dir
%         u = u(1:n,1:m,nt) u-velocity
%         v = v(1:n,1:m,nt) v-velocity

%% Output: 
%         D = Modal energy
%         FourCoef = temporal coefficients
%         uPOD,vPOD = u-vel POD and v-vel POD modes
        

%% Compute covariance matrix
C = zeros(nt,nt);         % init
for i = 1:nt
     i
    for j = i:nt
        C(i,j) = simpson_rule(u(:,:,i).*u(:,:,j) + ...
            v(:,:,i).*v(:,:,j),dx,dy,n-1,m-1);
        C(j,i) = C(i,j);          % use symmetry
    end
end
C = C./ nt;

%% Compute eigenvalues and eigenvectors of C
[V,D] = eig(C);
[D,I] = sort(diag(D),'descend'); 
V = V(:,I);

%% Scale Fourier Coefficients: <a_i a_i> = lambda_i
FourCoef = zeros(nt,nModes);
for i = 1:nModes
    i
    for j = 1:nt
        ModeFactor = sqrt(nt*D(i));
        FourCoef(j,i) = V(j,i)*ModeFactor;
    end
end

%% Compute POD modes
uPOD = zeros(n,m,nModes); 
vPOD = zeros(n,m,nModes); 
for i = 1:nModes
    i
    for j = 1:nt
        uPOD(:,:,i) = uPOD(:,:,i) + V(j,i)*u(:,:,j);
        vPOD(:,:,i) = vPOD(:,:,i) + V(j,i)*v(:,:,j);
    end
    % Normalize
    modeFactor = 1 ./ sqrt(nt*D(i));
    uPOD(:,:,i) = uPOD(:,:,i)*modeFactor;
    vPOD(:,:,i) = vPOD(:,:,i)*modeFactor;
end


function out = simpson_rule(U,dx,dy,NX,NY)
s1 = (U(1,1)+U(1,NY+1)+U(NX+1,1)+U(NX+1,NY+1));
ixo = 2:2:NX;
ixe = 3:2:NX-1;
iyo = 2:2:NY;
iye = 3:2:NY-1;
s2 = 2*(sum(U(1,iye))+sum(U(NX+1,iye))+sum(U(ixe,1))+sum(U(ixe,NY+1)) );
s3 = 4*(sum(U(1,iyo))+sum(U(NX+1,iyo))+sum(U(ixo,1))+sum(U(ixo,NY+1)) );
s4 = 16*sum(sum(U(ixo,iyo)))+ 4*sum(sum(U(ixe,iye)));
s5 = 8*sum(sum(U(ixe,iyo)))+ 8*sum(sum(U(ixo,iye)));
out = s1 + s2 + s3 + s4 + s5;
out = out*dx*dy/9.0;
