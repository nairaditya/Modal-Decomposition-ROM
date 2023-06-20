function [ ProjectedModes,DEv,ExactModes,Norm ] = DMD( X,Y,varargin )
% Dynamic Mode Decomposition as presented by 
% "On Dynamic Mode Decomposition: theory and applications" by Tu et al.,
% 10.3934/jcd.2014.1.391



% inputs : 
% Data Sets X and Y- should have the same size 
% each column of X is a set of measurements done at an instant 
% each column of Y is the image of the corresponding columns from X

% or 

% ( X,Y,Tol ) - with Tol (optional) being the threshold for filtering thru SVD - the
% default value is 1e-10



% outputs: 
% 1 - Projected Dynamic Modes
% 2 - Exact Dynamic Modes
% 3 - Dynamic Eigenvalues
% 4 - Norms - Euclidean (vector) norm of each mode, used to sort the data

% setting SVD hard threshold
if isempty(varargin)
    Tol=1e-10;
else
    Tol=varargin{1};
end


disp(['Tolerance used for filtering in DMD:',num2str(Tol)])

[U,S,V]=svd(X,'econ');

k = find(diag(S)>Tol,1,'last');
disp(['DMD subspace dimension:',num2str(k)])

U = U(:,1:k); V=V(:,1:k); S=S(1:k,1:k);

Atilde = ((U'*Y) *V )* diag((1./diag(S)));


[W,DEv]=eig(Atilde);
DEv = diag(DEv);

ProjectedModes = U*W;

ExactModes = bsxfun(@times,1./DEv.',((Y*V)*S^(-1))*W);



if nargout>3
   b = pinv(ProjectedModes)*X(:,end);     % terminal coordinates in the Koopman subspace
    [Norm,Index]=sort(abs(b),'descend');
    DEv = DEv(Index);
    ProjectedModes = ProjectedModes(:,Index);
    disp('modes sorted based on energy contribution to the last snapshot')
end
end


%=========================================================================%
% Hassan Arbabi - 08-17-2015
% Mezic research group
% UC Santa Barbara
% arbabiha@gmail.com
%=========================================================================%