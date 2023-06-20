function a_GP = Galerkin_Projection_ROM(X,Y,...
    uMean,vMean,uPOD,vPOD,a_POD,n,m,dx,dy,Re,t)
%% Galerkin Projection model for 2D incompressible flow
% -- Based on "Reduced-order modeling for flow control", Noack et al.
%
% Author: Aditya Nair, University of Nevada Reno
%
% AN provides no guarantees for this code.  Use as-is and for academic
% research use only; no commercial use allowed without permission.  For
% citations, please use the reference below:
%
% Ref: A. G. Nair, S. L. Brunton & ``````````````````````````````````````````````````````````K. Taira,
% "Networked-oscillator-based modeling and control of unsteady wake flows,"
% Physical Review E, vol 97, 2018
% (https://doi.org/10.1103/PhysRevE.97.063107)
%
% The code is written for educational clarity and not for speed.


%% Construct Galerkin Projection model terms

[Am1i,Am2i,BDmij,BAm1ij,BAm2ij,Cmijk] = ...
    Galerkin_projection(X,Y,uMean,vMean,uPOD,vPOD,n,m,dx,dy,Re);

%% Run Galerkin Projection model

initial_condition = a_POD(1,:); % set initial condition
options = odeset('RelTol',1e-8,'AbsTol',...
    1e-8*ones(1,size(initial_condition,2))); % set tolerance for integrator
[~,a_GP] = ode45(@Galerkin_proj_function,t,initial_condition,...
    options,Am1i,Am2i,BDmij,BAm1ij,BAm2ij,Cmijk);% run model

end


function [Am1i,Am2i,BDmij,BAm1ij,BAm2ij,Cmijk] = ...
    Galerkin_projection(X,Y,uMean,vMean,uPOD,vPOD,n,m,dx,dy,Re)
nModes = size(uPOD,3);         % number of POD modes for GP
Am1i = zeros(nModes,1);        % (viscous term, mean projection) - l_{i0}
Am2i = zeros(nModes,1);        % (inertial term, mean projection) - q_{i00}
BDmij = zeros(nModes,nModes);  % (viscous term, mode projection) - l_{ij}
BAm1ij = zeros(nModes,nModes); % (inertial term, modes and mean projection) - q_{ij0}
BAm2ij = zeros(nModes,nModes); % (inertial term, mean and modes projection) - q_{i0j}
Cmijk  = zeros(nModes,nModes,nModes); % (inertial term, mode projection) - q_{ijk}

% Calculate derivatives and laplacian of mean flow
duMdx = diffxy(X,uMean,2);duMdy = diffxy(Y,uMean,1);
dvMdx = diffxy(X,vMean,2);dvMdy = diffxy(Y,vMean,1);
duM2 = 4*del2(uMean,dx,dy);
dvM2 = 4*del2(vMean,dx,dy);

for i = 1:nModes
    % Compute Ai
    Am1i(i) = (1/Re)*simpson_rule((duM2(2:n-2,2:m-2).*uPOD(2:n-2,2:m-2,i) + ...
        dvM2(2:n-2,2:m-2).*vPOD(2:n-2,2:m-2,i)),dx,dy,n-4,m-4);

    Am2i(i) = simpson_rule(-((uMean(2:n-2,2:m-2).*duMdx(2:n-2,2:m-2) ...
        + vMean(2:n-2,2:m-2).*duMdy(2:n-2,2:m-2)).*uPOD(2:n-2,2:m-2,i) ...
        + (uMean(2:n-2,2:m-2).*dvMdx(2:n-2,2:m-2) ...
        + vMean(2:n-2,2:m-2).*dvMdy(2:n-2,2:m-2)).*vPOD(2:n-2,2:m-2,i)),dx,dy,n-4,m-4);

    for j = 1:nModes
        % Compute Bij
        % Calculate derivatives and laplacian of modes
        du2 = 4*del2(uPOD(:,:,j),dx,dy);dv2 = 4*del2(vPOD(:,:,j),dx,dy);
        dudx = diffxy(X,uPOD(:,:,j),2);dudy = diffxy(Y,uPOD(:,:,j),1);
        dvdx = diffxy(X,vPOD(:,:,j),2);dvdy = diffxy(Y,vPOD(:,:,j),1);

        BDmij(i,j) = (1/Re)*simpson_rule((du2(2:n-2,2:m-2).*uPOD(2:n-2,2:m-2,i) + ...
            dv2(2:n-2,2:m-2).*vPOD(2:n-2,2:m-2,i)),dx,dy,n-4,m-4);

        BAm1ij(i,j) = simpson_rule(-((uMean(2:n-2,2:m-2).*dudx(2:n-2,2:m-2) ...
            + vMean(2:n-2,2:m-2).*dudy(2:n-2,2:m-2)).*uPOD(2:n-2,2:m-2,i) + ...
            (uMean(2:n-2,2:m-2).*dvdx(2:n-2,2:m-2) ...
            + vMean(2:n-2,2:m-2).*dvdy(2:n-2,2:m-2)).*vPOD(2:n-2,2:m-2,i)),dx,dy,n-4,m-4);

        BAm2ij(i,j) = simpson_rule(-((uPOD(2:n-2,2:m-2,j).*duMdx(2:n-2,2:m-2) ...
            + vPOD(2:n-2,2:m-2,j).*duMdy(2:n-2,2:m-2)).*uPOD(2:n-2,2:m-2,i) + ...
            (uPOD(2:n-2,2:m-2,j).*dvMdx(2:n-2,2:m-2) ...
            + vPOD(2:n-2,2:m-2,j).*dvMdy(2:n-2,2:m-2)).*vPOD(2:n-2,2:m-2,i)),dx,dy,n-4,m-4);

        for k = 1:nModes
            % Compute Cijk
            dukdx = diffxy(X,uPOD(:,:,k),2);dukdy = diffxy(Y,uPOD(:,:,k),1);
            dvkdx = diffxy(X,vPOD(:,:,k),2);dvkdy = diffxy(Y,vPOD(:,:,k),1);

            Cmijk(i,j,k) = simpson_rule(-((uPOD(2:n-2,2:m-2,j).*dukdx(2:n-2,2:m-2) ...
                + vPOD(2:n-2,2:m-2,j).*dukdy(2:n-2,2:m-2)).*uPOD(2:n-2,2:m-2,i) + ...
                (uPOD(2:n-2,2:m-2,j).*dvkdx(2:n-2,2:m-2) ...
                + vPOD(2:n-2,2:m-2,j).*dvkdy(2:n-2,2:m-2)).*vPOD(2:n-2,2:m-2,i)),dx,dy,n-4,m-4);
        end
    end
end
end


function [dadt] = Galerkin_proj_function(t,a,A1i,A2i,BDij,BA1ij,BA2ij,Cijk)
%% Global variables
c = a;nModes = length(c);
dadt = zeros(nModes,1);
for iMode = 1:nModes
    dadt(iMode) = dadt(iMode) + A1i(iMode) + A2i(iMode);
    for j = 1:nModes
        dadt(iMode) = dadt(iMode) + BDij(iMode,j)*c(j) + ...
            BA2ij(iMode,j)*c(j)  + BA1ij(iMode,j)*c(j);
        for k = 1:nModes
            dadt(iMode) = dadt(iMode) + Cijk(iMode,j,k)*c(j)*c(k);
        end
    end
end
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
end


function dy = diffxy(x,y,varargin)
% DIFFXY - accurate numerical derivative/differentiation of Y w.r.t X.
%
%   DY = DIFFXY(X,Y) returns the derivative of Y with respect to X using a
%        pseudo second-order accurate method. DY has the same size as Y.
%   DY = DIFFXY(X,Y,DIM) returns the derivative along the DIM-th dimension
%        of Y. The default is differentiation along the first
%        non-singleton dimension of Y.
%   DY = DIFFXY(X,Y,DIM,N) returns the N-th derivative of Y w.r.t. X.
%        The default is 1.
%
%   Y may be an array of any dimension.
%   X can be any of the following:
%       - array X with size(X) equal to size(Y)
%       - vector X with length(X) equal to size(Y,DIM)
%       - scalar X denotes the spacing increment
%   DIM and N are both integers, with 1<=DIM<=ndims(Y)
%
%   DIFFXY has been developed especially to handle unequally spaced data,
%   and features accurate treatment for end-points.
%
%   Example:
%   % Data with equal spacing
%     x = linspace(-1,2,20);
%     y = exp(x);
%
%     dy = diffxy(x,y);
%     dy2 = diffxy(x,dy);  % Or, could use >> dy2 = diffxy(x,y,[],2);
%     figure('Color','white')
%     plot(x,(y-dy)./y,'b*',x,(y-dy2)./y,'b^')
%
%     Dy = gradient(y)./gradient(x);
%     Dy2 = gradient(Dy)./gradient(x);
%     hold on
%     plot(x,(y-Dy)./y,'r*',x,(y-Dy2)./y,'r^')
%     title('Relative error in derivative approximation')
%     legend('diffxy: dy/dx','diffxy: d^2y/dx^2',...
%            'gradient: dy/dx','gradient: d^2y/dx^2')
%
%   Example:
%   % Data with unequal spacing.
%     x = 3*sort(rand(20,1))-1;
%     % Run the example above from y = exp(x)
%
%   See also DIFF, GRADIENT
%        and DERIVATIVE on the File Exchange

% for Matlab (should work for most versions)
% version 1.0 (Nov 2010)
% (c) Darren Rowland
% email: darrenjrowland@hotmail.com
%
% Keywords: derivative, differentiation

[h,dy,N,perm] = parse_inputs(x,y,varargin);
if isempty(dy)
    return
end
n = size(h,1);
i1 = 1:n-1;
i2 = 2:n;

for iter = 1:N
    v = diff(dy)./h;
    if n>1
        dy(i2,:) = (h(i1,:).*v(i2,:)+h(i2,:).*v(i1,:))./(h(i1,:)+h(i2,:));
        dy(1,:) = 2*v(1,:) - dy(2,:);
        dy(n+1,:) = 2*v(n,:) - dy(n,:);
    else
        dy(1,:) = v(1,:);
        dy(n+1,:) = dy(1,:);
    end
end

% Un-permute the derivative array to match y
dy = ipermute(dy,perm);

end

%%% Begin local functions %%%
function [h,dy,N,perm] = parse_inputs(x,y,v)

numvarargs = length(v);
if numvarargs > 2
    error('diffxy:TooManyInputs', ...
        'requires at most 2 optional inputs');
end

h = [];
N = [];
perm = [];

% derivative along first non-singleton dimension by default
dim = find(size(y)>1);
% Return if dim is empty
if isempty(dim)
    dy = [];
    return
end
dim = dim(1);

% Set defaults for optional arguments
optargs = {dim 1};
newVals = ~cellfun('isempty', v);
optargs(newVals) = v(newVals);
[dim, N] = optargs{:};

% Error check on inputs
if dim<1 || dim>ndims(y) || dim~=fix(dim) || ~isreal(dim)
    error('diffxy:InvalidOptionalArg',...
        'dim must be specified as a non-negative integer')
end
if N~=fix(N) || ~isreal(N)
    error('diffxy:InvalidOptionalArg',...
        'N must be an integer')
end

% permutation which will bring the target dimension to the front
perm = 1:length(size(y));
perm(dim) = [];
perm = [dim perm];
dy = permute(y,perm);


if length(x)==1  % Scalar expansion to match size of diff(dy,[],1)
    sizeh = size(dy);
    sizeh(1) = sizeh(1) - 1;
    h = repmat(x,sizeh);
elseif ndims(x)==2 && any(size(x)==1) % Vector x expansion
    if length(x)~=size(dy,1)
        error('diffxy:MismatchedXandY',...
            'length of vector x must match size(y,dim)')
    end
    x = x(:);
    sizeh = size(dy);
    sizeh(1) = 1;
    h = repmat(diff(x),sizeh);
else
    if size(y) ~= size(x)
        error('diffxy:MismatchedXandY',...
            'mismatched sizes of arrays x and y');
    end
    % Permute x as for y, then diff
    h = diff(permute(x,perm),[],1);
end
end