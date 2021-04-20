% Generates a periodic triangulation for the unit square with a circular
% hole. This function uses the distmesh package: http://persson.berkeley.edu/distmesh/
% Input:
%  * r = radius of the hole
%  * h = mesh size
% Output:
%  * [c4n, n4e] = coordinate arrays of nodes and connectivity information,
%  * s = FEM stiffness matrix
%  * [left_bdry, right_bdry] left and right boundary nodes that match for
%  periodicity
%  * [lower_bdry, upper_bdry] lower and upper boundary nodes that match for
%  periodicity
function [c4n, n4e, s, left_bdry, right_bdry, lower_bdry, upper_bdry] = generate_mesh(r,h)
d=2;
npts = round(1/h)+1;
fx = linspace(0,1,npts)';
f0 = zeros(npts,1);
f1 = ones(npts,1);

%% --------Mesh-Generation---------
holes_hand=@(p,r) dcircle(p,0.5,0.5,r);
fd=@(p) ddiff(drectangle(p,0,1,0,1),holes_hand(p,r));
fh=@(p) sqrt(holes_hand(p,0.5*r));

box=[0,0;1,1];
fix=uniquetol([f0,fx; fx,f0; fx,f1; f1,fx],h/4,'ByRows',true);
try
    disp('running distmesh')
    [c4n,n4e]=distmesh2d(fd,@huniform,h,box,fix, false); % last argument = plotting
    disp('distmesh done')
    close all
    [c4n,n4e]=fixmesh(c4n,n4e);
    nC = length(c4n);
    nE = length(n4e);
catch
    error('Distmesh function not found. See https://github.com/ionhandshaker/distmesh')
end

%% Construct periodic version of n4e:
elemP = n4e;

right_bdry = find(c4n(:,1)==1);
left_bdry  = find(c4n(:,1)==0);

for i=1:length(left_bdry)
    idx = right_bdry(i);
    lidx = left_bdry(i);
    elemP(elemP==idx) = lidx;
end

upper_bdry = find(c4n(:,2)==1);
lower_bdry  = find(c4n(:,2)==0);

for i=1:length(left_bdry)
    idx = upper_bdry(i);
    lidx = lower_bdry(i);
    elemP(elemP==idx) = lidx;
end

%% Build stiffness matrix for periodic problem:
ctr = 0; 
ctr_max = (d+1)^2*length(n4e);
I = zeros(ctr_max,1); 
J = zeros(ctr_max,1); 
X = zeros(ctr_max,1);
for t = 1:length(n4e)
    X_T = [ones(1,d+1);c4n(n4e(t,:),:)'];
    grads_T = X_T\[zeros(1,d);eye(d)];
    vol_T = det(X_T)/factorial(d); % n4e is used to compute volumes
    for m = 1:d+1
        for n = 1:d+1
            ctr = ctr+1; 
            I(ctr) = elemP(t,m);   % elemP is used for indexing
            J(ctr) = elemP(t,n);
            X(ctr) = vol_T*grads_T(m,:)*grads_T(n,:)';
        end
    end
end
s = sparse(I,J,X,nC,nC);

end





