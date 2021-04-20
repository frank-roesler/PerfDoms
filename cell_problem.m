% Solves the cell problem on unit cell with circular hole. If a mesh is
% available, it is used, otherwise the function generate_mesh2(r,h) is
% called.
% Input:
%  * r = radius of hole,
%  * h = mesh size,
%  * xi = 2-vector that determines homogenized matrix,
%  * plotting = bool; if true, the solution is plotted.
% Output:
%  * [I1,I2] = integral of grad(u)+xi 
%    (this is equal to A*xi, where A is the homogenized matrix)
function [I1,I2] = cell_problem(r,h,xi, plotting)
try
    load(['mesh_r=',num2str(r),'.mat'], 'c4n','n4e','s','left_bdry','right_bdry','lower_bdry','upper_bdry');
catch
    disp('found no saved mesh. Generating new one...')
    [c4n, n4e, s, left_bdry, right_bdry, lower_bdry, upper_bdry] = generate_mesh(r,h);
end
nC = size(c4n,1);
nE = size(n4e,1);
TR = triangulation(n4e,c4n);
bdry = freeBoundary(TR);

b1 = bdry(:,1);
b2 = bdry(:,2);
[r1, c1] = find(vecnorm(c4n(b1,:)-[0.5,0.5],2,2)<r+h);

Nb = bdry(r1,:);
nNb = size(Nb,1);
%% Solve periodic FEM problem:
b = zeros(nC,1); % right-hand side
u = zeros(nC,1); % solution
for j = 1:nNb
    vol_S = norm(c4n(Nb(j,1),:)-c4n(Nb(j,2),:));
    mp_S = sum(c4n(Nb(j,:),:),1)/2;
    for k = 1:2
        b(Nb(j,k)) = b(Nb(j,k))+(1/2)*vol_S*neumann_data(mp_S,r,xi);
    end
end

% make s less singular:
ur_bdry = [right_bdry; upper_bdry];
iNodes = setdiff(1:nC,ur_bdry); % independent nodes (not identified by periodicity)

try
    Aeq = ones(1,length(u(iNodes)))/length(u(iNodes));
    beq = 0; % sets mean(u) to 0.

    u(iNodes) = lsqlin(s(iNodes,iNodes),b(iNodes),[],[],Aeq,beq,[],[]);
    % minimizes |su - b| under the constraint mean(u)=0.
catch
    disp('need to use singular matrix...')
    u(iNodes) = s(iNodes,iNodes)\b(iNodes); % less precise (s singular)
end
u(right_bdry) = u(left_bdry) ;
u(upper_bdry) = u(lower_bdry) ;

if plotting
    trisurf(n4e, c4n(:,1),c4n(:,2), u);
    colormap(autumn)
    shading interp
    light
end

%% Compute grad(u):
Du = zeros(nE,2);
for j = 1:nE
    X_T = [ones(1,2+1);c4n(n4e(j,:),:)'];
    grads_T = X_T\[zeros(1,2);eye(2)];
    corners = n4e(j,:);
    Du(j,:) = sum(u(corners).*grads_T);
end

% Compute integral of xi+grad(u):
I1 = 0;
I2 = 0;
for j = 1:nE
    X_T = [ones(1,2+1);c4n(n4e(j,:),:)'];
    vol_T = det(X_T)/factorial(2);
    I1 = I1 + (Du(j,1)+xi(1))*vol_T;
    I2 = I2 + (Du(j,2)+xi(2))*vol_T;
end





end


%%
function val = neumann_data(x,r,xi)
    val = (x(1)*xi(1) + x(2)*xi(2) - 0.5)/r;
end