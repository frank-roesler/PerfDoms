% Computes homogenized matrix by minimizing variational term in unit cell.
% Functions are discretized via finite differences; the hole in the unit
% cell is defined by the function chi.m
% Input:
%  * r = radius of hole,
%  * N = number of grid points,
%  * xi = 2-vector that determines homogenized matrix,
%  * plotting = bool; if true, the minimizer is plotted.
% Output: 
%  * Minimum of integral(|xi+grad(u)|). 
%    (This is equal to xi'*A*xi, where A is the homogenized matrix)
function integral = variational_problem(r,N,xi,plotting)
    %% Setup:
    X = linspace(0,1,N);
    Y = linspace(0,1,N);
    [XX,YY] = meshgrid(X,Y);

    mask = chi(XX,YY,r);

    e = N*ones(N,1);
    X = spdiags([-e,e], [-1,0], N,N);
    X(1,1) = -N;
    X(end,end) = N;
    X(1,2) = N;
    X(end,end-1) = -N;
    DY=[];
    for i=1:N
        DY = blkdiag(DY,X);
    end

    DX = sparse(N^2);
    for j=1:N^2-N
        DX(j,j) = -N;
        DX(j,j+N) = N;
    end
    DX(N^2-N+1:N^2,N^2-N+1:N^2) = N*eye(N);
    DX(N^2-N+1:N^2,N^2-2*N+1:N^2-N) = -N*eye(N);

    %% Periodicity constraints and minimization:

    maskv = mask(:);
    maskI = spdiags(maskv,0,N^2,N^2);

    H = [maskI*DX; maskI*DY];
    Aeq = sparse(2*N,N^2);
    Aeq(1:N,N:N:N^2) = eye(N);
    Aeq(1:N,1:N:N^2) = -eye(N);
    Aeq(N+1:2*N,1:N) = -eye(N);
    Aeq(N+1:2*N,end-N+1:end) = +eye(N);
    beq = sparse([xi(2)*ones(N,1); xi(1)*ones(N,1)]);

    [vnew,minimum,residual,exitflag,output,lambda] = lsqlin(H,zeros(2*N^2,1),[],[],Aeq,beq,[],[]);
    unew = reshape(vnew,N,N);
    unew = unew-mean(unew(:));

    integral = minimum/N^2;
    
    if plotting
        unew(logical(1-mask))=NaN;
        surf(XX,YY,(unew-xi(1)*XX-xi(2)*YY+0.5))
        colormap(autumn)
        shading interp
        light
    end
end





function val = chi(x,y,r)
    val = (x-0.5).^2+(y-0.5).^2>r^2;
end