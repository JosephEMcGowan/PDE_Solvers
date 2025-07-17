function pN = SOR_2D_Poisson_Solver_for_staggered_grid(Nx,Ny,dx,dy,dt,rho,uS,vS,xIds,yIds)

%---------------------------
% Initialize storage matrix
%---------------------------
rhsMAT = zeros(Ny+1,Nx+1);

%--------------------------------------------------------
% Non-homogeneous term (divergence of \vec{uStar})
%        uStar \in R^{(Ny+1)x(Nx+2)}
%        vStar \in R^{(Ny+2)x(Nx+1)}
%        P     \in R^{(Ny+1)x(Nx+1)}
%
%        (yIds,xIds): Interior node Ids for pressure
%               xIds: 2:Nx
%               yIds: 2:Ny
%--------------------------------------------------------
uS_x = (1/(dx)) * (uS(yIds, xIds+1) - uS(yIds, xIds));
vS_y = (1/(dy)) * (vS(yIds+1, xIds) - vS(yIds, xIds));

rhsMAT(yIds, xIds) = (rho/dt) * (uS_x + vS_y);

%----------------------
% Initial Guess
%----------------------
pPrev = zeros(Ny+1,Nx+1);

%----------------------------
% Perform SOR Iteration
%----------------------------
err = 1;    %initialize error (to start iterative process)
tol = 1e-6; %error tolerance (accuracy threshold)
ct = 0;     %just a counter 
pN = zeros(size(pPrev));
omega = 1.25;
maxiter = 500;
%
while err > tol && ct <= maxiter

    %Increment counter
    ct = ct + 1;
    
    %----------------------------------------------------
    % Apply SOR on interior nodes for Pressure
    %----------------------------------------------------
    for i=2:Nx
        for j=2:Ny
            pN(j,i) =  (1 - omega ) * pPrev(j,i) + omega * -0.25 * ( (dx)^2 * rhsMAT(j,i) ...
                      - pPrev(j,i+1) - pN(j,i-1) - pPrev(j+1,i) - pN(j-1,i) );
        end
    end

    %----------------------------------------------------
    % Apply Neumann BCs on Pressure
    %----------------------------------------------------
    pN(:,1) = pPrev(:,2);
    pN(:,Nx+1) = pPrev(:,Nx);
    pN(1,:) = pPrev(2,:);
    pN(Ny+1,:) = pPrev(Ny,:);
    
    %Calculate new error --> L2-error
    %err = sqrt( sum(sum( (pN - pPrev).^2 ) ) );
    
    %Calculate new error --> L-infinity error
    err = max(max( abs(pN - pPrev) ) );

    %redefine previous guess
    pPrev = pN;
    
    if ct > 100000
        fprintf(' --> error = %.4f \n', err);
    end

end