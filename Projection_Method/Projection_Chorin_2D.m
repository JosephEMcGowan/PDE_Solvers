function Projection_Chorin_2D()

%---------------------------------------------
% Grid Parameters
%---------------------------------------------
Nx = 100;            % Spatial Step (x)
Ny = 100;            % Spatial Step (y)
Lx = 1;              % Domain size (x)
Ly = 1;              % Domain size (y)
dx = Lx / (Nx - 1);  % Grid resolution
dy = Ly / (Ny - 1);  % Grid resolution

%---------------------------------------------
% Simulation Parameters (including fluid params) 
%---------------------------------------------
uLid = 2.5;          % Maximum lid flow rate
rho = 1;             % Fluid density (keep as 1)
mu = 0.01;           % Dynamic Viscosity
nu = mu/rho;         % Kinematic Viscosity (nu=mu/rho)


%---------------------------------------------
% Temporal Parameters
%---------------------------------------------
dt = 0.0001;         % Time step
TFinal = 5;          % Final time
t = 0;               % Initialize time
dp = 1000;           % How often to print to screen (# of dt's)


%---------------------------------------------
% Initialize Storage
%---------------------------------------------
u = zeros(Ny, Nx); 
v = zeros(Ny, Nx); 
p = zeros(Ny, Nx);
u_new = u; 
v_new = v; 


%---------------------------------------------
% Initial condition (e.g., lid-driven cavity)
%---------------------------------------------
u(end, :) = uLid*tanh(t);  % Top lid moving to the right

%---------------------------------------------
% Interior Point Matrices
%---------------------------------------------
xIds = 2:Nx-1;
yIds = 2:Ny-1;

%---------------------------------------------
% Grid Matrices (for plotting purposes)
%---------------------------------------------
[xMAT,yMAT] = meshgrid(2:Nx-1, 2:Ny-1);  % for vorticity
[x, y] = meshgrid(0:dx:Lx, 0:dy:Ly);     % for velocity field

%---------------------------------------------
% Time stepping
%---------------------------------------------
fprintf(' Re: %.2f\n', Lx*uLid/nu);

ct=0;
while t < TFinal
    
    %-----------------------------------
    % Increment Counters and Time
    %-----------------------------------
    ct=ct+1;
    t = t + dt;

    %----------------------------------------------------------
    % Step 1: Compute tentative velocities u*, v* (explicit Euler)
    %----------------------------------------------------------
    
    tic
    
    % Re-allocate storage
    u_star = u;
    v_star = v;

    % Advection terms (1st-order upwind)
    u_x = ( u(yIds,xIds) - u(yIds,xIds-1) ) / dx;
    u_y = ( u(yIds,xIds) - u(yIds-1,xIds) ) / dy;
    v_x = ( v(yIds,xIds) - v(yIds,xIds-1) ) / dx;
    v_y = ( v(yIds,xIds) - v(yIds-1,xIds) ) / dy;

    % Advection terms (greater numerical stability?)
    %u_sqr_x = ( u(yIds,xIds+1).^2 - u(yIds,xIds-1).^2 ) / (2*dx);
    %v_sqr_y = ( v(yIds+1,xIds).^2 - v(yIds-1,xIds).^2 ) / (2*dy);
    %uv_x = ( u(yIds,xIds+1).*v(yIds,xIds+1) - u(yIds,xIds-1).*u(yIds,xIds-1) ) / (2*dx);
    %uv_y = ( u(yIds+1,xIds).*v(yIds+1,xIds) - u(yIds-1,xIds).*v(yIds-1,xIds) ) / (2*dy);
    
    % Diffusion terms (central difference)
    u_xx = (u(yIds,xIds+1) - 2*u(yIds,xIds) + u(yIds,xIds-1)) / dx^2;
    u_yy = (u(yIds+1,xIds) - 2*u(yIds,xIds) + u(yIds-1,xIds)) / dy^2;
    v_xx = (v(yIds,xIds+1) - 2*v(yIds,xIds) + v(yIds,xIds-1)) / dx^2;
    v_yy = (v(yIds+1,xIds) - 2*v(yIds,xIds) + v(yIds-1,xIds)) / dy^2;

    % Find u-Star explicitly (w/ 1-st order upwind for advection)
    u_star(yIds,xIds) = u(yIds,xIds) + dt * ( ...
        - u(yIds,xIds) .* u_x - v(yIds,xIds) .* u_y + ...
          nu * (u_xx + u_yy) );
 
    % Find v-Star explicitly (w/ 1-st order upwind for advection)      
    v_star(yIds,xIds) = v(yIds,xIds) + dt * ( ...
        - u(yIds,xIds) .* v_x - v(yIds,xIds) .* v_y + ...
          nu * (v_xx + v_yy) );


    % Find u-Star explicitly (w/ more stable approach)
    %u_star(yIds,xIds) = u(yIds,xIds) + dt * ( ...
    %      - v_sqr_y - uv_x  + nu * (u_xx + u_yy) );

    % Find v-Star explicitly (w/ more stable approach)     
    %v_star(yIds,xIds) = v(yIds,xIds) + dt * ( ...
    %      - u_sqr_x - uv_y + nu * (v_xx + v_yy) );

    % Enforce boundary conditions for u_star, v_star
    u_star(1,:) = 0; u_star(end,:) = uLid*tanh(t);
    u_star(:,1) = 0; u_star(:,end) = 0;
    v_star(1,:) = 0; v_star(end,:) = 0;
    v_star(:,1) = 0; v_star(:,end) = 0;

    
    t1=toc;
    
    %-------------------------------------------------
    % Step 2: Solve Poisson equation for pressure
    %       --> Here using basic Jacobi solver
    %-------------------------------------------------
    tic
    div_u_star = (u_star(yIds,xIds+1) - u_star(yIds,xIds-1))/(2*dx) + ...
                 (v_star(yIds+1,xIds) - v_star(yIds-1,xIds))/(2*dy);

    rhs = (rho/dt) * div_u_star;

    % Iterative solver: Jacobi (simple for small grids)
    for it = 1:100
        p_old = p;
        p(yIds,xIds) = 0.25 * ( ...
            p_old(yIds,xIds+1) + p_old(yIds,xIds-1) + ...
            p_old(yIds+1,xIds) + p_old(yIds-1,xIds) - ...
            dx^2 * rhs );
        
        % Pressure Neumann BCs (dp/dn = 0)
        p(:,1) = p(:,2);
        p(:,end) = p(:,end-1);
        p(1,:) = p(2,:);
        p(end,:) = p(end-1,:);
        
    end
    t2=toc;
   
    %-------------------------------------------
    % Step 3: Correct velocity
    %-------------------------------------------
    tic
    u_new(yIds,xIds) = u_star(yIds,xIds) - dt/rho * (p(yIds,xIds+1) - p(yIds,xIds-1)) / (2*dx);
    v_new(yIds,xIds) = v_star(yIds,xIds) - dt/rho * (p(yIds+1,xIds) - p(yIds-1,xIds)) / (2*dy);

    % Apply BCs
    u_new(1,:) = 0;                 % bottom
    u_new(end,:) = uLid*tanh(2*t);  % top --> lid driven cavity flow
    u_new(:,1) = 0;                 % left
    u_new(:,end) = 0;               % right
    %
    v_new(1,:) = 0;                 % bottom
    v_new(end,:) = 0;               % top
    v_new(:,1) = 0;                 % left
    v_new(:,end) = 0;               % right

    % Update fields
    u = u_new;
    v = v_new;
    
    t3=toc;
    

    %-------------------------------------------
    % PLOT (if selected time-step)
    %-------------------------------------------
    if mod(ct,dp)==0
        
        % Calculate Divergence of Velocity Field ( divU = du/dx + dv/dy )
        div_u = ( u(yIds,xIds+1) - u(yIds,xIds-1) ) / (2*dx) + ...
                ( v(yIds+1,xIds) - v(yIds-1,xIds) ) / (2*dy);
        
        % Calculate Vorticity ( omega = dv/dx - du/dy )
        dvdx = ( v(yIds,xIds+1) - v(yIds,xIds-1) )/ (2*dx);            
        dudy = ( u(yIds+1,xIds) - u(yIds-1,xIds) )/ (2*dy);
        vorticity = dvdx - dudy;        
             
        % Print useful information to screen     
        fprintf('\n time: %.4f (timestep: %d)\n',ct*dt,ct)
        fprintf('    -->  CFL(adv)=%.4f\n', dt/dx*( max(max(abs(u))) + max(max(abs(v)))  )  );
        fprintf('    --> CFL(diff)=%.4f\n', 2*nu*dt*(1/dx^2+1/dy^2)*( max(max(abs(u))) + max(max(abs(v)))  )  );
        fprintf('    -->  max(|u|)=%.4f\n',max( max( abs(u) ) ) );
        fprintf('    -->  max(|v|)=%.4f\n',max( max( abs(v) ) ) );
        fprintf('    -->  max(|p|)=%.4f\n',max( max( abs(p) ) ) );
        fprintf('    -->  max(|w|)=%.4f\n',max( max( abs(vorticity) ) ) );
        fprintf('    -->    div(u)=%.10f\n',mean(mean(div_u)) );
        
        % Plot velocity field
        f1=figure(1);
        subplot(1,2,1)
        quiver(x, y, u, v);
        axis equal tight;
        title('2-D Velocity Field');
        xlabel('x'); ylabel('y');

        % Plot Vorticity
        subplot(1,2,2)
        colormap('jet');
        surf(xMAT,yMAT, vorticity ,'EdgeColor','none');
        set(gca,'FontSize',18);
        xlabel('x');
        ylabel('y');
        title('2-D Vorticity');
        view(2)
        xlim([0 Nx]);
        ylim([0 Ny]);
        %
        h=colorbar;
        %maxVort=0.7*max(max(abs(vorticity)));
        maxVort = 25;
        caxis([-maxVort maxVort])
        h.Ticks = linspace( -maxVort, maxVort, 5 );

        f1.Position = [250 250 1000 400];
        
    end
    
end
