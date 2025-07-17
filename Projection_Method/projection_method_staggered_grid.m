function projection_method_staggered_grid()
clc
clf

%------------------------------------------------------------
%                TEMPORAL AND GRID PARAMETERS
%------------------------------------------------------------
dt = 0.0001;        % Time step
printDump = 1000;   % Printing frequency
Time1 = 0;          % Speed checks Ustar Vstar solver
Time2 = 0;          % Speed checks Pressure solver
Time3 = 0;          % Speed checks Projection Step
Time4 = 0;          % Speed checks Plotting and Saving

Nx = 128;  % # of cells (faces) in x --> doesnt count ghost boundaries
Ny = 128;  % # of cells (faces) in y --> doesnt count ghost boundaries
Lx = 2;
Ly = 2;
xVec = linspace(0,Lx,Nx+1); %Pressure grid xVec (cell center locations)
yVec = linspace(0,Ly,Ny+1); %Pressure grid yVec (cell center locations)
dx = xVec(2) - xVec(1);
dy = yVec(2) - yVec(1);

xIds = 2:Nx; % Interior cell center nodes in x
yIds = 2:Ny; % Interior cell center nodes in y

%----------------------------
% Initialize Matrices
%----------------------------
u0 = zeros(Ny+1,Nx+2);
v0 = zeros(Ny+2,Nx+1);
uS = zeros(Ny+1,Nx+2);
vS = zeros(Ny+2,Nx+1);

TFinal = 5;

%----------------------------
% Storage Initialization
%----------------------------
uP = u0;                      %previous u solution value
vP = v0;                      %previous v solution value

t = 0;                        %start at time 0
ct = 1;                       %index for storing data
tVec(ct) = t;                 %vector for time

uLeftMAX = 0;                 %left u  max boundary condition
uRightMAX = 0;                %right u max boundary condition
uTopMAX = 2.5;                %top u max boundary condition
uBottomMAX = -2.5;            %bottom u max boundary condition

vLeftMAX = 0;                 %left v max boundary condition
vRightMAX = 0;                %right v max boundary condition
vTopMAX = 0;                  %top v max boundary condition
vBottomMAX = 0;               %bottom v max boundary condition

[XGrid,YGrid] = meshgrid(xVec, yVec); %matrices for plotting

mu = 0.02;                    %Dynamic Viscosity mu
rho = 1;                      %Density of liquid

BC_Vec = [uLeftMAX uRightMAX uTopMAX uBottomMAX vLeftMAX vRightMAX vTopMAX vBottomMAX]; %Boundary condtion vector for computing reynolds number

Re = (rho*Lx*abs(max(abs(BC_Vec))))/mu      %Reynolds number

%perform time stepping----------------------------------------------

ctsave = 0;
save_vtk_files(dx, dy, t, ctsave, uP, 'u')
save_vtk_files(dx, dy, t, ctsave, vP, 'v')
save_vtk_files(dx, dy, t, ctsave, uP, 'vort') %Vort, Pressure, uVecMag arent made yet, but everything starts at 0
save_vtk_files(dx, dy, t, ctsave, vP, 'P')
save_vtk_files(dx, dy, t, ctsave, uP, 'uMag')


%savevtk_vector_field(dx, dy, t, ctsave, uP, vP, 'uVec')

while t < TFinal
    %----------------------------------------------%
    %   Step 0: Update the boundary conditions     %
    %----------------------------------------------%
    uLeft = uLeftMAX * tanh(10*t);     %left u boundary condition
    uRight = uRightMAX * tanh(10*t);   %right u boundary condition
    uTop = uTopMAX * tanh(10*t);       %top u boundary condition
    uBottom = uBottomMAX * tanh(10*t); %bottom u boundary condition
    vLeft = vLeftMAX  * tanh(10*t);    %left v boundary condition
    vRight = vRightMAX * tanh(10*t);   %right v boundary condition
    vTop = vTopMAX * tanh(10*t);       %top v boundary condition
    vBottom = vBottomMAX * tanh(10*t); %bottom v boundary condition

    %----------------------------------------------%
    %   Step 1: Solve for auxilary variables       %
    %----------------------------------------------%
    % Dirichlet Boundary conditions(tops) and      %  
    % Neumann Boundary conditions(bottoms)         %
    %----------------------------------------------%

    %-------------------Set BCs--------------------%
    %    -> ghost nodes at u(j,i) and u(j, Nx+2)   %
    %    -> ghost nodes at v(j,i) and v(Ny+2, i)   %
    %----------------------------------------------%
    uS(:,1) = 2*uLeft - uP(:,2);        % Left BC (ghost node)
    %uS(:,1) = uP(:,2);
    uS(:,Nx+1) = 2*uRight - uP(:,Nx+1); % Right BC (ghost node)
    %uS(:,Nx+1) = uP(:,Nx);
    uS(1,:) = uBottom;                  % Bottom BC (physical boundary)
    %uS(1,:) = uP(2,:);
    uS(Ny+1,:) = uTop;                  % Top BC (physical boundary)
    %uS(Ny+1,:) = uP(Ny,:);
    

    vS(:,1) = vLeft;                    % Left BC (physical boundary)
    %vS(:,1) = vP(:,2);
    vS(:,Nx+1) = vRight;                % Right BC (physical boundary)
    %vS(:,Nx+1) = vP(:,Nx);
    vS(1,:) = 2*vBottom - vP(2,:);      % Bottom BC (ghost node)
    %vS(1,:) = vP(2,:);
    vS(Ny+1,:) = 2*vTop - vP(Ny+1,:);   % Top BC (ghost node)
    %vS(Ny+1,:) = vP(Ny,:);

    %----------------Solve for uStar---------------%
    %----------------------------------------------%
    tic
    for i = 2:Nx+1     % interior points; dont include ghost nodes
        for j = 2:Ny   % interior points; no physical boundaries
            
            % Calculate Laplacian of u_{j,i+1/2}
            u_xx = ( uP(j,i+1) - 2*uP(j,i) + uP(j,i-1) ) / dx^2;
            u_yy = ( uP(j+1,i) - 2*uP(j,i) + uP(j-1,i) ) / dy^2;

            % Calculate (d/dx)u_({j,i+1/2})^2
            u_sqr_x = ( ( uP(j,i+1) + uP(j,i) )^2 - ( uP(j,i) + uP(j,i-1) )^2 ) /(4*dx);

            % Calculate (d/dy) (uv)_{j,i+1/2}
            uv_y = ( ( uP(j+1,i) + uP(j,i) ) * ( vP(j,i) + vP(j,i-1) ) ...
                   - ( uP(j,i) + uP(j-1,i) ) * ( vP(j-1,i) + vP(j-1,i-1) ) ) / (4*dy);

            % Calculate uStar_{j,i+1/2}
            uS(j,i) = uP(j,i) + dt*( (mu/rho) * (u_xx + u_yy) - u_sqr_x - uv_y );

        end
    end

    %----------------Solve for vStar---------------%
    %----------------------------------------------%
    for j = 2:Ny+1     % interior points; dont include ghost nodes
        for i = 2:Nx   % interior points; no physical boundaries
            
            % Calculate Laplacian of v_{j+1/2,i}
            v_xx = ( vP(j,i+1) - 2*vP(j,i) + vP(j,i-1) ) / dx^2;
            v_yy = ( vP(j+1,i) - 2*vP(j,i) + vP(j-1,i) ) / dy^2;

            % Calculate (d/dy)v_({j+1/2,i})^2
            v_sqr_y = ( ( vP(j+1,i) + vP(j,i) )^2 - ( vP(j,i) + vP(j-1,i) )^2 ) /(4*dy);

            % Calculate (d/dx) (uv)_{j+1/2,i}
            uv_x = ( ( uP(j,i) + uP(j-1,i) ) * ( vP(j,i+1) + vP(j,i) ) ...
                   - ( uP(j,i-1) + uP(j-1,i-1) ) * ( vP(j,i) + vP(j,i-1) ) ) / (4*dx);

            % Calculate uStar_{j,i+1/2}
            vS(j,i) = vP(j,i) + dt*( (mu/rho) * (v_xx + v_yy) - v_sqr_y - uv_x );

        end
    end
    
    t1 = toc;
    Time1 = Time1 + t1;

    %-----------------------------------------------------%
    %   Step 2: Solve Poisson problem for pressure        %
    %       ->  Laplacian(p) = rho/dt * div(\vec{uStar})  %
    %                   w/ Neumann BCs (dp/dn = 0)        %
    %       ->  p \in R^{(Ny+1)x(Nx+1)}                   %
    %-----------------------------------------------------%
    tic
    Pressure = SOR_2D_Poisson_Solver_for_staggered_grid(Nx,Ny,dx,dy,dt,rho,uS,vS,xIds,yIds);
    t2 = toc;
    Time2 = Time2 + t2;
    %---------------------------------------------------------%
    %   Step 3: Projection Step (velocity correction)         % 
    %                                                         %
    %      -> (uN-uStar)/dt = -1/rho * grad(Pressure)         %
    %                                                         %
    %      ->  u, and uStar \in R^{(Ny+1)x(Nx+2)}             %
    %      ->  v, and vStar \in R^{(Ny+2)x(Nx+1)}             %
    %      ->  P     \in R^{(Ny+1)x(Nx+1)}                    %
    %                                                         %
    %      ->  BCs passed from (uStar, vStar) --> (uN, vN)    %
    %       (Correction doesn't touch "outer-rim")            %
    %---------------------------------------------------------%
    tic

    dpdx = zeros(Ny+1,Nx+2);        % Initialize storage for x derivative of pressure
    dpdy = zeros(Ny+2,Nx+1);        % Initialize storage for y derivative of pressure

    dpdx(2:Ny, 2:Nx+1) = ( Pressure(2:Ny,2:Nx+1) - Pressure(2:Ny, 1:Nx)) / (dx);   %dp/dx at grid face (j,i+1/2)
    dpdy(2:Ny+1, 2:Nx) = ( Pressure(2:Ny+1,2:Nx) - Pressure(1:Ny, 2:Nx)) / (dy);   %dp/dy at grid face (j+1/2,i)
    
    uN = uS - (dt/rho) * dpdx;      %Velocity correction for u
    vN = vS - (dt/rho) * dpdy;      %Velocity correction for v

    %------------------Update BCs------------------%
    %    -> ghost nodes at u(j,i) and u(j, Nx+2)   %
    %    -> ghost nodes at v(j,i) and v(Ny+2, i)   %
    %----------------------------------------------%
    uN(:,1) = 2*uLeft - uN(:,2);        % Left BC (ghost node)
    %uN(:,1) = uN(:,2);
    uN(:,Nx+1) = 2*uRight - uN(:,Nx+1); % Right BC (ghost node)
    %uN(:,Nx+1) = uN(:,Nx);
    uN(1,:) = uBottom;                  % Bottom BC (physical boundary)
    %uN(1,:) = uN(2,:);
    uN(Ny+1,:) = uTop;                  % Top BC (physical boundary)
    %uN(Ny+1,:) = uN(Ny,:);
    

    vN(:,1) = vLeft;                    % Left BC (physical boundary)
    %vN(:,1) = vN(:,2);
    vN(:,Nx+1) = vRight;                % Right BC (physical boundary)
    %vN(:,Nx+1) = vN(:,Nx);
    vN(1,:) = 2*vBottom - vN(2,:);      % Bottom BC (ghost node)
    %vN(1,:) = vN(2,:);
    vN(Ny+1,:) = 2*vTop - vN(Ny+1,:);   % Top BC (ghost node)
    %vN(Ny+1,:) = vN(Ny,:);

    %update time
    t = t+dt;

    %store solution values
    ct = ct+1;
    tVec(ct) = t;

    %update uP -> uN
    uP = uN;
    vP = vN;
    
    t3 = toc;
    Time3 = Time3 + t3;
    %--------------------------plot solution---------------------------%
    %------------------------------------------------------------------%
    if mod(ct,printDump) == 0
       %fprintf('timestep: %d | time: %.4f \n', ct, t);
        
       %------------------------------------------------------------------%
       %Interpolate velocity grids to cell center (from cell faces)
       %        --> uP, uN, uStar \in R^{(Ny+1)x(Nx+2)}
       %                uAvg \in R^{(Ny+1)x(Nx+1)}
       %
       %        --> vP, vN, vStar \in R^{(Ny+2)x(Nx+1)}
       %                vAvg \in R^{(Ny+1)x(Nx+1)}
       %------------------------------------------------------------------%
       uAvg = 0.5*( uP(:,2:Nx+2) + uP(:,1:Nx+1) );
       vAvg = 0.5*( vP(2:Ny+2,:) + vP(1:Ny+1,:) );
       
       %-----------------------------------------------------------------------%
       % Calculate Divergence of Velocity Field: div(vec(u)) = du/dx + dv/dy
       %            xIds --> 2:Nx   and   yIds --> 2:Ny
       %            (interior grid center locations)
       %-----------------------------------------------------------------------%
       div_u = ( uAvg(yIds,xIds+1) - uAvg(yIds,xIds-1) ) / (2*dx) + ...
                ( vAvg(yIds+1,xIds) - vAvg(yIds-1,xIds) ) / (2*dy);
        
        %------------------------------------------------------------------%
        % Calculate Vorticity ( omega = dv/dx - du/dy )
        %            xIds --> 2:Nx   and   yIds --> 2:Ny
        %            (interior grid center locations)
        %------------------------------------------------------------------%
        dvdx = ( vAvg(yIds,xIds+1) - vAvg(yIds,xIds-1) ) / (2*dx);            
        dudy = ( uAvg(yIds+1,xIds) - uAvg(yIds-1,xIds) ) / (2*dy);
        vorticity = dvdx - dudy;                
        
        % Print useful information to screen     
        fprintf('\n time: %.4f (timestep: %d)\n',ct*dt,ct)
        fprintf('    -->  CFL(adv)=%.8f\n', dt/dx*( max(max(abs(uAvg))) + max(max(abs(vAvg)))  )  );
        fprintf('    -->  CFL(diff)=%.8f\n', 2*mu/rho*dt*(1/dx^2+1/dy^2)*( max(max(abs(uAvg))) + max(max(abs(vAvg)))  )  );
        fprintf('    -->  max(|u|)=%.8f\n',max( max( abs(uAvg) ) ) );
        fprintf('    -->  max(|v|)=%.8f\n',max( max( abs(vAvg) ) ) );
        fprintf('    -->  max(|p|)=%.8f\n',max( max( abs(Pressure) ) ) );
        fprintf('    -->  max(|w|)=%.8f\n',max( max( abs(vorticity) ) ) );
        fprintf('    -->    div(u)=%.10f\n\n',mean(mean(div_u)) );
    
        
        % Mesh Grid for Vorticity (Interior nodes only of cell centers
        xVortGrid = XGrid(2:Ny, 2:Nx);
        yVortGrid = YGrid(2:Ny, 2:Nx);
        
        % Magnitude of velocity field: |\vec{u}|
        uVecMag = sqrt( uAvg.^2 + vAvg.^2 );
        
        %------------------------------------------------------------------%
        % Plot Colormaps in MATLAB
        %------------------------------------------------------------------%
        tic
        f1=figure(1);
        subplot(1,3,1)
        colormap("parula")                  %setting type of color map
        s = surf(xVortGrid,yVortGrid, vorticity);
        s.EdgeColor = "none";               %gets rid of black coloring grid lines
        h = colorbar;                       %adds color bar to the side
        maxVort = 5;
        caxis([-maxVort maxVort])
        set(gca,'fontsize', 23)             %font size
        xlabel('x') 
        ylabel('y')
        title('Vorticity')
        % zlim([0,7]) 
        view(2)                             %changes view of graph 1st number is angle from x to y, second number is angle from top to bottom

        subplot(1,3,2)
        colormap("parula")                  %setting type of color map
        s = surf(XGrid,YGrid, Pressure);
        s.EdgeColor = "none";               %gets rid of black coloring grid lines
        h = colorbar;                       %adds color bar to the side
        caxis([-6 1])
        set(gca,'fontsize', 23)             %font size
        xlabel('x') 
        ylabel('y')
        title('Pressure')
        % zlim([0,7]) 
        view(2)                            
        

        subplot(1,3,3)
        colormap("parula")                  %setting type of color map
        s = surf(XGrid,YGrid, uVecMag);
        s.EdgeColor = "none";               %gets rid of black coloring grid lines
        h = colorbar;                       %adds color bar to the side
        caxis([0 0.8*max(abs(BC_Vec))])
        set(gca,'fontsize', 23)             %font size
        xlabel('x') 
        ylabel('y')
        title('|vec(u)|')
        % zlim([0,7]) 
        view(2)

        f1.Position = [450 250 1500 400];
        pause(0.01)
        hold off

        %------------------------------------------------------------------%
        % Increment Storage counter and save data to vtk format
        %------------------------------------------------------------------%
        ctsave = ctsave + 1;
        save_vtk_files(dx, dy, t, ctsave, uP, 'u')
        save_vtk_files(dx, dy, t, ctsave, vP, 'v')
        save_vtk_files(dx, dy, t, ctsave, vorticity, 'vort')
        save_vtk_files(dx, dy, t, ctsave, Pressure, 'P')
        save_vtk_files(dx, dy, t, ctsave, uVecMag, 'uMag')
        savevtk_vector_field(dx, dy, t, ctsave, uAvg, vAvg, 'uVec')

        t4 = toc;
        Time4 = Time4 + t4;

        fprintf(' Speed Checks!\n')
        fprintf('    --> Solving U and V Star: %.10f\n', Time1/printDump );
        fprintf('    --> Pressure Solver: %.10f\n', Time2/printDump );
        fprintf('    --> Projection Step: %.10f\n', Time3/printDump );
        fprintf('    --> Ploting and Saving: %.10f\n\n', Time4/printDump );
        fprintf('---------------------------------------------\n');
    end
end