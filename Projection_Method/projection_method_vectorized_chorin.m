function projection_method_vectorized_chorin()
clc
clf

%INFO BOX-----------------------------------------------------------
%Updated version of my own code for projection method
%-------------------------------------------------------------------

%choices------------------------------------------------------------
%ONLY THINGS THAT NEED TO BE EDITED TO MAKE A NEW EXAMPLE ARE dt(maybe),
%AND u0%
%-------------------------------------------------------------------
dt = 0.0001;

printDump = 1000;

Nx = 100;

Ny = 100;

Lx = 2;

Ly = 2;

xVec = linspace(-Lx/2,Lx/2,Nx+1);

yVec = linspace(-Ly/2,Ly/2,Ny+1);

dx = xVec(2) - xVec(1);

dy = yVec(2) - yVec(1);

xIds = 2:Nx;
yIds = 2:Ny;

%----------------------------
% Initialize Matrices
%----------------------------
u0 = zeros(Ny+1,Nx+1);
v0 = zeros(Ny+1,Nx+1);
uS = zeros(Ny+1,Nx+1);
vS = zeros(Ny+1,Nx+1);

TFinal = 5;

%storage initialization---------------------------------------------
%ONLY THINGS THAT NEED TO BE EDITED TO MAKE A NEW EXAMPLE ARE K AND
%BOUNDARY CONDITIONS%
%-------------------------------------------------------------------
uP = u0;                      %previous u solution value

vP = v0;                      %previous v solution value

t = 0;                        %start at time 0

ct = 1;                       %index for storing data

tVec(ct) = t;                 %vector for time

uLeftMAX = 0;                 %left u  max boundary condition

uRightMAX = 0;                %right u max boundary condition

uTopMAX = 0;                %top u max boundary condition

uBottomMAX = 0;               %bottom u max boundary condition

vLeftMAX = 0;              %left v max boundary condition

vRightMAX = -2.5;                %right v max boundary condition

vTopMAX = 0;                  %top v max boundary condition

vBottomMAX = 0;               %bottom v max boundary condition

[X,Y] = meshgrid(xVec, yVec); %matrices for plotting

mu = 0.02;                    %Dynamic Viscosity mu

rho = 1;                      %Density of liquid

Re = (rho*Lx*abs(vRightMAX))/mu      %Reynolds number

%perform time stepping----------------------------------------------

ctsave = 0;
save_vtk_files(dx, dy, t, ctsave, uP, 'u')
save_vtk_files(dx, dy, t, ctsave, vP, 'v')
save_vtk_files(dx, dy, t, ctsave, uP, 'vort') %Vort and Pressure arent made yet, but everything starts at 0
save_vtk_files(dx, dy, t, ctsave, vP, 'P')
savevtk_vector_field(dx, dy, t, ctsave, uP, vP, 'uVec')

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
    uS(:,1) = uLeft;
    uS(:,Nx+1) = uRight;
    uS(1,:) = uBottom;
    uS(Ny+1,:) = uTop;
    vS(:,1) = vLeft;
    vS(:,Nx+1) = vRight;
    vS(1,:) = vBottom;
    vS(Ny+1,:) = vTop;
    
    %----------------------------------------------------%
    % Creating the needed derivative terms for computing %
    % u star and v star, than computing them.            %
    %----------------------------------------------------%
    uu_x = uP(yIds, xIds) .* ((1/dx) .* (uP(yIds, xIds) - uP(yIds, xIds-1)));
    vu_y = vP(yIds, xIds) .* ((1/dy) .* (uP(yIds, xIds) - uP(yIds-1, xIds)));
    uv_x = uP(yIds, xIds) .* ((1/dx) .* (vP(yIds, xIds) - vP(yIds, xIds-1)));
    vv_y = vP(yIds, xIds) .* ((1/dy) .* (vP(yIds, xIds) - vP(yIds-1, xIds)));
    u_xx = (1/(dx)^2) * (uP(yIds, xIds+1) - 2 * uP(yIds, xIds) + uP(yIds, xIds-1));
    u_yy = (1/(dy)^2) * (uP(yIds+1, xIds) - 2 * uP(yIds, xIds) + uP(yIds-1, xIds));
    v_xx = (1/(dx)^2) * (vP(yIds, xIds+1) - 2 * vP(yIds, xIds) + vP(yIds, xIds-1));
    v_yy = (1/(dy)^2) * (vP(yIds+1, xIds) - 2 * vP(yIds, xIds) + vP(yIds-1, xIds));

    uS(yIds, xIds) = uP(yIds, xIds) - dt*uu_x - dt*vu_y + ((dt*mu)/ rho) * (u_xx + u_yy);
    vS(yIds, xIds) = vP(yIds, xIds) - dt*uv_x - dt*vv_y + ((dt*mu)/ rho) * (v_xx + v_yy);
    
    %----------------------------------------------%
    %   Step 2: Find the Pressure Function         %
    %----------------------------------------------%

    Pressure = SOR_2D_Poisson_Solver(Nx,Ny,dx,dy,dt,rho,uS,vS,xIds,yIds);

    %----------------------------------------------%
    %   Step 3: Find uN and vN                     %
    %----------------------------------------------%

    dpdx = zeros(Ny+1,Nx+1);        %Initialize x derivative of pressure
    dpdy = zeros(Ny+1,Nx+1);        %Initialize y derivative of pressure

    dpdx(2:Ny, 2:Nx) = ( Pressure(2:Ny,3:Nx+1) - Pressure(2:Ny, 1:Nx-1)) / (2*dx);   %Find x derivative of pressure
    dpdy(2:Ny, 2:Nx) = ( Pressure(3:Ny+1,2:Nx) - Pressure(1:Ny-1, 2:Nx)) / (2*dy);   %Find y derivative of pressure
    
    uN = uS - (dt/rho) * dpdx;      %Solution to uN
    vN = vS - (dt/rho) * dpdy;      %Solution to vN

    %update time
    t = t+dt;

    %store solution values
    ct = ct+1;
    tVec(ct) = t;

    %update uP -> uN
    uP = uN;
    vP = vN;

%     if mod(ct,printDump/10) == 0
%         fprintf(' --> timestep: %d | time: %.4f \n', ct, t);
%         fprintf('    --> cfl(adv): %.4f | cfl(diff): %.4f \n', (dt/dx)*(max(max(abs(uP))) + max(max(abs(vP)))), (dt/(dx)^2)*mu*(max(max(abs(uP))) + max(max(abs(vP)))));
%         fprintf('    --> max u = %.4f | max v = %.4f | max pressure = %.4f \n', max(max(abs(uN))), max(max(abs(vN))), max(max(abs(Pressure))));
%        
%         dudx = ( uN(2:Ny, 3:Nx+1) - uN(2:Ny, 1:Nx-1) ) / (2*dx);
%         dvdy = ( vN(3:Ny+1, 2:Nx) - vN(1:Ny-1, 2:Nx) ) / (2*dy);
% 
%         fprintf(' --> div U = %.10f \n', sum(sum(dudx + dvdy)));
%     end

    %plot PDE solution--------------------------------------------------
    %ONLY THINGS THAT NEEDED TO BE CHANGED FOR A NEW EXAMPLE ARE AXIS AND
    %PAUSE(maybe)%
    %-------------------------------------------------------------------
    if mod(ct,printDump) == 0
       %fprintf('timestep: %d | time: %.4f \n', ct, t);
        
       % Calculate Divergence of Velocity Field ( divU = du/dx + dv/dy )
        div_u = ( uP(yIds,xIds+1) - uP(yIds,xIds-1) ) / (2*dx) + ...
                ( vP(yIds+1,xIds) - vP(yIds-1,xIds) ) / (2*dy);
        
        % Calculate Vorticity ( omega = dv/dx - du/dy )
        dvdx = ( vP(yIds,xIds+1) - vP(yIds,xIds-1) ) / (2*dx);            
        dudy = ( uP(yIds+1,xIds) - uP(yIds-1,xIds) ) / (2*dy);
        vorticity = dvdx - dudy;                
        
        % Print useful information to screen     
        fprintf('\n time: %.4f (timestep: %d)\n',ct*dt,ct)
        fprintf('    -->  CFL(adv)=%.8f\n', dt/dx*( max(max(abs(uP))) + max(max(abs(vP)))  )  );
        fprintf('    -->  CFL(diff)=%.8f\n', 2*mu/rho*dt*(1/dx^2+1/dy^2)*( max(max(abs(uP))) + max(max(abs(vP)))  )  );
        fprintf('    -->  max(|u|)=%.8f\n',max( max( abs(uP) ) ) );
        fprintf('    -->  max(|v|)=%.8f\n',max( max( abs(vP) ) ) );
        fprintf('    -->  max(|p|)=%.8f\n',max( max( abs(Pressure) ) ) );
        fprintf('    -->  max(|w|)=%.8f\n',max( max( abs(vorticity) ) ) );
        fprintf('    -->    div(u)=%.10f\n',mean(mean(div_u)) );        
        
        dvdx = ( vN(2:Ny, 3:Nx+1) - vN(2:Ny, 1:Nx-1) ) / (2*dx);
        dudy = ( uN(3:Ny+1, 2:Nx) - uN(1:Ny-1, 2:Nx) ) / (2*dy);
        xVort = X(2:Ny, 2:Nx);
        yVort = Y(2:Ny, 2:Nx);
        vorticity = dvdx - dudy;

        f1=figure(1);
        colormap("parula")                  %setting type of color map
        s = surf(xVort,yVort, vorticity);
        s.EdgeColor = "none";               %gets rid of black coloring grid lines
        h = colorbar;                       %adds color bar to the side
        maxVort = 25;
        caxis([-maxVort maxVort])
        set(gca,'fontsize', 23)             %font size
        xlabel('x') 
        ylabel('y')
        % zlim([0,7]) 
        view(2)                             %changes view of graph 1st number is angle from x to y, second number is angle from top to bottom
        pause(0.01)
        f1.Position = [250 250 520 400];
        hold off

        ctsave = ctsave + 1;
        save_vtk_files(dx, dy, t, ctsave, uP, 'u')
        save_vtk_files(dx, dy, t, ctsave, vP, 'v')
        save_vtk_files(dx, dy, t, ctsave, vorticity, 'vort')
        save_vtk_files(dx, dy, t, ctsave, Pressure, 'P')
        savevtk_vector_field(dx, dy, t, ctsave, uP, vP, 'uVec')
    end
end