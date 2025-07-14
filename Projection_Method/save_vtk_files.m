
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: saves simulation data to .vtk format in folder vizData/ 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_vtk_files(dx, dy, current_time, ctsave, uMAT, strName)

%----------------------------------------------
% Create folder to save data
%----------------------------------------------
if ctsave == 0
    mkdir('vizData');
end

%----------------------------------------------
% Print all Eulerian Data vizData folder
%----------------------------------------------
cd('vizData'); %Go into vizData directory
    

    %----------------------------------------------
    % Need transpose of matrix: 
    %       original: uMAT(j,i) <<-- i=x, j=y
    %          wants: "uMAT(x,y)"
    %----------------------------------------------
    uMATName = [strName '.' num2str(ctsave) '.vtk'];
    savevtk_scalar(uMAT', uMATName, strName, dx, dy, current_time);

cd ../;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints scalar matrix to vtk formated file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function savevtk_scalar(array, filename, colorMap,dx,dy, time)
%  savevtk Save a 3-D scalar array in VTK format.
%  savevtk(array, filename) saves a 3-D array of any size to
%  filename in VTK format.
    [ny, nx, nz] = size(array);
    fid = fopen(filename, 'wt');
    fprintf(fid, '# vtk DataFile Version 2.0\n');
    fprintf(fid, 'Comment goes here\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, '\n');
    fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
    %
    fprintf(fid, 'FIELD FieldData 1\n');
    fprintf(fid, 'TIME 1 1 double\n');
    fprintf(fid, '%.8f\n',time);
    %
    fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', ny, nx, nz);
    fprintf(fid, '\n');
    fprintf(fid, 'ORIGIN    0.000   0.000   0.000\n');
    %fprintf(fid, 'SPACING   1.000   1.000   1.000\n'); if want [1,32]x[1,32] rather than [0,Lx]x[0,Ly]
    fprintf(fid, ['SPACING   ' num2str(dx) ' '   num2str(dy) '   1.000\n']);
    fprintf(fid, '\n');
    fprintf(fid, 'POINT_DATA   %d\n', nx*ny*nz);
    fprintf(fid, ['SCALARS ' colorMap ' double\n']);
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '\n');
    for a=1:nz
        for b=1:nx
            for c=1:ny
                fprintf(fid, '%d ', array(c,b,a));
            end
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
return



