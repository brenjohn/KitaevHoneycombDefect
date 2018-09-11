% This function will flip the appropriate signs in the delta matrices, used
% to create the hamiltonian, to encode the vortex locations listed in V and
% to encode the homology sector labeled by lx and ly

function [delta_x,delta_y] = printVortexConfig(delta_x,delta_y,Nx,Ny,lx,ly,V)

% flip the sign of the links appropriately to encode the the values of
% the X and Y operators, according to the vortex configuration V
N = size(V);
for i=1:N(1)

    % flip signs to create branch cut above each vortex
    for j= floor(V(i,2))+1:Ny
        delta_x(index(floor(V(i,1)),j,Nx,Ny),:) = -1*delta_x(index(floor(V(i,1)),j,Nx,Ny),:);
    end

    % flip signs to connect branch cuts at top of the lattice
    for j= floor(V(i,1))+1:Nx
        delta_y(index(j,Ny,Nx,Ny),:) = -1*delta_y(index(j,Ny,Nx,Ny),:);
    end

    if (V(i,2)-floor(V(i,2)))>0
        delta_x(index(floor(V(i,1)),floor(V(i,2))+1,Nx,Ny),:) = -2*(V(i,2)-floor(V(i,2))-0.5)*delta_x(index(floor(V(i,1)),floor(V(i,2))+1,Nx,Ny),:);
    end

    if (V(i,1)-floor(V(i,1)))>0
        delta_y(index(floor(V(i,1))+1,Ny,Nx,Ny),:) = -2*(V(i,1)-floor(V(i,1))-0.5)*delta_y(index(floor(V(i,1))+1,Ny,Nx,Ny),:);
        for j= floor(V(i,2))+1:Ny
            delta_x(index(floor(V(i,1)+1),j,Nx,Ny),:) = -2*(V(i,1)-floor(V(i,1))-0.5)*delta_x(index(floor(V(i,1)+1),j,Nx,Ny),:);
        end
    end
end

% flip signs to impliment boundary conditions determined by loop
% operators
for i=1:Ny
    delta_x(index(Nx,i,Nx,Ny),:) = (-1)^Nx*lx*delta_x(index(Nx,i,Nx,Ny),:); %- %%(-1)^Nx*
end

for i=1:Nx
    delta_y(index(i,Ny,Nx,Ny),index(i,1,Nx,Ny)) = ly*delta_y(index(i,Ny,Nx,Ny),index(i,1,Nx,Ny)); %- %%
end


end