function i = index( x,y,Nx,Ny )

% torus boundary conditions
if x>Nx
    x = x-Nx;
end
    
if y>Ny
    y = y-Ny;
end
    
i = x + Nx*(y-1);

end