% This function creates the delta matrices used to create a hamiltonian for
% a Nx by Ny lattice with a defect of length Ld at (a,b)

function [delta_x,delta_y,delta_z,P] =  HamiltonianDelta(Nx,Ny,a,b,Ld)
    
    % preallocate memory for possible links
    delta_x = spalloc(Nx*Ny,Nx*Ny,Nx*Ny);
    delta_y = spalloc(Nx*Ny,Nx*Ny,Nx*Ny);
    delta_z = speye(Nx*Ny,Nx*Ny);
    
    % preallocate memory to record index of redundant degrees of freedom
    redDeg  = zeros(1,floor(Ld));
    
    % turn on the x and y links
    for y=1:Ny
        for x=1:Nx
            delta_x(index(x,y,Nx,Ny),index(x+1,y,Nx,Ny))=1;
            delta_y(index(x,y,Nx,Ny),index(x,y+1,Nx,Ny))=1;
        end
    end
    
    
    
    
    % Recouple to create diagonal defect of length Ld beginning at (a,b)
    for i=0:floor(Ld)
        
        % turn off old x-links along defect
        delta_x(index(a+i,b+i,Nx,Ny),index(a+i+1,b+i,Nx,Ny))=0;
        delta_x(index(a+i,b+i+1,Nx,Ny),index(a+i+1,b+i+1,Nx,Ny))=0;
        
        % turn on new x-links along defect
        delta_x(index(a+i,b+i+1,Nx,Ny),index(a+i+1,b+i,Nx,Ny))=1;
    end
    
    for i=1:floor(Ld)
        % turn off old y-links along defect
        delta_y(index(a+i,b+i,Nx,Ny),index(a+i,b+i+1,Nx,Ny))=0;
        delta_y(index(a+i,b-1+i,Nx,Ny),index(a+i,b+i,Nx,Ny))=0;
        delta_z(index(a+i,b+i,Nx,Ny),index(a+i,b+i,Nx,Ny))=0;
        
        % turn on new y-links along defect
        delta_y(index(a+i,b+i-1,Nx,Ny),index(a+i,b+i+1,Nx,Ny))=1;
        
        % record index of redundant degree of freedom
        redDeg(i) = index(a+i,b+i,Nx,Ny);
    end
    
    if Ld-floor(Ld)>0
        i = Ld-floor(Ld);

        delta_x(index(a+floor(Ld),b+floor(Ld),Nx,Ny),index(a+floor(Ld)+1,b+floor(Ld),Nx,Ny))=1-i;
        delta_x(index(a+floor(Ld),b+1+floor(Ld),Nx,Ny),index(a+floor(Ld)+1,b+floor(Ld)+1,Nx,Ny))=1-i;
      
        delta_x(index(a+floor(Ld),b+floor(Ld)+1,Nx,Ny),index(a+floor(Ld)+1,b+floor(Ld),Nx,Ny))=i;
        
        delta_y(index(a+floor(Ld),b+floor(Ld),Nx,Ny),index(a+floor(Ld),b+floor(Ld)+1,Nx,Ny))=1-i;
        delta_y(index(a+floor(Ld),b-1+floor(Ld),Nx,Ny),index(a+floor(Ld),b+floor(Ld),Nx,Ny))=1-i;
        delta_z(index(a+floor(Ld),b+floor(Ld),Nx,Ny),index(a+floor(Ld),b+floor(Ld),Nx,Ny))=1-i;
        
        delta_y(index(a+floor(Ld),b+floor(Ld)-1,Nx,Ny),index(a+floor(Ld),b+floor(Ld)+1,Nx,Ny))=i;
    end
    
    % project out the redundent degrees of freedom created by the defect
    P = speye(Nx*Ny,Nx*Ny);
    P(:,redDeg) = [];   
end