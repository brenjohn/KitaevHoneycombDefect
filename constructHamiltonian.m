function H = constructHamiltonian(Jx,Jy,Jz,k,Nx,Ny,delta_x,delta_y,delta_z,a,b,Ld,P)

% create the BdG blocks of the Hamiltonian
Xi = 2*Jz*delta_z + Jx*( delta_x+delta_x' ) + Jy*(delta_y+delta_y');
Del = Jx*(delta_x-delta_x') + Jy*(delta_y-delta_y');
    
% The three body terms at each end of the defect appear twice, once for
% the regular plaquette beside the defect end point and once again for
% the endpoint itself
if floor(Ld)>=0
    delta_y(index(a+floor(Ld)+1,b+floor(Ld)+1,Nx,Ny),:)=2*delta_y(index(a+floor(Ld)+1,b+floor(Ld)+1,Nx,Ny),:);
    delta_y(index(a,b,Nx,Ny),:)=2*delta_y(index(a,b,Nx,Ny),:);
    if (Ld-floor(Ld))>0
        delta_y(index(a+floor(Ld)+1,b+floor(Ld)+2,Nx,Ny),:)=2*delta_y(index(a+floor(Ld)+1,b+floor(Ld)+2,Nx,Ny),:);
    end
end
    
% create the three body magnetic terms
del = -1i*(delta_x-delta_x') + 1i*(delta_y-delta_y') - 1i*(delta_x-delta_x') + 1i*(delta_y-delta_y');
    
xi  =  - 1i*(delta_x*delta_y'- (delta_x*delta_y')' ) + 1i*( delta_y'*delta_x - (delta_y'*delta_x)' );
del = del + 1i*(delta_x*delta_y'- (delta_x*delta_y')' ) + 1i*( delta_y'*delta_x - (delta_y'*delta_x)' );
    
% include the magnetic terms in the BdG blocks
Xi = Xi + k.*xi;
Del = Del + k.*del;
    
% project out redundant degrees of freedom
Xi = P'*Xi*P;
Del = P'*Del*P;
    
% assemble the Hamiltonian
H = 0.5*[ Xi,Del ; Del',-(Xi).' ];

end