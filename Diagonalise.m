% This script calculates the two eigenstates of the effective fermionic
% Hamiltonian of the Kitaev Honeycomb model for a specified 
% vortex/homology sector that are nearest to zero energy and then plots
% them

clear;

% set the parameters of the system
Nx = 30; Ny = 30;
Jx = 1; Jy = 1; Jz = 1;
lx = -1; ly = -1; k = 0.1;

% set the location of the defect (a,b) and length Ld
a = 10; b = 10; Ld = 10;

% Calculate the total number of site on the lattice
if Ld > 0 
    N = Nx*Ny-Ld;
else
    N = Nx*Ny;
end

% set vortex configuration
V = [];
%V = [[10,20]',[30,20]']';

% construct Hamiltonian building blocks and the Hamiltonian
[delta_x,delta_y,Del_z,P] = HamiltonianDelta(Nx,Ny,a,b,Ld);
[delta_x,delta_y] = printVortexConfig(delta_x,delta_y,Nx,Ny,lx,ly,V);
H = constructHamiltonian(Jx,Jy,Jz,k,Nx,Ny,delta_x,delta_y,Del_z,a,b,Ld,P);

% diagonalise the Hamiltonain
[U,EH] = eig(full(H));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting the first mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% replace the degrees of freedom that were projected out when the defect
% was introduced with zeros
NN = size(P);
Z = spalloc(NN(1),NN(2),1);
W = U(:,N);
W = [[P',Z']',[Z',P']']*W;

% reshape array to match dimensions of the model's lattice
Psi = reshape(W(1:Nx*Ny),[Nx,Ny]);

% set colors to be used to represent complex phase
Ac = [0 0 1]; Bc = [0 1 0]; Cc = [1 0 0]; Dc = [1 1 0];

Psi = Psi.'; % get correct orientation of matrix

%allocate memory to color map
NN = size(Psi);
C = zeros(NN(1),NN(2),3);


% Fill color map with colors determined by complex phase
for x = 1:NN(1)
    for y = 1:NN(2)
        u = cos(angle(Psi(x,y))); v = sin(angle(Psi(x,y)));
        if u >= 0
            c = u*Ac;
        else
            c = -u*Cc;
        end
        if v >= 0
            c = c + v*Bc;
        else
            c = c - v*Dc;
        end
        C(x,y,:) = c;
    end
end

% Plot the wave function
figure(1);
surface(abs(Psi),C);
shading interp;
view(20,35)
xlabel('X')
ylabel('Y')
zlim([0 0.3])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting the second mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% replace the degrees of freedom that were projected out when the defect
% was introduced with zeros
NN = size(P);
Z = spalloc(NN(1),NN(2),1);
W = U(:,N+1);
W = [[P',Z']',[Z',P']']*W;

% reshape array to match dimensions of the model's lattice
Psi = reshape(W(1:Nx*Ny),[Nx,Ny]);

% set colors to be used to represent complex phase
Ac = [0 0 1]; Bc = [0 1 0]; Cc = [1 0 0]; Dc = [1 1 0];

Psi = Psi.'; % get correct orientation of matrix

%allocate memory to color map
NN = size(Psi);
C = zeros(NN(1),NN(2),3);


% Fill color map with colors determined by complex phase
for x = 1:NN(1)
    for y = 1:NN(2)
        u = cos(angle(Psi(x,y))); v = sin(angle(Psi(x,y)));
        if u >= 0
            c = u*Ac;
        else
            c = -u*Cc;
        end
        if v >= 0
            c = c + v*Bc;
        else
            c = c - v*Dc;
        end
        C(x,y,:) = c;
    end
end

% Plot the wave function
figure(2);
surface(abs(Psi),C);
shading interp;
view(20,35)
xlabel('X')
ylabel('Y')
zlim([0 0.3])