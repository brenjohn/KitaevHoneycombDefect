% This function uses the following arguments to draw a graph of the Kitaev model.
% Del_x, Del_y, Del_z are matrices which encode the couplings between sites of the Kitaev model.
% Nx, Ny are the number of sites in the x and y directions respectively.
function drawLattice(Del_x,Del_y,Del_z,Nx,Ny)

% create meshgrid of vertices
[x,y] = meshgrid(1:Nx,1:Ny);

% draw the vertices of the graph
x = reshape(x',Nx*Ny,1);
y = reshape(y',Nx*Ny,1);
scatter(x,y,'black','filled')
xlim([0 Nx+1])
ylim([0 Ny+1])
hold

% loop through each vertex of the graph
for iy = 1:Ny
    for ix = 1:Nx
        
        % for each site, find the index of the site it connects to
        i = find(Del_x(index(ix,iy,Nx,Ny),:));
        j = find(Del_y(index(ix,iy,Nx,Ny),:));

        
        % draw connecting x-link, black if postitive coupling, red otherwise
        if sum(Del_x(index(ix,iy,Nx,Ny),:))>=0
            colour = 'black';
        else
            colour = 'red';
        end
        if ix < Nx
            plot([ix,x(i)],[iy,y(i)],colour)
        else
            plot([ix,ix+1],[iy,iy],colour)
        end
        
        
        % draw connecting y-link, black if postitive coupling, red otherwise
        if sum(Del_y(index(ix,iy,Nx,Ny),:))>=0
            colour = 'black';
        else
            colour = 'red';
        end
        if iy < Ny
            plot([ix,x(j)],[iy,y(j)],colour)
        else
            plot([ix,ix],[iy,iy+1],colour)
        end
        
        
        % if a site has been deleted, colour it red
        if Del_z(index(ix,iy,Nx,Ny),index(ix,iy,Nx,Ny)) == 0
            scatter(ix,iy,'red','filled')
        end
        
    end
end
