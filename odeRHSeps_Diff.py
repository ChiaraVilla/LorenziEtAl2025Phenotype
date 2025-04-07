######################################################################################
### ODE approximation of RHS of PDE upon spatial discretisation (finite volume scheme)
###
### \dt n = \dx [ D(y)*\dx n - A(y,\rho,S)*n ] + \bar{D} \d2yy n + n*R(y,\rho,S)
### \rho(t,x) = \int_0^1 n(t,x,y) dy
###
### Diffusion-driven motion: 
### D(y) = Diffx * y ;    A(y,\rho,S) = 0 ;    R(y,\rho,S) = r*(1-y-\rho)
### 
######################################################################################

def odeRHSeps(t,sol):
    global Diffx, Diffy, r, k, Nx, Ny, dx, dy, x, y
    
    # Reshape sol into n(x,y), Integrate n to obtain rho
    n = np.reshape(sol,(Nx,Ny));
    rho = dy*np.sum(n,axis=1);
    
    # Phenotype-dependent diffusion coefficient D(y) = Diffx*y
    Dxy = np.tile(Diffx*y,(Nx,1));
    Dxy_mid = 0.5*(Dxy[0:-1,:]+Dxy[1:,:]); # Require approximation at grid cell interfaces in x
    
    # Flux in physical (flx) and phenotype (fly) space, at grid cell interfaces
    flx = np.zeros((Nx+1,Ny)); # This gives also zeroflux BCs at x=0 and x=Xm
    flx[1:-1,:] = Dxy_mid*(n[1:,:]-n[0:-1,:]);
    fly = np.zeros((Nx,Ny+1)); # This gives also zeroflux BCs at y=0 and y=1
    fly[:,1:-1] = Diffy*(n[:,1:]-n[:,0:-1]);
    
    # Flux contributions from spatial movement and phenotypic changes
    dflux = (flx[1:,:]-flx[:-1,:])/dx + (fly[:,1:]-fly[:,:-1])/dy;
    print(dflux.shape)
    
    # Cell proliferation and death, here F=n*R with R=r*(1-y-rho)
    kinetics = r*n*( np.ones((Nx,Ny))-np.tile(y,(Nx,1))-np.tile(rho,(Ny,1)).T );
    print(kinetics.shape)
    
    # Store RHS of ODE system obtained via the MOL in dndt
    dndt = dflux + kinetics;
    
    # Return dndt as vector
    return np.reshape(dndt,(Nx*Ny,)); 