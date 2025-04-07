### Function to plot travelling wave solution

def plot_TW(n,t):
    global Nx, Ny, dx, dy, x, y, Xm
    fig, (ax1,ax2) = plt.subplots(2,1,figsize=(11, 4))
    im = ax1.pcolormesh(x,y,n.T,shading='auto') 
    fig.colorbar(im,ax=ax1,orientation='vertical',label='n(x,y)')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_title(f'$n(t,x,y)$ at time t = {t}')
    rho = dy*np.sum(n,axis=1);
    ax2.plot(x, rho, 'b-', linewidth=2)
    ax2.set_xlabel('x')
    ax1.set_xlim([0, Xm])
    ax2.set_title(f'$\\rho(t,x)$ at time t = {t}')
    fig.subplots_adjust(hspace=1) 
    plt.draw()
    plt.pause(0.1)