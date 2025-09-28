import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import matplotlib.cm as cm
from sigfig import round

class parameters():
    def __init__(self):
        # Mesh Information
        # domain from from z=a to z=b
        self.Q = 1e-9 # m^3/s
        self.nu = 1.787e-6 # m^2/s
        self.G = 50e-3 # W/m^2
        self.ubtaub = 0.08 # W/m^2, np.mean(heat)
        self.rhoi = 917 # kg/m^3
        self.rhow = 1000 # kg/m^3
        self.mu = self.rhow*self.nu # Pa*s
        self.Latent = 334000 # J/kg
        self.A = 3.5e-25 # ice softness
        self.omega = 1/1000 # transition to Reynolds number
        self.N0 = 200*1000 # Pa
        self.nglen = 3 # glen's law exponent
        self.ri = 0.1 # m
        self.ro = 10000 # m
        self.rt = (self.omega*(self.Q))/(2*np.pi*self.nu) # m
        # gap opening
        self.b_scale = (self.G + self.ubtaub)/((self.rhoi*self.Latent*self.A)*(self.N0**3)) # m
        
class nondimensional():
    def __init__(self):
        prms = parameters()
        # scale of inner radius
        self.inner_radius = prms.ri/prms.ro
        
        # transition to turbulence scale
        self.R = prms.rt/prms.ro

        # melt scale
        self.M = (12*prms.mu*(prms.Q))/(2*np.pi*(prms.b_scale**3)*prms.N0)

        # dissipation (new)
        self.D = ((prms.Q)*prms.N0)/(2*np.pi*(prms.G+prms.ubtaub)*(prms.ro**2))
        
        # flux parameter
        self.W = (2*np.pi*(prms.G+prms.ubtaub)*(prms.ro**2))/((prms.Q)*prms.rhow*prms.Latent)
        
def Nfun(r,N,nondims):
    a = (nondims.D/r)
    b = (nondims.M/r)*(1+(nondims.R/r))*(N**9)
    find_deriv = lambda x : ((1+(a*x))**3)*x - b
    jac = lambda x : ((1+a*x)**3) + 3*((1+a*x)**2)*a*x
    if r<nondims.R:
        approx = (((nondims.R*nondims.M*abs(r))/(nondims.D**3))**(1/4))*(N**(9/4))
    else:
        approx = (nondims.M/r)*(1+(nondims.R/r))*(N**9)
        
    dNdr = -fsolve(find_deriv,approx,fprime=jac)
    return dNdr

def Nfun_flux(r,N,q,nondims):
    a = nondims.D*q
    b = nondims.M*q*(1+nondims.R*q)*(N**9)
    find_deriv = lambda x : ((1+(a*x))**3)*x - b
    jac = lambda x : ((1+a*x)**3) + 3*((1+a*x)**2)*a*x
    if r<nondims.R:
        approx = (((nondims.R*nondims.M*abs(r))/(nondims.D**3))**(1/4))*(N**(9/4))
    else:
        approx = (nondims.M/r)*(1+(nondims.R/r))*(N**9)
        
    dNdr = -fsolve(find_deriv,approx,fprime=jac)
    return dNdr
        
def rhsfun(r,y,parameters,nondims):
    N = y[0]
    F = y[1]
    q = F/r

    dNdr = Nfun_flux(r,N,q,nondims)
    mdot = 1 - nondims.D*q*dNdr[0]
    dFdr = -r*nondims.W*mdot
    return [dNdr[0],dFdr]

def shootfun(Ninner,parameters,nondims):
    sol = solve_ivp(rhsfun,[nondims.inner_radius,1],[Ninner,1],args=(parameters,nondims),max_step=1e-2)
    res = sol.y[0,-1]-1
    return res

def shootfun_otrolado(qouter,parameters,nondims):
    sol = solve_ivp(rhsfun,[1,nondims.inner_radius],[1,qouter],args=(parameters,nondims),max_step=1e-3)
    res = sol.y[1,-1]-1
    return res

if __name__ == "__main__":
    prms = parameters()
    nd = nondimensional()
    
    Natri = 1.0985 # nondimensional effective pressure
    Fi = 1
    
    Nival = fsolve(shootfun,Natri,args=(prms,nd))
    print(Nival)
    sol = solve_ivp(rhsfun,[nd.inner_radius,1],[Nival,Fi],args=(prms,nd),max_step=1e-2,t_eval=np.logspace(np.log10(nd.inner_radius+1e-5),0,128))
    print(sol.y[0,-1]-1)
    
    cmap = cm.get_cmap('viridis', 5)
    
    solfa = solve_ivp(Nfun,[1,nd.inner_radius],[1],args=(nd,),max_step=1e-4,t_eval=np.logspace(0,np.log10(nd.inner_radius+1e-5),32),method='RK45') # 'RK45'
    
    f1 = plt.figure()
    plt.semilogx(sol.t*prms.ro,(sol.y[0,:]-1)*prms.N0/1000,'.-',color=cmap(0),label='vary $q$',linewidth=1) 
    plt.plot(solfa.t*prms.ro,(solfa.y[0,:]-1)*prms.N0/1000,'^-',color=cmap(2),label='fixed $q$',linewidth=1) #^
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.yticks(fontsize=14); plt.xticks(fontsize=14)
    plt.xlabel('nondimensional radial distance, $r$ (m)',fontsize=16)
    plt.ylabel('scaled effective pressure, $N-N_0$ (kPa)',fontsize=16)
    plt.legend(fontsize=14)
    plt.text(0.16,0,'$Q='+str(prms.Q)+'$ m$^3$/s, '+'$N_0='+str(prms.N0/1000)+'$ kPa',fontsize=14)
    f1.savefig("effp_varyq.pdf",bbox_inches='tight')

    f2 = plt.figure()
    plt.loglog(sol.t,abs(sol.y[1,:])/sol.t,'.-',color=cmap(0),label='vary $q$')
    plt.loglog(sol.t,1/(sol.t),'-',color=cmap(2),label='fixed\n $q = 1/r$',linewidth=2)
    plt.loglog(sol.t,0.5*nd.W*sol.t)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    # plt.ylim(1e-2,1e5)
    plt.yticks(fontsize=14); plt.xticks(fontsize=14)
    plt.xlabel('nondimensional radial distance, $r/r_d$',fontsize=18)
    plt.ylabel('nondimensional flux, $2\pi r_d|q|/Q$',fontsize=18)
    plt.legend(fontsize=14)
    plt.text(0.18e-4,0.02,'$M='+str(round(nd.M,sigfigs=2))+'$, '+'$D='+str(round(nd.D,sigfigs=2))+'$',fontsize=14)
    f2.savefig("q_vary.pdf",bbox_inches='tight')
