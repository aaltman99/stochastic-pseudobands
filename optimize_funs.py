import numpy as np
from scipy.optimize import minimize_scalar, fsolve, Bounds
import matplotlib.pyplot as plt


### optimization of pseudobands critical parameter β
### we use ansatz for widths of width_j = α * exp(j * β)
### the sum constraint \sum_j width_j = Emax - E0 fixes 
### α in terms of β. Thus, we optimize only over β.

### Using this ansatz allows us to separate dependencies, i.e.
### the ith term in the loss function sum depends only on i,
### and not all j <= i as would be without this ansatz.





C = 1 / (12 * np.pi**8) # atomic units
B = np.sqrt(2) / np.pi**2


def alpha(beta, E0, Emax, nspbps, nslice): 
    dE = Emax - E0
    return (dE * (np.exp(beta)-1)) / (np.exp(beta) * (np.exp(nslice * beta) - 1))

def w(beta, j, E0, Emax, nspbps, nslice): 
    return alpha(beta, E0, Emax, nspbps, nslice) * np.exp(j * beta)

def Ebar(beta, j, E0, Emax, nspbps, nslice): 
    a = alpha(beta, E0, Emax, nspbps, nslice)
    width = w(beta, j, E0, Emax, nspbps, nslice)
    return E0 + a * np.exp(beta) * (np.exp(j * beta) - 1) / (np.exp(beta) - 1) - width / 2


def Loss(beta, E0=1, Emax=10, nspbps=1, nslice=10):
    out = 0
    for j in range(1, nslice+1):
        Ei = Ebar(beta, j, E0, Emax, nspbps, nslice)
        wi = w(beta, j, E0, Emax, nspbps, nslice)

        assert Ei > 0, Ei
        assert wi > 0, wi
        
        # FEG approximation to dimension of the subspace breaks down, let dim = 1 in this case
        if B * wi * np.sqrt(Ei) < 1:
            out += C * wi**2 / Ei**2
            continue
            
        l = C * wi**2 / Ei**2 + 1 / (Ei**2 * nspbps) * (1 - 1 / (B * wi * np.sqrt(Ei)))
        out += l
        
    return out



def optimize(E0=1, Emax=10, nspbps=1, nslice=10):
    
    dE = Emax - E0
    
    min_dim = fsolve(lambda x: B * x * np.sqrt(E0 + x/2) - 1, 4)
    print(min_dim)
    
    low = max(fsolve(lambda x: alpha(x, E0, Emax, nspbps, nslice)*np.exp(x) - max(min_dim), .2*np.ones(3), maxfev=1000))
    high = max(fsolve(lambda x: B * w(x, nslice, E0, Emax, nspbps, nslice) * np.sqrt(Ebar(x, nslice, E0, Emax, nspbps, nslice)) - 1,  .2*np.ones(3), maxfev=1000))
    print(low)
    print(high)
   
    tol=1e-9
    # bnds = (max(tol,low), min(1-tol,high))
    
    # low and high are based on the FEG appx to the subspace dimension, which is bad
    bnds = (tol, 1-tol)
    print('bounds: ', bnds)
    
    
    # X = np.linspace(bnds[0], bnds[1],500)
    # toplot=[Loss(x, E0, Emax, nspbps, nslice) for x in X]
    # plt.scatter(X, toplot)
    # plt.savefig('loss.png')
    # print(np.argmin(toplot), min(toplot))
    # plt.scatter(X[100:480], np.gradient(toplot)[100:480])
    # plt.savefig('lossdiff.png')
             
             
    result = minimize_scalar(Loss, args=(E0, Emax, nspbps, nslice), bounds=bnds, method='bounded')
    print(f'minimization results: \n{result}')
    print(f'alpha: {alpha(result.x, E0, Emax, nspbps, nslice)}')
    return result
    
    
