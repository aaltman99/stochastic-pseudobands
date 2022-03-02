import numpy as np
import matplotlib.pyplot as plt
from optimize_funs import optimize, alpha, Loss



if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--E0', type=float, required=True)
    parser.add_argument('--Emax', type=float, required=True)
    parser.add_argument('--nspbps', type=int, required=True)
    parser.add_argument('--nslice', type=int, required=True)
    
    args = parser.parse_args()
    
    
    res = optimize(**vars(args))
    
    print(res)
    print(alpha(res.x, **vars(args)))
    print(Loss(res.x, **vars(args)))
    
    x = np.arange(1,vars(args)['nslice']+1)
    plt.plot(x, alpha(res.x, **vars(args))*np.exp(res.x*x))
    plt.xlim(0,x[-1]+1)
    plt.ylim(0, max(alpha(res.x, **vars(args))*np.exp(res.x*x))*1.2)
    plt.savefig('widths.png')
