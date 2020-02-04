import numpy as np

def f(x):
    return 1.0/(x**2+1.0);

def df(x):
    return -2.0*x*f(x)**2;

def ddf(x):
    return 8.0*x**2*f(x)**3-2.0*f(x)**2;

def h(x):
    return x*f(x);

def dh(x):
    return f(x) + x*df(x);

def ddh(x):
    return 2*df(x)+x*ddf(x);

def dp(y,p):
    dp = p*dh(y)
    y_ki = y[:,np.newaxis]-y
    dh_y_ki = dh(y_ki)
    dp += p*np.sum(dh_y_ki,axis=0) - p*dh(0) # pk -p*dh(0) ??
    dp += - np.sum(p*dh_y_ki,axis=1) + p*dh(0)
    return dp;

def dy(y,u):
    y_ki = y[:,np.newaxis]-y
    dy = np.sum(h(y_ki)) - h(y) - u
    return dy;

def u_sing(y,p):
    y_ki = y[:,np.newaxis]-y
    
    Dy = dy(y,0) # u later cancels out
    Dy_ki = Dy[:,np.newaxis]-Dy
    
    Dp = dp(y,p)
    
    a = Dp*( dh(y)+np.sum(dh(y_ki),axis=0)-dh(0) )
    
    a += p*ddh(y)*( np.sum(h(y_ki),axis=0)-h(y) )
    
    a += p*np.sum(Dy_ki*ddh(y_ki),axis=0)
    
    a += - np.sum(Dp*dh(y_ki),axis=0) + Dp*dh(0)
    
    a += - np.sum(p*ddh(y_ki)*Dy_ki,axis=0)
    
    u = np.sum(a)
    
    u /= np.sum(p*ddh(y),axis=0)
    
    return u;

def u_sing2(y,p):
    
    Dy = dy(y,0)
    
    Dp = dp(y,p)
    
    u = np.sum(p*ddh(y)*Dy, axis=0)
     
    u += np.sum(Dp*(dh(y)+dh(0)), axis=0)
    
    u /= np.sum(p*ddh(y), axis=0)
    
    return u;
































