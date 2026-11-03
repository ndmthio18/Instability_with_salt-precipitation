"""Descriptive growth and spatial evidence helpers, not physical onset labels."""
import numpy as np

def neumann_spectrum(xs,zs,anomaly):
    xs=np.asarray(xs,float);zs=np.asarray(zs,float);d=np.asarray(anomaly,float)
    if d.shape!=(len(zs),len(xs)):raise ValueError('field/grid shape mismatch')
    wx=np.diff(np.r_[xs[0],(xs[1:]+xs[:-1])/2,xs[-1]])
    wz=np.diff(np.r_[zs[0],(zs[1:]+zs[:-1])/2,zs[-1]])
    cm=np.cos(np.pi*np.arange(1,len(xs))[:,None]*(xs-xs[0])/(xs[-1]-xs[0]))
    coeff=2*(d*wx)@cm.T/(xs[-1]-xs[0]);energy=.5*coeff**2;energy[:,-1]*=.5
    return coeff,np.sum(energy*wz[:,None],axis=0)

def log_growth_fit(t,amplitude,start,end):
    t=np.asarray(t,float);a=np.asarray(amplitude,float)
    m=(t>=start-1e-7)&(t<=end+1e-7)&(a>0)&np.isfinite(a)
    x=t[m];y=np.log(a[m])
    if len(x)<3:return dict(status='insufficient_samples',slope_per_s=None)
    slope,intercept=np.polyfit(x-x[0],y,1);res=y-(intercept+slope*(x-x[0]));ss=sum((y-y.mean())**2)
    return dict(status='complete_window' if x[-1]>=end-1e-7 else 'partial_window',samples=len(x),first_actual_s=float(x[0]),last_actual_s=float(x[-1]),slope_per_s=float(slope),e_folding_s=float(1/slope) if slope>0 else None,R_squared=float(1-sum(res*res)/ss) if ss>0 else None,interpretation='descriptive log-amplitude slope; not LSA eigenvalue')

def energy_depth(z,sigma,q):
    z=np.asarray(z,float);s=np.asarray(sigma,float)
    if not 0<q<=1 or np.any(np.diff(z)<=0) or len(z)<2:raise ValueError('invalid depth grid or quantile')
    widths=np.diff(np.r_[z[0],(z[1:]+z[:-1])/2,z[-1]])
    e=s*s*widths
    if not np.isfinite(e).all() or sum(e)<=0:return None
    i=min(np.searchsorted(np.cumsum(e[::-1])/sum(e),q),len(z)-1)
    return float((z[-1]-z[::-1][i])*1000)
