import numpy as np

import yfinance as yf

start = "2014-01-01"
end = '2019-1-01'
aapl = yf.download('AAPL',start,end)

returns_aapl = list(aapl["Open"].pct_change()[1:])

def C(m):
    return 2*np.log(m)**(1/2)-(np.log(np.pi)+np.log(np.log(m)))/(2*(2*np.log(m))**(1/2))

def S(m):
    return 1/(2*np.log(m)**(1/2))

def IV(t,list_returns,n):
    coeff=(np.pi/2)*(n/(n-1))*(1/(n-1))
    res=0
    for i in range(n-2):
        res=res+np.abs(list_returns[t-i])*np.abs(list_returns[t-i-1])
    return res*coeff

def presence_jumps(list_returns,Am,m,n,eps):
    value=[]
    Jump=[]
    for i in range(Am):
        value.append((np.abs(list_returns[i]/IV(i,list_returns,n))-C(m))/S(m))
        Max=max(value)
        t=value.index(Max)
        psi=np.exp(np.exp(-t))
        if psi+eps<Max<=psi+eps:
            Jump.append(True)
        else :
            Jump.append(False)
    return Jump