import math
from scipy.optimize import fsolve


### Numerical Compute
class thresholdCompute(object):
    
    def __init__(self, p,G):
        self.p = p
        self.G = G
        self.eta = 2*p - 4/3*pow(p,2)
        
    def f(self, x):
        return x* math.log(x /0.75) + (1-x)*math.log((1-x)/0.25)  -2*(x*math.log(x/self.eta)  + (1-x)*math.log((1-x)/(1-self.eta))) 

    def findRoot(self):
        alpha= fsolve(self.f, 0.4)
        I = alpha * math.log(alpha / 0.75) + (1-alpha) * math.log((1-alpha)/0.25)
        liid = 2/I * math.log(self.G) 
        
        self.liidInt = (int) (math.ceil(liid))
        self.threshold = (int) (math.ceil(alpha*liid))
        return self.liidInt,self.threshold 
        

class Ncompute(object):
    def __init__(self, G,L,epsilon):
        self.G =G
        self.L = L
        self.epsilon = epsilon
    
    def f(self,x):
        return x - (self.G/ float(self.L))*math.log(x/self.epsilon) 
        
    def findRoot(self):
        initialGuess = (self.G/ float(self.L) )* math.log(self.G/self.epsilon)
        N = fsolve(self.f, initialGuess)
        return round(N)