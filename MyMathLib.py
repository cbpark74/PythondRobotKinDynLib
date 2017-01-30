'''
Created on Mar 7, 2010

@author: cbpark
'''
from numpy.linalg import *
from numpy import *

from UnitConversion import *

def arange2(start,end, interval):
    return arange(start, end+interval, interval)


def GetDictFromKeys(keys):
    
    dict = OrderedDict()
    idx = 0
    for key in keys:
        dict[key]= idx
        idx += 1
    
    return dict

def cross(u, v):
    n = zeros(3);
#    print u, v, n
    n[0] = u[1]*v[2] - u[2]*v[1];
    n[1] = u[2]*v[0] - u[0]*v[2];
    n[2] = u[0]*v[1] - u[1]*v[0];
    
    return n

def norm(v):
    return sqrt(dot(v,v))
    
import traceback
def dprint(__x):
    print traceback.extract_stack(limit=2)[0][3][7:][:-1]
    print __x

def nprint(__x):
    print traceback.extract_stack(limit=2)[0][3][7:][:-1]
    print __x

    
def afs(ss): # 1D & 2D array from string
    sl = ss.split(";")
    
    if (len(sl)==1):
        return fromstring(ss,sep=' ')
    
#    print sl
    vl =[]
    for s in sl:
        v = fromstring(s,sep=' ')
        vl.append(v)
    
    return array(vl)


def skew(a):
    
    A = array([\
               [0       ,-a[2]  , a[1]],
               [a[2]    ,0      , -a[0]],
               [-a[1]   ,a[0]   , 0]
               ])

    return A

def LAd(T, A1):
    R = T[0:3, 0:3]
    P = T[0:3, 3]
    
    Op = zeros((6,6))
    Op[0:3,0:3] = R
    Op[3:6,3:6] = R
    Op[3:6,0:3] = dot(skew(P),R)
    
    #dprint(P)
    #dprint(Op)
    #dprint(A1)
    A2 = dot(Op,A1)
    return A2

def dLAd(T, A1):
    R = T[0:3, 0:3]
    P = T[0:3, 3]
    
    Op = zeros((6,6))
    Op[0:3,0:3] = R
    Op[3:6,3:6] = R
    Op[3:6,0:3] = dot(skew(P),R)
    
    Op = Op.T
    
    A2 = dot(Op,A1)
    
    return A2

def sad(V1, V2):
    w1= V1[0:3]
    v1= V1[3:6]
    
    Op = zeros((6,6))
    Op[0:3,0:3] = skew(w1)
    Op[3:6,3:6] = skew(w1)
    Op[3:6,0:3] = skew(v1)
    
    V = dot(Op, V2)
    
    return V


def dsad(V1, V2):
    w1= V1[0:3]
    v1= V1[3:6]
    
    Op = zeros((6,6))
    Op[0:3,0:3] = skew(w1)
    Op[3:6,3:6] = skew(w1)
    Op[3:6,0:3] = skew(v1)
    
    Op = Op.T
    
    V = dot(Op, V2)
    
    return V


def ToSE3(R,P):
    T = eye(4)
    T[0:3,0:3] = R
    T[0:3,3]= P

    return T
def FromSE3(T):
    R = T[0:3,0:3];
    P = T[0:3,3];
    return R, P




def ExpOfso3(w):
    I = eye(3)
    eps = 1e-6
    
    norm_w = norm(w)
    skew_w = skew(w)
    
    if (norm_w<eps):
        R = I
    else:
        R = I + sin(norm_w)*skew_w/norm_w + (1-cos(norm_w))*dot(skew_w,skew_w) / norm_w**2
    
    return R
def DRotX(q):
    q =d2r(q)
    R = array([
              [1,0,0],
              [0,cos(q),-sin(q)],
              [0,sin(q),cos(q)],
              ])
    return R

def DRotY(q):
    q =d2r(q)
    R = array([
              [cos(q),0,sin(q)],
              [0,1,0],
              [-sin(q),0,cos(q)],
              ])
    return R

def DRotZ(q):
    q =d2r(q)
    R = array([
              [cos(q),-sin(q),0],
              [sin(q),cos(q),0],
              [0,0,1],
              ])
    return R
def RotX_rad(q):
    R = array([
              [1,0,0],
              [0,cos(q),-sin(q)],
              [0,sin(q),cos(q)],
              ])
    return R

def RotY_rad(q):
    R = array([
              [cos(q),0,sin(q)],
              [0,1,0],
              [-sin(q),0,cos(q)],
              ])
    return R

def RotZ_rad(q):
    R = array([
              [cos(q),-sin(q),0],
              [sin(q),cos(q),0],
              [0,0,1],
              ])
    return R
    
def ExpOfse3(V):
    I = eye(3)
    eps = 1e-6
    
#    dprint(V)
    
    w = V[0:3]
    v = V[3:6]
    norm_w = norm(w)
    
#    dprint(v)
    if (norm_w < eps):
        T = eye(4)
        T[0:3,3]= v
        
        return T
    
    skew_w = skew(w)
    
    R = ExpOfso3(w)
    
    A = I + (1-cos(norm_w))*skew_w/norm_w**2+(norm_w-sin(norm_w))*dot(skew_w,skew_w)/(norm_w**3)
    
    P = dot(A,v)
    T = ToSE3(R,P)
    
    return T
    
    
def POE(A, q):
    DOF = len(A)
    
    T = eye(4)
    
    for i in range(DOF):
        E = ExpOfse3(A[i]*q[i])
        T = dot(T,E)
        
    return T
    

def FwdKin (Robot, q):
    T = dot(POE(Robot.A, q),Robot.M0)
    return T


def FwdKinPOE (A, M0, q): # I am not sure if this function is used or not
    T = dot(POE(A,q), M0)
    return T

def To4Vec(p):
    p4 = append(p,1)
    return p4

def ComputeJointPosition( Robot, q):
    T_EE = dot(POE(Robot.A, q), Robot.M0)
    
    T = eye(4)
    JP = zeros((3,Robot.DOF+1))
    for i in range(Robot.DOF):
        p4 = append(Robot.p[i],1)
        P = dot(T, p4)
        JP[:,i] = P[0:3]
        T = dot(T, ExpOfse3(Robot.A[i]*q[i]) )
        
    JP[:,Robot.DOF] = T_EE[0:3,3]
    
    return JP, T_EE
def trace(A):
    t = sum(diag(A))
    return t

def LogOfSO3(R):
    eps = 1e-7 # ???
    
    a = (trace(R)-1)/2
    if a >1:
        a = 1
    elif a <-1:
        a = -1

    q = math.acos(a)

    if (abs(q)<eps):
        w = zeros(3)
    else:
        if (abs(q-pi)<eps):
            skew_w = q/(2*sin(q))*(R-R.T)
            w = array([skew_w[2,1], skew_w[0,2], skew_w[1,0]])
            if (abs(w[0])>abs(w[1]))and(abs(w[0])>abs(w[2])):
                w = array([1,0,0])*q*sign(w[0])
            if (abs(w[1])>abs(w[2]))and(abs(w[1])>abs(w[0])):
                w = array([0,1,0])*q*sign(w[1])
            if (abs(w[2])>abs(w[1]))and(abs(w[2])>abs(w[0])):
                w = array([0,0,1])*q*sign(w[2])
            
        else:
            skew_w = q/(2*sin(q))*(R-R.T)
            w = array([skew_w[2,1], skew_w[0,2], skew_w[1,0]])
        
    return w


def InvKinNum(Robot, T_tgt, q0):
    
    R_tgt, P_tgt = FromSE3(T_tgt)
    q_cur = q0
    
    
    
    ret = False
#    step = 0.5 # 0 ~ 1
    step = 0.3 # 0 ~ 1
    iter = 0
    MaxIter = 40
    eps = 1e-6
    
    while (ret == False):
        iter += 1
        T_cur = FwdKin(Robot, q_cur)
        R_cur, P_cur = FromSE3(T_cur)
        dR = dot(R_tgt,R_cur.T)
        
        w = LogOfSO3(dR)
        v = P_tgt - P_cur
        
        dX = hstack([w,v])
        
        J = ComputeJacobian( Robot, q_cur)
        dq = dot(pinv(J),dX)*step
        
        q_cur += dq

#        dprint(dR)
#
#        dprint(w)
#        dprint(v)
        
        norm_dPos = norm(v)
        norm_dOri = norm(w)
        norm_dX = max(norm_dPos , norm_dOri)
#        dprint(norm_dPos)
#        dprint(norm_dOri)
#        dprint(norm_dX)
        if (norm_dX < eps):
            q_tgt = q_cur
            ret = True
#            dprint(iter)
            return q_tgt, ret
        
        if (iter > MaxIter):
            q_tgt = q_cur
            ret = False
#            print 'iter > MaxIter'
#            dprint(iter)
#            dprint(eps)
#            dprint(norm_dX)
            return q_tgt, ret
        
    
    return q_tgt, ret


def ComputeJs(A, q):
    DOF = len(q)
    
    POA = eye(4)
    Js = zeros((6,DOF))
    
    for i in range(DOF):
        Js[:,i]= LAd(POA, A[i])
        #dprint(LAd(POA, A[i]))
        POA = dot(POA, ExpOfse3(A[i]*q[i]) ) 
    
    return Js
        
def ComputeJacobian(Robot, q):
    T = FwdKin(Robot, q)
    
    p = T[0:3,3]
    Js = ComputeJs(Robot.A, q)

    s2ee = eye(6)
    s2ee[3:6,0:3]=skew(-1*p)
    
    Jee_p = dot(s2ee,Js)

    Jee = eye(6)

    Jee[0:3,0:6] = Jee_p[3:6,0:6]    
    Jee[3:6,0:6] = Jee_p[0:3,0:6]    
    
    return Jee



def ComputeTorque(Robot, q, dq, ddq, F_ext, gravity ):
    '''
    prerequisite:
      
    '''
    if gravity == None:
        gravity = array([0, 0, -9.81])
    
    V_0 = array([0,0,0,0,0,0])
    dV_0 = zeros(6)
    dV_0[3:6] = -gravity
    
    f = range(Robot.DOF+1)
    V = range(Robot.DOF)
    dV = range(Robot.DOF)
    for i in range(Robot.DOF):
        f[i] = dot(Robot.M[i], ExpOfse3(Robot.S[i]*q[i]))
        if (i==0):
            V_prev = V_0
            dV_prev = dV_0
        else:
            V_prev = V[i-1]
            dV_prev = dV[i-1]
        
        V[i] = (LAd(inv(f[i]),V_prev)+Robot.S[i]*dq[i])
        dV[i] = Robot.S[i] * ddq[i] + LAd(inv(f[i]),dV_prev) + sad(V[i], Robot.S[i]*dq[i])
    
    f[Robot.DOF] = Robot.M[Robot.DOF]
    
    m = Robot.m
    r = Robot.r
    I = Robot.I

    J = range(Robot.DOF+1)
    tau = range(Robot.DOF)
    F = range(Robot.DOF)
    for i in range(Robot.DOF-1,-1,-1):
        if (i==Robot.DOF-1):
            F_prev = zeros(6)
        else:
            F_prev = F[i+1]
        
        J[i] = zeros((6,6))
        J[i][0:3,0:3]= I[i]-m[i]*dot(skew(r[i]),skew(r[i]))
        J[i][0:3,3:6]= m[i]*skew(r[i])
        J[i][3:6,0:3]= -m[i]*skew(r[i])
        J[i][3:6,3:6]= m[i]*eye(3)
        
        if (i==Robot.DOF-1):
            F_ext_i = F_ext
        else:
            F_ext_i = zeros(6)
  
            
        F[i] = dLAd(inv(f[i+1]),F_prev) + dot(J[i],dV[i])-dsad(V[i], dot(J[i],V[i])) - dLAd(inv(f[i+1]), F_ext_i)
        
        tau[i] = dot(Robot.S[i], F[i])
        
    return tau

def ComputeDeflection(Robot, q, tau):
    """
    def ComputeDeflection(Robot, q, tau):
    
    Prerequisite:
      ComputeTorque
    """    
    dq = dot(-Robot.InvK, tau)
    J = ComputeJacobian(Robot, q)
    dX = dot(J, dq)
    
    return dX, dq

## {{{ Recipe 576693 (r5): Ordered Dictionary for Py2.4 
from UserDict import DictMixin

class OrderedDict(dict, DictMixin):

    def __init__(self, *args, **kwds):
        if len(args) > 1:
            raise TypeError('expected at most 1 arguments, got %d' % len(args))
        try:
            self.__end
        except AttributeError:
            self.clear()
        self.update(*args, **kwds)

    def clear(self):
        self.__end = end = []
        end += [None, end, end]         # sentinel node for doubly linked list
        self.__map = {}                 # key --> [key, prev, next]
        dict.clear(self)

    def __setitem__(self, key, value):
        if key not in self:
            end = self.__end
            curr = end[1]
            curr[2] = end[1] = self.__map[key] = [key, curr, end]
        dict.__setitem__(self, key, value)

    def __delitem__(self, key):
        dict.__delitem__(self, key)
        key, prev, next = self.__map.pop(key)
        prev[2] = next
        next[1] = prev

    def __iter__(self):
        end = self.__end
        curr = end[2]
        while curr is not end:
            yield curr[0]
            curr = curr[2]

    def __reversed__(self):
        end = self.__end
        curr = end[1]
        while curr is not end:
            yield curr[0]
            curr = curr[1]

    def popitem(self, last=True):
        if not self:
            raise KeyError('dictionary is empty')
        if last:
            key = reversed(self).next()
        else:
            key = iter(self).next()
        value = self.pop(key)
        return key, value

    def __reduce__(self):
        items = [[k, self[k]] for k in self]
        tmp = self.__map, self.__end
        del self.__map, self.__end
        inst_dict = vars(self).copy()
        self.__map, self.__end = tmp
        if inst_dict:
            return (self.__class__, (items,), inst_dict)
        return self.__class__, (items,)

    def keys(self):
        return list(self)

    setdefault = DictMixin.setdefault
    update = DictMixin.update
    pop = DictMixin.pop
    values = DictMixin.values
    items = DictMixin.items
    iterkeys = DictMixin.iterkeys
    itervalues = DictMixin.itervalues
    iteritems = DictMixin.iteritems

    def __repr__(self):
        if not self:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, self.items())

    def copy(self):
        return self.__class__(self)

    @classmethod
    def fromkeys(cls, iterable, value=None):
        d = cls()
        for key in iterable:
            d[key] = value
        return d

    def __eq__(self, other):
        if isinstance(other, OrderedDict):
            return len(self)==len(other) and \
                   min(p==q for p, q in  zip(self.items(), other.items()))
        return dict.__eq__(self, other)

    def __ne__(self, other):
        return not self == other
## End of recipe 576693 }}}


#if __name__ == "__main__":
#    
#    from RobotKinDyn_Schunk_LWA3 import Schunk_LWA3_Robot
#    from RobotKinDyn_Schunk_LWA3 import Schunk_LWA3_Robot
#    
#    Robot = Schunk_LWA3_Robot()
#    q_zero = array([0, 0, 0, 0, 0, 0, 0]) * pi / 180
#    q_arb = array([10, 10, 10, 90, 0, 90, 0]) * pi / 180
#    q_home = array([0, 0, 0, 90, 0, 90, 0]) * pi / 180
#    q = q_arb
#    dq = zeros(7)
#    ddq = zeros(7)
#    Jee = ComputeJacobian(Robot, q)
#    Force = array([50, 0, 0])
#    Torque = array([0, 0, 0])
#    F_ext = hstack([Torque, Force])
#    
#    gravity = array([0, 0, -9.81])
##    gravity = array([0,0,0])
#    
#    tau = ComputeTorque(Robot, q, zeros(7), ddq, F_ext, gravity)
#
#    dX, dq = ComputeDeflection(Robot, q, tau)
#    
#    dX_Pos_mm = 1000*dX[3:6]
##    nprint(Jee)
##    nprint(tau)
##    nprint(dX)
##    nprint(dX_Pos_mm)
#    
#    T = FwdKin(Robot, q)
#    nprint(T)
#
#    q2, ret = InvKinNum(Robot, T, q_home)
#    nprint(q2)
#    
#    q2_deg = floor(rad2deg(q2)*10)/10
#    nprint(q2_deg)
#    nprint(rad2deg(q))
#    
#    T2 = FwdKin(Robot, q2)
#    
#    nprint(T)
#    nprint(T2)
#    nprint(T2-T)
#    
#    a = array([1,2])
#    a = append(a, 3)
#    nprint(a)
#    
#
#    JP, EE = ComputeJointPosition(Robot, q)
#    
#    nprint(JP)

def PrintMatrixRaw(m):
    for i in range(len(m)):
        for j in range(len(m[0])):
            txt = "%g\t" % (m[i][j])
            print txt,
        print
    

def Differ(X, dT, N = 5):
    DOF = len(X[0])
    dX = zeros((len(X),DOF))
    for i in range(N,len(X)):
        if (i<len(X)-N):
            dX[i,:] = (X[i+N,:]-X[i-N,:])/(2*N*dT)
        else:
            dX[i,:] = dX[i-1,:] 

    for i in range(N):
        dX[i,:] = dX[N,:]
        
    return dX

def Differ1D(X, dT, N=5):
    dX = zeros((len(X)))
    for i in range(N,len(X)):
        if (i<len(X)-N):
            dX[i] = (X[i+N]-X[i-N])/(2*N*dT)
        else:
            dX[i] = dX[i-1] 

    for i in range(N):
        dX[i] = dX[N]

    return dX

def StrListToArray(list):
    size = len(list)
    arr = zeros((size,1))
    
    i=0
    for str in list:
        arr[i]=float(str)
        i+=1
        
    return arr    


def Tx(d):
    R = eye(3)
    P = array([d, 0, 0])
    T = ToSE3(R,P)
    return T

def Ty(d):
    R = eye(3)
    P = array([0, d, 0])
    T = ToSE3(R,P)
    return T

def Tz(d):
    R = eye(3)
    P = array([0, 0, d])
    T = ToSE3(R,P)
    return T
    
def Rx(deg):
    R = DRotX(deg)
    P = array([0, 0, 0])
    T = ToSE3(R,P)
    return T

def Ry(deg):
    R = DRotY(deg)
    P = array([0, 0, 0])
    T = ToSE3(R,P)
    return T

def Rz(deg):
    R = DRotZ(deg)
    P = array([0, 0, 0])
    T = ToSE3(R,P)
    return T

def Trim(M):
    m = len(M)
    n = len(M[0])
    Eps = 1e-10
    for i in range(m):
        for j in range(n):
            if (abs(M[i,j]) <Eps):
                M[i,j]=0
                
    return M
