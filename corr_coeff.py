from numpy import ones,mean,sqrt,diag
from numpy.matlib import repmat

def corr_coef(X):
    n=len(X)
    Y=X-repmat(mean(X),n,1)
    one=ones(1,n).reshape(n)
    w=diag((1/n)*one)
    X_s=sqrt(np.diag(np.matmul(np.matmul(np.transpose(Y),w),Y)))
    D=np.diag(1/X_s)
    C=np.matmul(np.matmul(np.transpose(Y),w),Y)
    Cf = np.matmul(np.matmul(D, C), D)
    return Cf
