import numpy as np
import mat73
import matplotlib.pyplot as plt
import pandas as pd

def plot_result(X1,X2,X3):
    plt.plot(np.arange(0,rowC-1),X1,'ro', label='crs')
    plt.plot(np.arange(rowC,rowC+rowF-1),X2,'bo', label='fn')
    plt.plot(np.arange(rowC+rowF,rowC+rowF+rowS-1),X3,'go', label='sed')
    plt.title("AFD 2D")
    plt.legend()
    plt.show()

def read_mat(filePath,variable):
    data_dict=mat73.loadmat(filePath,'r')
    data=np.array(data_dict[variable])
    return data

folderPath='E:\\ETUDES\\M2\\Analyses des données\\TPs\\TP1'
crs=read_mat(folderPath+'\\crs2.mat','ign_crs')
fn=read_mat(folderPath+'\\fn.mat','ign_fn')
sed=read_mat(folderPath+'\\sed.mat','sed_crs')


rowC,colC=crs.shape
rowF,colF=fn.shape
rowS,colS=sed.shape
X=np.vstack((crs,fn,sed))
g=X.mean(0)
g1=X[0:rowC-1,:].mean(0)-g
g2=X[rowC:rowC+rowF-1,:].mean(0)-g
g3=X[rowF:rowF+rowS-1,:].mean(0)-g
C1=np.cov(np.transpose(X[0:rowC-1,:]))
C2=np.cov(np.transpose(X[rowC:rowC+rowF-1,:]))
C3=np.cov(np.transpose(X[rowF+rowF:rowC+rowF+rowS-1,:]))
G=np.vstack((g1,g2,g3))

B=np.matmul(np.transpose(G),G)
C=np.cov(np.transpose(X))
K=C1+C2+C3
P=np.matmul(B,np.linalg.inv(C))
D,V=numpy.linalg.eig(P)
DIn=pd.index(D)
D,I=DIn.sort_values(return_indexer=True,ascending=False)
for i in range(1, len(I)):
    Vzz[i]=V[I[i]]

Xf11=np.matmul(X[0:rowC-1,:],V[:,1])
Xf12=np.matmul(X[0:rowC-1,:],V[:,2])

Xf21=np.matmul(X[rowC:rowC+rowF-1,:],V[:,1])
Xf22=np.matmul(X[rowC:rowC+rowF-1,:],V[:,2])

Xf31=np.matmul(X[rowC+rowF:rowC+rowF+rowS-1,:],V[:,1])
Xf32=np.matmul(X[rowC+rowF:rowC+rowF+rowS-1,:],V[:,2])

def plot_result(X1,X2,X3):
    plt.plot(np.arange(0,rowC-1),X1,'ro', label='crs')
    plt.plot(np.arange(rowC,rowC+rowF-1),X2,'bo', label='fn')
    plt.plot(np.arange(rowC+rowF,rowC+rowF+rowS-1),X3,'go', label='sed')
    plt.title("AFD")
    plt.legend()
    plt.show()

"""""
[V D]=eig(P);%V :les axes descriminents
[D I]=sort(diag(abs(D)),'descend');
V=V(:,I);




load_igcp_1;
load_ign_crs;
load_ign_fn;
load_sed_crs;
crs=ign_crs;
fn=ign_fn;
sed=sed_crs;
close all
clc
X=[crs;fn;sed];
g=mean(X);
X=X;
g1=mean(X(1:34,:))-g;
g2=mean(X(35:67,:))-g;
g3=mean(X(68:end,:))-g;
C1=cov(X(1:34,:));
C2=cov(X(35:67,:));
C3=cov(X(68:end,:));

G=[g1;g2;g3];
B=G'*G;
C=cov(X);%matrice de covariance (X non centré)
K=C1+C2+C3;
P=B*inv(C);
[V D]=eig(P);%V :les axes descriminents
[D I]=sort(diag(abs(D)),'descend');
V=V(:,I);

Xf11=X(1:34,:)*V(:,1);
Xf12=X(1:34,:)*V(:,2);

Xf21=X(35:67,:)*V(:,1);
Xf22=X(35:67,:)*V(:,2);

Xf31=X(68:end,:)*V(:,1);
Xf32=X(68:end,:)*V(:,2);
"""
