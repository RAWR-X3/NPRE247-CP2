import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# inputting files through pandas
input2Gscatter = pd.read_excel('kep2groupscatterinput.xlsx')
input2GXs = pd.read_excel('kep2groupXSinput.xlsx')
input8Gscatter = pd.read_excel('kep8groupscatterinput.xlsx')
input8GXs = pd.read_excel('kep8groupXSinput.xlsx')

def pit(B,l):
    k = 0.69 # nice
    flux = np.ones(B.shape)
    for i in range(l):
        flux = (B@flux)/np.linalg.norm(B@flux)
        k = (((B@flux).T)@flux)/(flux.T@flux)
    return [k,flux]

# 2 group equation
scgain2G = np.array(input2Gscatter)
scloss2G = np.diag(np.sum(scgain2G,axis = 0))
migration2G = np.diag(input2GXs['Sigma_a']) + scloss2G - scgain2G
fissionmat2G = np.array([input2GXs['X']]).T@np.array([input2GXs['Sigma_f']])
bmat2G = np.linalg.inv(migration2G)@fissionmat2G
eigval2G,eigvec2G = np.linalg.eig(bmat2G)
eigvec2G[:,np.argmax(eigval2G)]
pit2G = pit(bmat2G,10)

# 8 group equation
scgain8G = np.array(input8Gscatter)
scloss8G = np.diag(np.sum(scgain8G,axis = 0))
migration8G = np.diag(input8GXs['Sigma_a']) + scloss8G - scgain8G
fissionmat8G = np.array([input8GXs['X']]).T@np.array([input8GXs['Sigma_f']])
bmat8G = np.linalg.inv(migration8G)@fissionmat8G
eigval8G,eigvec8G = np.linalg.eig(bmat8G)
eigvec8G[:,np.argmax(eigval8G)]
pit8G = pit(bmat8G,10)