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
# scattering loss and gain terms
scgain2G = np.array(input2Gscatter)
scloss2G = np.diag(np.sum(scgain2G,axis = 0))
# migration and fission terms
migration2G = np.diag(input2GXs['Sigma_a']) + scloss2G - scgain2G
fissionmat2G = np.array([input2GXs['X']]).T@np.array([input2GXs['Sigma_f']])
bmat2G = np.linalg.inv(migration2G)@fissionmat2G
# final calculation
eigval2G,eigvec2G = np.linalg.eig(bmat2G)
eigvec2G[:,np.argmax(eigval2G)]
# power iteration calculation
pit2G = pit(bmat2G,10)

# 2 group output
# guess what. all formatting.
# power iteration frames
pit2Geigval = pd.DataFrame(np.array([pit2G[0][0,0]]), columns = ["power eigenvalues"])
pit2Geigvec = pd.DataFrame(pit2G[1][:,0], columns = ["power eigenvector"])
pit2Gframe = pd.concat([pit2Geigval,pit2Geigvec],axis = 1)
# computation frames
eigvec2Gframe = pd.DataFrame(eigvec2G.T, columns = ["flux1","flux2"])
eigval2Gframe = pd.DataFrame(eigval2G, columns = ["eigenvalues"])
compout2Gframe = pd.concat([eigval2Gframe,eigvec2Gframe], axis = 1)
# combining frames and outputting to file
input2Gcombined = pd.concat([input2Gscatter,input2GXs], axis = 1)
Output2Ga = pd.concat([input2Gcombined,compout2Gframe], axis = 1)
Output2G = pd.concat([Output2Ga,pit2Gframe], axis = 1)
Output2G.to_excel("kep2groupoutput.xlsx")

# 8 group equation
# scattering loss and gain terms
scgain8G = np.array(input8Gscatter)
scloss8G = np.diag(np.sum(scgain8G,axis = 0))
# migration and fission terms
migration8G = np.diag(input8GXs['Sigma_a']) + scloss8G - scgain8G
fissionmat8G = np.array([input8GXs['X']]).T@np.array([input8GXs['Sigma_f']])
bmat8G = np.linalg.inv(migration8G)@fissionmat8G
# final calculation
eigval8G,eigvec8G = np.linalg.eig(bmat8G)
eigvec8G[:,np.argmax(eigval8G)]
# power iteration calculation
pit8G = pit(bmat8G,10)

# 8 group output
# all formatting... TWO
# power iteration frames
pit8Geigval = pd.DataFrame(np.array([pit8G[0][0,0]]), columns = ["power eigenvalues"])
pit8Geigvec = pd.DataFrame(pit8G[1][:,0], columns = ["power eigenvector"])
pit8Gframe = pd.concat([pit8Geigval,pit8Geigvec],axis = 1)
# computation frames
eigvec8Gframe = pd.DataFrame(eigvec8G.T, columns = ["flux1","flux2","flux3", "flux4", "flux5", "flux6", "flux7", "flux8"])
eigval8Gframe = pd.DataFrame(eigval8G, columns = ["eigenvalues"])
compout8Gframe = pd.concat([eigval8Gframe,eigvec8Gframe], axis = 1)
# combining frames and outputting to file
input8Gcombined = pd.concat([input8Gscatter,input8GXs], axis = 1)
Output8Ga = pd.concat([input8Gcombined,compout8Gframe], axis = 1)
Output8G = pd.concat([Output8Ga,pit8Gframe], axis = 1)
Output8G.to_excel("kep8Groupoutput.xlsx")