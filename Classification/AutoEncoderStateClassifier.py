import numpy as np
import scipy.io
import time
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances
import matplotlib.pyplot as plt
import torch
from torch import nn
import h5py


''' Function for loading matlab files '''
def loadmatfile(fname):
    try:
        f = h5py.File(fname,'r')
    except:
        f = scipy.io.loadmat(fname)
    data = {}
    for k, v in f.items():
        data[k] = np.array(v)
    return data


'''  Majority filter '''
def majorityFilt(x, window):
    x = x.astype(int)
    length = int(np.shape(x)[0])
    # New array to store data in
    x_filt = np.zeros(length)
    # Loop through array
    for i in range(length):
        # Get a window of the array
        left = max(0, i-int(window/2))
        right = min(length, i+int(window/2)+1)
        # Counts of each value 
        counts = np.bincount(x[left:right])
        # Indices of maximum values
        maxcount = np.max(counts)
        temp = np.where(counts == maxcount)[0]

        # If only one index is maximum, use that as the new value
        if np.shape(temp) == 1:
            x_filt[i] = temp[0]
        else:
            # If there are ties, the original value wins if it's a part of it.
            # Otherwise, just choose the first one
            if np.any(x[i] == temp):
                x_filt[i] = x[i]
            else:
                x_filt[i] = temp[0]

    return np.squeeze(x_filt)


''' File names and load spectra '''
file = 'R:/Yun/Kronk/Neurochip/Kronk_20191015_01/Spectra.mat'
movefile = 'R:/Yun/Kronk/Neurochip/Kronk_20191015_01/Movement.mat'
idxfile = 'R:/Yun/Kronk/Neurochip/Kronk_20191015_01/SortedIdx.mat'
modelfile = 'R:/Yun/Kronk/Neurochip/Kronk_20191015_01/Autoencoder.pt'

data = loadmatfile(file)
bins = data['bins']
dwnfs = data['dwnfs']
f = data['f']
finish = data['finish']
spectra_raw = data['spectra_raw']
start = data['start']
    
del data
    

''' Preprocessing '''
spectra = np.log(spectra_raw)
spectra = np.mean(spectra, axis=0)

# Normalized spectra for later plotting and testing
spectra_norm = spectra - np.min(spectra, axis=0)
spectra_norm = spectra_norm / np.sum(spectra_norm, axis=0)

# Normalize spectra for training. Multiply by 1000 to rescale.
spectra = spectra - np.min(spectra, axis=0)
spectra = spectra / np.sum(spectra, axis=0)
spectra = np.transpose(spectra)*1000

# # Another normalization, does not work as well
# spectra = np.transpose(spectra) #- np.min(spectra, axis=1)
# avg = np.mean(spectra, axis=0)
# spectra = (spectra-avg) / np.std(spectra, axis=0)

# Convert into tensor
spectra = torch.from_numpy(spectra)
spectra = spectra.type(torch.FloatTensor)

samples = np.size(spectra, axis=0)


''' Autoencoder '''
# Define model
r = 32
epochs = 500  # number of epochs
Losses = np.zeros((epochs))
Losses_weighted = np.zeros((epochs))
avgLoss = 0 # exponentially weighted average
beta = 0.98

d = np.size(spectra,axis=1)

class autoencode(nn.Module):                    
    def __init__(self, d, r):                         
        super(autoencode,self).__init__()      
        h1 = 256 # number of neurons in hidden layers    
        h2 = 128
        h3 = 64
        self.encoder=nn.Sequential(             
            nn.Linear(d, h1),
            nn.BatchNorm1d(h1),
            nn.ReLU(),   
            nn.Linear(h1, h2),
            nn.BatchNorm1d(h2),
            nn.ReLU(),   
            nn.Linear(h2, h3),
            nn.BatchNorm1d(h3),
            nn.ReLU(),   
            nn.Linear(h3, r),
            nn.ReLU(), 
        )                                       
        self.decoder=nn.Sequential(    
            nn.Linear(r, h3),
            nn.BatchNorm1d(h3),
            nn.ReLU(), 
            nn.Linear(h3, h2),
            nn.BatchNorm1d(h2),
            nn.ReLU(),   
            nn.Linear(h2, h1),
            nn.BatchNorm1d(h1),
            nn.ReLU(),   
            nn.Linear(h1, d),
            nn.ReLU(),
        )                                       
                                                
    def forward(self, x):                       
        encoded1 = self.encoder(x)                 
        decoded1 = self.decoder(encoded1)    
        encoded2 = self.encoder(decoded1)
        decoded2 = self.decoder(encoded2)
        return encoded2, decoded2    
    
model = autoencode(d, r)
model.train()
lr = 1e-3
l1_lambda = 1e-6
optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=0)
loss_fn = nn.MSELoss(reduction='mean')


''' Train autoencoder with mini-batches across e epochs '''
batch = 64  # size of batch
loops = int(samples/batch) 
prev_change = False

for e in range(epochs):
    print(e)
    start = time.time()
    # mini-batch
    perm = torch.randperm(samples)
    for b in range(loops):
        ind = perm[b * batch:((b + 1) * batch)]
        tempx = spectra[ind, :]
        _, x_hat = model(tempx)
       
        # MSE loss
        loss = loss_fn(x_hat, tempx)
        
        # Calculate L1 norm         
        l1_norm = sum(p.abs().sum() for p in model.parameters())
        loss += l1_lambda * l1_norm

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
    # Save losses per epoch
    _, x_hat = model(spectra)
    loss = nn.functional.mse_loss(x_hat, spectra)
    Losses[e] = loss.item()
    
    print(Losses[e])
        
    prevloss = avgLoss
    avgLoss = (beta*prevloss + (1-beta) * loss.item())
    Losses_weighted[e] = avgLoss / (1-beta ** (e+1))
    change = Losses_weighted[e-1] -  Losses_weighted[e]

    print(change)
    
    # # If gains are minimal, break from training
    # if change > 0 and change < 1e-5:
    #     if(prev_change):
    #         break
    #     else:
    #         prev_change = True
    # else:
    #     prev_change = False
        
    print(time.time()-start)


''' Save trained autoencoder '''
# Save so we don't have to retrain 
torch.save(model.state_dict(), modelfile)
print('Saved')

# # To load
# model = autoencode(d, r)
# model.load_state_dict(torch.load(modelfile))


''' Plot losses '''
plt.figure()
plt.plot(np.transpose(Losses))
plt.plot(np.transpose(Losses_weighted))
plt.ylabel('Loss')
plt.xlabel('Epochs')
plt.legend(('Loss', 'Exponential Average'))
plt.show()


''' Plot original spectra with reconstructions '''
[reduced, x_hat] = model(spectra)
spectra_np = spectra.numpy()
reconst = x_hat.detach().numpy()
inds = np.floor(np.random.rand((6))*samples)
inds = inds.astype(int)
fig, axes = plt.subplots(2, 3)
for i in range(6):
    axes[np.floor(i/3).astype(int),np.mod(i,3)].plot(spectra_np[inds[i],:])
    axes[np.floor(i/3).astype(int),np.mod(i,3)].plot(reconst[inds[i],:])
plt.show()


''' Dimensionality Reduction '''
# Get reduced representation after encoding
reduced_np = reduced.detach().numpy()


''' k-means '''
# Load movement data and normalize
data = loadmatfile(movefile)
move = data['move']
move = np.log(move)
move = (move-np.mean(move))/np.std(move)
std = np.std(reduced_np, axis=0)
move = move * np.max(std)

reduced_np = np.hstack([reduced_np, move.T])

dist = pairwise_distances(reduced_np)
avgdist = np.mean(dist, axis=0)
lim = np.percentile(avgdist, 90)
coredata = reduced_np[avgdist<lim, :]

kmeans = KMeans(n_clusters=4, max_iter=1000).fit(coredata)
labels = kmeans.predict(reduced_np)


''' Rearrange labels to move, rest, REM, NREM '''
avgmove = np.zeros(4)
avgspect = np.zeros((np.shape(spectra_norm)[0], 4))
for i in range(4):
    avgmove[i] = np.mean(move[0,labels == i])
    avgspect[:,i] = np.mean(spectra_norm[:, labels==i], axis=1)

# New label 
fixedlabels = np.empty(4)
fixedlabels[:] = np.nan

# Most movement is awake and moving
fixedlabels[np.argmax(avgmove)] = 0
remlabel = np.isnan(fixedlabels)

# Highest delta is NREM
deltalim, _ = np.where(f>5)
delta = np.sum(avgspect[0:deltalim[0], :], axis=0)
delta = delta * remlabel
fixedlabels[np.argmax(delta)] = 3
remlabel = np.isnan(fixedlabels)

# Most movement of the remaining labels is awake and at rest
tempmove = avgmove * remlabel
tempmove[tempmove==0] = -np.inf
fixedlabels[np.argmax(tempmove)] = 1

# The remaining label is REM 
fixedlabels[np.isnan(fixedlabels)] = 2

prev_labels = np.copy(labels)
for i in range(4):
    labels[prev_labels==i] = fixedlabels[i]


''' Plot states over time '''
fig, axes = plt.subplots(1, 2)
x = range(samples)
for i in range(4):
    axes[0].scatter([i for i, x in enumerate(labels==i) if x], labels[labels==i])

data = loadmatfile(idxfile)
idx = data['idx']
idx = idx.flatten()
idx = idx-1

for i in range(4):
    axes[1].scatter([i for i, x in enumerate(idx==i) if x], idx[idx==i])
plt.show()

np.sum(labels!=idx)

smooth_labels = majorityFilt(labels, 4)
smoothidx = data['smoothidx']
smoothidx = smoothidx.flatten()
smoothidx = smoothidx-1
fig, axes = plt.subplots(1, 2)
x = range(samples)
for i in range(4):
    axes[0].scatter([i for i, x in enumerate(smooth_labels==i) if x], smooth_labels[smooth_labels==i])
for i in range(4):
    axes[1].scatter([i for i, x in enumerate(smoothidx==i) if x], smoothidx[smoothidx==i])
plt.show()

np.sum(smooth_labels!=smoothidx)


''' Plot averaged spectra '''
plt.figure()
for i in range(4):
    plt.plot(f, np.mean(spectra_norm[:,labels==i], axis=1))
plt.legend(('0','1','2','3'))
plt.show()


''' Maximum values for each dimension of low dimensional representation '''
plt.figure()
for i in range(32):
    temp = np.argsort(reduced_np[:,i])
    plt.plot(f, np.mean(spectra_norm[:, temp[0:100]], axis=1))





