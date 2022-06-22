#%%
from spclustering import SPC
import numpy as np
#%%
def test_SPC():
#%%
    rng = np.random.RandomState(0)
    center1= np.array([10, 10])
    center2= np.array([-10, -10])
    nc1 = 1500
    nc2 = 1000
    
    cl1 = rng.multivariate_normal(center1, [[2.5,0.3],[0.3,2.5]], size=nc1)
    cl2 = rng.multivariate_normal(center2, [[3,0.9],[0.9,4]], size=nc2)
    data = np.concatenate([cl1,cl2])

    clustering = SPC(mintemp=0.0,maxtemp=0.03)
    results, sizes = clustering.run(data, return_sizes=True)
    assert sum(results[0,:]==0) == sizes[0,0]
    assert sum(results[1,:]==1) == sizes[1,1]
    
    err = 0.1
    gt = np.ones(nc1+nc2)
    gt[:nc1]=0
    assert (1-sum(results[1,:]==gt)/len(gt))<err

    #import matplotlib.pyplot as plt
    #plt.plot(*data[results[0,:]==0,:].T,'.r')
    #plt.plot(*data[results[0,:]==1,:].T,'.b')


#%%
if __name__ == '__main__':
    test_SPC()