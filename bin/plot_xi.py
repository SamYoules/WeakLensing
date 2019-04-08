

modelfile = ['ximodel-noiseless-cnstwe.txt', 'ximodel_gaussian-v3.txt']
modeltitle = ['Noiseless Mocks: R_t * dxi/drt', 'Noisy Mocks: R_t * dxi/drt']
nbins=50

ion()
for j,k in enumerate(modelfile):
    data_rp, data_rt, xi_dist = loadtxt(k, unpack=1)
    xi_dist_2d = xi_dist.reshape(50,50)
    rt_min = data_rt.reshape(50, 50)[:, 0].max()
    rt_max = data_rt.reshape(50, 50)[:, -1].min()
    rt = linspace(rt_min, rt_max, nbins)
    dxidrt=[]
    for i in range(rt.size - 1):
        dxidrt.append(rt[i]*(xi_dist_2d[i+1]-xi_dist_2d[i])/(rt[i+1]-rt[i]))
    figure()
    pcolormesh(dxidrt, vmin=-0.0002, vmax=0.0002)
    colorbar()
    title(modeltitle[j])

#pcolormesh(dxidrt, vmin=-0.00001, vmax=0.00001)
#############################

modelfile = ['ximodel-noiseless-cnstwe.txt', 'ximodel_gaussian-v3.txt']
modeltitle = ['Noiseless Mocks: xi', 'Noisy Mocks: xi']
nbins=50
ion()
for j,k in enumerate(modelfile):
    data_rp, data_rt, xi_dist = N.loadtxt(k, unpack=1)
    xi_dist_2d = xi_dist.reshape(50,50)
    figure()
    pcolormesh(xi_dist_2d, vmin=-0.0002, vmax=0.0002)
    colorbar()
    title(modeltitle[j])
#####################################

