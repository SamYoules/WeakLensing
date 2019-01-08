import scipy as sp

def cov(cl,we):

    mcl = (cl*we).sum(axis=0)  ## sum all the (data points * weights)
    swe = we.sum(axis=0)       ## sum all the weights
    w = swe>0.                 ## mask any negative or zero weights
    mcl[w] /= swe[w]           ## apply mask to both and divide sum of
			       ##       weighted data by sum of weights

    wcl = we*(cl-mcl)          ## wcl vector is weight * variation from the (mean)

    print("Computing cov...")

    co = wcl.T.dot(wcl)
    sswe = swe*swe[:,None]
    w = sswe>0.
    co[w] /= sswe[w]

    return co

def smooth_cov(cl,we,rp,rt,drt=4,drp=4):

    co = cov(cl,we)

    ncl = cl.shape[1]
    var = sp.diagonal(co)
    if sp.any(var==0.):
        print('WARNING: data has some empty bins, impossible to smooth')
        print('WARNING: returning the unsmoothed covariance')
        return co

    cor = co/sp.sqrt(var*var[:,None])

    cor_smooth = sp.zeros([ncl,ncl])

    dcor={}
    dncor={}

    for i in range(ncl):
        sys.stderr.write("\rsmoothing {}".format(i))
        for j in range(i+1,ncl):
            idrp = round(abs(rp[j]-rp[i])/drp)
            idrt = round(abs(rt[i]-rt[j])/drt)
            if not (idrp,idrt) in dcor:
                dcor[(idrp,idrt)]=0.
                dncor[(idrp,idrt)]=0

            dcor[(idrp,idrt)] +=cor[i,j]
            dncor[(idrp,idrt)] +=1

    for i in range(ncl):
        cor_smooth[i,i]=1.
        for j in range(i+1,ncl):
            idrp = round(abs(rp[j]-rp[i])/drp)
            idrt = round(abs(rt[i]-rt[j])/drt)
            cor_smooth[i,j]=dcor[(idrp,idrt)]/dncor[(idrp,idrt)]
            cor_smooth[j,i]=cor_smooth[i,j]


    sys.stderr.write("\n")
    co_smooth = cor_smooth * sp.sqrt(var*var[:,None])
    return co_smooth


