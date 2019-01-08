import picca
from picca import constants
from picca.wedgize import wedge
from picca.data import delta
from picca.fitter.parameters import parameters
from picca.fitter.cosmo import model
from picca.fitter import metals
from astropy.io import fits
import numpy as N
import pylab as P
import scipy as sp
import configargparse


### Get the parser
param = parameters()
parser = configargparse.ArgParser()

parser.add('-c', '--config_file', required=False, \
           is_config_file=True, help='config file path')

for i in param.dic_init_float:
    parser.add('--'+i, type=float, required=False, \
               help=param.help[i], default=param.dic_init_float[i])
for i in param.dic_init_int:
    parser.add('--'+i, type=int,   required=False,  \
               help=param.help[i], default=param.dic_init_int[i])
for i in param.dic_init_bool:
    parser.add('--'+i, action='store_true',  required=False,  \
               help=param.help[i], default=param.dic_init_bool[i])
for i in param.dic_init_string:
    parser.add('--'+i, type=str, required=False, \
               help=param.help[i], default=param.dic_init_string[i])
for i in param.dic_init_list_string:
    parser.add('--'+i, type=str, required=False, \
               help=param.help[i], \
               default=param.dic_init_list_string[i],nargs="*")

args, unknown = parser.parse_known_args()
param.set_parameters_from_parser(args=args,unknown=unknown)
param.test_init_is_valid()
dic_init = param.dic_init

m = model(dic_init)
m.add_auto(dic_init)

met = None
if dic_init['metals'] is not None:
    met=metals.model(dic_init) 
    met.templates = not dic_init['metal_dmat']
    met.grid = not dic_init['metal_xdmat']
    met.add_auto()

data = fits.open(dic_init['data_auto'])[1].data 

xi = m.valueAuto( data.RP, data.RT, data.Z, dic_init)
if met: 
    xi += met.valueAuto(dic_init)

xid = data.DM.dot(xi)

d0 = fits.open('cf-unlensed-exp.out.gz')[1].data

for mus in [[0., 0.5], [0.5, 0.8], [0.8, 0.95], [0.95, 1.]]:
        w = picca.wedgize.wedge(rpmin=0.0, rpmax=200.0, nrp=50, \
                rtmin=0.0, rtmax=200.0, nrt=50, \
                rmin=0.0, rmax=200.0, nr=50, \
                mumin=mus[0], mumax=mus[1], ss=10)
        r, wed0, wedcov0 = w.wedge(d0.DA, d0.CO)
        r, wedm, wedcov0 = w.wedge(xid, data.CO)

        dwed0 = N.sqrt(N.diag(wedcov0))

        P.figure()
        P.rcParams.update({'font.size':22})
        P.errorbar(r, wed0*r**2, dwed0*r**2, fmt='o', color='#3498db')
        P.plot(r, wedm*r**2, color='#a569bd')
        P.legend()
        P.ylim(-1.2, 0.5)
        P.title(r'$%.1f < \mu < %.1f$'%(mus[0], mus[1]))
        P.xlabel(r'$r [Mpc/h]$')
        P.ylabel(r'$r^2 \xi(r, \mu) [(Mpc/h)^2]$')
P.show()

