import gammalib
import ctools
import cscripts
import numpy as np
from ebltable.tau_from_model import OptDepth
from random import uniform
import xml_generator as xml
from astropy.io import fits
import argparse

tau =  OptDepth.readmodel(model = 'dominguez')

parser = argparse.ArgumentParser()
parser.add_argument('-alert', action='store', dest='alertfile',
                        default='3e-9_all.out.alert', help='File with alerts')
parser.add_argument('--irf', action='store', dest='irf',
                        default='North_z20_average_30m', help='IRF')
parser.add_argument('--obs', action='store', dest='tobs',
                        type=float, default=600.,
                        help='Observation duration time in [s]')
parser.add_argument('--inter', action='store', dest='interaction',
                        default='no',
                        help='Interaction type: pp (proton-proton), pph (proton-photon), txs (TXS-like sources), no (no scaling)')
parser.add_argument('--trans', action='store', dest='trans',
                        default=100.,
                        help='Transient duration [s]')
options = parser.parse_args()

input_model= options.alertfile

imin = 0

gam = 2.19

ep = 100.

ttrans = argparse.trans

tobscta = argparse.tobs

debug = True
edisp = True

caldb='prod3b-v1'
irf=argparse.irf

hdr = fits.Header()
hdr['EXTNAME'] = 'Time profile'
hdr['MJDREFI'] = '59000'
hdr['MJDREFF'] = '5.0000000000E-01'
hdr['TIMEUNIT'] = 's'
hdr['TIMESYS'] = 'TT'
hdr['TIMEREF'] = 'LOCAL'

declination,redshift,A = np.loadtxt(input_model, skiprows=11, unpack=True)
imax = len(redshift)

# flux scaling according to intearction type pp, p-gamma or no scaling
if options.interaction == 'no':
    A_prefix = 1.0
if options.interaction == 'pp':
    A_prefix = np.pow(2.,-gam-1)
if options.interaction == 'pph':
    A_prefix = np.pow(2.,-gam)

realsrc=open('nu_src_ts_'+str(int(ttrans))+'s_'+irf+'_'+str(int(tobscta))+'s_'+str(imin+1)+'-'+str(imax)+'.dat', 'w')
lowrealsrc=open('nu_src_low_ts_'+str(int(ttrans))+'s_'+irf+'_'+str(int(tobscta))+'s_'+str(imin+1)+'-'+str(imax)+'.dat', 'w')
fakesrc=open('nu_src_fake_'+str(int(ttrans))+'s_'+irf+'_'+str(int(tobscta))+'s_'+str(imin+1)+'-'+str(imax)+'.dat', 'w')

for i in xrange(imin, imax):
    z = redshift[i]
    if z < 4.:
        delem = uniform(0.,ttrans)
        delal = uniform(20.,80.)
        delrp = uniform(20.,50.)
        tobs = ttrans * (1. + z)
        delobs = delem * (1. + z)
        tsigstart = delobs + delal + delrp
        if tsigstart < tobs:
            lib,doc = xml.CreateLib()
            LCfile = 'LC_nu_'+str(i+1)+'.fits'
            ta=np.empty(500)
            na=np.empty(500)
            for j in xrange(0, 500):
                ta[j]=j
                if ta[j] < tsigstart:
                    na[j] = 0.
                elif tsigstart < ta[j] < tobs:
                    na[j] = 1.
                else:
                    na[j] = 0.
            time = fits.Column(name='TIME', array=ta, format='1D', unit='s')
            norm = fits.Column(name='NORM', array=na, format='1D')
            t = fits.BinTableHDU.from_columns([time, norm],header=hdr)
            t.writeto(LCfile,overwrite=True) 
            tsig = tobs - tsigstart
            ra = uniform(0.,360.)
            dec = declination[i]
            ETeV = np.logspace(-2,2.5,45)
            EMeV = ETeV * 1e6
            if z < 0.01:
                atten = 1.
            else:
                atten = np.exp(-1. * tau.opt_depth(z,ETeV))
            if options.interaction == 'txs': # reference: https://arxiv.org/abs/1811.07439
                prefac = A[i] * 1e-13
                spec = prefac * (ETeV / ep) ** (-2) * exp(-0.1*(z+1)/ETeV - ETeV/(20.*(z+1)))
            else:
                prefac = A[i] * A_prefix * 1e-13
                spec = prefac * (ETeV / ep) ** (-gam)
            specebl = spec * atten
            sourcename = 'nu'+str(i+1)
            Filefunction = 'spec_nu_ebl_'+str(i+1)+'.dat'
            np.savetxt(Filefunction, np.column_stack([EMeV,specebl + 1.e-300]))
            speci = xml.addFileFunction(lib, sourcename, type = "PointSource", filefun=Filefunction, flux_free=1, flux_value=1., flux_scale=1., flux_max=100000000.0, flux_min=0.0)
            spatial = xml.AddPointLike(doc,ra,dec)
            temporal = xml.AddLCTrans(doc, LCfile, 1.)
            speci.appendChild(spatial)
            speci.appendChild(temporal)
            lib.appendChild(speci)
    
            bkg = xml.addCTAIrfBackground(lib)
            lib.appendChild(bkg)

            open('nu_sources_'+str(i+1)+'.xml', 'w').write(doc.toprettyxml('  '))
            
            sim = ctools.ctobssim()
            sim['inmodel']   = 'nu_sources_'+str(i+1)+'.xml'
            sim['caldb']     = caldb
            sim['irf']       = irf
            sim['ra']        = ra
            sim['dec']       = dec
            sim['rad']       = 5.0
            sim['tmin']      = '2020-05-31T12:00:00'
            sim['tmax']      = '2020-05-31T12:10:00'
            sim['emin']      = 0.02
            sim['emax']      = 199.0
            sim['maxrate']   = 1.0e9
            sim['debug']     = debug
            sim['edisp']     = edisp
            sim.run()

            like = ctools.ctlike(sim.obs())
            like['debug']    = debug
            like['edisp']    = edisp
            like.run()
            
            nuts = like.obs().models()[sourcename].ts()
            nunormsp = like.obs().models()[sourcename].spectral()['Normalization'].value()
            nunormsp_error = like.obs().models()[sourcename].spectral()['Normalization'].error()
            
            if nuts >= 25.:
                if nunormsp > 2. or nunormsp < 0.5:
                    fake = str(i+1)+' '+str(nuts)+' '+str(nunormsp)+' '+str(nunormsp_error)+' '+str(ra)+' '+str(dec)+' '+str(tsig)+'\n'
                    fakesrc.write(fake)
                else:
                    real_nu = str(i+1)+' '+str(nuts)+' '+str(nunormsp)+' '+str(nunormsp_error)+' '+str(ra)+' '+str(dec)+' '+str(tsig)+'\n'
                    realsrc.write(real_nu)
            else:
                lowreal_nu = str(i+1)+' '+str(nuts)+' '+str(nunormsp)+' '+str(nunormsp_error)+' '+str(ra)+' '+str(dec)+' '+str(tsig)+'\n'
                lowrealsrc.write(lowreal_nu)
                    
realsrc.close()
lowrealsrc.close()
fakesrc.close()
