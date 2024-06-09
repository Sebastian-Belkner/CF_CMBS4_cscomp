"""
MF estimator helper on simulated CMB-S4 component separation comparison project observed CMB polarization data.

QE and iterative reconstruction uses anisotropic filters. 
"""

import numpy as np
import os
from os.path import join as opj
import healpy as hp

import delensalot
from delensalot import utils
from delensalot.utility.utils_hp import gauss_beam
import delensalot.core.power.pospace as pospace
from delensalot.config.config_helper import LEREPI_Constants as lc
from delensalot.config.metamodel.dlensalot_mm import *

fg = '13'
desc_flag = 'Lmin11'

mask_fn = opj('/global/homes/s/sebibel/data/mask/cmbs4_s08b_mask.fits')
rhits_fn = opj(os.environ['SCRATCH'], 'data/cmbs4/reanalysis/mapphi_intermediate/s08b/masks/08b_rhits_positive_nonan.fits')

data_dir = '/pscratch/sd/f/fbianc/CMBS4/lowellbb/'
data_dir_noise = opj(os.environ['CFS'], 'cmbs4/awg/lowellbb/reanalysis/foreground_cleaned_maps/')

def func(data):
    return data * 1e6 #np.nan_to_num(utils.cli(hp.read_map(rhits_fn)))

dlensalot_model = DLENSALOT_Model(
    defaults_to = 'default_CMBS4_maskedsky_polarization',
    validate_model = False,
    job = DLENSALOT_Job(
        jobs = ["QE_lensrec"]
    ),
    computing = DLENSALOT_Computing(
        OMP_NUM_THREADS = 8
    ),                              
    analysis = DLENSALOT_Analysis(
        key = 'p_p',
        version = '',
        simidxs = np.arange(0,96),
        simidxs_mf = np.arange(0,96),
        TEMP_suffix = 'cscomp_08b_fg{}_{}'.format(fg, desc_flag),
        Lmin = 11, 
        lm_max_ivf = (4000, 4000),
        lmin_teb = (30, 30, 200),
        # zbounds = ('mr_relative', 10.),
        # zbounds_len = ('extend', 5.),
        zbounds = (-1,1),
        zbounds_len = (-1,1),
        beam = 2.3,
        mask = mask_fn
    ),
    simulationdata = DLENSALOT_Simulation(
        space = 'cl', 
        flavour = 'unl',
        lmax = 4096,
        phi_lmax = 5120,
        nlev = {'P': 0.5, 'T': np.sqrt(1)},
        transfunction = gauss_beam(2.3/180/60 * np.pi, lmax=4096),
        geominfo = ('healpix', {'nside': 2048}),
    ),
    noisemodel = DLENSALOT_Noisemodel(
        OBD = 'trunc',
        sky_coverage = 'masked',
        spectrum_type = 'white',
        nlev = {'P': 0.5, 'T': np.sqrt(1)},
        geominfo = ('healpix',{'nside': 2048}),
        rhits_normalised = (rhits_fn, np.inf)
    ),
    # obd = DLENSALOT_OBD(
    #     libdir = '/global/cfs/cdirs/cmbs4/awg/lowellbb/reanalysis/mapphi_intermediate/s08b',
    #     nlev_dep = 1e4,
    #     rescale = (0.5/0.3505)**2
    # ),
    qerec = DLENSALOT_Qerec(
        tasks = ["calc_phi", "calc_blt"],
        filter_directional = 'anisotropic',
        lm_max_qlm = (4000, 4000),
        cg_tol = 1e-4
    ),
    itrec = DLENSALOT_Itrec(
        tasks = ["calc_phi", "calc_blt"],
        filter_directional = 'anisotropic',
        itmax = 10,
        cg_tol = 1e-4,
        lm_max_unl = (4200, 4200),
        lm_max_qlm = (4000, 4000),
        stepper = DLENSALOT_Stepper(
            typ = 'harmonicbump',
            lmax_qlm = 4000,
            mmax_qlm = 4000,
            a = 0.5,
            b = 0.499,
            xa = 400,
            xb = 1500
        ),
        lenjob_geominfo= ('thingauss',{'lmax': 4500, 'smax': 3}),
        lenjob_pbdgeominfo= ('pbd', (0., 1*np.pi)),
        # mfvar = '/pscratch/sd/s/sebibel/analysis/cscomp_08b_fg13_Lmin11_lminB200/QE/mf_allsims.npy'
    ),
    madel = DLENSALOT_Mapdelensing(
        data_from_CFS = False,
        edges = lc.AoA_edges,
        iterations = [10,12,14,16,18],
        masks_fn = ['/global/homes/s/sebibel/data/mask/cmbs4_s08b_mask.fits'],
        lmax = 1024,
        Cl_fid = 'obs', #this doesn't make sense right now..
        basemap = 'obs',
        libdir_it = None,
        binning = 'binned',
        spectrum_calculator = pospace,
    ),
)
    