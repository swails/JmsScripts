# All of the variables in the pb namelist
class pb:

   def __init__(self):
      self.sander = {'epsin' : 1.0, 'epsout' : 80.0, 'smoothopt' : 0, 'istrng' : 0.0, 'pbtemp' : 300.0,    
                     'radiopt' : 1, 'dprob' : 0.0, 'iprob' : 2.0, 'npbopt' : 0, 'solvopt' : 1, 'accept' : 0.001, 'maxitn' : 100, 
                     'fillratio' : 2.0, 'space' : 0.5, 'nbuffer' : 0, 'nfocus' : 2, 'fscale' : 8, 'npbgrid' : 1,     
                     'arcres' : 0.0625,'dbfopt' : -1,'bcopt' : 5,'scalec' : 0,'eneopt' : -1,'frcopt' : 0,'cutres' : 99.0,'cutfd' : 5.0,  
                     'cutnb' : 0.0, 'nsnbr' : 1, 'nsnba' : 1,'phiout' : 0, 'phiform' : 0, 'npbverb' : 0, 'npopt' : 2,    
                     'decompopt' : 1, 'use_rmin' : 0, 'sprob' : 1.6, 'vprob' : 1.28, 'rhow_effect' : 1.0, 'use_sav' : 0,
                     'cavity_surften' : 0.04356, 'cavity_offset' : -1.008, 'maxsph' : 400, 'maxarc' : 256,          
                     'cutsa' : 9.0, 'ndofd' : 1, 'ndosas' : 1, 'fmiccg' : -0.3, 'ivalence' : 1.0, 'laccept' : 0.1, 'wsor' : 1.9,  
                     'lwsor' : 1.95, 'pbkappa' : 0, 'radinc' : 0.8, 'expthresh' : 0.2, 'offx' : 0.0, 'offy' : 0.0, 'offz' : 0.0,    
                     'sepbuf' : 4.0, 'mpopt' : 0, 'lmax' : 80, 'maxarcdot' : 0 }

      self.sanderAPBS = {'dime' : '0,0,0', 'pdime' : '0,0,0', 'cglen' : '0,0,0', 'fglen' : '0,0,0', 'grid' : '0,0,0', 'nlev' : 4,
                         'nonlin' : 0, 'bcfl' : 1, 'nion' : 0, 'pdie' : 2.0, 'sdie' : 78.4, 'srfm' : 2, 'chgm' : 1, 
                         'temp' : 298.15, 'gamma' : 0.105, 'sdens' : 10.0, 'calc_type' : 1, 'calcnpenergy' : 1, 
                         'cmeth' : 1, 'ccmeth' : 1, 'fcmeth' : 1, 'ionq' : 'NO_DEFAULT', 'ionc' : 'NO_DEFAULT', 'ionrr' : 'NO_DEFAULT',
                         'smvolume' : 0, 'smsize' : 0, 'srad' : 1.4, 'swin' : 0.3,
                         'calcenergy' : 1, 'calcforce' : 0, 'apbs_debug' : 0, 'sp_apbs' : '.false.', 'apbs_print' : 1,
                         'wpot' : 0, 'wchg' : 0, 'wsmol' : 0, 'ispara' : 0, 'radiopt' : 0, 'geom_upd_limit' : 0,
                         'evdw_upd_limit' : 0, 'pqr' : 0, 'dime_updates' : 0,
                         'wkappa' : 0, 'wdiel' : 0, 'rchg' : 0, 'rkappa' : 0, 'rdiel' : 0 }
