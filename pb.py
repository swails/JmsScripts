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
