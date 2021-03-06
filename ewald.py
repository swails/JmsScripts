# All of the variables in the &ewald namelist
class ewald:

   def __init__(self):
      self.sander = {'dsum_tol' : 1.0e-5,'ew_coeff' : 0.0,
         'skinnb' : 0.0,'diptol' : 1.0e-4,'dipmass' : 0.33,'diptau' : 11.0,
         'nfft1' : 0,'nfft2' : 0,'nfft3' : 0,'order' : 4,'opt_infl' : 1,
         'verbose' : 0,'nbflag' : 0,'nbtell' : 0,'netfrc' : 12344321,
         'ew_type' : 0,'vdwmeth' : 1,'eedmeth' : 0,'ee_type' : 1,
         'eedtbdns' : 5000.0,'rsum_tol' : 5.0e-5,'maxexp' : 0.0,'mlimit' : "0,0,0",'use_pme' : 0,
         'maxiter' : 20,'indmeth' : 3,'irstdip' : 0,'nquench' : 0,
         'frameon' : 1,'chngmask' : 1,'scaldip' : 1,
         'gridpointers' : 1,'column_fft' : 1 }

      self.pmemd = {'nfft1' : 0, 'nfft2' : 0, 'nfft3' : 0, 'order' : 4, 'verbose' : 0, 'ew_type' : 0,
            'dsum_tol' : 1.0e-5, 'rsum_tol' : 5.0e-5, 'mlimit' : "0,0,0", 'ew_coeff' : 0.0, 'nbflag' : 1,
            'skinnb' : 2, 'nbtell' : 0, 'netfrc' : -1, 'use_pme' : 1, 'vdwmeth' : 1,
            'eedmeth' : 1, 'eedtbdns' : 5000.0, 'ee_type' : 0, 'frameon' : 1,
            'chngmask' : 1, 'alpha' : 0.0, 'beta' : 0.0, 'gamma' : 0.0,
            'a' : 0.0, 'b' : 0.0, 'c' : 0.0, 'use_axis_opt' : -1, 'fft_grids_per_ang' : 1.0, 
            'block_fft' : -1, 'fft_blk_y_divisor' : -1, 'excl_recip' : -1, 'excl_master' : -1,
            'atm_redist_freq' : -1 }
