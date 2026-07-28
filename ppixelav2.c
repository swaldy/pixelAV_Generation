/* Subroutine */ int deposit(float *vect, float thick, int new_drde, int nmax, vect_or_f xeh[2][NEHSTORE], int *neh)
{
    /* Initialized data */
	
    static float twome = (float)1.022e6;
    static int fcall = -1;
	
    /* System generated locals */
    int i__1;
    double r__1, r__2, r__3;

    /* Local variables */
    static float path, rvec[5], phip, costd, phid;
    static vect_or_f pxvec, pvvec, dxvec, dvvec;
    static vect_or_f xin, avec, vzero, y1, y2;
    vect_or_f d,x1,x2,vsplat,dvvec2,vb;
#ifdef __POWERPC__
    static vect_or_c mask1, mask2;
#endif
    static int ierr;
    static float cost, step, sint, e;
    static int i__, npair;
    static float ehnum, pdelta, ppion, etar, rumfp;
    static float costp, rnorm, twopi, r1, r2, te, deltae;
    static int ndelta;
    static float pathmf;
    static float phi;
    static int ipr;
	
	
	/* ********************************************************************** */
	/* * This routine propagates a pion (position and direction described   * */
	/* * by VECT) through a slab of silicon of thickness THICK (between     * */
	/* * z=0 and z=THICK um).  During the passage, a total of NEH e- and    * */
	/* * holes are created.  The positions and momenta of e- h are logged   * */
	/* * in VEH.                                                            * */
	/* * Parameters: VECT[6] - VECT[0-2] are the x,y,z coordinates of a     * */
	/* *                         point on the track in microns,             * */
	/* *                         VECT[3-5] are the velocity in um/ps        * */
	/* *                 THICK - the thickness of the silicon in um         * */
	/* *                  NMAX - the maximum number of e-hole pairs         * */
	/* *     XEH[2][NMAX].f[4] - the positions of each e and h.  The 2nd    * */
	/* *                         index = 0 for e and 1 for h.               * */
	/* *                   NEH - the actual number of generated pairs       * */
	/* ********************************************************************** */
	
	
	/*  Twice the electron mass in eV */
	
	
    /* Function Body */
	
	/*  The flag for first call */
	
	
    if (fcall) {

		twopi = acos((float)-1.) * (float)2.;
		
		/*  Define some numerical constant vectors */
		
		for (i__ = 0; i__ < 4; ++i__) {
			vzero.f[i__] = (float)0.;
		}
		
        pxvec.v=vzero.v;
        pvvec.v=vzero.v;
        dxvec.v=vzero.v;
        dvvec.v=vzero.v;
        xin.v = vzero.v;
		
#ifdef __POWERPC__
        mask1.c[0]=0x08;
        mask1.c[1]=0x09;
        mask1.c[2]=0x0a;
        mask1.c[3]=0x0b;
        mask1.c[4]=0x00;
        mask1.c[5]=0x01;
        mask1.c[6]=0x02;
        mask1.c[7]=0x03;
        mask1.c[8]=0x04;
        mask1.c[9]=0x05;
        mask1.c[10]=0x06;
        mask1.c[11]=0x07;
        mask1.c[12]=0x0c;
        mask1.c[13]=0x0d;
        mask1.c[14]=0x0e;
        mask1.c[15]=0x0f;
        mask2.c[0]=0x04;
        mask2.c[1]=0x05;
        mask2.c[2]=0x06;
        mask2.c[3]=0x07;
        mask2.c[4]=0x08;
        mask2.c[5]=0x09;
        mask2.c[6]=0x0a;
        mask2.c[7]=0x0b;
        mask2.c[8]=0x00;
        mask2.c[9]=0x01;
        mask2.c[10]=0x02;
        mask2.c[11]=0x03;
        mask2.c[12]=0x0c;
        mask2.c[13]=0x0d;
        mask2.c[14]=0x0e;
        mask2.c[15]=0x0f;
        
        y1.v=vec_perm(bfield.v,vzero.v,mask1.v);
        y2.v=vec_perm(bfield.v,vzero.v,mask2.v);		
#else
		
		y1.v=_mm_set_ps(0,bfield.f[1],bfield.f[0],bfield.f[2]);
		y2.v=_mm_set_ps(0,bfield.f[0],bfield.f[2],bfield.f[1]);
#endif
		
		fcall = 0;
    }
	
/* Allow the field direction to change event to event */

#ifdef __POWERPC__
	
	y1.v=vec_perm(bfield.v,vzero.v,mask1.v);
	y2.v=vec_perm(bfield.v,vzero.v,mask2.v);		
#else
	
	y1.v=_mm_set_ps(0,bfield.f[1],bfield.f[0],bfield.f[2]);
	y2.v=_mm_set_ps(0,bfield.f[0],bfield.f[2],bfield.f[1]);
#endif
	
	
	/* Computing 2nd power */
	r__1 = vect[3];
	/* Computing 2nd power */
	r__2 = vect[4];
	/* Computing 2nd power */
	r__3 = vect[5];
	ppion = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
	
	 e = eloss(ppion, 0.5, &path);
	 pathmf = path;
	
/*  Copy the pion vector and insure normalization */

    for (i__ = 0; i__ < 3; ++i__) {
	pxvec.f[i__] = vect[i__];
    }
/* Computing 2nd power */
    r__1 = vect[3];
/* Computing 2nd power */
    r__2 = vect[4];
/* Computing 2nd power */
    r__3 = vect[5];
    rnorm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
    for (i__ = 0; i__ < 3; ++i__) {
	pvvec.f[i__] = vect[i__ + 3] / rnorm;
    }

	
	/*  The pion direction in polar angles */
	
    costp = pvvec.f[2];
    phip = atan2((double)pvvec.f[1], (double)pvvec.f[0]);
	
	/*  Calculate the entry and exit points of the pion */
	
    if (costp > (float)0.) {
		xin.f[0] = pxvec.f[0] + pvvec.f[0] / pvvec.f[2] * (-pxvec.f[2]);
		xin.f[1] = pxvec.f[1] + pvvec.f[1] / pvvec.f[2] * (-pxvec.f[2]);
		xin.f[2] = (float)0.;
    } else {
		xin.f[0] = pxvec.f[0] + pvvec.f[0] / pvvec.f[2] * (thick - pxvec.f[2]);
		xin.f[1] = pxvec.f[1] + pvvec.f[1] / pvvec.f[2] * (thick - pxvec.f[2]);
		xin.f[2] = thick;
    }
	
	/*  Initialize the pion position vector */
	
    pxvec.v = xin.v;
	
	/*  Now begin stepping through the silicon and generating e-h pairs */
	
    *neh = -1;
L100:
    ranlux_(rvec, &c__3);
    r1 = (float)1. - rvec[0];
    if (r1 < (float)1e-25) {
		r1 = (float)1e-25;
    }
    vsplat.f[0] = -pathmf * log(r1);
#ifdef __POWERPC__
    avec.v = vec_splat(vsplat.v,0);
    pxvec.v=vec_madd(avec.v,pvvec.v,pxvec.v);
#else
    avec.v = _mm_load_ps1(&vsplat.f[0]);
    pxvec.v=_mm_add_ps(_mm_mul_ps(avec.v,pvvec.v),pxvec.v);
#endif
    if (pxvec.f[2] < 0. || pxvec.f[2] > thick) {
		goto L200;
    }
	
	/*  If still in the Si, generate e- hole pairs */
	
	/*  Calculate the kinetic energy of the primary electron */
	
    r2 = rvec[1];
    te = eloss(ppion, r2, &path);
	
	/*  Calculate the mean number of e-hole pairs */
	
    ehnum = te / (float)3.68;
	
	/*  Choose the actual number from a Poisson distribution */
	
    rnpssn_(&ehnum, &npair, &ierr);
	
	/*  If the number is zero, take another step */
	
    if (npair <= 0) {
		goto L100;
    }
	
	/*  Create an e-h pair at the current location of the pion */
	
    if (*neh >= nmax-1) {
		goto L200;
    }
    ++(*neh);
    xeh[0][*neh].v = pxvec.v;
    xeh[1][*neh].v = pxvec.v;
	
	/*  See if there are additional e-hole pairs to create */
	
    ndelta = npair - 1;
    if (ndelta <= 0) {
		goto L100;
    }
	
	/*  Determine the direction of the primary (Delta Ray) wrt the pion */
	
    cost = sqrt(te / (twome + te));
	/* Computing 2nd power */
    r__1 = cost;
    sint = sqrt((float)1. - r__1 * r__1);
    phi = twopi * rvec[2];
	
	/*  Construct the six vector for the delta ray (first wrt to the z axis) */
	
    dvvec.f[2] = cost;
    dvvec.f[1] = sint * sin(phi);
    dvvec.f[0] = sint * cos(phi);
	
	/*  now rotate the z-axis to the pion direction */
	
    newdir(&dvvec, costp, phip);
	
	/*  One end begins at the current pion position */
	
    dxvec.v = pxvec.v;
	
	/*  Calculate the energy/e-h pair */
	
    deltae = te / npair;
	
	/*  Loop over each pair and step along the delta ray track */
	
    i__1 = ndelta;
    for (ipr = 1; ipr <= i__1; ++ipr) {
		
		/*  Calculate the step size along the delta track */
		
	    step = drde(te, new_drde) * deltae;
        vsplat.f[0] = step;
#ifdef __POWERPC__
        avec.v = vec_splat(vsplat.v,0);
        dxvec.v=vec_madd(avec.v,dvvec.v,dxvec.v);
#else
        avec.v = _mm_load_ps1(&vsplat.f[0]);
        dxvec.v=_mm_add_ps(_mm_mul_ps(avec.v,dvvec.v),dxvec.v);
#endif
		
		/* Now add magnetic deflection to the delta track */
		
        pdelta = sqrt(te*te+twome*te);
        if (pdelta > 32000.) {
			
#ifdef __POWERPC__		       
			/*  Calculate the cross product of the velocity and the b-field */
			
			x1.v=vec_perm(dvvec.v,vzero.v,mask1.v);
			x2.v=vec_perm(dvvec.v,vzero.v,mask2.v);
			d.v=vec_madd(x2.v,y1.v,vzero.v);
			vb.v=vec_nmsub(x1.v,y2.v,d.v);
			
			/* update the direction vector  */
			
			vsplat.f[0] = -300.*step/pdelta;
			avec.v = vec_splat(vsplat.v,0);
			dvvec.v=vec_madd(avec.v,vb.v,dvvec.v);
			
			/* update the position vector  */
			
			vsplat.f[0] = vsplat.f[0] * step/2.;
			avec.v = vec_splat(vsplat.v,0);
			dxvec.v=vec_madd(avec.v,vb.v,dxvec.v);
			
			/* renorm the direction vector  */
			
			dvvec2.v = vec_madd(dvvec.v,dvvec.v,vzero.v);
			vsplat.f[0] = 1./sqrt(dvvec2.f[0]+dvvec2.f[1]+dvvec2.f[2]);
			avec.v = vec_splat(vsplat.v,0);
			dvvec.v=vec_madd(avec.v,dvvec.v,vzero.v);
#else
			/*  Calculate the cross product of the velocity and the b-field */
			
			x1.v=_mm_set_ps(0,dvvec.f[1],dvvec.f[0],dvvec.f[2]);
			x2.v=_mm_set_ps(0,dvvec.f[0],dvvec.f[2],dvvec.f[1]);
			d.v=_mm_mul_ps(x2.v,y1.v);
			vb.v=_mm_sub_ps(d.v,_mm_mul_ps(x1.v,y2.v));
			
			/* update the direction vector  */
			
			vsplat.f[0] = -300.*step/pdelta;
			avec.v = _mm_load_ps1(&vsplat.f[0]);
			dvvec.v=_mm_add_ps(_mm_mul_ps(avec.v,vb.v),dvvec.v);
			
			/* update the position vector  */
			
			vsplat.f[0] = vsplat.f[0] * step/2.;
			avec.v = _mm_load_ps1(&vsplat.f[0]);
			dxvec.v=_mm_add_ps(_mm_mul_ps(avec.v,vb.v),dxvec.v);;
			
			/* renorm the direction vector  */
			
			dvvec2.v = _mm_mul_ps(dvvec.v,dvvec.v);
			vsplat.f[0] = 1./sqrt(dvvec2.f[0]+dvvec2.f[1]+dvvec2.f[2]);
			avec.v = _mm_load_ps1(&vsplat.f[0]);
			dvvec.v=_mm_mul_ps(avec.v,dvvec.v);
#endif
			
			/* Add multiple scattering if NIST Estar is chosen */ 
			
			if(new_drde) {
				rumfp = (float)rutherford(te, &etar);
				ranlux_(rvec, &c__1);
				if(rvec[0] < step/rumfp) {
					
					/* Rutherford scatter the delta ray */	
					/*  First get the current direction */
					
					costd = dvvec.f[2];
					phid = atan2((double)dvvec.f[1], (double)dvvec.f[0]);
					
					/* Next, get two random numbers */
					
					ranlux_(rvec, &c__2);
					
					/*  Determine the direction of the scattered electron wrt the incident one */
					
					cost = (rvec[0]-etar+2.*etar*rvec[0])/(rvec[0]+etar);
					/* Computing 2nd power */
					r__1 = cost;
					sint = sqrt((float)1. - r__1 * r__1);
					phi = twopi * rvec[1];
					
					/*  Construct the direction vector for the scattered delta ray (first wrt to the z axis) */
					
					dvvec.f[2] = cost;
					dvvec.f[1] = sint * sin(phi);
					dvvec.f[0] = sint * cos(phi);
					
					/*  now rotate the z-axis to the incident direction */
					
					newdir(&dvvec, costd, phid);
					
				}
				
			}
			
		}
        
		if (dxvec.f[2] <= (float)0. || dxvec.f[2] >= thick) {goto L100;}
		
		/*  Create an e-h pair at the current location of the delta ray */
		
		if (*neh >= nmax-1) {goto L200;}
		++(*neh);
		xeh[0][*neh].v = dxvec.v;
		xeh[1][*neh].v = dxvec.v;
		
		/*  Decrement the energy and loop back */
		
		te -= deltae;
		if(te < 1.) {te = 1. ;}
    }
    goto L100;
L200: (*neh) = (*neh)+1;
    return 0;
} /* deposit */

/* Subroutine */ int propag(float thick, float xsize, float ysize, float temp, float flux[2], float rhe, float rhh, float max_drift_time, float stimstp, int ehole, int neh, vect_or_f xeh[2][NEHSTORE], int indeh[2][NEHSTORE], vect_or_f xhisteh[NCRRC][2][NEHSTORE])
{
	/* Initialized data */
	
	static float rh[2] = { (float)1.12,(float).9 };
	static float qeh[2] = { (float)-1.,(float)1. };
	static float trapc[2] = { (float)5.65e-5,(float)7.70e-5 };
	static double trape[2] = { -0.86, -1.52 };    
	/*  The initial time step (in ps) */
	static float timstp = (float)100.;
	static float timmax = (float)50000.;
	/*  The target step precision in um */
	static double deltax = 0.0025;
	static double deltx2, deltmax, deltinit;
	static double sexp = 0.10;
	static int fcall = -1;
	
	/* System generated locals */
	int i__1;
	float r__1, r__2, r__3;
	double d__1;
	
	/* Local variables */
	static vect_or_f avec, bvec, cvec, tvec;
	static vect_or_f vtwo, vhalf, v[6], x[7];
	static vect_or_f a[6], b1[1], b2[2], b3[3], b4[4], b5[5], c[6], dc[6];
	static int nind;
	static float etot, zbound;
	static int i__, j, j__;
	static float dconv;
	static int nstep, ehcomp;
	static vect_or_f vzero, vzone, bf, bf2, ef, exb, edb, dx[7], difvec, y1, y2;
	static float qe[2], diflen, mu, mu2, zerom2, thickp2, difcon, rnorm, trptim[2];
	float timexp, zproj, rvec[3];
	static int nindeh[2];
	static float bconst, cconst, dconst, econst;
	static int nshift;
	static int ntpart;
	static vect_or_f xgauss, vsixth, xdummy;
	static int ieh, jeh, ind, ntmp;
	static float qrh[2], qrh2[2], rhbf2[2];
	vect_or_f d,x1,x2,vsplat;
#ifdef __POWERPC__		       
	static vect_or_c mask1, mask2;
#endif
	float timeh[2][NEHSTORE], timemx[2][NEHSTORE];
	int itimstp[2][NEHSTORE];
	static double eps, omeps;
	double scalef, dx2;
	
	
	
	/* ********************************************************************** */
	/* * This routine propagates the electron-hole pairs in the applied E   * */
	/* * and B fields until all of them have reached the boundaries of the  * */
	/* * slab (at z=0 and z=THICK um).                                      * */                                 
	/* * Parameters:     THICK - the thickness of the silicon in um         * */
	/* *                 XSIZE - xsize of the pixel                         * */
	/* *                 YSIZE - ysize of the pixel                         * */
	/* *                  TEMP - the temperature of the silicon in K        * */
	/* *                  FLUX - the neq fluence in units of 10^{14}cm^{-2} * */
	/* *                   RHE - the electron Hall factor (if<=0, use 1.12) * */
	/* *                   NEH - the number of electron hole pairs to track * */
	/* *      XEH[2][NEH].f[4] - the positions of each e and h.             * */
	/* *                         The 2nd index = 0 for e and 1 for h.       * */
	/* *         INDEH[2][NEH] - a dummy array used to index particles      * */
	/* ********************************************************************** */
	
	
	/*  The e/h: effective masses (in me units), Hall factors, and charges */
	
	/* Function Body */
	
	
	/*  The flag for first call */
	
	if (fcall) {
		
	   printf("flux_e = %f, flux_h = %f \n", flux[0], flux[1]);
	   if(rhe > 0.) { rh[0] = rhe; }
	   if(rhh > 0.) { rh[1] = rhh; }
	   if(ehole) { ehcomp = 0;} else { ehcomp = 1;}
		if(max_drift_time > 0.) {timmax = max_drift_time;}
		
	   for (i__ = 0; i__ < 2; ++i__) {
			
			/*  Choose units of um and ps */
			
	      qe[i__] = qeh[i__] * (float)1e-6;
			qrh[i__] = rh[i__] * (float)1e-6;
	      qrh2[i__] = qeh[i__] * rh[i__] * rh[i__] * (float)1e-6;
	      rhbf2[i__] = rh[i__] * rh[i__] *(bfield.f[0]*bfield.f[0] 
														+ bfield.f[1]*bfield.f[1] +  bfield.f[2]*bfield.f[2]);
	      d__1 = (double) (temp / 263.);
	      r__1 = (float) pow(d__1, trape[i__]) * flux[i__] * trapc[i__];
			if (r__1 > 0.) {
	         trptim[i__] = (float) (-1.)/r__1;
			} else {
				trptim[i__] = (float) -1.e10;
			}
	   }
		
		/*  square of target precision */
		
	   deltx2= deltax * deltax;
	   deltmax = 4. * deltx2;
	   deltinit = 8. * deltx2;
		
		/*  limits of exponential */
		
	   eps = 1.e-10;
		
		/*  step adjustment boundaries */
		
	   zerom2 = -0.1;
	   thickp2 = thick + 0.1;
		
		
		/*  DCONV = 2kT/e is the constant needed to convert mobility to 2 x diffusion */
		/*  constant DIFCON.  Units are chosen to convert mu [cm**2/(V*s)] into */
		/*  DIFCON [um**2/ps] */
		/*  mu is now delivered in [um**2/(V*ps)] so increase dconv by 10**4 (12/11/03) */
		
		dconv = 2. *  temp * (float)8.617e-5;
		
		/*  Define some numerical constant vectors */
		
	   for (i__ = 0; i__ < 4; ++i__) {
	      vzero.f[i__] = (float)0.;
	      vtwo.f[i__] = (float)2.;
	      vhalf.f[i__] = (float).5;
	      vsixth.f[i__] = (float).16666666666666666;
			/* Constants for adaptive step sizing using Cash-Karp embedded 5th-order technique */
	      b1[0].f[i__] = 1./5.;
	      b2[0].f[i__] = 3./40.;
	      b2[1].f[i__] = 9./40.;
	      b3[0].f[i__] = 3./10.;
	      b3[1].f[i__] =-9./10.;
	      b3[2].f[i__] = 6./5.;
	      b4[0].f[i__] = -11./54.;
	      b4[1].f[i__] = 5./2.;
	      b4[2].f[i__] = -70./27.;
	      b4[3].f[i__] = 35./27.;
	      b5[0].f[i__] = 1631./55296.;
	      b5[1].f[i__] = 175./512.;
	      b5[2].f[i__] = 575./13824.;
	      b5[3].f[i__] = 44275./110592.;
	      b5[4].f[i__] = 253./4096.;
	      c[0].f[i__] = 37./378.;
	      c[1].f[i__] = 0.;
	      c[2].f[i__] = 250./621.;
	      c[3].f[i__] = 125./594.;
	      c[4].f[i__] = 0.;
	      c[5].f[i__] = 512./1771.;
	      dc[0].f[i__] = c[0].f[i__] - 2825./27648.;
	      dc[1].f[i__] = 0.;
	      dc[2].f[i__] = c[2].f[i__] - 18575./48384.;
	      dc[3].f[i__] = c[3].f[i__] - 13525./55296.;
	      dc[4].f[i__] = c[4].f[i__] - 277./14336.;
	      dc[5].f[i__] = c[5].f[i__] - 1./4.;
	   }
		vzone.v = vzero.v;
		vzone.f[3] = 1.;
		xdummy.f[0] = -100. * xsize;
		xdummy.f[1] = -100. * ysize;
		xdummy.f[2] = -100. * thick;
		xdummy.f[3] = 0.;
		xgauss.f[3] = 0.;
		
#ifdef __POWERPC__		              
		mask1.c[0]=0x08;
		mask1.c[1]=0x09;
		mask1.c[2]=0x0a;
		mask1.c[3]=0x0b;
		mask1.c[4]=0x00;
		mask1.c[5]=0x01;
		mask1.c[6]=0x02;
		mask1.c[7]=0x03;
		mask1.c[8]=0x04;
		mask1.c[9]=0x05;
		mask1.c[10]=0x06;
		mask1.c[11]=0x07;
		mask1.c[12]=0x0c;
		mask1.c[13]=0x0d;
		mask1.c[14]=0x0e;
		mask1.c[15]=0x0f;
		mask2.c[0]=0x04;
		mask2.c[1]=0x05;
		mask2.c[2]=0x06;
		mask2.c[3]=0x07;
		mask2.c[4]=0x08;
		mask2.c[5]=0x09;
		mask2.c[6]=0x0a;
		mask2.c[7]=0x0b;
		mask2.c[8]=0x00;
		mask2.c[9]=0x01;
		mask2.c[10]=0x02;
		mask2.c[11]=0x03;
		mask2.c[12]=0x0c;
		mask2.c[13]=0x0d;
		mask2.c[14]=0x0e;
		mask2.c[15]=0x0f;
#endif
		
		
	   fcall = 0;
	   printf("time step = %f, exit propag initialization \n", timstp);
	}
	
#ifdef __POWERPC__		              
	
	y1.v=vec_perm(bfield.v,vzero.v,mask1.v);
	y2.v=vec_perm(bfield.v,vzero.v,mask2.v);
#else
	y1.v=_mm_set_ps(0,bfield.f[1],bfield.f[0],bfield.f[2]);
	y2.v=_mm_set_ps(0,bfield.f[0],bfield.f[2],bfield.f[1]);
#endif
	
	
	/*  Index all particles */
	
	ntmp = -1;
	for (ieh = 0; ieh < neh; ++ieh) {
		if((ieh/Nscale)*Nscale == ieh) {
			
			/* If particle is one of the prescaled, index it (don't index holes/electron if there is no flux/trapping) */
			
			++ntmp;
			indeh[ehole][ntmp] = ieh;
			
			/* also initialize the time steps, and the maximum time */
			
			timeh[ehole][ieh] = timstp;
			ranlux_(rvec, &c__2);
			d__1 = (double) rvec[0];
			d__1 = fmax(d__1, eps);
			timexp = (float) (trptim[ehole] * log(d__1));
			timemx[ehole][ieh] = fminf(timmax, timexp);
			itimstp[ehole][ieh] = 1;
			
			indeh[ehcomp][ntmp] = ieh;
			timeh[ehcomp][ieh] = timstp;			   
			d__1 = (double) rvec[1];
			d__1 = fmax(d__1, eps);
			timexp = (float) (trptim[ehcomp] * log(d__1));
			timemx[ehcomp][ieh] = fminf(timmax, timexp);
			itimstp[ehcomp][ieh] = 1;
		} else {
			
			/* If particle is not one of the prescaled, remove it to a dummy location */
			
			xeh[ehole][ieh].v = xdummy.v;
			xeh[ehcomp][ieh].v = xdummy.v;
			for(i__ = 0; i__ < NCRRC; ++i__) { 
				xhisteh[i__][ehole][ieh].v = xdummy.v; 
				xhisteh[i__][ehcomp][ieh].v = xdummy.v;
			}
		}
	}
	nindeh[ehole] = ntmp+1;
	nindeh[ehcomp] = ntmp+1; 
	
	/*  Use the Cash-Karp Runge-Kutta method to solve the equations of motion for the */
	/*  position of each electron/hole */
	
	/*  Begin the first step in time */
	
	nstep = 0;
L100:
	++nstep;
	
	/*  Loop over the electron/hole index */
	
	for (jeh = 0; jeh < 2; ++jeh) {
		
		/*  See if there are any particles left to propagate */
		
		nind = nindeh[jeh];
		if (nind > 0) {
			
			/*  Loop over the number of particles */
			
			for (ind = 0; ind < nind; ++ind) {
				ieh = indeh[jeh][ind];
				dx2 = deltinit;
				
				/* Calculate the step while the current precision is large */
				
				while(dx2 > deltmax) {
					
					/* Check to see if we must reduce the step size to account for trapping */
					
					if ((xeh[jeh][ieh].f[3]+timeh[jeh][ieh]) > (timemx[jeh][ieh]+0.05)) { 
						timeh[jeh][ieh] = timemx[jeh][ieh] - xeh[jeh][ieh].f[3] + 0.05; 
					} 
					
					/* Check to see if we must reduce the step size to hit the next history point */
					
					if ((xeh[jeh][ieh].f[3]+timeh[jeh][ieh]) > (itimstp[jeh][ieh]*stimstp+0.1)) { 
						timeh[jeh][ieh] = itimstp[jeh][ieh]*stimstp - xeh[jeh][ieh].f[3] + 0.1; 
					} 
					
					/* Load the position into local vector */
					
					x[0].v = xeh[jeh][ieh].v;
					
					/*  Evaluate the E and B fields at the current particle location */
					/*  (which is the first step of six in the improved Runge-Kutta method) */
					
					fieldsav(&x[0], &etot, &ef, &bf, thick, xsize, ysize);
					
					/*  Calculate the mobility */
					
					mu = mobil(jeh, etot, temp);   
					mu2 = mu*mu;             
					
					/* Calculate useful constants */
					
				   rnorm = 1. + rhbf2[jeh] * mu2;
#ifdef __POWERPC__		              
					edb.v = vec_madd(ef.v, bf.v, vzero.v);
					
					/*  Define a useful vectors which are constant over this field/mu step */
					
					vsplat.f[0] = qe[jeh] * mu / rnorm;
					avec.v = vec_splat(vsplat.v,0);
					vsplat.f[0] = qrh[jeh] * mu2 / rnorm;
					bvec.v = vec_splat(vsplat.v,0);
					vsplat.f[0] = qrh2[jeh] * mu * mu2 *(edb.f[0]+edb.f[1]+edb.f[2])/ rnorm;
					cvec.v = vec_splat(vsplat.v,0);
					
					/*  Calculate the cross product of the e-field and the b-field */
					
					x1.v=vec_perm(ef.v,vzero.v,mask1.v);
					x2.v=vec_perm(ef.v,vzero.v,mask2.v);
					d.v=vec_madd(x2.v,y1.v,vzero.v);
					exb.v=vec_nmsub(x1.v,y2.v,d.v);
					
					/*  Calculate the first Runge-Kutta step for x */
					
					v[0].v = vec_madd(avec.v, ef.v, vzone.v);
					v[0].v = vec_madd(bvec.v, exb.v, v[0].v);
					v[0].v = vec_madd(cvec.v, bf.v, v[0].v);
					
#else
					edb.v = _mm_mul_ps(ef.v, bf.v);
					
					/*  Define a useful vectors which are constant over this field/mu step */
					
					vsplat.f[0] = qe[jeh] * mu / rnorm;
					avec.v = _mm_load_ps1(&vsplat.f[0]);
					vsplat.f[0] = qrh[jeh] * mu2 / rnorm;
					bvec.v =  _mm_load_ps1(&vsplat.f[0]);
					vsplat.f[0] = qrh2[jeh] * mu * mu2 *(edb.f[0]+edb.f[1]+edb.f[2])/ rnorm;
					cvec.v =  _mm_load_ps1(&vsplat.f[0]);
					
					/*  Calculate the cross product of the e-field and the b-field */
					
					
					x1.v=_mm_set_ps(0,ef.f[1],ef.f[0],ef.f[2]);
					x2.v=_mm_set_ps(0,ef.f[0],ef.f[2],ef.f[1]);
					d.v=_mm_mul_ps(x2.v,y1.v);
					exb.v=_mm_sub_ps(d.v,_mm_mul_ps(x1.v,y2.v));
					
					/*  Calculate the first Runge-Kutta step for x */
					
					v[0].v = _mm_add_ps(_mm_mul_ps(avec.v, ef.v), vzone.v);
					v[0].v = _mm_add_ps(_mm_mul_ps(bvec.v, exb.v), v[0].v);
					v[0].v = _mm_add_ps(_mm_mul_ps(cvec.v, bf.v), v[0].v);
					
#endif
					
					/*  Test to see if this step will hit a boundary and reduce the step size if it will */
					
				   zproj = x[0].f[2]+timeh[jeh][ieh]*v[0].f[2];
				   if (zproj < zerom2 && v[0].f[2] < 0.) {
				      timeh[jeh][ieh] = (zerom2-x[0].f[2])/v[0].f[2];
					} 
				   if (zproj  > thickp2 && v[0].f[2] > 0.) { 
				      timeh[jeh][ieh] = (thickp2-x[0].f[2])/v[0].f[2]; 
					} 
					
					/* Load the time step into a vector */
					
					vsplat.f[0] = timeh[jeh][ieh];
#ifdef __POWERPC__		              
					tvec.v = vec_splat(vsplat.v,0);
					dx[0].v = vec_madd(v[0].v,tvec.v,vzero.v);
					x[1].v =  vec_madd(dx[0].v,b1[0].v,x[0].v);
#else
					tvec.v = _mm_load_ps1(&vsplat.f[0]);
					dx[0].v = _mm_mul_ps(v[0].v,tvec.v);
					x[1].v =  _mm_add_ps(_mm_mul_ps(dx[0].v,b1[0].v),x[0].v);
#endif
					
					/*  2nd step, evaluate the E and B fields at the current particle location */
					
					fieldsav(&x[1], &etot, &ef, &bf, thick, xsize, ysize);
					
					/*  Calculate the mobility */
					
					mu = mobil(jeh, etot, temp);   
					mu2 = mu*mu;             
					
					/* Calculate useful constants */
					
					rnorm = 1. + rhbf2[jeh] * mu2;
#ifdef __POWERPC__		              
					edb.v = vec_madd(ef.v, bf.v, vzero.v);
					
					/*  Define a useful vectors which are constant over this field/mu step */
					
					vsplat.f[0] = qe[jeh] * mu / rnorm;
					avec.v = vec_splat(vsplat.v,0);
					vsplat.f[0] = qrh[jeh] * mu2 / rnorm;
					bvec.v = vec_splat(vsplat.v,0);
					vsplat.f[0] = qrh2[jeh] * mu * mu2 *(edb.f[0]+edb.f[1]+edb.f[2])/ rnorm;
					cvec.v = vec_splat(vsplat.v,0);
					
					/*  Calculate the cross product of the e-field and the b-field */
					
					x1.v=vec_perm(ef.v,vzero.v,mask1.v);
					x2.v=vec_perm(ef.v,vzero.v,mask2.v);
					d.v=vec_madd(x2.v,y1.v,vzero.v);
					exb.v=vec_nmsub(x1.v,y2.v,d.v);
					
					/*  Calculate the 2nd Runge-Kutta step for x */
					
					v[1].v = vec_madd(avec.v, ef.v, vzone.v);
					v[1].v = vec_madd(bvec.v, exb.v, v[1].v);
					v[1].v = vec_madd(cvec.v, bf.v, v[1].v);
					dx[1].v = vec_madd(v[1].v,tvec.v,vzero.v);
					x[2].v =  vec_madd(dx[0].v,b2[0].v,x[0].v);
					x[2].v =  vec_madd(dx[1].v,b2[1].v,x[2].v);
#else
					edb.v =_mm_mul_ps(ef.v, bf.v);
					
					/*  Define a useful vectors which are constant over this field/mu step */
					
					vsplat.f[0] = qe[jeh] * mu / rnorm;
					avec.v = _mm_load_ps1(&vsplat.f[0]);
					vsplat.f[0] = qrh[jeh] * mu2 / rnorm;
					bvec.v = _mm_load_ps1(&vsplat.f[0]);
					vsplat.f[0] = qrh2[jeh] * mu * mu2 *(edb.f[0]+edb.f[1]+edb.f[2])/ rnorm;
					cvec.v = _mm_load_ps1(&vsplat.f[0]);
					
					/*  Calculate the cross product of the e-field and the b-field */
					
					x1.v=_mm_set_ps(0,ef.f[1],ef.f[0],ef.f[2]);
					x2.v=_mm_set_ps(0,ef.f[0],ef.f[2],ef.f[1]);
					d.v=_mm_mul_ps(x2.v,y1.v);
					exb.v=_mm_sub_ps(d.v,_mm_mul_ps(x1.v,y2.v));
					
					/*  Calculate the 2nd Runge-Kutta step for x */
					
					v[1].v = _mm_add_ps(_mm_mul_ps(avec.v, ef.v), vzone.v);
					v[1].v = _mm_add_ps(_mm_mul_ps(bvec.v, exb.v), v[1].v);
					v[1].v = _mm_add_ps(_mm_mul_ps(cvec.v, bf.v), v[1].v);
					dx[1].v = _mm_mul_ps(v[1].v,tvec.v);
					x[2].v =  _mm_add_ps(_mm_mul_ps(dx[0].v,b2[0].v),x[0].v);
					x[2].v =  _mm_add_ps(_mm_mul_ps(dx[1].v,b2[1].v),x[2].v);
#endif
					
					/*  3rd step, evaluate the E and B fields at the current particle location */
					
					fieldsav(&x[2], &etot, &ef, &bf, thick, xsize, ysize);
					
					/*  Calculate the mobility */
					
					mu = mobil(jeh, etot, temp);   
					mu2 = mu*mu;             
					
					/* Calculate useful constants */
					
					rnorm = 1. + rhbf2[jeh] * mu2;
#ifdef __POWERPC__		              
					edb.v = vec_madd(ef.v, bf.v, vzero.v);
					
					/*  Define a useful vectors which are constant over this field/mu step */
					
					vsplat.f[0] = qe[jeh] * mu / rnorm;
					avec.v = vec_splat(vsplat.v,0);
					vsplat.f[0] = qrh[jeh] * mu2 / rnorm;
					bvec.v = vec_splat(vsplat.v,0);
					vsplat.f[0] = qrh2[jeh] * mu * mu2 *(edb.f[0]+edb.f[1]+edb.f[2])/ rnorm;
					cvec.v = vec_splat(vsplat.v,0);
					
					/*  Calculate the cross product of the e-field and the b-field */
					
					x1.v=vec_perm(ef.v,vzero.v,mask1.v);
					x2.v=vec_perm(ef.v,vzero.v,mask2.v);
					d.v=vec_madd(x2.v,y1.v,vzero.v);
					exb.v=vec_nmsub(x1.v,y2.v,d.v);
					
					/*  Calculate the 3rd Runge-Kutta step for x */
					
					v[2].v = vec_madd(avec.v, ef.v, vzone.v);
					v[2].v = vec_madd(bvec.v, exb.v, v[2].v);
					v[2].v = vec_madd(cvec.v, bf.v, v[2].v);
					dx[2].v = vec_madd(v[2].v,tvec.v,vzero.v);
					x[3].v =  vec_madd(dx[0].v,b3[0].v,x[0].v);
					x[3].v =  vec_madd(dx[1].v,b3[1].v,x[3].v);
					x[3].v =  vec_madd(dx[2].v,b3[2].v,x[3].v);
#else
					edb.v = _mm_mul_ps(ef.v, bf.v);
					
					/*  Define a useful vectors which are constant over this field/mu step */
					
					vsplat.f[0] = qe[jeh] * mu / rnorm;
					avec.v = _mm_load_ps1(&vsplat.f[0]);
					vsplat.f[0] = qrh[jeh] * mu2 / rnorm;
					bvec.v = _mm_load_ps1(&vsplat.f[0]);
					vsplat.f[0] = qrh2[jeh] * mu * mu2 *(edb.f[0]+edb.f[1]+edb.f[2])/ rnorm;
					cvec.v = _mm_load_ps1(&vsplat.f[0]);
					
					/*  Calculate the cross product of the e-field and the b-field */
					
					x1.v=_mm_set_ps(0,ef.f[1],ef.f[0],ef.f[2]);
					x2.v=_mm_set_ps(0,ef.f[0],ef.f[2],ef.f[1]);
					d.v=_mm_mul_ps(x2.v,y1.v);
					exb.v=_mm_sub_ps(d.v,_mm_mul_ps(x1.v,y2.v));
					
					/*  Calculate the 3rd Runge-Kutta step for x */
					
					v[2].v = _mm_add_ps(_mm_mul_ps(avec.v, ef.v), vzone.v);
					v[2].v = _mm_add_ps(_mm_mul_ps(bvec.v, exb.v), v[2].v);
					v[2].v = _mm_add_ps(_mm_mul_ps(cvec.v, bf.v), v[2].v);
					dx[2].v = _mm_add_ps(_mm_mul_ps(v[2].v,tvec.v),vzero.v);
					x[3].v =  _mm_add_ps(_mm_mul_ps(dx[0].v,b3[0].v),x[0].v);
					x[3].v =  _mm_add_ps(_mm_mul_ps(dx[1].v,b3[1].v),x[3].v);
					x[3].v =  _mm_add_ps(_mm_mul_ps(dx[2].v,b3[2].v),x[3].v);
#endif
					
					/*  4th step, evaluate the E and B fields at the current particle location */
					
				   fieldsav(&x[3], &etot, &ef, &bf, thick, xsize, ysize);
					
					/*  Calculate the mobility */
					
					mu = mobil(jeh, etot, temp);   
					mu2 = mu*mu;             
					
					/* Calculate useful constants */
					
				   rnorm = 1. + rhbf2[jeh] * mu2;
#ifdef __POWERPC__		              
					edb.v = vec_madd(ef.v, bf.v, vzero.v);
					
					/*  Define a useful vectors which are constant over this field/mu step */
					
					vsplat.f[0] = qe[jeh] * mu / rnorm;
					avec.v = vec_splat(vsplat.v,0);
					vsplat.f[0] = qrh[jeh] * mu2 / rnorm;
					bvec.v = vec_splat(vsplat.v,0);
					vsplat.f[0] = qrh2[jeh] * mu * mu2 *(edb.f[0]+edb.f[1]+edb.f[2])/ rnorm;
					cvec.v = vec_splat(vsplat.v,0);
					
					/*  Calculate the cross product of the e-field and the b-field */
					
					x1.v=vec_perm(ef.v,vzero.v,mask1.v);
					x2.v=vec_perm(ef.v,vzero.v,mask2.v);
					d.v=vec_madd(x2.v,y1.v,vzero.v);
					exb.v=vec_nmsub(x1.v,y2.v,d.v);
					
					/*  Calculate the 4th Runge-Kutta step for x */
					
					v[3].v = vec_madd(avec.v, ef.v, vzone.v);
					v[3].v = vec_madd(bvec.v, exb.v, v[3].v);
					v[3].v = vec_madd(cvec.v, bf.v, v[3].v);
					dx[3].v = vec_madd(v[3].v,tvec.v,vzero.v);
					x[4].v =  vec_madd(dx[0].v,b4[0].v,x[0].v);
					x[4].v =  vec_madd(dx[1].v,b4[1].v,x[4].v);
					x[4].v =  vec_madd(dx[2].v,b4[2].v,x[4].v);
					x[4].v =  vec_madd(dx[3].v,b4[3].v,x[4].v);
#else
					edb.v = _mm_mul_ps(ef.v, bf.v);
					
					/*  Define a useful vectors which are constant over this field/mu step */
					
					vsplat.f[0] = qe[jeh] * mu / rnorm;
					avec.v = _mm_load_ps1(&vsplat.f[0]);
					vsplat.f[0] = qrh[jeh] * mu2 / rnorm;
					bvec.v = _mm_load_ps1(&vsplat.f[0]);
					vsplat.f[0] = qrh2[jeh] * mu * mu2 *(edb.f[0]+edb.f[1]+edb.f[2])/ rnorm;
					cvec.v = _mm_load_ps1(&vsplat.f[0]);
					
					/*  Calculate the cross product of the e-field and the b-field */
					
					x1.v=_mm_set_ps(0,ef.f[1],ef.f[0],ef.f[2]);
					x2.v=_mm_set_ps(0,ef.f[0],ef.f[2],ef.f[1]);
					d.v=_mm_mul_ps(x2.v,y1.v);
					exb.v=_mm_sub_ps(d.v,_mm_mul_ps(x1.v,y2.v));
					
					/*  Calculate the 4th Runge-Kutta step for x */
					
					v[3].v = _mm_add_ps(_mm_mul_ps(avec.v, ef.v), vzone.v);
					v[3].v = _mm_add_ps(_mm_mul_ps(bvec.v, exb.v), v[3].v);
					v[3].v = _mm_add_ps(_mm_mul_ps(cvec.v, bf.v), v[3].v);
					dx[3].v = _mm_mul_ps(v[3].v,tvec.v);
					x[4].v =  _mm_add_ps(_mm_mul_ps(dx[0].v,b4[0].v),x[0].v);
					x[4].v =  _mm_add_ps(_mm_mul_ps(dx[1].v,b4[1].v),x[4].v);
					x[4].v =  _mm_add_ps(_mm_mul_ps(dx[2].v,b4[2].v),x[4].v);
					x[4].v =  _mm_add_ps(_mm_mul_ps(dx[3].v,b4[3].v),x[4].v);
#endif
					
					/*  5th step, evaluate the E and B fields at the current particle location */
					
				   fieldsav(&x[4], &etot, &ef, &bf, thick, xsize, ysize);
					
					/*  Calculate the mobility */
					
					mu = mobil(jeh, etot, temp);   
					mu2 = mu*mu;             
					
					/* Calculate useful constants */
					
					rnorm = 1. + rhbf2[jeh] * mu2;
#ifdef __POWERPC__		              
					edb.v = vec_madd(ef.v, bf.v, vzero.v);
					
					/*  Define a useful vectors which are constant over this field/mu step */
					
					vsplat.f[0] = qe[jeh] * mu / rnorm;
					avec.v = vec_splat(vsplat.v,0);
					vsplat.f[0] = qrh[jeh] * mu2 / rnorm;
					bvec.v = vec_splat(vsplat.v,0);
					vsplat.f[0] = qrh2[jeh] * mu * mu2 *(edb.f[0]+edb.f[1]+edb.f[2])/ rnorm;
					cvec.v = vec_splat(vsplat.v,0);
					
					/*  Calculate the cross product of the e-field and the b-field */
					
					x1.v=vec_perm(ef.v,vzero.v,mask1.v);
					x2.v=vec_perm(ef.v,vzero.v,mask2.v);
					d.v=vec_madd(x2.v,y1.v,vzero.v);
					exb.v=vec_nmsub(x1.v,y2.v,d.v);
					
					/*  Calculate the 5th Runge-Kutta step for x */
					
					v[4].v = vec_madd(avec.v, ef.v, vzone.v);
					v[4].v = vec_madd(bvec.v, exb.v, v[4].v);
					v[4].v = vec_madd(cvec.v, bf.v, v[4].v);
					dx[4].v = vec_madd(v[4].v,tvec.v,vzero.v);
					x[5].v =  vec_madd(dx[0].v,b5[0].v,x[0].v);
					x[5].v =  vec_madd(dx[1].v,b5[1].v,x[5].v);
					x[5].v =  vec_madd(dx[2].v,b5[2].v,x[5].v);
					x[5].v =  vec_madd(dx[3].v,b5[3].v,x[5].v);
					x[5].v =  vec_madd(dx[4].v,b5[4].v,x[5].v);
#else
					edb.v = _mm_mul_ps(ef.v, bf.v);
					
					/*  Define a useful vectors which are constant over this field/mu step */
					
					vsplat.f[0] = qe[jeh] * mu / rnorm;
					avec.v = _mm_load_ps1(&vsplat.f[0]);
					vsplat.f[0] = qrh[jeh] * mu2 / rnorm;
					bvec.v = _mm_load_ps1(&vsplat.f[0]);
					vsplat.f[0] = qrh2[jeh] * mu * mu2 *(edb.f[0]+edb.f[1]+edb.f[2])/ rnorm;
					cvec.v = _mm_load_ps1(&vsplat.f[0]);
					
					/*  Calculate the cross product of the e-field and the b-field */
					
					x1.v=_mm_set_ps(0,ef.f[1],ef.f[0],ef.f[2]);
					x2.v=_mm_set_ps(0,ef.f[0],ef.f[2],ef.f[1]);
					d.v=_mm_mul_ps(x2.v,y1.v);
					exb.v=_mm_sub_ps(d.v,_mm_mul_ps(x1.v,y2.v));
					
					/*  Calculate the 5th Runge-Kutta step for x */
					
					v[4].v = _mm_add_ps(_mm_mul_ps(avec.v, ef.v), vzone.v);
					v[4].v = _mm_add_ps(_mm_mul_ps(bvec.v, exb.v), v[4].v);
					v[4].v = _mm_add_ps(_mm_mul_ps(cvec.v, bf.v), v[4].v);
					dx[4].v = _mm_mul_ps(v[4].v,tvec.v);
					x[5].v =  _mm_add_ps(_mm_mul_ps(dx[0].v,b5[0].v),x[0].v);
					x[5].v =  _mm_add_ps(_mm_mul_ps(dx[1].v,b5[1].v),x[5].v);
					x[5].v =  _mm_add_ps(_mm_mul_ps(dx[2].v,b5[2].v),x[5].v);
					x[5].v =  _mm_add_ps(_mm_mul_ps(dx[3].v,b5[3].v),x[5].v);
					x[5].v =  _mm_add_ps(_mm_mul_ps(dx[4].v,b5[4].v),x[5].v);
#endif
					
					/*  6th step, evaluate the E and B fields at the current particle location */
					
				   fieldsav(&x[5], &etot, &ef, &bf, thick, xsize, ysize);
					
					/*  Calculate the mobility */
					
					mu = mobil(jeh, etot, temp);   
					mu2 = mu*mu;             
					
					/* Calculate useful constants */
					
					rnorm = 1. + rhbf2[jeh] * mu2;
#ifdef __POWERPC__		              
					edb.v = vec_madd(ef.v, bf.v, vzero.v);
					
					/*  Define a useful vectors which are constant over this field/mu step */
					
					vsplat.f[0] = qe[jeh] * mu / rnorm;
					avec.v = vec_splat(vsplat.v,0);
					vsplat.f[0] = qrh[jeh] * mu2 / rnorm;
					bvec.v = vec_splat(vsplat.v,0);
					vsplat.f[0] = qrh2[jeh] * mu * mu2 *(edb.f[0]+edb.f[1]+edb.f[2])/ rnorm;
					cvec.v = vec_splat(vsplat.v,0);
					
					/*  Calculate the cross product of the e-field and the b-field */
					
					x1.v=vec_perm(ef.v,vzero.v,mask1.v);
					x2.v=vec_perm(ef.v,vzero.v,mask2.v);
					d.v=vec_madd(x2.v,y1.v,vzero.v);
					exb.v=vec_nmsub(x1.v,y2.v,d.v);
					
					/*  Calculate the 6th Runge-Kutta step for x */
					
					v[5].v = vec_madd(avec.v, ef.v, vzone.v);
					v[5].v = vec_madd(bvec.v, exb.v, v[5].v);
					v[5].v = vec_madd(cvec.v, bf.v, v[5].v);
					dx[5].v = vec_madd(v[5].v,tvec.v,vzero.v);
					
					/*  Make the final step */
					
					x[6].v = vec_madd(dx[0].v,c[0].v,x[0].v);
					x[6].v = vec_madd(dx[2].v,c[2].v,x[6].v);
					x[6].v = vec_madd(dx[3].v,c[3].v,x[6].v);
					x[0].v = vec_madd(dx[5].v,c[5].v,x[6].v);
					dx[6].v = vec_madd(dx[0].v,dc[0].v,vzero.v);
					dx[6].v = vec_madd(dx[2].v,dc[2].v,dx[6].v);
					dx[6].v = vec_madd(dx[3].v,dc[3].v,dx[6].v);
					dx[6].v = vec_madd(dx[4].v,dc[4].v,dx[6].v);
					dx[6].v = vec_madd(dx[5].v,dc[5].v,dx[6].v);                  
#else
					edb.v = _mm_mul_ps(ef.v, bf.v);
					
					/*  Define a useful vectors which are constant over this field/mu step */
					
					vsplat.f[0] = qe[jeh] * mu / rnorm;
					avec.v =_mm_load_ps1(&vsplat.f[0]);
					vsplat.f[0] = qrh[jeh] * mu2 / rnorm;
					bvec.v = _mm_load_ps1(&vsplat.f[0]);
					vsplat.f[0] = qrh2[jeh] * mu * mu2 *(edb.f[0]+edb.f[1]+edb.f[2])/ rnorm;
					cvec.v = _mm_load_ps1(&vsplat.f[0]);
					
					/*  Calculate the cross product of the e-field and the b-field */
					
					x1.v=_mm_set_ps(0,ef.f[1],ef.f[0],ef.f[2]);
					x2.v=_mm_set_ps(0,ef.f[0],ef.f[2],ef.f[1]);
					d.v=_mm_mul_ps(x2.v,y1.v);
					exb.v=_mm_sub_ps(d.v,_mm_mul_ps(x1.v,y2.v));
				   
					/*  Calculate the 6th Runge-Kutta step for x */
					
					v[5].v = _mm_add_ps(_mm_mul_ps(avec.v, ef.v), vzone.v);
					v[5].v = _mm_add_ps(_mm_mul_ps(bvec.v, exb.v), v[5].v);
					v[5].v = _mm_add_ps(_mm_mul_ps(cvec.v, bf.v), v[5].v);
					dx[5].v = _mm_mul_ps(v[5].v,tvec.v);
					
					/*  Make the final step */
					
					x[6].v = _mm_add_ps(_mm_mul_ps(dx[0].v,c[0].v),x[0].v);
					x[6].v = _mm_add_ps(_mm_mul_ps(dx[2].v,c[2].v),x[6].v);
					x[6].v = _mm_add_ps(_mm_mul_ps(dx[3].v,c[3].v),x[6].v);
					x[0].v = _mm_add_ps(_mm_mul_ps(dx[5].v,c[5].v),x[6].v);
					dx[6].v = _mm_mul_ps(dx[0].v,dc[0].v);
					dx[6].v = _mm_add_ps(_mm_mul_ps(dx[2].v,dc[2].v),dx[6].v);
					dx[6].v = _mm_add_ps(_mm_mul_ps(dx[3].v,dc[3].v),dx[6].v);
					dx[6].v = _mm_add_ps(_mm_mul_ps(dx[4].v,dc[4].v),dx[6].v);
					dx[6].v = _mm_add_ps(_mm_mul_ps(dx[5].v,dc[5].v),dx[6].v);                  
#endif
					
					/*  dx2 is the estimated square uncertainty on this step */
					
					dx2 = (double) (dx[6].f[0]*dx[6].f[0] + dx[6].f[1]*dx[6].f[1] + dx[6].f[2]*dx[6].f[2]);
					
					/*   Recale the time step to achieve the correct precision on next step */
					
				   if (dx2 > 0.) {
						scalef = pow((deltx2/dx2), sexp);
						scalef = fmin(scalef, 5.);
						scalef = fmax(scalef, 0.2);
				   } else {
						scalef = 5.0;
					}
					
				   timeh[jeh][ieh] *= scalef;
				   
				   if(timeh[jeh][ieh] < 0.1) {
						timeh[jeh][ieh] = 0.1;
						break;
				   } 
				   if(timeh[jeh][ieh] > 250.) {
						timeh[jeh][ieh] = 250.;
						break;
					}
				}
				
				/*  Calculate the diffusion constant and diffusion length for the chosen */
				/*  time step */
				
				difcon = dconv * mu;
				vsplat.f[0] = sqrt(difcon * tvec.f[0]);
#ifdef __POWERPC__		              
				difvec.v = vec_splat(vsplat.v,0);
				
				/*  Add diffusion offsets to all three coordinates */
				
				triplg(&xgauss);
				xeh[jeh][ieh].v = vec_madd(xgauss.v,difvec.v,x[0].v);
#else
				difvec.v = _mm_load_ps1(&vsplat.f[0]);
				
				/*  Add diffusion offsets to all three coordinates */
				
				triplg(&xgauss);
				xeh[jeh][ieh].v = _mm_add_ps(_mm_mul_ps(xgauss.v,difvec.v),x[0].v);
#endif
				
				/*  decide if this particle has reached its maximum drift time, it traps */
				
				if(xeh[jeh][ieh].f[3] >= timemx[jeh][ieh]) {
					indeh[jeh][ind] = -1;
					xhisteh[itimstp[jeh][ieh]-1][jeh][ieh].v = xeh[jeh][ieh].v;
					for(i__ = itimstp[jeh][ieh]; i__ < NCRRC; ++i__) { xhisteh[i__][jeh][ieh].v = xeh[jeh][ieh].v;}
				}
				
				
				/*   Test to see if this particle has left the slab */
				
	        	if (xeh[jeh][ieh].f[2] <= 0. || xeh[jeh][ieh].f[2] >= thick) {
					indeh[jeh][ind] = -1;
					xhisteh[itimstp[jeh][ieh]-1][jeh][ieh].v = xeh[jeh][ieh].v;
					for(i__ = itimstp[jeh][ieh]; i__ < NCRRC; ++i__) { xhisteh[i__][jeh][ieh].v = xeh[jeh][ieh].v;}
				}
				
				/*  If we've reached a time boundary, update the history */
				
				if(((int)(xeh[jeh][ieh].f[3]/stimstp)) == itimstp[jeh][ieh]) {
					xhisteh[itimstp[jeh][ieh]-1][jeh][ieh].v = xeh[jeh][ieh].v;
					if(itimstp[jeh][ieh] < NCRRC) {++itimstp[jeh][ieh];}
				}
				
				
				/*		if (ieh == 100) {            */                       
				/*                    printf("Step %d NIND = %d\n",nstep,nind);         */
				/*                    printf("Etot = %f, e(3) = %f %f %f\n",etot,ef.f[0],ef.f[1],ef.f[2]);         */       
				/*                    printf("%12.5e %12.5e %12.5e %12.5e\n",v[0].f[0],v[0].f[1],v[0].f[2], v[0].f[3]);   */   
				/*                    printf("%12.5e %12.5e %12.5e %12.5e\n",dx[0].f[0],dx[0].f[1],dx[0].f[2], dx[0].f[3]);  */
				/*                    printf("%12.5e %12.5e %12.5e %12.5e\n",v[1].f[0],v[1].f[1],v[1].f[2], v[1].f[3]);    */  
				/*                    printf("%12.5e %12.5e %12.5e %12.5e\n",dx[1].f[0],dx[1].f[1],dx[1].f[2], dx[1].f[3]);  */
				/*                    printf("%12.5e %12.5e %12.5e %12.5e\n",v[2].f[0],v[2].f[1],v[2].f[2], v[2].f[3]);   */   
				/*                    printf("%12.5e %12.5e %12.5e %12.5e\n",dx[2].f[0],dx[2].f[1],dx[2].f[2], dx[2].f[3]);  */
				/*                    printf("%12.5e %12.5e %12.5e %12.5e\n",v[3].f[0],v[3].f[1],v[3].f[2], v[3].f[3]);    */  
				/*                    printf("%12.5e %12.5e %12.5e %12.5e\n",dx[3].f[0],dx[3].f[1],dx[3].f[2], dx[3].f[3]);  */
				/*                    printf("X0 = %12.5e %12.5e %12.5e %12.5e\n",x[0].f[0],           */
				/*                    x[0].f[1],x[0].f[2],x[0].f[3]);                */
				/*                    printf("XEH = %12.5e %12.5e %12.5e %12.5e\n",xeh[jeh][ieh].f[0],        */
				/*                    xeh[jeh][ieh].f[1],xeh[jeh][ieh].f[2],xeh[jeh][ieh].f[3]);  */
				/*					printf("dx2 = %12.5e, scalef = %12.5e, timeh = %12.5e\n", dx2, scalef, timeh[jeh][ieh]); */
				/*					printf("ieh = %d, jeh = %d, timemx = %12.5e\n", ieh, jeh, timemx[jeh][ieh]);      */                                         
				/*		}       */
			}                                                                                                                                       
		}
		
		/*  Check which particles have left the slab and should be deleted */
		/*  from the list */
		
		nindeh[jeh] = 0;
		ind = 0;
	L500:
		if (ind > nind-1) {
			goto L1000;
		}
		if (indeh[jeh][ind] > -1) {
			++nindeh[jeh];
		} else {
			nshift = nind - ind - 1;
			if (nshift > 0) {
				for (i__ = 0; i__ < nshift; ++i__) {
					indeh[jeh][ind + i__] = indeh[jeh][ind + i__ + 1];
				}
			}
			--nind;
			goto L500;
		}
		++ind;
		goto L500;
	L1000:
		;
	}
	
	/*  Determine the total number of particles remaining in the slab */
	
	ntpart = nindeh[0] + nindeh[1];
	
	/*  If any left, continue stepping */
	
	if (ntpart > 0) {
		goto L100;
	}
	return 0;
} /* propag_ */

double mobil(int j, float et, float t)
{
   /* Initialized data */
   
   static double cvm[2] = { 1.53e9,1.62e8 };
   static double evm[2] = { -.87,-.52 };
   static double cec[2] = { 1.01,1.24 };
   static double eec[2] = { 1.55,1.68 };
   static double cbeta[2] = { .0257,.46 };
   static double ebeta[2] = { .66,.17 };
   static int fcall = -1;
   static int jold = -1;
   
   /* System generated locals */
   double ret_val;
   double d__1, d__2, d__3, d__4;
   
   /* Local variables */
   static double beta[2];
   static int i__;
   static double ibeta[2], betaj, ec[2], ibetaj, vm[2], mu0[2], ecj, vmj, mu0j;
   
   
   /* ***************************************************************** */
   /* * This function calculates the mobility of electrons and holes. * */
   /* * Parameters:  J - int index = 0 for e and 1 for h          * */
   /* *             ET - the total electric field in V/m              * */
   /* *              T - the temperature in K                         * */
   /* * Change units to m**2/(V*s) from cm**2/(V*s)                   * */
   /* ***************************************************************** */
   
   
   /*  Initalize the parameters */
   
   if (fcall) {
      for (i__ = 0; i__ < 2; ++i__) {
         d__1 = (double) t;
         vm[i__] = cvm[i__] * pow(d__1, evm[i__]);
         ec[i__] = cec[i__] * pow(d__1, eec[i__]);
         beta[i__] = cbeta[i__] * pow(d__1, ebeta[i__]);
         mu0[i__] = vm[i__] / ec[i__];
         ibeta[i__] = 1. / beta[i__];
      }
      fcall = 0;
   }
   
   /*  Check that J is in bounds */
   
   if (j < 0 || j > 1) {
      ret_val = (float)0.;
      return ret_val;
   }
   
   /*  This routine is called frequently, avoid unnecessary lookups */
   
   if (j != jold) {
      vmj = vm[j];
      ecj = ec[j];
      betaj = beta[j];
      mu0j = mu0[j];
      ibetaj = ibeta[j];
      jold = j;
   }
   
   /*  Calculate the mobility */
   
   d__2 = (double) (et /(100. * ecj));
   d__1 = pow(d__2, betaj) + 1.;
   ret_val = 1.e-4 * mu0j / pow(d__1, ibetaj);
   return ret_val;
} /* mobil */

double eloss(float ppion, float prob, float *path)
{
    /* Initialized data */

    static double emin = 1.8;
    static double emax = 1e6;
    static float ppold = -1.;

    /* System generated locals */
    int i__1;
    double ret_val;

    /* Local variables */
    static int ilow;
    static double xmin, xmax;
    static int i__;
    static double p, coeff;
	static float x;
    static int ihigh;
    static double delog;
    static int ibins, nb;
    static double xb;
    static double xlebin, rinteg[301], xlemin, xlemax, xinteg;


/* *********************************************************************** */
/* * This function returns an energy loss (in eV) for a probability PROB * */
/* * between 0 and 1.  On the first call, the program also calculates    * */
/* * the total mean free path PATH in microns for single scattering      * */
/* * energy loss events.                                                 * */
/* *********************************************************************** */


/*  The minimum and maximum energies */


/*  The flag for first call */


    if (fabs((double)(ppion-ppold)) > 0.005) {

/*  Integrate the correct distribution for each of NBINS upper limits */

	xlemin = log(emin);
	xlemax = log(emax);
	xlebin = (xlemax - xlemin) / 300;
	rinteg[0] = 0.;
	for (ibins = 1; ibins <= 300; ++ibins) {

/*  Work out parameters of each (Simpson's Rule) integral */

	    xinteg = 0.;
	    xmin = xlemin + (ibins - 1) * xlebin;
	    xmax = xmin + xlebin;
	    nb = 20;
	    xb = xlebin / nb;
	    ++nb;
	    i__1 = nb;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		x = (float)(xmin + (i__ - 1) * xb);

/*  DETERMINE THE COEFFICIENT OF THE INTEGRAND (ACCORDING TO SIMPSON'S */
/*  RULE) */
		if (i__ == 1 || i__ == nb) {
		    coeff = 1.;
		} else {
		    coeff = 4.;
		    if (i__ / 2 << 1 != i__) {
			coeff = 2.;
		    }
		}

/*  SUM THE INTEGRAND */
		xinteg += coeff * dsdlne(ppion,x);
	    }

/*  NOW MULTIPLY BY THE NORMALIZATION AND THE BIN WIDTH */

	    xinteg = xb / 3. * xinteg;

/*  NOW MULTIPLY BY THE NORMALIZATION AND THE BIN WIDTH */

	    rinteg[ibins] = rinteg[ibins - 1] + xinteg;
	}

/*  Express the mean free path in microns */

	*path = (float)1. / rinteg[300] * (float)1e4;

/*  Sum the running integral */

	for (ibins = 1; ibins <= 300; ++ibins) {
	    rinteg[ibins] /= rinteg[300];
	}
	ppold = ppion;
    }

/*  Check that the argument is in bounds */

    p = (double)prob;
    if (p < 0.) {
	p = 0.;
    }
    if (p > 1.) {
	p = 1.;
    }

/*  Interpolate the answer */

    for (ibins = 1; ibins <= 300; ++ibins) {
	if (p <= rinteg[ibins]) {
	    ilow = ibins - 1;
	    ihigh = ibins;
	    goto L100;
	}
    }
    ilow = 299;
    ihigh = 300;
L100:
    delog = xlemin + ilow * xlebin + xlebin * (p - rinteg[ilow]) / (rinteg[
	    ihigh] - rinteg[ilow]);

/*  Now return the quantity energy instead of log(energy) */

    ret_val = exp(delog);
    return ret_val;
} /* eloss */

double dsdlne(float ppion, float einlog)
{
    /* Initialized data */

    static float mpion = 0.13957;
	static float twome = 510999.;
	static float kruthr = 178327.;
    static int fcall = -1;

    /* System generated locals */
    double ret_val;

    /* Local variables */
	static double elog;
    static float gamma2, gb, beta2, emax, xruthr, scale, ppold;
    static int ilow, i, j, k;
    static float e;
    static int ihigh, ibins, jbins;
    static float dslog;
    static float engy[1251], gambet[50], xs[1251][50];
	static float xsect[1250];
	static int nbins[50];
    static char extfile[80], inp1[80], inp2[80], inp3[80], inp4[80];
    FILE *extfp;


/* ************************************************************************ */
/* * This routine calculates the product N_A*dsigma/d(ln E) in cm**-1     * */
/* * The function shown in Fig 9 of H. Bichsel, RMP 60, 663 (1988) is     * */
/* * used (valid for 45 GeV pi+- on Si).                                  * */
/* * Parameters: ppion - momentum of incident pion in GeV/c               * */
/* *             EINLOG - the natural logarithm of the lost energy (eV)   * */
/* ************************************************************************ */


/*  The flag for first call */


    if (fcall) {

/* Make sure file is available */

    extfp = fopen("SIRUTH.SPR", "r");
    if (extfp==NULL) {
       printf("cant find SIRUTH.SPR");
       return 1;
    }
    
/* read-in d\sigma/dE/(d\sigma_R/dE) */

       for(i=0; i<46; i+=5){
         fscanf(extfp, "%s %s %s", inp1, inp2, inp3);
		 fscanf(extfp, "%s %f %f %f %f %f", inp4, &gambet[i+0], &gambet[i+1], &gambet[i+2], &gambet[i+3], &gambet[i+4]);
		 for(j=0; j<1251; ++j) {
			fscanf(extfp,"%d %f %f %f %f %f %f", &k, &engy[j], &xs[j][i+0], &xs[j][i+1], &xs[j][i+2],
					&xs[j][i+3], &xs[j][i+4]);
			if(k==99999) break;
		 }
		 nbins[i+0]=j; nbins[i+1]=j; nbins[i+2]=j; nbins[i+3]=j; nbins[i+4]=j;
	   }
        
/* close xsection file */
        
       fclose(extfp);
       ppold = -1.;
       fcall = 0;
    }
        
    if (fabs((double)(ppion-ppold)) > 0.005) {        
        
/* calculate some useful constants */

       gb=ppion/mpion;
	   gamma2=1.+gb*gb;
	   beta2=(gamma2-1.)/gamma2;
	   emax=twome*gb*gb;
        
/* determine which distributions to interpolate */

       if(gb <= gambet[0]) {
         ilow=0;
	     ihigh=0;
	     goto L50;
	     } else if (gb > gambet[49]) {
	     ilow=49;
	     ihigh=49;
	     goto L50;
	  } else {
	     for(i=1; i<50; ++i) {
	       if(gb < gambet[i]) {
		     ihigh = i;
		     ilow = i-1;
	         goto L50;
           }
	     }
	   }
L50:   xruthr = kruthr/beta2;
       for(j=0; j<nbins[ilow]; ++j) {
	     if(ilow == ihigh) {
	        scale = xs[j][ilow];
	     } else {
	        scale = xs[j][ilow]+(xs[j][ihigh]-xs[j][ilow])*(gb-gambet[ilow])/(gambet[ihigh]-gambet[ilow]);
	     }
	  
/* multiply correction by Rutherford x-section ds/d(lnE) = E * ds/dE */

	     xsect[j] = xruthr*(1.-beta2*engy[j]/emax)/engy[j]*scale;
	   }
	   ibins = nbins[ilow];
	   jbins = ibins - 1;
	   ppold = ppion;
	}

/*  Check that the argument is in bounds */

    elog = einlog;
	e = exp(elog);
    if (e < engy[0]) {
	   e = engy[0];
    }
    if (e > engy[jbins]) {
	  e = engy[jbins];
    }

/*  Interpolate the answer */

	for(i=1; i<ibins; ++i) {
	   if(e < engy[i]) {
		  ihigh = i;
		  ilow = i-1;
	      goto L100;
        }
	  }
	
L100: dslog = xsect[ilow]+(xsect[ihigh]-xsect[ilow])*(e - engy[ilow])/(engy[ihigh]-engy[ilow]);

/*  Now return the quantity N_A*dsigma/dlogE */

    ret_val = dslog;
    return ret_val;
} /* dsdlne */


double drde(float e, int new)
{
   /* System generated locals */
   double ret_val, loge, loge2, logdr;
   double d__1;
   
   /* Local variables */
   static float d__;
   
   /* ********************************************************************** */
   /* * This function returns the energy derivative of RANGE in microns/eV * */
   /* * for incident energy E in eV. new > 0 selects the NIST Etar dR/dE.  * */
   /* ********************************************************************** */
   
   d__1 = fabs((double) e);
   
   if(new <= 0) {
      if (d__1 < 5236.) {
         d__ = pow(d__1, 1.216) * 3.118e-9;
      } else if (d__1 < 35320.) {
         d__ = pow(d__1, .912) * 3.635e-8;
      } else {
         d__ = pow(d__1, .509) * 1.951e-6;
      }
      ret_val = d__;
   } else {
      if(d__1 < 1000.) {d__1 = 1000.;}
      if(d__1 > 1000000.) {d__1 = 1000000.;}
      loge = log(d__1);
      loge2 = loge*loge;
      /*	logdr = -18.099 + 1.4344*loge - 0.038823*loge2; */
      logdr = -2.7601 - 3.1528*loge + 0.40799*loge2 - 0.014195*loge*loge2;
      ret_val = exp(logdr);
   }
   return ret_val;
} /* drde */

double rutherford(float e, float *etar)
{
   
   /* ************************************************************************** */
   /* * This function returns Rutherford scattering the mean free path (in um) * */ 
   /* * for an electron of kinetic energy e (in eV). The screening parameter   * */
   /* * etar is also returned for cos(theta_s) selection.                      * */
   /* ************************************************************************** */
   
   static double pi, x, ctf, tau, tau2, alphap, beta, beta2, eta, eta0, rho, third, twothird, const1, const2, const3, xr;
   static double me = 510998.9;
   static double NA = 6.0221415e23;
   static double r02 = 7.9407877e-26;
   static double alpha = 7.297352568e-3;
   static double rhomass = 2.33;
   static double A = 28.0855;
   static double Z = 14.;
   static int fcall = -1;
   if(fcall) {
      pi = acos((float)-1.);
      third = 1./3.;
      twothird = 2.*third;
      x = 9.*pi*pi/128.;
      ctf = pow(x,third);
      rho = 2.33/28.0855*NA;
      const1 = alpha*alpha*pow(Z, twothird)/(4*ctf*ctf);
      const2 = Z*alpha;
      /* include conversion from cm to microns */
      const3 = 1.e4/(rhomass/A*NA*pi*r02*Z*Z);
      fcall = 0;
   }
   tau = e/me;
   tau2 = tau*(2.+tau);
   eta0 = const1/tau2;
   beta2 = 1. - 1./((1.+tau)*(1.+tau));
   beta = sqrt(beta2);
   alphap = const2/beta;
   eta = eta0*(1.13+3.76*alphap*alphap);
   xr = const3*beta2*tau2*eta*(1.+eta);
   *etar = (float)eta;
   return xr;
} /* rutherford */

/* Subroutine */ int newdir(vect_or_f *v, float cost, float phi)
{
   /* Initialized data */
   
   static double coslst = 0.;
   static double philst = 0.;
   static double cphi = 1.;
   static double sphi = 0.;
   
   /* System generated locals */
   double dcost, dphi;
   
   /* Local variables */
   static float sint;
   static int i__;
   static vect_or_f rotcol[3], tmp, vzero;
   vect_or_f vsplat, vt;
   static int fcall = -1;
   
   /* ************************************************** */
   /*  THIS ROUTINE ROTATES THE Z AXIS TO THE DIRECTION */
   /*  ACOS(COST),PHI. */
   /* ************************************************** */
   
   /*  The flag for first call */
   
   
   if (fcall) {
      
      /*  Define some numerical constant vectors */
      
      for (i__ = 0; i__ < 4; ++i__) {
         vzero.f[i__] = (float)0.;
      }
      rotcol[0].f[3] = (float)0.;
      rotcol[1].f[3] = (float)0.;
      rotcol[2].f[3] = (float)0.;
      fcall = 0;
   }
   
   /* Function Body */
   dcost = (double)cost;
   dphi = (double)phi;
   if (dcost != coslst || philst != dphi) {
      coslst = dcost;
      sint = sqrt(1. - dcost * dcost);
      philst = dphi;
      cphi = cos(dphi);
      sphi = sin(dphi);
      rotcol[0].f[0] = (float)(dcost * cphi);
      rotcol[0].f[1] = (float)(dcost * sphi);
      rotcol[0].f[2] = -(float)(sint);
      rotcol[1].f[0] = -(float)(sphi);
      rotcol[1].f[1] = (float)(cphi);
      rotcol[1].f[2] = (float)0.;
      rotcol[2].f[0] = (float)(sint * cphi);
      rotcol[2].f[1] = (float)(sint * sphi);
      rotcol[2].f[2] = cost;
   }
#ifdef __POWERPC__		              
   vt.v = vec_splat(v[0].v,0);
   tmp.v = vec_madd(rotcol[0].v,vt.v,vzero.v);
   vt.v = vec_splat(v[0].v,1);
   tmp.v = vec_madd(rotcol[1].v,vt.v,tmp.v);
   vt.v = vec_splat(v[0].v,2);
   v[0].v = vec_madd(rotcol[2].v,vt.v,tmp.v);
#else
   vt.v = _mm_load_ps1(&v[0].f[0]);
   tmp.v = _mm_mul_ps(rotcol[0].v,vt.v);
   vt.v = _mm_load_ps1(&v[0].f[1]);
   tmp.v = _mm_add_ps(_mm_mul_ps(rotcol[1].v,vt.v),tmp.v);
   vt.v = _mm_load_ps1(&v[0].f[2]);
   v[0].v = _mm_add_ps(_mm_mul_ps(rotcol[2].v,vt.v),tmp.v);
#endif
   return 0;
} /* newdir */

/* Subroutine */ int crrc(float peaktim, float samptim, float *stimstp, float respfcn[NCRRC])
{
	/* Initialized data */
	
	double Q1, Q2, P, t, Dt, h1, h0;
	int i;
   P = peaktim; 
	t = samptim;
	Dt = t/(double)NCRRC;
	for (i = 0; i < NCRRC; ++i) {
		Q1 = Dt*i;
		Q2 = Q1 + Dt;
		if(Q2 <= t) {
//			h1 = exp(-t/P)*((P*exp(t/P) - exp(Q1/P)*(P - Q1 +t))
//			 - (P*exp(t/P) - exp(Q2/P)*(P - Q2 +t)));
            h1 = 1.;
		} else {
			h1 = 0.;
		}
		if(i == 0) {h0 = h1;}
		respfcn[i] = h1/h0;
		printf("respfcn[%lf] = %f \n", Q1, respfcn[i]);
	}
	*stimstp = (float)Dt;
	return 0;
} /* crrc */

/* Subroutine */ int triplg(vect_or_f *x)
{
   /* Initialized data */
   
   static int fcall = -1;
   
   /* Local variables */
   static float rvec[120];
   static int i__;
   static float r__;
   static int ibase;
   static float rbuff[120], twopi;
   static float arg, phi;
   
   
   /* **************************************************************** */
   /* * This routine calculates 3 gaussianly-distributed random nums * */
   /* * Parameters:  X(3) - 3 random numbers (rms = 1.)              * */
   /* **************************************************************** */
   
   /* Function Body */
   
   /*  Initalize the parameters */
   
   if (fcall) {
      twopi = 2.*acos((float)-1.);
      ibase = 120;
      fcall = 0;
   }
   
   /*  If all random numbers used up, generate 120 more */
   
   if (ibase == 120) {
      ranlux_(rvec, &c__120);
      for (i__ = 1; i__ <= 119; i__ += 2) {
         arg = 1. - rvec[i__ - 1];
         if (arg < (float)1e-30) {
            arg = (float)1e-30;
         }
         r__ = sqrt(log(arg) * (float)-2.);
         phi = twopi * rvec[i__];
         rbuff[i__ - 1] = r__ * cos(phi);
         rbuff[i__] = r__ * sin(phi);
      }
      ibase = 0;
   }
   for (i__ = 1; i__ <= 3; ++i__) {
      (*x).f[i__-1] = rbuff[ibase + i__ - 1];
   }
   ibase += 3;
   return 0;
} /* triplg */

/* Subroutine */ int detect(float xsize, float ysize, float thick, int ehole, vect_or_f xhisteh[NCRRC][2][NEHSTORE], int neh, float respfcn[NCRRC], float pixel[TXSIZE][TYSIZE], float pixhist[NCRRC+1][TXSIZE][TYSIZE])
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int j, ix, iy, ieh, jhist, ixmax, ixmin, iymax, iymin, i, k, ixm, iym;
    static float qmax;
    static float qeh[2];
	// DS - Changed from T_SIZE/2 to the following to have offset at center of pixel even for even number of pixels.
	static float xoffset = (float)((TXSIZE+1)/2)-0.5;
	static float yoffset = (float)((TYSIZE+1)/2)-0.5;
	float wgtmat[3][3];
	vect_or_f xp;



/* ********************************************************************** */
/* * This routine determines the signal in a 21X7 array of pixels where * */
/* * the pixel (11,4) is centered on the coordinates (0,0).             * */
/* * Parameters: XSIZE - x-dimension (in um) of the pixels              * */
/* *             YSIZE - y-dimension (in um) of the pixels              * */
/* *                 THICK - the thickness of the silicon in um         * */
/* *                 ehole - 0 for electrons, 1 for holes               * */
/* *	 XEH[2][NMAX].f[4] - the positions and velocities of each e and * */
/* *                         h after propagation to the pixel ends      * */
/* *                   NEH - the actual number of generated pairs       * */
/* *	   PIXEL[2][21][13] - the number of detected electrons or holes * */
/* ********************************************************************** */


/*  Zero pixel array */

    /* Function Body */
	for (jhist = 0; jhist < NCRRC+1; ++jhist) {
	   for (iy = 0; iy < TYSIZE; ++iy) {
		   for (ix = 0; ix < TXSIZE; ++ix) {
			   pixhist[jhist][ix][iy] = 0.;
			}
		}
	}
	
/* Set-up charges and anti-charges, define positive signal for collected carrier */
 
	if(ehole) { qeh[0] = -1.; qeh[1] = 1.;} else { qeh[0] = 1.; qeh[1] = -1.;}
	
/*  Loop over the histories of each particle */
	
	for (jhist = 0; jhist < NCRRC; ++jhist) {
		
/*  Calculate the coordinates of each electron/hole and see if it was collected */

       for (ieh = 0; ieh < neh; ++ieh) {
           if (xhisteh[jhist][ehole][ieh].f[2] <= 0. || xhisteh[jhist][ehole][ieh].f[2] >= thick) {
              ix = (int) (xhisteh[jhist][ehole][ieh].f[0]/xsize + xoffset);
              iy = (int) (xhisteh[jhist][ehole][ieh].f[1]/ysize + yoffset);
              if (ix > -1 && ix < TXSIZE && iy > -1 && iy < TYSIZE) {
                 ++pixhist[jhist+1][ix][iy];
              }
           }
       }    

/*  Loop over all trapped charge and calculate the charge induced on the n+ implants */

       for (j=0; j<2; ++j) {        
          for (ieh = 0; ieh < neh; ++ieh) {   
             if (xhisteh[jhist][j][ieh].f[2] > 0. && xhisteh[jhist][j][ieh].f[2] < thick) { 
                ix = (int) (xhisteh[jhist][j][ieh].f[0]/xsize + xoffset);  
                iy = (int) (xhisteh[jhist][j][ieh].f[1]/ysize + yoffset);  
                if (ix > -1 && ix < TXSIZE && iy > -1 && iy < TYSIZE) { 
                 xp.v = xhisteh[jhist][j][ieh].v;
/* hole collection uses the same e- weighting potential but with complementary z position */
                 if(ehole) {xp.f[2] = thick - xp.f[2];}
                 wgtpotsav(&xp, wgtmat, thick, xsize, ysize);
                 ixmin = ix-1;
                 ixmax = ix+1;
                 iymin = iy-1;
                 iymax = iy+1;
                 for (i=ixmin; i<=ixmax; ++i) {
                   for(k=iymin; k<=iymax; ++k) {
                      if (i > -1 && i < TXSIZE && k > -1 && k < TYSIZE) { pixhist[jhist+1][i][k] += wgtmat[k-iymin][i-ixmin]*qeh[j];}
                   }
                 }
               } 
            }    
          }       
        }
		/* End this time slice */
	}	
	
	/* Now, include the preamp response function */
	
	ixm = -1;
	iym = -1;
	qmax = 0.f;
	for (iy = 0; iy < TYSIZE; ++iy) {
		for (ix = 0; ix < TXSIZE; ++ix) {
		   pixel[ix][iy] = 0.;
		   for (jhist = 0; jhist < NCRRC; ++jhist) {
			   pixel[ix][iy] += respfcn[jhist]*(pixhist[jhist+1][ix][iy] - pixhist[jhist][ix][iy]);
		   }
/*		   if(pixel[ix][iy] > qmax) {
		     qmax = pixel[ix][iy];
		     ixm = ix;
		     iym = iy;
		   } */
		}
	}
	
/*    printf("\n Time vs current for pixel with max signal = %f \n", qmax);
    for (jhist = 0; jhist < NCRRC; ++jhist) {
		printf(" step %d, current %f \n", jhist, (pixhist[jhist+1][ixm][iym] - pixhist[jhist][ixm][iym]));
	}
	
    printf("\n Time vs current for pixel with max neighbor signal \n");
    for (jhist = 0; jhist < NCRRC; ++jhist) {
		printf(" step %d, current %f \n", jhist, (pixhist[jhist+1][ixm][iym+1] - pixhist[jhist][ixm][iym+1]));
	} */
	
    return 0;
} /* detect */

/* Subroutine */ int fieldsav(vect_or_f *xp, float *etot, vect_or_f *ef, vect_or_f *bf, float thick, float xsize, float ysize)
{
   /* Initialized data */
   
   static int fcall = -1;
   
   /* Local variables */
   static int i,j;
#ifndef __POWERPC__		              
   static int mm__;
#endif
   static unsigned int ix[3][4];
   static float tmfive;
   static vect_or_f igrid, vzero, vhalf, pixel, pgrid, vfloor, xshift;
   static vect_or_f x, ax, xind, x0ind, dxind, adxind;
   static vect_or_i ind0;
   static vect_or_f eft, ef0, def1, def2, def3, alpha1, alpha2, alpha3;
   
   /* ***************************************************************** */
   /* * This routine calculates the electric and magnetic fields at   * */
   /* * position X.                                                   * */
   /* * Parameters: XP.f[4] - position of the particle (in microns)   * */
   /* *             EF.f[4] - the electric field in V/cm              * */
   /* *             BF.f[4] - the magnetic field in Tesla             * */
   /* ***************************************************************** */
   
   /* Function Body */
   
   /*  Initalize the parameters */
   
   if (fcall) {
      
      /*  Set up vector with inverse grid pitch */
      
      igrid.f[0] = (float)(npixx-1)*2./xsize;
      igrid.f[1] = (float)(npixy-1)*2./ysize;
      igrid.f[2] = (float)(npixz-1)/thick;
      igrid.f[3] = (float)0.;
      
      /*  Set up vector with inverse pixel dimensions */
      
      pgrid.f[0] = 1./xsize;
      pgrid.f[1] = 1./ysize;
      pgrid.f[2] = 0.;
      pgrid.f[3] = 0.;
      
      /*  Set up vector with pixel dimensions */
      
      pixel.f[0] = -xsize;
      pixel.f[1] = -ysize;
      pixel.f[2] = 0.;
      pixel.f[3] = 0.;
      
      /*  Useful constant vectors */
      
      vzero.f[0] = 0.;
      vzero.f[1] = 0.;
      vzero.f[2] = 0.;
      vzero.f[3] = 0.;
      
      vhalf.f[0] = 0.5;
      vhalf.f[1] = 0.5;
      vhalf.f[2] = 0.5;
      vhalf.f[3] = 0.5;
      
      /*  Useful constant */
      
      tmfive = thick - 5.;
      
      fcall = 0;
   }
   
   /* Load the constant B-field */
   
   (*bf).v = bfield.v;
   
   /* Make a local copy and check the z-coordinate */
   
   x.v = (*xp).v;
   x.f[2] = fmaxf(x.f[2], 5.);
   x.f[2] = fminf(x.f[2], tmfive);
   
   /* xshift contains the coordinates modulo 1 pixel */
   
#ifdef __POWERPC__		              
   vfloor.v = vec_floor(vec_madd(x.v,pgrid.v,vhalf.v));
   xshift.v = vec_madd(vfloor.v,pixel.v,x.v);
   
   /* Find the floating int coordinates of the point x */
   /* First take the abs of x to fold it into the 1st quadrant */
   
   ax.v=vec_abs(xshift.v);
   
   /* Find x in lattice spacing units */
   
   xind.v=vec_madd(ax.v,igrid.v,vzero.v);
   
   /* Find coordinates of nearest point */
   
   x0ind.v=vec_round(xind.v);
   
   /* Find the displacement vector from nearest point */
   
   dxind.v=vec_sub(xind.v,x0ind.v);
   adxind.v=vec_abs(dxind.v);
   
   /* Convert nearest point to int coordinates */
   
   ind0.v=vec_ctu(x0ind.v,0);
#else
   vfloor.v= _mm_add_ps(_mm_mul_ps(x.v,pgrid.v),vhalf.v);	
   for(mm__=0;mm__<4;mm__++){
      vfloor.f[mm__]=floorf(vfloor.f[mm__]);
   }
   xshift.v = _mm_add_ps(_mm_mul_ps(vfloor.v,pixel.v),x.v);
   
   
   /* Find the floating int coordinates of the point x */
   /* First take the abs of x to fold it into the 1st quadrant */
   
   for(mm__=0;mm__<4;mm__++){
      ax.f[mm__]=fabs(xshift.f[mm__]);
   }
   
   /* Find x in lattice spacing units */
   
   xind.v= _mm_mul_ps(ax.v,igrid.v);
   
   /* Find coordinates of nearest point */
   for(mm__=0;mm__<4;mm__++){
      x0ind.f[mm__]=round(xind.f[mm__]);
   }
   
   /* Find the displacement vector from nearest point */
   
   dxind.v=_mm_sub_ps(xind.v,x0ind.v);
   for(mm__=0;mm__<4;mm__++){
      adxind.f[mm__]=fabs(dxind.f[mm__]);
   }
   
   /* Convert nearest point to int coordinates */
   
   for(mm__=0;mm__<4;mm__++){
      if(x0ind.f[mm__]>0)			
         ind0.i[mm__]=(unsigned int) x0ind.f[mm__];
      else ind0.i[mm__]=0;
   }
#endif
   
   /* unload coordinates and work-out coordinates of three next nearest points */
   
   for (j = 0; j < 4; ++j) {
      ix[0][j]=ind0.i[0];
      ix[1][j]=ind0.i[1];
      ix[2][j]=ind0.i[2];
   }
   for (j = 1; j < 4; ++j) {
      if(dxind.f[j-1] > 0.) {
         if(ix[j-1][j] < mnode[j-1]) ++ix[j-1][j];
      } else {
         if(ix[j-1][j] > 0) --ix[j-1][j];
      }
   }
   
   /* evaluate field at nearest and three next nearest points */
   
   ef0.v = efield[ix[0][0]][ix[1][0]][ix[2][0]].v;
   
#ifdef __POWERPC__		              
   def1.v = vec_sub(efield[ix[0][1]][ix[1][1]][ix[2][1]].v,ef0.v);
   def2.v = vec_sub(efield[ix[0][2]][ix[1][2]][ix[2][2]].v,ef0.v);
   def3.v = vec_sub(efield[ix[0][3]][ix[1][3]][ix[2][3]].v,ef0.v);
   
   /* prepare the coefficients of each term */
   
   alpha1.v = vec_splat(adxind.v,0);
   alpha2.v = vec_splat(adxind.v,1);
   alpha3.v = vec_splat(adxind.v,2);
   
   /* do the arithmetic */
   
   eft.v = vec_madd(def1.v,alpha1.v,ef0.v);
   eft.v = vec_madd(def2.v,alpha2.v,eft.v);
   eft.v = vec_madd(def3.v,alpha3.v,eft.v);
#else
   def1.v = _mm_sub_ps(efield[ix[0][1]][ix[1][1]][ix[2][1]].v,ef0.v);
   def2.v = _mm_sub_ps(efield[ix[0][2]][ix[1][2]][ix[2][2]].v,ef0.v);
   def3.v = _mm_sub_ps(efield[ix[0][3]][ix[1][3]][ix[2][3]].v,ef0.v);
   
   /* prepare the coefficients of each term */
   
   alpha1.v = _mm_load_ps1(&adxind.f[0]);
   alpha2.v = _mm_load_ps1(&adxind.f[1]);
   alpha3.v = _mm_load_ps1(&adxind.f[2]);
   
   /* do the arithmetic */
   
   eft.v = _mm_add_ps(_mm_mul_ps(def1.v,alpha1.v),ef0.v);
   eft.v = _mm_add_ps(_mm_mul_ps(def2.v,alpha2.v),eft.v);
   eft.v = _mm_add_ps(_mm_mul_ps(def3.v,alpha3.v),eft.v);
#endif
   
   /* do the necessary reflections in x and y */
   
   if(xshift.f[0] < 0.){
      eft.f[0] = -eft.f[0];
   }
   if(xshift.f[1] < 0.){
      eft.f[1] = -eft.f[1];
   }
   
   /* copy answer to output address */
   
   (*ef).v = eft.v;
   
   (*etot) = sqrt(eft.f[0]*eft.f[0] + eft.f[1]*eft.f[1] + eft.f[2]*eft.f[2]);
   
   return 0;
} /* fieldsav */

/* Subroutine */ int wgtpotsav(vect_or_f *xp, float wgtmat[3][3], float thick, float xsize, float ysize)
{
   /* Initialized data */
   
   static int fcall = -1;
   
   /* Local variables */
   static int i,j,k,l;
#ifndef __POWERPC__		              
   static int mm__;
#endif
   static unsigned int ix[3][4];
   static vect_or_f igrid, vzero, vhalf, pixel, pgrid, vfloor, xshift;
   static vect_or_f x, ax, xind, x0ind, dxind, adxind;
   static vect_or_i ind0;
   static vect_or_f wgt[3], wg0[3], def1[3], def2[2], def3[3], alpha1, alpha2, alpha3;
   
   
   /* ***************************************************************** */
   /* * This routine calculates the weighting potential at position x * */
   /* * Parameters: XP.f[4] - position of the particle (in microns)   * */
   /* *             wgtmat.f[4] - the 3x3 weighting potential matrix  * */
   /* *             thick - the sensor thickness in um                * */
   /* *             xsize - the sensor size in x in um                * */
   /* *             ysize - the sensor size in y in um                * */
   /* ***************************************************************** */
   
   /* Function Body */
   
   /*  Initalize the parameters */
   
   if (fcall) {
      
      /*  Set up vector with inverse grid pitch */
      
      igrid.f[0] = (float)(npixx-1)*2./xsize;
      igrid.f[1] = (float)(npixy-1)*2./ysize;
      igrid.f[2] = (float)(npixz-1)/thick;
      igrid.f[3] = (float)0.;
      
      /*  Set up vector with inverse pixel dimensions */
      
      pgrid.f[0] = 1./xsize;
      pgrid.f[1] = 1./ysize;
      pgrid.f[2] = 0.;
      pgrid.f[3] = 0.;
      
      /*  Set up vector with pixel dimensions */
      
      pixel.f[0] = -xsize;
      pixel.f[1] = -ysize;
      pixel.f[2] = 0.;
      pixel.f[3] = 0.;
      
      /*  Useful constant vectors */
      
      vzero.f[0] = 0.;
      vzero.f[1] = 0.;
      vzero.f[2] = 0.;
      vzero.f[3] = 0.;
      
      vhalf.f[0] = 0.5;
      vhalf.f[1] = 0.5;
      vhalf.f[2] = 0.5;
      vhalf.f[3] = 0.5;
      
      fcall = 0;
   }
   
   
   /* Make a local copy */
   
   x.v = (*xp).v;
   
   /* xshift contains the coordinates modulo 1 pixel */
   
#ifdef __POWERPC__		              
   vfloor.v = vec_floor(vec_madd(x.v,pgrid.v,vhalf.v));
   xshift.v = vec_madd(vfloor.v,pixel.v,x.v);
   
   /* Find the floating int coordinates of the point x */
   /* First take the abs of x to fold it into the 1st quadrant */
   
   ax.v=vec_abs(xshift.v);
   
   /* Find x in lattice spacing units */
   
   xind.v=vec_madd(ax.v,igrid.v,vzero.v);
   
   /* Find coordinates of nearest point */
   
   x0ind.v=vec_round(xind.v);
   
   /* Find the displacement vector from nearest point */
   
   dxind.v=vec_sub(xind.v,x0ind.v);
   adxind.v=vec_abs(dxind.v);
   
   /* Convert nearest point to int coordinates */
   
   ind0.v=vec_ctu(x0ind.v,0);
#else
   vfloor.v= _mm_add_ps(_mm_mul_ps(x.v,pgrid.v),vhalf.v);	
   for(mm__=0;mm__<4;mm__++){
      vfloor.f[mm__]=floorf(vfloor.f[mm__]);
   }
   xshift.v = _mm_add_ps(_mm_mul_ps(vfloor.v,pixel.v),x.v);
   
   
   /* Find the floating int coordinates of the point x */
   /* First take the abs of x to fold it into the 1st quadrant */
   
   for(mm__=0;mm__<4;mm__++){
      ax.f[mm__]=fabs(xshift.f[mm__]);
   }
   
   /* Find x in lattice spacing units */
   
   xind.v= _mm_mul_ps(ax.v,igrid.v);
   
   /* Find coordinates of nearest point */
   for(mm__=0;mm__<4;mm__++){
      x0ind.f[mm__]=round(xind.f[mm__]);
   }
   
   /* Find the displacement vector from nearest point */
   
   dxind.v=_mm_sub_ps(xind.v,x0ind.v);
   for(mm__=0;mm__<4;mm__++){
      adxind.f[mm__]=fabs(dxind.f[mm__]);
   }
   
   /* Convert nearest point to int coordinates */
   
   for(mm__=0;mm__<4;mm__++){
      if(x0ind.f[mm__]>0)			
         ind0.i[mm__]=(unsigned int) x0ind.f[mm__];
      else ind0.i[mm__]=0;
   }
#endif
   
   /* unload coordinates and work-out coordinates of three next nearest points */
   
   for (j = 0; j < 4; ++j) {
      ix[0][j]=ind0.i[0];
      ix[1][j]=ind0.i[1];
      ix[2][j]=ind0.i[2];
   }
   for (j = 1; j < 4; ++j) {
      if(dxind.f[j-1] > 0.) {
         if(ix[j-1][j] < mnode[j-1]) ++ix[j-1][j];
      } else {
         if(ix[j-1][j] > 0) --ix[j-1][j];
      }
   }
   
   /* evaluate field at nearest and three next nearest points */
   
   for(k=0; k<3; ++k) {wg0[k].v = wgtpot[ix[0][0]][ix[1][0]][ix[2][0]][k].v;}
   
#ifdef __POWERPC__	
   for(k=0; k<3; ++k) {	              
      def1[k].v = vec_sub(wgtpot[ix[0][1]][ix[1][1]][ix[2][1]][k].v,wg0[k].v);
      def2[k].v = vec_sub(wgtpot[ix[0][2]][ix[1][2]][ix[2][2]][k].v,wg0[k].v);
      def3[k].v = vec_sub(wgtpot[ix[0][3]][ix[1][3]][ix[2][3]][k].v,wg0[k].v);
   }
   
   /* prepare the coefficients of each term */
   
   alpha1.v = vec_splat(adxind.v,0);
   alpha2.v = vec_splat(adxind.v,1);
   alpha3.v = vec_splat(adxind.v,2);
   
   /* do the arithmetic */
   
   for(k=0; k<3; ++k) {
      wgt[k].v = vec_madd(def1[k].v,alpha1.v,wg0[k].v);
      wgt[k].v = vec_madd(def2[k].v,alpha2.v,wgt[k].v);
      wgt[k].v = vec_madd(def3[k].v,alpha3.v,wgt[k].v);
   }
#else
   for(k=0; k<3; ++k) {
      def1[k].v = _mm_sub_ps(wgtpot[ix[0][1]][ix[1][1]][ix[2][1]][k].v,wg0[k].v);
      def2[k].v = _mm_sub_ps(wgtpot[ix[0][2]][ix[1][2]][ix[2][2]][k].v,wg0[k].v);
      def3[k].v = _mm_sub_ps(wgtpot[ix[0][3]][ix[1][3]][ix[2][3]][k].v,wg0[k].v);
   }
   
   /* prepare the coefficients of each term */
   
   alpha1.v = _mm_load_ps1(&adxind.f[0]);
   alpha2.v = _mm_load_ps1(&adxind.f[1]);
   alpha3.v = _mm_load_ps1(&adxind.f[2]);
   
   /* do the arithmetic */
   
   for(k=0; k<3; ++k) {
      wgt[k].v = _mm_add_ps(_mm_mul_ps(def1[k].v,alpha1.v),wg0[k].v);
      wgt[k].v = _mm_add_ps(_mm_mul_ps(def2[k].v,alpha2.v),wgt[k].v);
      wgt[k].v = _mm_add_ps(_mm_mul_ps(def3[k].v,alpha3.v),wgt[k].v);
   }
#endif
   
   /* do quadrant based unfolding into the 3x3 weighting potential matrix */
   
   for(i=0; i<3; ++i) {
      for(j=0; j<3; ++j) {
         if(xshift.f[0] >= 0.) {l=j;} else {l=2-j;}
         if(xshift.f[1] >= 0.) {k=2-i;} else {k=i;}
         wgtmat[i][j] = wgt[k].f[l];
      }
   }
   
   return 0;
} /* wgtpotsav */


/* Subroutine */ int pixinit(float *pimom, float *thick, float *xsize, float *ysize, float *temp, float flux[2], float *rhe, float *rhh, float *peaktim, float *samptim, int *ehole, int *new_drde, int *filebase)
{
	/* Initialize everything from external file */
	
    /* Local variables */
    static int  i, j, ix, iy, iz, idumx, idumy, idumz;
    static float tdummy, xdummy, ydummy, wdum[3][3];
	char c;
	FILE *ifp;
	static vect_or_f vzero, vhund;
	
	/* **************************************************************** */
	/* * This routine initializes several parameters and the electric * */
	/* * field array from an external file pixel.init                 * */
	/* * Parameters: thick - thickness of the detector (in microns)   * */
	/* *             xsize - x-size of each pixel (in microns)        * */
	/* *             ysize - y-size of each pixel (in microns)        * */
	/* *              temp - the silicon temperature (in deg K)       * */
	/* *              flux - flux of hadrons (in 10**14 h/cm**2)      * */
	/* *               rhe - electron Hall factor (if <=0, use 1.12)  * */
	/* *               rhh - electron Hall factor (if <=0, use 0.90)  * */
	/* *           peaktim - CRRC peaking time                        * */
	/* *           samptim - full drift sampling time                 * */
	/* *             ehole - flag to select e (0) or (h) signal calc  * */
	/* *          new_drde - flag to select NIST Estar drde (>0)      * */
	/* *          filebase - the base index of all output files       * */
	/* **************************************************************** */
	
	/* Function Body */
	
	/*  Useful constant vectors */
	
	vzero.f[0] = 0.;
	vzero.f[1] = 0.;
	vzero.f[2] = 0.;
	vzero.f[3] = 0.;
	
	vhund.f[0] = 100.;
	vhund.f[1] = 100.;
	vhund.f[2] = 100.;
	vhund.f[3] = 100.;
	
	/*  Initalize the parameters */
	
	ifp = fopen("ppixel2.init", "r");
	if (ifp==NULL) {
      printf("no ppixel2.init initialization file found/n");
      return 0;
	}
	
	/* Read-in a header string first and print it*/    
	
	for (i=0; (c=getc(ifp)) != '\n'; ++i) {
		header[i] = c;
	}
	printf("%s\n", header);
	
	/* next, the random number seed and the hour of program shutdown */    
	
	fscanf(ifp,"%f %d", pimom, filebase);
	
	printf("pion momentum (< 1.1 defauts to 45 GeV) = %f, file base index = %d\n", *pimom, *filebase);
	
	/* next, the magnetic field */    
	
	fscanf(ifp,"%f %f %f", &bfield.f[0], &bfield.f[1], &bfield.f[2]);
	bfield.f[3] = 0.;
	bfield_z = bfield.f[2];
	
	printf("bfield = %f, %f, %f, %f\n", bfield.f[0], bfield.f[1], bfield.f[2], bfield.f[3]);
	
	/* next a bunch of parameters */    
	
	fscanf(ifp,"%f %f %f %f %f %f %f %f %f %f %d %d %d %d %d", thick, xsize, ysize, temp, &flux[0], &flux[1], rhe, rhh, peaktim, samptim, ehole, new_drde, &npixx, &npixy, &npixz);
	if((*ehole) !=0) { *ehole = 1;}
	
	printf("thickness = %f, x/y sizes = %f/%f, temp = %f, flux_e = %f, flux_h = %f, rhe = %f, rhh = %f, peaktim = %f, samptim = %f, ehole = %d, new_drde = %d \n", 
	       *thick, *xsize, *ysize, *temp, flux[0], flux[1], *rhe, *rhh, *peaktim, *samptim, *ehole, *new_drde);
	printf("x/y/z field array dimensions = %d/%d/%d\n", npixx, npixy, npixz);
   if(npixx > NARRAYX || npixy > NARRAYY || npixz > NARRAYZ) {printf("Efield array too large, resize arrays \n"); return 0;}
	/* Set up boundary limit array */
	mnode[0] = npixx-1; mnode[1] = npixy-1; mnode[2] = npixz-1;
	
	/* Now do the E-field */
	
	for (iz = 0; iz < npixz; ++iz) {
      for (iy = 0; iy < npixy; ++iy) {
			for (ix = 0; ix < npixx; ++ix) {
				fscanf(ifp,"%d %d %d %f %f %f", &idumx, &idumy, &idumz,
						 &efield[ix][iy][iz].f[0], &efield[ix][iy][iz].f[1], &efield[ix][iy][iz].f[2]);
				efield[ix][iy][iz].f[3]=0.;
				
				/* Force boundary conditions */			  
				
				if(ix == 0) efield[ix][iy][iz].f[0]=1.e-6;
				if(ix == mnode[0]) efield[ix][iy][iz].f[0]=1.e-6;
				if(iy == 0) efield[ix][iy][iz].f[1]=1.e-6;
				if(iy == mnode[1]) efield[ix][iy][iz].f[1]=1.e-6;
				
				/* Convert from V/cm to V/m */          
				
#ifdef __POWERPC__		              
				efield[ix][iy][iz].v = vec_madd(vhund.v,efield[ix][iy][iz].v,vzero.v);
#else
				efield[ix][iy][iz].v = _mm_mul_ps(vhund.v,efield[ix][iy][iz].v);
#endif
			}
      }
	}
	printf("central p-side field = %f %f %f\n", efield[0][0][1].f[0], efield[0][0][1].f[1], efield[0][0][1].f[2]);
	printf("central .5 way field = %f %f %f\n", efield[0][0][npixz/2].f[0], efield[0][0][npixz/2].f[1],
			 efield[0][0][npixz/2].f[2]);
	printf("central n-side field = %f %f %f\n", efield[0][0][npixz-2].f[0], efield[0][0][npixz-2].f[1],
			 efield[0][0][npixz-2].f[2]);
	
	fclose(ifp);
    
/* Read in weighting potential file */

    ifp = fopen("wgt_pot.init", "r");
    if (ifp==NULL) {
		printf("no wgt_pot.init initialization file found/n");
		return 0;
    }
    fscanf(ifp,"%f %f %f %d %d %d", &tdummy, &xdummy, &ydummy, &idumx, &idumy, &idumz);
    if(tdummy != *thick || xdummy != *xsize || ydummy != *ysize || idumx != npixx || idumy != npixy || idumz != npixz) {
       printf("thickness/xsize/ysize = %f/%f/%f, npixx/npixy/npixz = %d/%d/%d match problem\n",
              tdummy, xdummy, ydummy, idumx, idumy, idumz);
       return 0;
    }
    
	/* Now do the weighting potential */
	
    for (iz = 0; iz < npixz; ++iz) {
		for (iy = 0; iy < npixy; ++iy) {
			for (ix = 0; ix < npixx; ++ix) {
				fscanf(ifp,"%d %d %d %f %f %f %f %f %f %f %f %f",&idumx,&idumy,&idumz,&wdum[0][0],&wdum[0][1],
				       &wdum[0][2],&wdum[1][0],&wdum[1][1],&wdum[1][2],&wdum[2][0],&wdum[2][1],&wdum[2][2]);
				for(i=0;i<3;++i) {
				   for(j=0;j<3;++j) {
/* large meshes can sometimes produce small negative values ... zero them */
				      if(wdum[i][j] > 0.) {wgtpot[ix][iy][iz][i].f[j] = wdum[i][j];} else {wgtpot[ix][iy][iz][i].f[j] = 0.;}
				   }
				   wgtpot[ix][iy][iz][i].f[3]=0.;	
				}		
			}
		}
    }
    
    printf("central p-side weighting pot = %f %f %f\n", wgtpot[0][0][1][1].f[0], wgtpot[0][0][1][1].f[1], wgtpot[0][0][1][1].f[2]);
    printf("central .5 way weighting pot = %f %f %f\n", wgtpot[0][0][npixz/2][1].f[0], wgtpot[0][0][npixz/2][1].f[1],
		   wgtpot[0][0][npixz/2][1].f[2]);
    printf("central n-side weighting pot = %f %f %f\n", wgtpot[0][0][npixz-2][1].f[0], wgtpot[0][0][npixz-2][1].f[1],
		   wgtpot[0][0][npixz-2][1].f[2]);
	
    fclose(ifp);
	return 0;
} /* pixinit_ */

/* Subroutine */ int ranlux_0_(int n__, float *rvec, int *lenv, int *
							   isdext, int *lout, int *inout, int *k1, int *k2, 
							   int *lux, int *ins)
{
    /* Initialized data */
	
    static int notyet = 1;
    static int luxlev = 3;
    static int in24 = 0;
    static int kount = 0;
    static int mkount = 0;
    static int i24 = 24;
    static int j24 = 10;
    static float carry = 0.f;
    static int ndskip[5] = { 0,24,73,199,365 };
	
    /* System generated locals */
    int i__1, i__2;
	
    /* Local variables */
    static int i__, k, lp, isd, isk;
    static float uni;
    static int ilx, ivec, izip, next[24], izip2, jseed;
    static float seeds[24];
    static int inner, nskip;
    static float twom12, twom24;
    static int inseed, iseeds[24], iouter;
	
	
	/*         Subtract-and-borrow random number generator proposed by */
	/*         Marsaglia and Zaman, implemented by F. James with the name */
	/*         RCARRY in 1991, and later improved by Martin Luescher */
	/*         in 1993 to produce "Luxury Pseudorandom Numbers". */
	/*     Fortran 77 coded by F. James, 1993 */
	
	/*       references: */
	/*  M. Luscher, Computer Physics Communications  79 (1994) 100 */
	/*  F. James, Computer Physics Communications 79 (1994) 111 */
	
	/*   LUXURY LEVELS. */
	/*   ------ ------      The available luxury levels are: */
	
	/*  level 0  (p=24): equivalent to the original RCARRY of Marsaglia */
	/*           and Zaman, very long period, but fails many tests. */
	/*  level 1  (p=48): considerable improvement in quality over level 0, */
	/*           now passes the gap test, but still fails spectral test. */
	/*  level 2  (p=97): passes all known tests, but theoretically still */
	/*           defective. */
	/*  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible */
	/*           correlations have very small chance of being observed. */
	/*  level 4  (p=389): highest possible luxury, all 24 bits chaotic. */
	
	/* !!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
	/* !!!  Calling sequences for RANLUX:                                  ++ */
	/* !!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++ */
	/* !!!                   32-bit random floating point numbers between  ++ */
	/* !!!                   zero (not included) and one (also not incl.). ++ */
	/* !!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++ */
	/* !!!               one 32-bit int INT and sets Luxury Level LUX  ++ */
	/* !!!               which is int between zero and MAXLEV, or if   ++ */
	/* !!!               LUX .GT. 24, it sets p=LUX directly.  K1 and K2   ++ */
	/* !!!               should be set to zero unless restarting at a break++ */
	/* !!!               point given by output of RLUXAT (see RLUXAT).     ++ */
	/* !!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four ints++ */
	/* !!!               which can be used to restart the RANLUX generator ++ */
	/* !!!               at the current point by calling RLUXGO.  K1 and K2++ */
	/* !!!               specify how many numbers were generated since the ++ */
	/* !!!               initialization with LUX and INT.  The restarting  ++ */
	/* !!!               skips over  K1+K2*E9   numbers, so it can be long.++ */
	/* !!!   A more efficient but less convenient way of restarting is by: ++ */
	/* !!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++ */
	/* !!!                   ISVEC of 25 32-bit ints (see RLUXUT)      ++ */
	/* !!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++ */
	/* !!!                 32-bit int seeds, to be used for restarting ++ */
	/* !!!      ISVEC must be dimensioned 25 in the calling program        ++ */
	/* !!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
    /* Parameter adjustments */
    if (rvec) {
		--rvec;
	}
    if (isdext) {
		--isdext;
	}
	
    /* Function Body */
    switch(n__) {
		case 1: goto L_rluxin;
		case 2: goto L_rluxut;
		case 3: goto L_rluxat;
		case 4: goto L_rluxgo;
	}
	
	/*                               default */
	/*  Luxury Level   0     1     2   *3*    4 */
	/* orresponds to p=24    48    97   223   389 */
	/*     time factor 1     2     3     6    10   on slow workstation */
	/*                 1    1.5    2     3     5   on fast mainframe */
	
	/*  NOTYET is .TRUE. if no initialization has been performed yet. */
	/*              Default Initialization by Multiplicative Congruential */
    if (notyet) {
		notyet = 0;
		jseed = 314159265;
		inseed = jseed;
		printf("Ranlux default initialization %d \n",jseed);
		luxlev = 3;
		nskip = ndskip[luxlev];
		lp = nskip + 24;
		in24 = 0;
		kount = 0;
		mkount = 0;
		printf("Ranlux default luxury level = %d, p = %d \n", luxlev, lp);
		twom24 = 1.f;
		for (i__ = 1; i__ <= 24; ++i__) {
			twom24 *= .5f;
			k = jseed / 53668;
			jseed = (jseed - k * 53668) * 40014 - k * 12211;
			if (jseed < 0) {
				jseed += 2147483563;
			}
			iseeds[i__ - 1] = jseed % 16777216;
			/* L25: */
		}
		twom12 = twom24 * 4096.f;
		for (i__ = 1; i__ <= 24; ++i__) {
			seeds[i__ - 1] = (float) iseeds[i__ - 1] * twom24;
			next[i__ - 1] = i__ - 1;
			/* L50: */
		}
		next[0] = 24;
		i24 = 24;
		j24 = 10;
		carry = 0.f;
		if (seeds[23] == 0.f) {
			carry = twom24;
		}
    }
	
	/*          The Generator proper: "Subtract-with-borrow", */
	/*          as proposed by Marsaglia and Zaman, */
	/*          Florida State University, March, 1989 */
	
    i__1 = *lenv;
    for (ivec = 1; ivec <= i__1; ++ivec) {
		uni = seeds[j24 - 1] - seeds[i24 - 1] - carry;
		if (uni < 0.f) {
			uni += 1.f;
			carry = twom24;
		} else {
			carry = 0.f;
		}
		seeds[i24 - 1] = uni;
		i24 = next[i24 - 1];
		j24 = next[j24 - 1];
		rvec[ivec] = uni;
		/*  small numbers (with less than 12 "significant" bits) are "padded". */
		if (uni < twom12) {
			rvec[ivec] += twom24 * seeds[j24 - 1];
			/*        and zero is forbidden in case someone takes a logarithm */
			if (rvec[ivec] == 0.f) {
				rvec[ivec] = twom24 * twom24;
			}
		}
		/*        Skipping to luxury.  As proposed by Martin Luscher. */
		++in24;
		if (in24 == 24) {
			in24 = 0;
			kount += nskip;
			i__2 = nskip;
			for (isk = 1; isk <= i__2; ++isk) {
				uni = seeds[j24 - 1] - seeds[i24 - 1] - carry;
				if (uni < 0.f) {
					uni += 1.f;
					carry = twom24;
				} else {
					carry = 0.f;
				}
				seeds[i24 - 1] = uni;
				i24 = next[i24 - 1];
				j24 = next[j24 - 1];
				/* L90: */
			}
		}
		/* L100: */
    }
    kount += *lenv;
    if (kount >= 1000000000) {
		++mkount;
		kount += -1000000000;
    }
    return 0;
	
	/*           Entry to input and float int seeds from previous run */
	
L_rluxin:
    twom24 = 1.f;
    for (i__ = 1; i__ <= 24; ++i__) {
		next[i__ - 1] = i__ - 1;
		/* L195: */
		twom24 *= .5f;
    }
    next[0] = 24;
    twom12 = twom24 * 4096.f;
    printf(" Full initialization of Ranlux with 25 integers: \n");
	for(i__ = 0; i__< 5; ++i__) {
		printf("%d  %d  %d  %d  %d \n", isdext[1+5*i__], isdext[2+5*i__], isdext[3+5*i__], isdext[4+5*i__], isdext[5+5*i__]);
    }
    for (i__ = 1; i__ <= 24; ++i__) {
		seeds[i__ - 1] = (float) isdext[i__] * twom24;
		/* L200: */
    }
    carry = 0.f;
    if (isdext[25] < 0) {
		carry = twom24;
    }
    isd = abs(isdext[25]);
    i24 = isd % 100;
    isd /= 100;
    j24 = isd % 100;
    isd /= 100;
    in24 = isd % 100;
    isd /= 100;
    luxlev = isd;
    if (luxlev <= 4) {
		nskip = ndskip[luxlev];
		printf("ranlux luxury level set by rluxin to: %d\n", luxlev);
    } else if (luxlev >= 24) {
		nskip = luxlev - 24;
		printf("ranlux p-value set by rluxin to: %d\n", luxlev);
    } else {
		nskip = ndskip[4];
		printf("ranlux ILLEGAL LUXURY rluxin: %d\n", luxlev);
		luxlev = 4;
    }
    inseed = -1;
	notyet = 0;
    return 0;
	
	/*                    Entry to ouput seeds as ints */
	
L_rluxut:
    for (i__ = 1; i__ <= 24; ++i__) {
		isdext[i__] = (int) (seeds[i__ - 1] * 4096.f * 4096.f);
		/* L300: */
    }
    isdext[25] = i24 + j24 * 100 + in24 * 10000 + luxlev * 1000000;
    if (carry > 0.f) {
		isdext[25] = -isdext[25];
    }
    return 0;
	
	/*                    Entry to output the "convenient" restart point */
	
L_rluxat:
    *lout = luxlev;
    *inout = inseed;
    *k1 = kount;
    *k2 = mkount;
    return 0;
	
	/*                    Entry to initialize from one or three ints */
	
L_rluxgo:
    if (*lux < 0) {
		luxlev = 3;
    } else if (*lux <= 4) {
		luxlev = *lux;
    } else if (*lux < 24 || *lux > 2000) {
		luxlev = 4;
		printf("ranlux ILLEGAL LUXURY rluxgo: %d\n", *lux);
    } else {
		luxlev = *lux;
		for (ilx = 0; ilx <= 4; ++ilx) {
			if (*lux == ndskip[ilx] + 24) {
				luxlev = ilx;
			}
			/* L310: */
		}
    }
    if (luxlev <= 4) {
		nskip = ndskip[luxlev];
		i__1 = nskip + 24;
		printf("ranlux luxury level set by rluxgo: %d, p = %d\n", luxlev, i__1);
    } else {
		nskip = luxlev - 24;
		printf("ranlux luxury level set by rluxgo to: %d\n", luxlev);
    }
    in24 = 0;
    if (*ins < 0) {
		printf("Illegal initialization by RLUXGO, negative input seed\n");
    }
    if (*ins > 0) {
		jseed = *ins;
		printf("ranlux initialized by rluxgo from seeds: %d, %d, %d \n", jseed, *k1, *k2);
    } else {
		jseed = 314159265;
		printf("ranlux initialized by rluxgo from default seed \n");
    }
    inseed = jseed;
    notyet = 0;
    twom24 = 1.f;
    for (i__ = 1; i__ <= 24; ++i__) {
		twom24 *= .5f;
		k = jseed / 53668;
		jseed = (jseed - k * 53668) * 40014 - k * 12211;
		if (jseed < 0) {
			jseed += 2147483563;
		}
		iseeds[i__ - 1] = jseed % 16777216;
		/* L325: */
    }
    twom12 = twom24 * 4096.f;
    for (i__ = 1; i__ <= 24; ++i__) {
		seeds[i__ - 1] = (float) iseeds[i__ - 1] * twom24;
		next[i__ - 1] = i__ - 1;
		/* L350: */
    }
    next[0] = 24;
    i24 = 24;
    j24 = 10;
    carry = 0.f;
    if (seeds[23] == 0.f) {
		carry = twom24;
    }
	/*        If restarting at a break point, skip K1 + IGIGA*K2 */
	/*        Note that this is the number of numbers delivered to */
	/*        the user PLUS the number skipped (if luxury .GT. 0). */
    kount = *k1;
    mkount = *k2;
    if (*k1 + *k2 != 0) {
		i__1 = *k2 + 1;
		for (iouter = 1; iouter <= i__1; ++iouter) {
			inner = 1000000000;
			if (iouter == *k2 + 1) {
				inner = *k1;
			}
			i__2 = inner;
			for (isk = 1; isk <= i__2; ++isk) {
				uni = seeds[j24 - 1] - seeds[i24 - 1] - carry;
				if (uni < 0.f) {
					uni += 1.f;
					carry = twom24;
				} else {
					carry = 0.f;
				}
				seeds[i24 - 1] = uni;
				i24 = next[i24 - 1];
				j24 = next[j24 - 1];
				/* L450: */
			}
			/* L500: */
		}
		/*         Get the right value of IN24 by direct calculation */
		in24 = kount % (nskip + 24);
		if (mkount > 0) {
			izip = 1000000000 % (nskip + 24);
			izip2 = mkount * izip + in24;
			in24 = izip2 % (nskip + 24);
		}
		/*       Now IN24 had better be between zero and 23 inclusive */
		if (in24 > 23) {
			printf("Error in RESTARTING with RLUXGO: \n");
			printf("The values %d, %d, %d cannot occur at luxury level %d\n", *ins, *k1, *k2, luxlev);
			in24 = 0;
		}
    }
    return 0;
} /* ranlux_ */

/* Subroutine */ int ranlux_(float *rvec, int *lenv)
{
    return ranlux_0_(0, rvec, lenv, (int *)0, (int *)0, (int *)0, 
					 (int *)0, (int *)0, (int *)0, (int *)0);
}

/* Subroutine */ int rluxin_(int *isdext)
{
    return ranlux_0_(1, (float *)0, (int *)0, isdext, (int *)0, (
																 int *)0, (int *)0, (int *)0, (int *)0, (int *)
					 0);
}

/* Subroutine */ int rluxut_(int *isdext)
{
    return ranlux_0_(2, (float *)0, (int *)0, isdext, (int *)0, (
																 int *)0, (int *)0, (int *)0, (int *)0, (int *)
					 0);
}

/* Subroutine */ int rluxat_(int *lout, int *inout, int *k1, 
							 int *k2)
{
    return ranlux_0_(3, (float *)0, (int *)0, (int *)0, lout, inout, 
					 k1, k2, (int *)0, (int *)0);
}

/* Subroutine */ int rluxgo_(int *lux, int *ins, int *k1, int 
							 *k2)
{
    return ranlux_0_(4, (float *)0, (int *)0, (int *)0, (int *)0, (
																   int *)0, k1, k2, lux, ins);
}

/* Subroutine */ int rnpssn_0_(int n__, float *amu, int *n, int *ierr, float *amx)
{
    /* Initialized data */
	
    static float amu0 = -12345.67f;
    static float amax = 88.f;
	
	static int c__1 = 1;
	
    /* Builtin functions */
    double sqrt(double), exp(double);
	
    /* Local variables */
    static int j;
    static float p, r__, emu;
    extern /* Subroutine */ int ranlux_(float *, int *), rnormx_(float *, 
																 int *, int (*)(float *, int *));
	
    switch(n__) {
		case 1: goto L_rnpset;
	}
	
    *ierr = 0;
    if (*amu <= 0.f) {
		*ierr = 1;
		j = 0;
    } else if (*amu > amax) {
		rnormx_(&r__, &c__1, ranlux_);
		j = r__ * sqrt(*amu) + *amu + .5f;
    } else {
		if (*amu != amu0) {
			amu0 = *amu;
			emu = exp(-(*amu));
		}
		p = 1.f;
		j = -1;
	L1:
		++j;
		ranlux_(&r__, &c__1);
		p *= r__;
		if (p > emu) {
			goto L1;
		}
    }
	/* PN */
    if (j < 0) {
		printf("RNPSSN: Warning: J<0; J = %d\n", j);
		printf("Correction: J=0\n");
		printf("Increase AMAX value!");
	    j = 0;
    }
	/* PN */
    *n = j;
    return 0;
	
L_rnpset:
    amax = fminf(*amx,88.f);
	printf("+++++ CERN V136 RNPSSN: SWITCH TO NORMAL APPROXIMATION FOR AMU > %f \n", amax);
    return 0;
} /* rnpssn_ */

/* Subroutine */ int rnpssn_(float *amu, int *n, int *ierr)
{
    return rnpssn_0_(0, amu, n, ierr, (float *)0);
}

/* Subroutine */ int rnpset_(float *amx)
{
    return rnpssn_0_(1, (float *)0, (int *)0, (int *)0, amx);
}

/* Subroutine */ int rnormx_(float *devias, int *ndev, int (*routin)(float *, int *))
{
    /* Initialized data */
	
    static float s = .449871f;
    static float t = -.386595f;
    static float a = .196f;
    static float b = .25472f;
    static float r1 = .27597f;
    static float r2 = .27846f;
	
	static int c__2 = 2;
	
    /* System generated locals */
    int i__1;
    float r__1, r__2;
	
    /* Builtin functions */
    double log(double);
	float fabsf(float);
	
    /* Local variables */
    static float q, u[2], v, x, y;
    static int idev;
    static float deviat;
	
	/*        Generator of a vector of independent Gaussian-distributed */
	/*        (pseudo-)random numbers, of mean zero and variance one, */
	/*        making use of a uniform pseudo-random generator (RANMAR). */
	/*        The algorithm for converting uniform numbers to Gaussian */
	/*        is that of "Ratio of Uniforms with Quadratic Bounds."  The */
	/*        method is in principle exact (apart from rounding errors), */
	/*        and is based on the variant published by Joseph Leva in */
	/*        ACM TOMS vol. 18(1992), page 449 for the method and 454 for */
	/*        the Fortran algorithm (ACM No. 712). */
	/*        It requires at least 2 and on average 2.74 uniform deviates */
	/*        per Gaussian (normal) deviate. */
	/*   WARNING -- The uniform generator should not produce exact zeroes, */
	/*   since the pair (0.0, 0.5) provokes a floating point exception. */
    /* Parameter adjustments */
    --devias;
	
    /* Function Body */
	/*         generate pair of uniform deviates */
    i__1 = *ndev;
    for (idev = 1; idev <= i__1; ++idev) {
	L50:
		(*routin)(u, &c__2);
		v = (u[1] - .5f) * 1.7156f;
		x = u[0] - s;
		y = fabsf(v) - t;
		/* Computing 2nd power */
		r__1 = x;
		q = r__1 * r__1 + y * (a * y - b * x);
		/*           accept P if inside inner ellipse */
		if (q < r1) {
			goto L100;
		}
		/*           reject P if outside outer ellipse */
		if (q > r2) {
			goto L50;
		}
		/*           reject P if outside acceptance region */
		/* Computing 2nd power */
		r__1 = v;
		/* Computing 2nd power */
		r__2 = u[0];
		if (r__1 * r__1 > log(u[0]) * -4.f * (r__2 * r__2)) {
			goto L50;
		}
		/*           ratio of P's coordinates is normal deviate */
	L100:
		deviat = v / u[0];
		/* L200: */
		devias[idev] = deviat;
    }
    return 0;
} /* rnormx_ */

