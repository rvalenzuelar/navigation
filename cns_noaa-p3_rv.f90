!******************************************************************
!**** THIS PROGRAM CALCULATES NAVIGATIONAL ERRORS:
!****
!****    D_TILT_aft
!****    D_TILT_fore
!****    D_ROLL_aft
!****    D_ROLL_fore
!****    D_PITCH
!****    D_HEADING
!****    RANGE_DELAY_aft
!****    RANGE_DELAY_fore
!****    D_Xwe
!****    D_Ysn
!****    D_Z
!****    D_VH_acft
!****
!**** FROM COMPARISONS BETWEEN:
!****
!****    - RADAR-DERIVED SURFACE AND DIGITAL TERRAIN MAP 
!****     (or CONSTANT GROUND LEVEL);
!****    - DOPPLER VELOCITY AT SURFACE LEVEL AND ZERO;
!****    - DOPPLER  VELOCITY AT LOW-ELEVATION CLOSE TO THE AIRCRAFT 
!****      AND THE PROJECTION OF THE FLIGHT-LEVEL (IN SITU) WIND;
!****
!**** DIFFERENT OPTIONS ARE AVAILABLE (see DATA_cns_run FILE)
!****
!**** THIS PROGRAM ALSO PRODUCES A "RADAR-DERIVED" SURFACE MAP
!**** (FILE "SURF_EL_*" HAS THE SAME STRUCTURE AS FILE "SURF_DTM_*")
!****
!****
!******************************************************************
!**** Author: Frank ROUX (LA, UPS-CNRS, Toulouse), March 2000  ****
!******************************************************************
!     Modified by Huaqing Cai at NCAR, April, 2010
!******************************************************************
!		Modified by Raul Valenzuela at NOAA, 2014
!******************************************************************
!		Algorithm based on Georgis et al (2000, GRH), 
!		Testud et al (1995, THL) and Lee et al (1994)
!		(RV)
!********************************************************************


		PROGRAM cns_eldo

		!IMPLICIT NONE

		INTEGER, PARAMETER :: nvar=12						! dTa-dTf-dRa-dRf-dP-dH-RDa-RDf-dX-dY-dZ-dV
		INTEGER, PARAMETER :: nxysurfmax=721 	! max number of cells in DTM (RV)
		!
		!     include '/home/users/rouf/SOURCES/ELDO/mes_parametres'
		! CAI-START: Inlcude the parameter file mes_parametres directlybelow
		!      instead of using the inlcude function above
		!******************************************************************
		!**** FICHIER 'mes_parametres'
		!******************************************************************

		INTEGER, PARAMETER :: MAXRAD=2 				! number of radars? (RV)
		INTEGER, PARAMETER :: maxport=2000		! buffer size for fields? (RV)
		INTEGER, PARAMETER :: MAXPORAD=MAXRAD*maxport

		! Variable for reading text files
		INTEGER ntimes		! number of beams (RV)
		INTEGER nranges		! number of gates per beam (RV)
		INTEGER nsweep		! Sweep number in each netcdf file; Used to identify different sweep
		INTEGER*4 counter  ! ray ID for the entire period (RV)
		INTEGER start_year
		INTEGER start_mon
		INTEGER start_day
		INTEGER start_hour
		INTEGER start_min
		INTEGER start_sec


		!  Scaler variable for coccrection factors
		REAL azimuth_correction
		REAL elevation_correction
		REAL rangevalue_correction
		REAL longitude_correction
		REAL latitude_correction
		REAL pressure_altitude_correction
		REAL radar_altitude_correction
		REAL ew_gound_speed_correction
		REAL ns_ground_speed_correction
		REAL vertical_velocity_correction
		REAL heading_correction
		REAL roll_correction
		REAL pitch_correction
		REAL drift_correction
		REAL rotation_correction
		REAL tilt_correction

		! Variable for cfac files
		REAL tilt_corr_aft
		REAL tilt_corr_fore
		REAL rot_angle_corr_aft
		REAL rot_angle_corr_fore
		REAL pitch_corr_cfac
		REAL drift_corr_cfac
		REAL rangevalue_delay_corr_aft
		REAL rangevalue_delay_corr_fore
		REAL pressure_alt_corr
		REAL ew_gndspd_corr

		! Scaler variable for each ray
		INTEGER sweep_number
		REAL*8  time
		REAL azimuth			! beam azimuth angle
		REAL elevation			! beam elevation angle
		REAL*8 latitude			! radar latitude
		REAL*8 longitude		! radar longitude
		REAL*8 altitude			! radar altitude (pressure)
		REAL altitude_agl		! radar altitude (above ground level)
		REAL heading			! platform heading
		REAL roll					! platform roll
		REAL pitch				! platform pitch
		REAL drift					! platform drift (dif between heading and track)
		REAL rotation			! radar rotation (what's the difference with azimuth? values are close)
		REAL tilt					! radar tilt (same as elevation)
		REAL acftspd_we		! aircraft speed west-eas
		REAL acftspd_sn		! aircraft speed south-north
		REAL acftspd_nz		! aircraft speed vertical
		REAL wind_we			! wind speed west-east
		REAL wind_sn			! wind speed south-north
		REAL wind_nz			! wind speed vertical

		! One dimensional array of DBZ, VR, SW, NCP, etc
		REAL rangevalue(maxport)	! beam range
		REAL ZE(maxport)					! reflectivity
		REAL NCP(maxport)				! normalized coherent power
		REAL VR(maxport)				! radial velocity
		REAL SW(maxport)				! spectral width
		REAL VS(maxport)				! block commented by a previous editor
		REAL VL(maxport)					! block commented by a previous editor

		! Variables for input file list
		CHARACTER(len=80) infilename
		CHARACTER(len=80) filename		!CfRadial file name
		INTEGER nfile										! total number of netcdf text file & file number
		INTEGER ifile
		INTEGER lastfile     			
		INTEGER iopen										! flag for file status (0=closed, 1=open)

		! Variables declarations previous in Franks' common block, which has been deleted  

		! From COMMON /CSPD_OU_CELV/
		INTEGER(4)	::	nb_portes
		REAL 					d_porte(MAXPORAD)

		! From COMMON /CFAC/ 
		REAL corr_azest(MAXRAD),			&
					 corr_elhor(MAXRAD),		&
					 corr_dist(MAXRAD),			&
					 corr_lon(MAXRAD),			&
					 corr_lat(MAXRAD),			&
					 corr_p_alt(MAXRAD),		&
					 corr_r_alt(MAXRAD),		&
					 corr_vwe_av(MAXRAD),	&
					 corr_vsn_av(MAXRAD),	&
					 corr_vnz_av(MAXRAD),	&
					 corr_cap(MAXRAD),			&
					 corr_roul(MAXRAD),			&
					 corr_tang(MAXRAD),		&
					 corr_derv(MAXRAD),		&
					 corr_rota(MAXRAD),			&
					 corr_incl(MAXRAD)

		! From COMMON /RYIB/ 
		INTEGER(2)	::	ih_rdl					! beam hour
		INTEGER 			im_rdl					! beam minute
		INTEGER 			is_rdl					! beam second
		INTEGER 			ims_rdl				! beam milisecond
		INTEGER(4)	::	num_swp

		! CAI-STOP
		REAL(4) ::	dgate_corr(maxport),				&									! range gate
								dgate_true(maxport), 		&
								vdop_corr(maxport), 		&
								rota_start(2) 					= 	[-999.,-999.], 	&
								rota_end(2) 					= 	[-999.,-999.],	&
								rota_prev(2) 					=	[-999.,-999.],	&								
								xp(2),						 								& ! counter of rays (float)
								ssc(2)							= 	[0.,0.],				& ! sum of shdg (sine of aircraft heading)
								scc(2)							= 	[0.,0.], 				&	! sum of chdg (cosine of aircraft heading)
								sxa(2)							= 	[0.,0.],				&	! sum of x_acft
								sya(2)							= 	[0.,0.],				&	! sum of y_acft
								sza(2)							= 	[0.,0.], 				&	! sum of z_acft
								sacfthspd(2)					= 	[0.,0.],				&	! sum of acftspd_hor (aircraft horizontal speed)
								stime(2)							= 	[0.,0.], 				&	! sum of time_ks
								xp_acft(2)						= 	[0.,0.],				& ! counter of xp_acft
								su_acft(2)						= 	[0.,0.],				& ! sum of acftspd_we
								sv_acft(2)						= 	[0.,0.],				&	! sum of acftspd_sn
								SW_acft(2)					= 	[0.,0.], 				&	! sum of acftspd_nz
								su_wind(2)						= 	[0.,0.],				&	! sum of wind_we
								sv_wind(2)						= 	[0.,0.],				&	! sum of wind_sn
								SW_wind(2)					= 	[0.,0.],				&	! sum of wind_nz
								xp_wind(2)						= 	[0.,0.], 				&	! counter of xp_wind
								stilt(2) 							= 	[0.,0.],				& ! sum of tilts
								stilt2(2) 							= 	[0.,0.],				&	! sum of tilts squared
								xsweeps(2)					= 	[0.,0.], 				& ! counter of xsweeps
								swdzsurf_sweep(2)		= 	[0.,0.],				&	! sum of wghtsurf_ray for DZ_surf
								dzsurfsweep_mean(2)	= 	[0.,0.], 				&	! sum of wghtsurf_ray*d_hsurf
								dzsurfsweep_rms(2)		= 	[0.,0.], 				&	! sum of wghtsurf_ray*d_hsurf*d_hsurf
								swvsurf_sweep(2)			= 	[0.,0.],				&	! sum of wghtsurf_ray for VDOP_surf
								vsurfsweep_mean(2)		= 	[0.,0.],				&	! sum of wghtsurf_ray*vdopsurf_ray
								vsurfsweep_rms(2)			= 	[0.,0.], 				&	! sum of wghtsurf_ray*vdopsurf_ray*vdopsurf_ray
								swinsitu_sweep(2)			= 	[0.,0.],				&	! sum of wghtinsitu_ig
								dvinsitusweep_mean(2)	= 	[0.,0.],				&	! sum of wghtinsitu_ig*dv_dopinsitu
								dvinsitusweep_rms(2)		= 	[0.,0.]						! sum of wghtinsitu_ig*dv_dopinsitu*dv_dopinsitu
					
		REAL(4) ::	var(nvar),															&
								xmat(nvar,nvar),											&
								vect(nvar), 												&
								xinv(nvar,nvar),											&
								res(nvar), 													&
								vect_dzsurf(nvar),										&
								xmat_dzsurf(nvar,nvar), 								&
								vect_vsurf(nvar),										&
								xmat_vsurf(nvar,nvar), 								&
								vect_dvinsitu(nvar),									&
								xmat_dvinsitu(nvar,nvar),							&
								alt_dtm(nxysurfmax,nxysurfmax), 				&	! aster DEM altitude
								swdzsurf_wri(nxysurfmax,nxysurfmax), 		&
								SW_or_altsurf_wri(nxysurfmax,nxysurfmax)
							
		REAL(4)	::	zs_rot(2,500),zs_el(2,500),zs_az(2,500), 				&
								zs_dsurf(2,500),zs_dhor(2,500), 						&
								zs_zsurf(2,500),	zs_hsurf(2,500),						&
								vs_dhor(2,500),vs_vdopsurf(2,500), 				&
								vi_dhor(2,500),vi_vdop(2,500),vi_vinsitu(2,500),&
								rms_var_zsurf(nvar),rms_var_vsurf(nvar), 		&
								rms_var_vinsitu(nvar),corr_var(nvar,nvar),		&
								s_vpv(2,2),sv_vpv(2,2),svv_vpv(2,2), 				&
								x_vpv(2,2),xv_vpv(2,2),xvv_vpv(2,2)

		INTEGER(2) ::  iyymmdd(3),	&
										ig_dismiss(15)

		INTEGER(4) ::	nb_ray(2)								=	[0,0],	&
										nb_sweep(2)					=	[0,0],	&
										n_dzsurf(2)					=	[0,0],	&
										n_vsurf(2)						=	[0,0],	&
										n_dvinsitu(2)					=	[0,0],	&
										nsurf_wri(2) 					=	[0,0],	&
										ndismiss_vhacft(2)			=	[0,0],	&
										ndismiss_vdopcorr(2)		=	[0,0],	&
										ndismiss_vdopsurf(2)		=	[0,0],	&
										swp(2),										&
										swp_prev(2),								&
										ndop_ok(2),								&
										nref_ok(2),									&
										istart_sweep(2),							&
										itab(nxysurfmax),						&
										ihms_dtm(6),								&
										ialtsurf_wri(nxysurfmax)   

		REAL 	vg(15),	&
						vu(15)

		CHARACTER directory*60,dir_READ*60,									&
										fich_sis*30,												&
										fich_cornav*30,											&
										fich_log*30,												&
										fich_surf*30,												&										
										dtm_file*50,												&
										wrisurfile*50,												&
										yymmdd_dtm*12,suff_dtm*20,					&
										yymmdd*12,c_hms_min*7,c_hms_max*7,	&
										projname*12, 											&
						  				argu*30 ! command line arguments

		CHARACTER(len=24) :: cfac_text01="azimuth_corr"
		CHARACTER(len=24) :: cfac_text02="elevation_corr"
		CHARACTER(len=24) :: cfac_text03="rangevalue_delay_corr"
		CHARACTER(len=24) :: cfac_text04="longitude_corr"
		CHARACTER(len=24) :: cfac_text05="latitude_corr"
		CHARACTER(len=24) :: cfac_text06="pressure_alt_corr"
		CHARACTER(len=24) :: cfac_text07="radar_alt_corr"
		CHARACTER(len=24) :: cfac_text08="ew_gndspd_corr"
		CHARACTER(len=24) :: cfac_text09="ns_gndspd_corr"
		CHARACTER(len=24) :: cfac_text10="vert_vel_corr"
		CHARACTER(len=24) :: cfac_text11="heading_corr"
		CHARACTER(len=24) :: cfac_text12="roll_corr"
		CHARACTER(len=24) :: cfac_text13="pitch_corr"
		CHARACTER(len=24) :: cfac_text14="drift_corr"
		CHARACTER(len=24) :: cfac_text15="rot_angle_corr"
		CHARACTER(len=24) :: cfac_text16="tilt_corr"

		COMMON	/cosinang/	crr,srr,cti,sti, 				&
														chdg,shdg,cdri,			&
														sdri,cpit,spit, 				&
														caze,saze,celh,selh

		!******************************************************************
		!**** CONSTANT PARAMETRES
		!******************************************************************
		deg_lon0					=	111.32			! [km], = 1 deg of longitude at equator 
		deg_lat						=	111.13			! [km], = 1 deg of latitude at 45deg of lat
		conv							=	3.14159/180.	! degrees to radians
		rayter							=	6370.  			! [km], Earth radius 
		xncp_min					=	0.25				! min	NCP
		SW_max					=	5.					! max SW
		vdop_max				=	200.				! max doppler velocity
		selh_surf					=	-0.15				! threshold for selh=sin(conv*elhor_ray)
		zacftmin_surf			=	0.5				! threshold for z_acft (km?) (1.5 in original code)
		igstart_surf				=	5					! first gate to dismiss ?? (should be read from DATA_cns)
		refsurf_min0			=	20.				! equivalent reflectivity (Georgis et al 2000)
		gradrefsurf_min0	=	50.				! gradient of equivalent reflectivity (Georgis et al 2000)
		dhsurf_max				=	999.				! ???
		vdopsurf_max			=	999.				! apparently not necessary
		dmax_insitu				=	5.					! ???
		xpmin_contray		=	3.					! used in continuity of the ray, not sure about its interpretation
		dvdop_max				=	10.				! threshold for dv (?)
		dvdopinsitu_max	=	999.				! theshold for dv_dopinsitu (?)
		selhinsitu_max		=	0.1				! threshold for selh=sin(conv*elhor_ray)
		ssurfins_min			=	1.					! threshold for ssurfins

		!********************************************************************
		!**** READ INPUT PARAMETERS ON FILE "DATA_cns"
		!*********************************************************************

		! READ in command line arguments
		CALL GETARG(1, argu)
		OPEN(99,file=argu,status='unknown')

		!==== DIRECTORY WITH DTM AND FOR OUTPUT FILES====
		PRINT *,' '
		READ(99,*)directory
		! length of output path directory
		ndir=1
		DO WHILE(directory(ndir:ndir).NE.' ')
				ndir=ndir+1
		ENDDO
		ndir=ndir-1 
		PRINT *,' DIRECTORY FOR THE OUTPUT FILES :',directory(1:ndir)

		!==== INTPUT DIRECTORY ====
		PRINT *,' '
		READ(99,*)dir_READ
		! length of input path directory
		ndirr=1
		DO WHILE(dir_READ(ndirr:ndirr).NE.' ')
				ndirr=ndirr+1
		ENDDO
		ndirr=ndirr-1
		PRINT *,' DIRECTORY FOR READ FILES :',dir_READ(1:ndirr)

		!==== NUMBER OF INPUT FILES (CfRadial in text format) ====
		PRINT *,' '
		READ(99,*)nfile
		PRINT *,' Total Number of Sweep Files :',nfile

		!==== FIRST INPUT FILES (CfRadial in text format) ====
		READ(99,*)ifile

		!=== INPUT DATE ====
		READ(99,*)iyymmdd
		WRITE(yymmdd,"(i12)")10000*iyymmdd(3)+100*iyymmdd(2)+iyymmdd(1)
		PRINT *,' YYYYMMDD :',yymmdd

		!==== INPUT TIME====
		READ(99,*)ihms_min,ihms_max
		PRINT *,' HHMMSS (min,max) :',ihms_min,ihms_max
		ih_min=ihms_min/10000
		im_min=(ihms_min-10000*ih_min)/100
		is_min=ihms_min-10000*ih_min-100*im_min
		ih_max=ihms_max/10000
		im_max=(ihms_max-10000*ih_max)/100
		is_max=ihms_max-10000*ih_max-100*im_max

		!==== LAT/LON ORIGIN ====
		READ(99,*)orig_lat,orig_lon
		PRINT *,' ORIGIN_LATITUDE,_LONGITUDE :',orig_lat,orig_lon

		!==== GATES TO DISMISS ====
		PRINT *,' '
		READ(99,*)ig_dismiss
		PRINT *,' 15 GATES TO DISMISS (0 IF not) :',ig_dismiss

		!==== MIN/MAX DISTANCE FROM RADAR ====
		PRINT *,' '
		READ(99,*)dmin,dmax0
		PRINT *,' DMIN,DMAX FROM RADAR [km]:',dmin,dmax0

		!==== ALTITUDE MODE ==== (pressure or agl?)
		PRINT *,' '
		READ(99,*)ipr_alt
		PRINT *,' ALTITUDE (1:pressure,2:radar) :',ipr_alt

		!==== REFERENCE DBZ (?) ====
		PRINT *,' '
		READ(99,*)ref_min0,ref_max
		PRINT *,' REF_min(at 10km),REF_max [dBZ]:',ref_min0,ref_max

		!==== DOPPLER VELOCITY MODE (?) ====
		PRINT *,' '
		READ(99,*)ichoice_vdop
		PRINT *,' WHICH VDOP (1:VR,2:VG,3:VU) :', ichoice_vdop
		IF(ichoice_vdop.EQ.1.OR.ichoice_vdop.EQ.2)THEN
		ictrl_contray=0
		PRINT *,' -> WILL NOT CONTROL CONTINUITY ALONG THE RAY !!!!'
		ELSE
		ictrl_contray=0
		ENDIF

		PRINT *,' '
		PRINT *,' CORRECTION OF NAVIGATIONAL ERRORS'
		PRINT *,' '
		PRINT *,' FIELDS TAKEN INTO ACCOUNT:'
		READ(99,*)

		!==== FLAGS FOR GRH J FUNCTIONS ====
		kdzsurf=1		! J1
		kvsurf=1		! J2
		kdvinsitu=1	! J3


		!==== RELATIVE WEIGHTS dZsurf, Vsurf, dVinsitu ====
		READ(99,*)rw_dzsurf,rw_vsurf,rw_dvinsitu
		PRINT *,'   -> REL.WGHT_dZsurf,Vsurf,dVinsitu (1/0) :',rw_dzsurf,rw_vsurf,rw_dvinsitu
		PRINT *,' '

		!==== FLAGS FOR RETRIEVALS ====
		PRINT *,' CORRECTIONS TO CALCULATE:'
		READ(99,*)
		READ(99,*)idtiltaft,idtiltfore
		PRINT *,'   -> D_TILT_AFT,D_TILT_FORE (1/0) :',idtiltaft,idtiltfore
		READ(99,*)idrotaaft,idrotafore
		PRINT *,'   -> D_ROTA_AFT,D_ROTA_FORE (1/0) :',idrotaaft,idrotafore
		READ(99,*)idpitch,idhdg
		PRINT *,'   -> D_PITCH,D_HEADING (1/0) :',idpitch,idhdg
		READ(99,*)irdaft,irdfore
		PRINT *,'   -> RANGE_DELAY_AFT,RANGE_DELAY_FORE (1/0) :',irdaft,irdfore
		READ(99,*)idxwe,idysn,idzacft
		PRINT *,'   -> D_XWE,D_YSN,D_ZACFT (1/0) :',idxwe,idysn,idzacft
		READ(99,*)idvh
		PRINT *,'   -> D_VHACFT (1/0) :',idvh
		PRINT *,' '
		READ(99,*)

		!==== FLAG FOR SIMULATION WITH GUESS====
		READ(99,*)isim
		PRINT *,' SIMULATION AVEC +dXXX_GUESS INITIAUX (1/0) :',isim
		PRINT *,' '
		READ(99,*)

		!==== GUESS FOR RETRIEVALS ====
		READ(99,*)dtiltaft_guess,dtiltfore_guess
		PRINT *,' D_TILT_AFT,D_TILT_FORE (deg) guess :',dtiltaft_guess,dtiltfore_guess
		READ(99,*)drotaaft_guess,drotafore_guess
		PRINT *,' D_ROTA_AFT,D_ROTA_FORE (deg) guess :',drotaaft_guess,drotafore_guess
		READ(99,*)dpitch_guess,dhdg_guess
		PRINT *,' D_PITCH,D_HEADING (deg) guess :',droll_guess,dpitch_guess,dhdg_guess
		READ(99,*)rdaft_guess,rdfore_guess
		PRINT *,' RANGE_DELAY_AFT,RANGE_DELAY_FORE (km) guess :',rdaft_guess,rdfore_guess
		READ(99,*)dxwe_guess,dysn_guess,dzacft_guess
		PRINT *,' D_XWE,D_YSN,D_ZACFT (km) guess :',dxwe_guess,dysn_guess,dzacft_guess
		READ(99,*)dvh_guess
		PRINT *,' D_VHACFT (m/s) guess :',dvh_guess
		READ(99,*)
		PRINT *,'  '

		!==== DIGITAL TERRAIN MODEL ====
		READ(99,*)idtmfile,dtm_file,zsurf_cst
		ndtmfile=0
		IF(idtmfile.NE.0)THEN
				DO WHILE(dtm_file(ndtmfile+1:ndtmfile+1).NE.' ')
						ndtmfile=ndtmfile+1
				ENDDO
		ENDIF
		IF(idtmfile.EQ.0)PRINT *,' NO "SURF_DTM_*" FILE WILL BE READ ', &
														'-> ZSURF_CST (km) =',zsurf_cst
		IF(idtmfile.EQ.1)PRINT *,' WILL READ "SURF_DTM_*" FILE :', dtm_file(1:ndtmfile)

		!==== SIMULATE SURFACE DTM MAP (?) ====
		READ(99,*)iwrisurfile,wrisurfile,xywidth_wrisurf,hxy_wrisurf
		nsf=0
		IF(iwrisurfile.EQ.1)THEN
				DO WHILE(wrisurfile(nsf+1:nsf+1).NE.' ')
						nsf=nsf+1
				ENDDO
				PRINT *,' WILL WRITE "SURF_EL_*" FILE : ', wrisurfile(1:nsf)
				xmin_wrisurf=-xywidth_wrisurf/2.
				xmax_wrisurf=+xywidth_wrisurf/2.
				ymin_wrisurf=-xywidth_wrisurf/2.
				ymax_wrisurf=+xywidth_wrisurf/2.
				PRINT *,' -> Xmin,max_wrisurf:',xmin_wrisurf,xmax_wrisurf
				PRINT *,'    Ymin,max_wrisurf:',ymin_wrisurf,ymax_wrisurf
				PRINT *,'    Hx,y_wrisurf:',hxy_wrisurf        
				nx_wrisurf=((xmax_wrisurf-xmin_wrisurf)/hxy_wrisurf+1.)
				ny_wrisurf=((ymax_wrisurf-ymin_wrisurf)/hxy_wrisurf+1.)
				PRINT *,'    Nx,Ny_wrisurf:',nx_wrisurf,ny_wrisurf
				IF(nx_wrisurf.GT.nxysurfmax.OR.ny_wrisurf.GT.nxysurfmax)THEN
						PRINT *,' !!!! Nx,Ny_wrisurf :',nx_wrisurf,ny_wrisurf, &
												' > NxySURFmax !!!!'
						PRINT *,' !!!! MODIFY l.30 AND RECOMPILE THE PROGRAM !!!!'
						GOTO 3	! exit program
				ENDIF       
				!***************************************************************************************
				!**** OPEN "SURF_EL_*" FILE #30 FOR WRITING (IF IWRISURFILE=1)
				!***************************************************************************************
				PRINT *,' OPEN "SURF_EL_*" FILE #30 FOR WRITING :', &
									directory(1:ndir) // '/' // wrisurfile(1:nsf)
				OPEN(30 ,file=directory(1:ndir) // '/' // wrisurfile(1:nsf), &
									form='formatted',status='unknown')
				iolat_wrisurf=(1000.*orig_lat)
				iolon_wrisurf=(1000.*orig_lon)
				ixmin_wrisurf=(1000.*xmin_wrisurf)
				iymin_wrisurf=(1000.*ymin_wrisurf)
				ihxy_wrisurf=(1000.*hxy_wrisurf)
				WRITE(30,111)yymmdd,'NOAA-P3',		&
								iolat_wrisurf,iolon_wrisurf,			&
								0,0,0,0,0, 											&
								ih_min,im_min,is_min, 				&
								ih_max,im_max,is_max, 			&
								ixmin_wrisurf,iymin_wrisurf,0, 	&
								nx_wrisurf,ny_wrisurf,1, 				&
								ihxy_wrisurf,ihxy_wrisurf,0
		ELSE
				PRINT *,' NO "SURF_EL_*" FILE WILL BE WRITTEN'
		ENDIF
		
		READ(99,*)
		!====WRITE OUTPUT LOG FILE FOR EACH RAY===
		READ(99,*)ioutputlog
		READ(99,*)iprintscreen

		!====CLOSE FILE DATA_cns====			
		CLOSE(99)

		IF(no_lect.GT.900)	THEN
				PRINT *, 'no_lect=',no_lect
				GOTO 3	! exit program
		END IF

		!******************************************************************
		!**** GENERATE THE SURFACE ARRAYS
		!******************************************************************
		PRINT *,' '
		PRINT *,' GENERATE THE SURFACE ARRAYS'
		PRINT *,' '

		IF(idtmfile.EQ.1)THEN
				!------------------------------------------------------------------
				!---- FROM THE INPUT "SURF_DTM_*" FILE
				!------------------------------------------------------------------
				 PRINT *,' IFIDTM=1 -> READ THE "SURF_DTM_*" FILE #20 :', &
									directory(1:ndir)//'/'//dtm_file(1:ndtmfile)
				 OPEN(20,file=directory(1:ndir) // '/' //dtm_file(1:ndtmfile), &
									form='formatted',status='unknown') ! open DTM ascii file

				!--------------I modified the previous-------------
				READ(20,111)projname, yymmdd_dtm, 	& 
											iolat_dtm,iolon_dtm, 			&
											ixmin_dtm,iymin_dtm, 		&
											nx_dtm,ny_dtm, 					&
											ihx_dtm,ihy_dtm       
				111   format(a12,a12,22i7)
				!---------------------------------------------------------(RV)

				!CAI START
				WRITE(*,*)'nx_dtm,ny_dtm,ixmin_dtm',nx_dtm,ny_dtm,ixmin_dtm
				! CAI STOP
				olat_dtm=float(iolat_dtm)/1000.
				olon_dtm=float(iolon_dtm)/1000.
				xlatm_surf=(olat_dtm+orig_lat)/2.
				deg_lon=deg_lon0*cos(conv*xlatm_surf)	!	[km/deg]
				dx_dtm=(olon_dtm-orig_lon)*deg_lon		!	[km]
				dy_dtm=(olat_dtm-orig_lat)*deg_lat			!	[km]
				xmin_dtm=float(ixmin_dtm)/1000.+dx_dtm
				ymin_dtm=float(iymin_dtm)/1000.+dy_dtm
				hx_dtm=float(ihx_dtm)/1000.					!	[km]
				hy_dtm=float(ihy_dtm)/1000.					!	[km]
				xmax_dtm=xmin_dtm+float(nx_dtm-1)*hx_dtm
				ymax_dtm=ymin_dtm+float(ny_dtm-1)*hy_dtm
				PRINT *,'    X_DTM_min,max :',xmin_dtm,xmax_dtm
				PRINT *,'    Y_DTM_min,max :',ymin_dtm,ymax_dtm
				PRINT *,'    Hx,y_DTM :',hx_dtm,hy_dtm
				PRINT *,'    Nx,y_DTM:',nx_dtm,ny_dtm
				saltdtm=0.
				altdtm_mean=0.
				altdtm_min=+999.
				altdtm_max=-999.
				DO jdtm=1,ny_dtm
						READ(20,333)(itab(idtm),idtm=1,nx_dtm)
						333      format(721i6) 	! 721 integer values made of 6 characters (RV)
						DO idtm=1,nx_dtm
								IF(itab(idtm).GT.-1000)THEN
										h_dtm=float(itab(idtm))/1000. ! dtm altitude in km
										alt_dtm(idtm,jdtm)=h_dtm
										saltdtm=saltdtm+1.											! counter 
										altdtm_mean=altdtm_mean+h_dtm
										altdtm_min=amin1(altdtm_min,h_dtm)		! minimum dtm altitude
										altdtm_max=amax1(altdtm_max,h_dtm)	! maximum dtm altitude
								ELSE
										alt_dtm(idtm,jdtm)=-999.
								ENDIF
						ENDDO
				ENDDO
				altdtm_mean=altdtm_mean/amax1(1.,saltdtm)
				CLOSE(20) ! close DTM ascii file
		ELSEIF(idtmfile.EQ.0)THEN
				!------------------------------------------------------------------
				!---- FROM ZSURF_CST (READ in DATA_cns)
				!------------------------------------------------------------------
				PRINT *,' IFIDTM=0 -> ALT_SURF(i,j)=cst',' (',zsurf_cst,' )'
				xmin_dtm=xmin_wrisurf
				ymin_dtm=ymin_wrisurf
				hx_dtm=hxy_wrisurf
				hy_dtm=hxy_wrisurf
				nx_dtm=nx_wrisurf
				ny_dtm=ny_wrisurf
				xmax_dtm=xmin_dtm+float(nx_dtm-1)*hx_dtm
				ymax_dtm=ymin_dtm+float(ny_dtm-1)*hy_dtm
				DO jdtm=1,ny_dtm
						DO idtm=1,nx_dtm
								alt_dtm(idtm,jdtm)=zsurf_cst
						ENDDO
				ENDDO
				altdtm_mean=zsurf_cst
				altdtm_min=zsurf_cst
				altdtm_max=zsurf_cst
		ENDIF

		PRINT *,'     -> NPTS:',int(saltdtm), 											&
								 ' ALTSURF_mean,min,max:',altdtm_mean, 	&
									altdtm_min,altdtm_max
		zsurfrad_min=altdtm_min-1.			! what's this for?
		zsurfrad_max=altdtm_max+1.		!
		PRINT *,' '
		PRINT *,' -> AUTHORIZED ZSURF_RAD_min,max :', &
								zsurfrad_min,zsurfrad_max

		!******************************************************************
		!**** MIN AND MAX TIMES
		!******************************************************************
		tmin=3.6*float(ih_min)+0.06*float(im_min) +0.001*float(is_min)
		tmax=3.6*float(ih_max)+0.06*float(im_max)+0.001*float(is_max)
		WRITE(c_hms_min,"(i7)")1000000+ihms_min
		WRITE(c_hms_max,"(i7)")1000000+ihms_max
		WRITE(fich_cornav,"('CORNAV_E_',a6,'_',a6)") c_hms_min(2:7),c_hms_max(2:7)
		WRITE(fich_sis,"('SIS_E_',a6,'_',a6)")  c_hms_min(2:7),c_hms_max(2:7)
		WRITE(fich_log,"('OUTPUT_LOG_',a6,'_',a6)") c_hms_min(2:7),c_hms_max(2:7)
		WRITE(fich_surf,"('OUTPUT_SURF_',a6,'_',a6)") c_hms_min(2:7),c_hms_max(2:7)		     

		!******************************************************************
		!**** OPEN THE OUPUT "CORNAV_EL_*" FILE #10
		!******************************************************************
		PRINT *,' '
		PRINT *,' OPEN "CORNAV_EL_*" FILE #10 :', &
							directory(1:ndir)//'/'//fich_cornav
		OPEN(10,file=directory(1:ndir)//'/'//fich_cornav, &
							form='formatted',status='unknown')
		WRITE(10,"(' YYYYMMDD : ',a12)")yymmdd
		WRITE(10,"(' HHMMSS_min HHMMSS_max : ',a6,3x,a6,/)") &
				 			 c_hms_min(2:7),c_hms_max(2:7)
		WRITE(10,"( ' FIELDS TAKEN INTO ACCOUNT',// &
							'  -> REL.WGHT_dZsurf,Vsurf,dVinsitu : ',3f6.3,/)") &
								rw_dzsurf,rw_vsurf,rw_dvinsitu
		WRITE(10,"( ' VARIABLES TAKEN INTO ACCOUNT',// &
							'  -> D_TILT_AFT,D_TILT_FORE (1/0) : ',2i2,// &
							'  -> D_ROTA_AFT,D_ROTA_FORE (1/0) : ',2i2,// &
							'  -> D_PITCH,D_HEADING (1/0) : ',2i2,// &
							'  -> RANGE_DELAY_AFT,RANGE_DELAY_FORE (1/0) : ',2i2,// &
							'  -> D_XWE,D_YSN,D_ZACFT (1/0) : ',3i2,// &
							'  -> D_VHACFT (1/0) : ',i2)") &    
								idtiltaft,idtiltfore, 																							&
								idrotaaft,idrotafore, 																						&
								idpitch,idhdg, 																									&
								irdaft,irdfore, 																									&
								idxwe,idysn,idzacft, 																					&
								idvh	
     
		IF(idtmfile.EQ.1)THEN
		  WRITE(10,"(' READS THE SURF_DTM_* FILE :',a50)") &
								 directory(1:ndir)//'/'//dtm_file(1:ndtmfile)
		ELSE
		  WRITE(10,"( ' NO SURF_DTM_* FILE TO READ ',// &
								'-> ALT_SURF(x,y)=CST (',f6.3,')')") &
								zsurf_cst
		ENDIF

		IF(iwrisurfile.EQ.1)THEN
		  WRITE(10,"(' WRITES THE SURF_EL_* FILE :',a50,//)") &
								 directory(1:ndir)//'/'//wrisurfile(1:nsf)
		ELSE
		  WRITE(10,"(' NO SURF_EL_* FILE TO WRITE ',//)")
		ENDIF

		!******************************************************************
		!**** OPEN THE OUTPUT "SIS_EL_*" FILE #50
		!******************************************************************
		PRINT *,' '
		PRINT *,' OPEN "SIS_EL_*" FILE #50 :', directory(1:ndir)//'/'//fich_sis
		OPEN(50,file=directory(1:ndir)//'/'//fich_sis,form='unformatted',status='unknown')
		WRITE(50)iyymmdd,orig_lon,orig_lat,ihms_min,ihms_max


		!******************************************************************
		!**** OPEN  OUTPUT LOG FILE #60
		!******************************************************************
		IF (ioutputlog.EQ.1) THEN
				PRINT *,' '
				PRINT *,' OPEN "OUTPUT_LOG_*" FILE #60 :', directory(1:ndir)//'/'//fich_log
				OPEN(60,file=directory(1:ndir)//'/'//fich_log,form='formatted',status='unknown')
		ENDIF

		!******************************************************************
		!**** OPEN  OUTPUT SURF_POINTS FILE #70
		!******************************************************************
		PRINT *,' '
		PRINT *,' OPEN "OUTPUT_SURFPOINTS_*" FILE #70 :', directory(1:ndir)//'/'//fich_surf
		OPEN(70,file=directory(1:ndir)//'/'//fich_surf,form='formatted',status='unknown')
		
		!******************************************************************
		!**** INITIALIZATIONS
		!******************************************************************
		time_prev=-999
		ihms_prev=-999
		DO iradar=1,2
				istart_sweep(iradar)=0
				nsurf_wri(iradar)=0
				DO jgd=1,2
						s_vpv(iradar,jgd)=0.
						sv_vpv(iradar,jgd)=0.
						svv_vpv(iradar,jgd)=0.
						x_vpv(iradar,jgd)=0.
						xv_vpv(iradar,jgd)=0.
						xvv_vpv(iradar,jgd)=0.
				ENDDO
		ENDDO     

		DO i=1,nvar
				rms_var_zsurf(i)=0.
				rms_var_vsurf(i)=0.
				rms_var_vinsitu(i)=0.				
				var(i)=0.				
				vect_dzsurf(i)=0.
				vect_vsurf(i)=0.
				vect_dvinsitu(i)=0.
				vect(i)=0.
				res(i)=0.				
				DO j=1,nvar
						corr_var(i,j)=0.
						xmat_dzsurf(i,j)=0.
						xmat_vsurf(i,j)=0.
						xmat_dvinsitu(i,j)=0.
						xmat(i,j)=0.						
				ENDDO				
		ENDDO
		
		DO i_wrisurf=1,nxysurfmax
				DO j_wrisurf=1,nxysurfmax
						swdzsurf_wri(i_wrisurf,j_wrisurf)=0.
						SW_or_altsurf_wri(i_wrisurf,j_wrisurf)=0.
				ENDDO
		ENDDO

		DO ig=1,maxport
				vdop_corr(ig)=-999.
		ENDDO

		DO i=1,2
				DO n=1,500
						zs_rot(i,n)=0.
						zs_el(i,n)=0.
						zs_az(i,n)=0.
						zs_dsurf(i,n)=0.
						zs_dhor(i,n)=0.
						zs_zsurf(i,n)=0.
						zs_hsurf(i,n)=0.
						vs_dhor(i,n)=0.
						vs_vdopsurf(i,n)=0.
						vi_dhor(i,n)=0.
						vi_vdop(i,n)=0.
						vi_vinsitu(i,n)=0.
				ENDDO
		ENDDO
		
		ssurfins=0.	! Olivier (réel)
		swdzsurf_tot=0.
		swdzmsurf_tot=0.
		swdz2surf_tot=0.
		swadzsurf_tot=0.		
		swvsurf_tot=0.
		swvmsurf_tot=0.
		swv2surf_tot=0.
		swavsurf_tot=0.
		swdvinsitu_tot=0.
		swdvminsitu_tot=0.
		swdv2insitu_tot=0.
		swadvinsitu_tot=0.
		nb1=0
		nb2=0
		nb3=0
		nb4=0
		nb5=0
		nb6=0
		nb7=0
		nb8=0
		nsup=0
		nbtotals=0	! Olivier
		nbon=0			! Olivier
		nmauvais=0	! Olivier

      IF(no_lect.EQ.999) &
      	GOTO 3 !	exit program

		!******************************************************************
		!*** READ THE NOAA-P3 DATA FROM CFRADIAL TEXT FILES 
		!*** CREATED WITH readnetcdf_DBZ_VR.f90 
		!******************************************************************
      PRINT *,' '
      PRINT *,'********************************************************************'
      PRINT *,' READ NOAA-P3 DATA FROM CFRADIAL TEXT FILES'
      PRINT *,'*********************************************************************'
      PRINT *,' '
      iopen = 0		! flag for CfRadial file state (0=closed, 1=open)
      !ifile = 1		! CfRadial text file number
      					! ifile is now an input in the DATA_cns file
      filecounter=1
      lastfile = 0
      
1	   IF(iopen .EQ. 0) THEN
				WRITE(infilename,'(I0.3)') ifile ! integer of the form 001,002,etc
				infilename = dir_READ(1:ndirr) // '/' // TRIM(ADJUSTL(infilename)) // '.txt'
				OPEN(55,file=infilename, status='old')
				iopen = 1
		ENDIF
      
      ! ==== READ DATA FROM A SINGLE BEAM ====
		READ(55,101,END=5)counter, 																	&
						nsweep,NTIMES,NRANGES, 													&
						start_year,start_mon,start_day, 													&
						start_hour,start_min,start_sec,time, 											&
						azimuth,elevation,latitude,longitude,											&
						altitude, altitude_agl,																	&
						heading,roll,pitch,drift,rotation,tilt,												&
						acftspd_we,acftspd_sn,acftspd_nz,												&
						wind_we,wind_sn,wind_nz,															&
						azimuth_correction,elevation_correction, 										&
						rangevalue_correction,longitude_correction,latitude_correction, 		&
						pressure_altitude_correction,radar_altitude_correction, 			        &
						ew_gound_speed_correction,ns_ground_speed_correction, 	        &
						vertical_velocity_correction,heading_correction, 						    &
						roll_correction,pitch_correction,drift_correction, 							&
						rotation_correction,tilt_correction
		READ(55,102,END=5)counter,(rangevalue(J),J=1,nranges)
		READ(55,102,END=5)counter,(ZE(J),J=1,nranges)
		READ(55,102,END=5)counter,(NCP(J),J=1,nranges) ! <---fake -999 field
		READ(55,102,END=5)counter,(VR(J),J=1,nranges)
		READ(55,102,END=5)counter,(SW(J),J=1,nranges) ! <---fake -999 field
		101  format(I10,2x,50x,3I10,I5,5I3,d20.8,2f10.4,3d20.8,29f10.4) 
		102  format(I10,800f10.4) 

		! When READ reads successfully
		GOTO 6
		
        ! When READ reaches end-of-file
5      CLOSE(55)
        iopen = 0
        IF(filecounter .EQ. nfile) THEN
		        lastfile=1 
		        GOTO 7 
        ELSE
		        ifile = ifile +1
		        filecounter=filecounter+1
		        GOTO 1
        ENDIF

6      CONTINUE

		! ************ Get the ray time *************
		ih_rdl1 = start_hour
		im_rdl1 = start_min
		is_rdl1 = start_sec
		ims_rdl1 = (time-INT(time))*1000

		! add to the start seconds by time, which is the elpased time after start time
		is_rdl1 = is_rdl1+INT(time)
		
		! adjusting hh,mm,ss for passing 60 mark, assign to Frank's ray time variables
		ims_rdl = ims_rdl1
		is_rdl = MOD(is_rdl1,60)
		im_rdl1 = im_rdl1+is_rdl1/60
		im_rdl = MOD(im_rdl1, 60)
		ih_rdl1 = ih_rdl1+im_rdl1/60
		ih_rdl = MOD(ih_rdl1, 60)

		! Assign  the total number of gates AND rangevalue of each gates,
		!  The aft/fore radar are different
		nb_portes = nranges
		IF (tilt .LT. 0) THEN !AFT,iradar_ray=1,iaftfore= -1
				DO ig = 1, nranges
						d_porte(ig) = rangevalue(ig)	! In SOLO, these ranges are consistent with values in
																				! "Examine->cell values" but inconsistent with values
																			 	!	in "Data Widget" (RV)
				ENDDO
				! Assign the correction factors to Frank's variable
				! NOTE: Here the correction factors are arrays with two elements
				! This is different from any other variables				
				corr_azest(1) = azimuth_correction
				corr_elhor(1) = elevation_correction
				corr_dist(1) = rangevalue_correction
				corr_lon(1) = longitude_correction
				corr_lat(1) = latitude_correction
				corr_p_alt(1) = pressure_altitude_correction
				corr_r_alt(1) = radar_altitude_correction
				corr_vwe_av(1) = ew_gound_speed_correction
				corr_vsn_av(1) = ns_ground_speed_correction
				corr_vnz_av(1) = vertical_velocity_correction
				corr_cap(1) = heading_correction
				corr_roul(1) = roll_correction
				corr_tang(1) = pitch_correction
				corr_derv(1) = drift_correction
				corr_rota(1) = rotation_correction
				corr_incl(1) = tilt_correction				
		ELSEIF(tilt .GT. 0) THEN !FORE,iradar_ray=2,iaftfore= +1
				DO ig = 1, nranges
						d_porte(maxport+ig) = rangevalue(ig)  ! This change fixed icorrupted infilename
				ENDDO
				corr_azest(2) = azimuth_correction
				corr_elhor(2) = elevation_correction
				corr_dist(2) = rangevalue_correction
				corr_lon(2) = longitude_correction
				corr_lat(2) = latitude_correction
				corr_p_alt(2) = pressure_altitude_correction
				corr_r_alt(2) = radar_altitude_correction
				corr_vwe_av(2) = ew_gound_speed_correction
				corr_vsn_av(2) = ns_ground_speed_correction
				corr_vnz_av(2) = vertical_velocity_correction
				corr_cap(2) = heading_correction
				corr_roul(2) = roll_correction
				corr_tang(2) = pitch_correction
				corr_derv(2) = drift_correction
				corr_rota(2) = rotation_correction
				corr_incl(2) = tilt_correction				
		ENDIF
		
		! Assign the swp number READ from text file to num_swp
       num_swp = nsweep 


		!******************************************************************
		!**** CONTROL FOR THE LAST CFRADIAL TEXT FILE
		!******************************************************************
 7   	iend=0
		IF(filecounter.EQ.nfile .AND. lastfile .EQ.1)THEN   
				iend=2
				PRINT *,' '
				PRINT *,'**** END OF READING ALL CFRADIAL TEXT FILES ****'
		ENDIF
      
		!******************************************************************
		!**** CONTROL OF CURRENT TIME
		!******************************************************************
		ih_ray=ih_rdl
		im_ray=im_rdl
		is_ray=is_rdl
		ims_ray=ims_rdl
		ihhmmss=10000*ih_ray+100*im_ray+is_ray
		
		IF(ihhmmss.LE.0) &
				GOTO 1	! if iopen=1 read next CfRadial line
		time_ks=3.6*float(ih_ray)+0.06*float(im_ray) +0.001*float(is_ray)+1.e-6*float(ims_ray)
		IF( time_ks-time_prev.LT.-80..OR.time_ks-tmin.LT.-80.)THEN
				time_ks=time_ks+86.4
				ihhmmss=ihhmmss+240000
		ENDIF
		
		time_prev=time_ks
		
		IF(time_ks.lt.tmin)THEN
				IF(ihhmmss/10.GT.ihms_prev)THEN
						PRINT *,' HHMMSS:',ihhmmss,' < HHMMSS_min:',ihms_min
						ihms_prev=ihhmmss/10
				ENDIF
				IF(iend .NE. 2) &
						GOTO 1   ! if iopen=1 read next CfRadial line
		ENDIF

		IF(time_ks.GT.tmax)THEN
				iend=2
				PRINT *,' '
				PRINT *,' HHMMSSms:',100*ihhmmss+ims_rdl,' > HHMMSSms_max:',100*ihms_max
		ENDIF
		
		IF(iend.EQ.2) &
				GOTO 2

		!******************************************************************
		!**** CONTROL OF LAT, LON, P_ALT AND R_ALT
		!******************************************************************
		IF( abs(latitude).LT.0.001 .OR.abs(longitude).LT.0.001.OR. &
			(abs(altitude).LT.0.001 .AND.abs(altitude_agl).LT.0.001)) &
				GOTO 1	! if iopen=1 read next CfRadial line


		!******************************************************************
		!**** RADAR IDENTIFICATION THROUGH TILT_RAY (=INCL_RDL)
		!****  -> AFT : IRADAR_RAY=1, IAFTFORE=-1
		!****  -> FORE : IRADAR_RAY=2, IAFTFORE=1
		!******************************************************************
		tilt_ray=tilt
		IF(abs(tilt_ray).LT.15.)THEN
				GOTO 1	! if iopen=1 read next CfRadial line
		ELSEIF(abs(tilt_ray).LT.30.)THEN
				IF(tilt_ray.LT.-15.)THEN
						iradar_ray=1
						iaftfore=-1
						swp(iradar_ray)=num_swp
				ENDIF
				IF(tilt_ray.GT.+15.)THEN
						iradar_ray=2
						iaftfore=+1
						swp(iradar_ray)=num_swp
				ENDIF
		ELSE
				GOTO 1	! if iopen=1 read next CfRadial line
		ENDIF

		!******************************************************************
		!**** NYQUIST VELOCITY
		!******************************************************************
		! CAI-START----- Oliver's modification is for P3, for ELDORA, we don't
		!                need these NYQUIST velocity stuff!!!!
		! CAI-STOP

		!     vnyq=vit_nonamb(iradar_ray)	! Mod Oliv
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!!!! THIS VALUE IS CORRECT FOR IPP1/IPP2=4/5 ONLY !!!!
		!      vnyq_el=vnyq/20.
		!      vnyq_s=5.*vnyq_el	!Olivier
		!      vnyq_l=4.*vnyq_el	!Olivier
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!******************************************************************
		!**** CONTROL FOR AN END OF SWEEP
		!******************************************************************
		IF(nb_ray(iradar_ray).EQ.1)THEN
				tandrot=0.
		ELSE
				tandrot=tan(conv*(rotation-rota_prev(iradar_ray)))
		ENDIF

		IF(nb_ray(iradar_ray).GT.1 .AND.( (swp(iradar_ray).NE.swp_prev(iradar_ray)).OR. &
				(abs(tandrot).GT.0.2) ) ) &
				iend=1
				

		!******************************************************************
		!****    END OF A SWEEP ( IEND = 1 )                                 
		!**** or END OF THE TAPE or END OF CONSIDERED PERIOD ( IEND = 2 )    
		!******************************************************************
 2   CONTINUE
 
		IF(iend.GE.1)THEN 
				rota_end(iradar_ray)=rota_prev(iradar_ray)

				!******************************************************************
				!**** MEAN VALUES FOR THE PAST SWEEP
				!******************************************************************
				IF(nb_ray(iradar_ray).GT.1)THEN
						xp(iradar_ray)=float(nb_ray(iradar_ray))
						tilt_mean=stilt(iradar_ray)/xp(iradar_ray)
						tilt_rms=sqrt(amax1(0.,(xp(iradar_ray)*stilt2(iradar_ray) 	&
												-stilt(iradar_ray)*stilt(iradar_ray)) 					&
												/(xp(iradar_ray)*(xp(iradar_ray)-1))))
						nb_sweep(iradar_ray)=nb_sweep(iradar_ray)+1
						xacft_mean=sxa(iradar_ray)/xp(iradar_ray)
						yacft_mean=sya(iradar_ray)/xp(iradar_ray)
						zacft_mean=sza(iradar_ray)/xp(iradar_ray)
						acfthspd_mean=sacfthspd(iradar_ray)/xp(iradar_ray)
						time_ks_mean=stime(iradar_ray)/xp(iradar_ray)
						ihmean=time_ks_mean/3.6
						immean=(time_ks_mean-3.6*float(ihmean))/0.06
						ismean=(time_ks_mean-3.6*float(ihmean) -0.06*float(immean))/0.001
						ihms=10000*ihmean+100*immean+ismean
						hdg_mean=atan2( ssc(iradar_ray)/xp(iradar_ray), &
						scc(iradar_ray)/xp(iradar_ray))/conv
						uacft_mean=su_acft(iradar_ray)/xp_acft(iradar_ray)
						vacft_mean=sv_acft(iradar_ray)/xp_acft(iradar_ray)
						wacft_mean=SW_acft(iradar_ray)/xp_acft(iradar_ray)
						uwind_mean=su_wind(iradar_ray)/xp_wind(iradar_ray)
						vwind_mean=sv_wind(iradar_ray)/xp_wind(iradar_ray)
						wwind_mean=SW_wind(iradar_ray)/xp_wind(iradar_ray)

						!******************************************************************
						!**** CONTROL PRINTS FOR THE PAST SWEEP
						!******************************************************************
						PRINT *,' '
						PRINT *,' '
						PRINT *,' *******************************************'
						PRINT *,' **** CONTROL PRINTS FOR THE PAST SWEEP ****'
						PRINT *,' *******************************************'
						PRINT *,' '
						PRINT *,' HHMMSS :',ihms
						PRINT *,' RADAR(aft=1,fore=2) :',iradar_ray
						PRINT *,' SWEEP(aft=-1,fore=+1) :',iaftfore, &
											' NO_SWEEP(this program) :',nb_sweep(iradar_ray), &
											'  [ on tape :',swp(iradar_ray),' ]'
						PRINT *,' X_we/OLON,Y_sn/OLAT,Z_acft :', &
											xacft_mean,yacft_mean,zacft_mean
						PRINT *,' HEADING :',hdg_mean
						PRINT *,' U,V,W_acft :',uacft_mean,vacft_mean,wacft_mean
						PRINT *,' U,V,W_insitu :',uwind_mean,vwind_mean,wwind_mean
						PRINT *,' -> NB_RAYS_THIS_SWEEP :',nb_ray(iradar_ray)
						PRINT *,' -> TILT_mean,rms :',tilt_mean,tilt_rms
						PRINT *,' -> ROTA_start,end :',rota_start(iradar_ray), rota_end(iradar_ray)
						PRINT *,' '
						PRINT *,' -> NREF_OK:',nref_ok(iradar_ray),'  NDOP_OK:',ndop_ok(iradar_ray)
						PRINT *,' '

						!======WRITE TO OUTPUT LOG FILE=====
						IF (ioutputlog.EQ.1) THEN
								WRITE(60,800) ' '
								WRITE(60,800) ' '
								WRITE(60,800) '*******************************************'
								WRITE(60,800) '**** CONTROL PRINTS FOR THE PAST SWEEP ****'
								WRITE(60,800) '*******************************************'
								WRITE(60,800) ' '            
								WRITE(60,801) 'HHMMSS :',ihms
								WRITE(60,801) 'RADAR(aft=1,fore=2) :',iradar_ray
								WRITE(60,803) 'SWEEP(aft=-1,fore=+1) :',iaftfore, &
																' NO_SWEEP(this program) :',nb_sweep(iradar_ray), &
																'  [ on tape :',swp(iradar_ray),' ]'
								WRITE(60,804) ' X_we/OLON,Y_sn/OLAT,Z_acft :', &
																	xacft_mean,yacft_mean,zacft_mean
								WRITE(60,805) ' HEADING :',hdg_mean
								WRITE(60,804) ' U,V,W_acft :',uacft_mean,vacft_mean,wacft_mean
								WRITE(60,804) ' U,V,W_insitu :',uwind_mean,vwind_mean,wwind_mean
								WRITE(60,801) ' -> NB_RAYS_THIS_SWEEP :',nb_ray(iradar_ray)
								WRITE(60,804) ' -> TILT_mean,rms :',tilt_mean,tilt_rms
								WRITE(60,804) ' -> ROTA_start,end :',rota_start(iradar_ray), rota_end(iradar_ray)
								WRITE(60,800) ' '
								WRITE(60,802) ' -> NREF_OK:',nref_ok(iradar_ray),' NDOP_OK:',ndop_ok(iradar_ray)
								WRITE(60,800) ' '

								800   format (X,A)
								801   format (A,I10)
								802   format (A,I10,A,I10)
								803   format (A,I10,A,I10,A,I10,A)
								804   format (A,F10.2,F10.2,F10.2)
								805   format (A,F10.2)
						ENDIF
						!===========================================================(RV)

						IF(kdzsurf.EQ.1)THEN
								IF(swdzsurf_sweep(iradar_ray).GT.0.)THEN
										bias_dzsurf=dzsurfsweep_mean(iradar_ray) /swdzsurf_sweep(iradar_ray)
										stdv_dzsurf=sqrt( swdzsurf_sweep(iradar_ray) 			&
																			*dzsurfsweep_rms(iradar_ray) 		&
																			- dzsurfsweep_mean(iradar_ray)	&
																			*dzsurfsweep_mean(iradar_ray))	&
																			/swdzsurf_sweep(iradar_ray)
										PRINT *,' -> dZHSURF_npts,swghts,mean,stdv :', 				&
															n_dzsurf(iradar_ray),swdzsurf_sweep(iradar_ray), &
															bias_dzsurf,stdv_dzsurf
										IF(iwrisurfile.EQ.1) &
										PRINT *,'     [ NPTS_SURF FOR SURF_EL_*:', &
																nsurf_wri(iradar_ray),' ]'
								ELSE
										PRINT *,' !!!! NPTS_dZHSURF :',n_dzsurf(iradar_ray),' !!!!'
								ENDIF
						ENDIF

						IF(kvsurf.EQ.1)THEN
								IF(swvsurf_sweep(iradar_ray).GT.0.)THEN
										bias_vsurf=vsurfsweep_mean(iradar_ray)/swvsurf_sweep(iradar_ray)
										stdv_vsurf=sqrt(swvsurf_sweep(iradar_ray) 		&
																		*vsurfsweep_rms(iradar_ray) 	&
																		- vsurfsweep_mean(iradar_ray) &
																		*vsurfsweep_mean(iradar_ray)) &
																		/swvsurf_sweep(iradar_ray)
										PRINT *,' -> VSURF_npts,swghts,mean,stdv :', &
																	n_vsurf(iradar_ray),swvsurf_sweep(iradar_ray), &
																	bias_vsurf,stdv_vsurf
								ELSE
										PRINT *,' !!!! NPTS_VSURF :',n_vsurf(iradar_ray),' !!!!'
										PRINT *,' !!!! Ndismissed_VACFT,VDOPCORR,VDOPSURF:', &
																	ndismiss_vhacft(iradar_ray), &
																	ndismiss_vdopcorr(iradar_ray), &
																	ndismiss_vdopsurf(iradar_ray),' !!!!'
								ENDIF
						ENDIF

						IF(kdvinsitu.EQ.1)THEN
								IF(swinsitu_sweep(iradar_ray).GT.0.)THEN
										bias_dvinsitu=dvinsitusweep_mean(iradar_ray)/swinsitu_sweep(iradar_ray)
										stdv_dvinsitu=sqrt(swinsitu_sweep(iradar_ray) 				&
																			*dvinsitusweep_rms(iradar_ray)		&
																			- dvinsitusweep_mean(iradar_ray)	&
																			*dvinsitusweep_mean(iradar_ray))	&
																			/swinsitu_sweep(iradar_ray)
										PRINT *,' -> dVINSITU_npts,swghts,mean,stdv :', &
																n_dvinsitu(iradar_ray),swinsitu_sweep(iradar_ray), &
																bias_dvinsitu,stdv_dvinsitu
										PRINT *,'     -> LEFT_swghts,mean,stdv:', &
																	s_vpv(iradar_ray,1), &
																	sv_vpv(iradar_ray,1)/amax1(0.001,s_vpv(iradar_ray,1)), &
																	sqrt( s_vpv(iradar_ray,1)*svv_vpv(iradar_ray,1) &
																	-sv_vpv(iradar_ray,1)*sv_vpv(iradar_ray,1)) &
																	/amax1(0.001,s_vpv(iradar_ray,1))
										PRINT *,'     -> RIGHT_swghts,mean,stdv:', &
																	s_vpv(iradar_ray,2), &
																	sv_vpv(iradar_ray,2)/amax1(0.001,s_vpv(iradar_ray,2)), &
																	sqrt( s_vpv(iradar_ray,2)*svv_vpv(iradar_ray,2) &
																	-sv_vpv(iradar_ray,2)*sv_vpv(iradar_ray,2)) &
																	/amax1(0.001,s_vpv(iradar_ray,2))
								ELSE
										PRINT *,' !!!! NPTS_VINSITU :',n_dvinsitu(iradar_ray),' !!!!'
								ENDIF
						ENDIF

						PRINT *,' '
						PRINT *,' *******************************************'
						PRINT *,' '
						PRINT *,' '

						!******************************************************************
						!**** WRITE THE RESULTS FOR THE PAST SWEEP 
						!**** ON THE "SIS_EL_*" FILE #50
						!******************************************************************

						!**** SWEEP HEADER
						WRITE(50)iaftfore,nb_sweep(iradar_ray), &
												xacft_mean,yacft_mean,zacft_mean, &
												time_ks_mean,hdg_mean, &
												u_mean,v_mean,w_mean

						!******************************************************************
						!**** SWEEP DATA: DZ_surf
						!******************************************************************
						PRINT *,' '
						IF(kdzsurf.EQ.1.AND.n_dzsurf(iradar_ray).GT.0)THEN
								WRITE(50)n_dzsurf(iradar_ray)
								WRITE(50)( zs_rot(iradar_ray,n), &
														zs_el(iradar_ray,n),zs_az(iradar_ray,n), &
														zs_dsurf(iradar_ray,n),zs_dhor(iradar_ray,n), &
														zs_zsurf(iradar_ray,n),zs_hsurf(iradar_ray,n), &
														n=1,n_dzsurf(iradar_ray))			
								!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
								!!!!      PRINT *,' '
								!!!!      PRINT *,' SIS_* -> NPTS_Zsurf:',n_dzsurf(iradar_ray)
								!!!!      PRINT *,' [ ROT - DH - Z_surf - H_surf ]'
								!!!!      DO n=1,n_dzsurf(iradar_ray)
								!!!!         PRINT *,zs_rot(iradar_ray,n),zs_dhor(iradar_ray,n)
								!!!!     &          ,zs_zsurf(iradar_ray,n),zs_hsurf(iradar_ray,n)
								!!!!      ENDDO
								!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						ELSE
								WRITE(50)0
						ENDIF

						!******************************************************************
						!**** SWEEP DATA: VDOP_surf
						!******************************************************************
						IF(kvsurf.EQ.1.AND.n_vsurf(iradar_ray).GT.0)THEN
								WRITE(50)n_vsurf(iradar_ray)
								WRITE(50)(vs_dhor(iradar_ray,n),vs_vdopsurf(iradar_ray,n), &
								n=1,n_vsurf(iradar_ray))
						ELSE
								WRITE(50)0
						ENDIF

						!******************************************************************
						!**** SWEEP DATA: DVDOP_insitu
						!******************************************************************
						IF(kdvinsitu.EQ.1.AND.n_dvinsitu(iradar_ray).GT.0)THEN
								WRITE(50)n_dvinsitu(iradar_ray)
								WRITE(50)( vi_dhor(iradar_ray,n), &
								vi_vdop(iradar_ray,n),vi_vinsitu(iradar_ray,n), &
								n=1,n_dvinsitu(iradar_ray))
						ELSE
								WRITE(50)0
						ENDIF
				ENDIF    !!  of  !! IF(nb_ray(iradar_ray).GT.1)THEN  !!
		
				!******************************************************************
				!**** END OF THE TAPE or END OF CONSIDERED PERIOD ( IEND = 2 )
				!******************************************************************
				IF(iend.EQ.2)THEN 
		
						PRINT *,' '
						PRINT *,' ****************************************************'
						PRINT *,'   HHMMSS :',ih_ray,im_ray,is_ray, &
										'   -> END OF CONSIDERED PERIOD'
						PRINT *,'   NB_SWEEPS_READ FOR AFT AND FORE RADARS :', nb_sweep
						PRINT *,' ****************************************************'
						PRINT *,' '
						PRINT *,' '
						WRITE(10,"(' NB_SWEEPS FOR THE AFT AND FORE RADARS: ',2i5,/)"), &
											int(xsweeps(1)),int(xsweeps(2))

						!*****************************************************************************************
						!**  SUM OF INDIVIDUAL WEIGHTS, MEAN, RMS VALUES OF
						!**    - DIFFERENCE OF SURFACE ALTITUDE (RADAR-DTM);
						!**    - DOPPLER VELOCITY OF THE GROUND CLUTTER;
						!**    - DIFFERENCE OF DOPPLER VELOCITY (RADAR-FLIGHT_LEVEL);
						!*****************************************************************************************
						PRINT *,' '
						PRINT *,' **********************************************'
						PRINT *,' ************ MEAN AND RMS ERRORS *************'
						PRINT *,' **********************************************'
						PRINT *,' '

						IF(kdzsurf.EQ.1)THEN
								bias_dzsurf=swdzmsurf_tot/amax1(1.,swdzsurf_tot)
								stdv_dzsurf=sqrt( swdzsurf_tot*swdz2surf_tot &
												 -swdzmsurf_tot*swdzmsurf_tot) &
										  /amax1(1.,swdzsurf_tot)
								PRINT *,' '
								PRINT *,' dZ_surf (km) sum_wghts,mean,stdv :', &
												swdzsurf_tot,bias_dzsurf,stdv_dzsurf
								WRITE(10,"(//,' dZ_surf (km) sum_wghts,mean,stdv :',f10.1,2f8.3,/)") &
								swdzsurf_tot,bias_dzsurf,stdv_dzsurf
						ENDIF

						IF(kvsurf.EQ.1)THEN
								bias_vsurf=swvmsurf_tot/amax1(1.,swvsurf_tot)
								stdv_vsurf=sqrt( swvsurf_tot*swv2surf_tot &
													-swvmsurf_tot*swvmsurf_tot) &
											 /amax1(1.,swvsurf_tot)
								PRINT *,' '
								PRINT *,' VDOP_surf (m/s) sum_wghts,mean,stdv :', &
												swvsurf_tot,bias_vsurf,stdv_vsurf
								WRITE(10,"(' VDOP_surf (m/s) sum_wghts,mean,stdv :',f10.1,2f8.3,/)") &
									  swvsurf_tot,bias_vsurf,stdv_vsurf
						ENDIF

						IF(kdvinsitu.EQ.1)THEN
								bias_dvinsitu=swdvminsitu_tot/amax1(1.,swdvinsitu_tot)
								stdv_dvinsitu=sqrt( swdvinsitu_tot*swdv2insitu_tot &
																		-swdvminsitu_tot*swdvminsitu_tot) &
																		/amax1(1.,swdvinsitu_tot)
								PRINT *,' '
								PRINT *,' dVDOP_insitu (m/s) sum_wghts,mean,stdv :', &
													swdvinsitu_tot,bias_dvinsitu,stdv_dvinsitu
								WRITE(10,"(' dVDOP_insitu (m/s) sum_wghts,mean,stdv :',f10.1,2f8.3)") &
															swdvinsitu_tot,bias_dvinsitu,stdv_dvinsitu
								DO iradar=1,2
										PRINT *,'   IRADAR (AR=1,AV=2) :',iradar
										bias_dvinsitu_ir_g=xv_vpv(iradar,1)/amax1(1.,x_vpv(iradar,1))
										stdv_dvinsitu_ir_g=sqrt( x_vpv(iradar,1)*xvv_vpv(iradar,1) &
																					-xv_vpv(iradar,1)*xv_vpv(iradar,1)) &
																					/amax1(1.,x_vpv(iradar,1))
										PRINT *,'     -> VDOP-PROJWIND_LEFT_npts,mean,stdv:', &
																	x_vpv(iradar,1), &
																	bias_dvinsitu_ir_g,stdv_dvinsitu_ir_g
										WRITE(10,"('   IRADAR (AR=1,AV=2) :',i1,// &
																	'    -> VDOP-PROJWIND_LEFT_npts,mean,stdv:', &
																	f10.1,2f8.3)") &
																	iradar,x_vpv(iradar,1), &
																	bias_dvinsitu_ir_g,stdv_dvinsitu_ir_g
										bias_dvinsitu_ir_d=xv_vpv(iradar,2)/amax1(1.,x_vpv(iradar,2))
										stdv_dvinsitu_ir_d=sqrt( x_vpv(iradar,2)*xvv_vpv(iradar,2) &
																						-xv_vpv(iradar,2)*xv_vpv(iradar,2)) &
																						/amax1(1.,x_vpv(iradar,2))
										PRINT *,'     -> VDOP-PROJWIND_RIGHT_npts,mean,stdv:', &
																	x_vpv(iradar,2), &
																	bias_dvinsitu_ir_d,stdv_dvinsitu_ir_d
										WRITE(10,"('    -> VDOP-PROJWIND_RIGHT_npts,mean,stdv:', &
																	f10.1,2f8.3)") &
																	x_vpv(iradar,2), &
																	bias_dvinsitu_ir_d,stdv_dvinsitu_ir_d
								ENDDO
								PRINT *,' '
								WRITE(10,"(//)")
						ENDIF

						PRINT *,' '
						PRINT *,' **********************************************'
						PRINT *,' '
						PRINT *,' '
						PRINT *,' '

						!******************************************************************
						!****  (IF SUM_WGHTS_surf+insitu > SUM_WGHTS_min)
						!****   -> NAVIGATIONAL ERROS CAN BE CALCULATED
						!******************************************************************
						IF(ssurfins.GT.ssurfins_min)THEN		

								!******************************************************************
								!**** RMS VALUES OF THE NORMALIZED VARIABLES 
								!******************************************************************
								PRINT *,' '
								PRINT *,' **********************************************'
								PRINT *,' *** RMS VALUES OF THE NORMALIZED VARIABLES ***'
								PRINT *,' **********************************************'
								PRINT *,' '

								IF(swdzsurf_tot.GT.1.)THEN
										PRINT *,' DZ_surf -> sWGHTs:',swdzsurf_tot
										PRINT *,'          rms_VAR(dTaft,dTfore):', &
													(sqrt(rms_var_zsurf(i)/swadzsurf_tot),i=1,2)
										PRINT *,'          -------(dRaft,dRfore):', &
													(sqrt(rms_var_zsurf(i)/swadzsurf_tot),i=3,4)
										PRINT *,'          -------(dPitch,dHdg):', &
													(sqrt(rms_var_zsurf(i)/swadzsurf_tot),i=5,6)
										PRINT *,'          -------(RDaft,RDfore):', &
													(sqrt(rms_var_zsurf(i)/swadzsurf_tot),i=7,8)
										PRINT *,'          -------(dXwe,dYsn,dZacft):', &
													(sqrt(rms_var_zsurf(i)/swadzsurf_tot),i=9,11)
										PRINT *,'          -------(dVHacft):', &
													sqrt(rms_var_zsurf(12)/swadzsurf_tot)
								ELSE
										PRINT *,' !!!! DZ_surf -> sWGHTs:',swdzsurf_tot,' !!!!' 
								ENDIF

								IF(swvsurf_tot.GT.1.)THEN
										PRINT *,' VDOP_surf -> sWGHTs:',swvsurf_tot
										PRINT *,'          rms_VAR(dTaft,dTfore):', &
													(sqrt(rms_var_vsurf(i)/swavsurf_tot),i=1,2)
										PRINT *,'          -------(dRaft,dRfore):', &
													(sqrt(rms_var_vsurf(i)/swavsurf_tot),i=3,4)
										PRINT *,'          -------(dPitch,dHdg):', &
													(sqrt(rms_var_vsurf(i)/swavsurf_tot),i=5,6)
										PRINT *,'          -------(RDaft,RDfore):', &
													(sqrt(rms_var_vsurf(i)/swavsurf_tot),i=7,8)
										PRINT *,'          -------(dXwe,dYsn,dZacft):', &
													(sqrt(rms_var_vsurf(i)/swavsurf_tot),i=9,11)
										PRINT *,'          -------(dVHacft):', &
													sqrt(rms_var_vsurf(12)/swavsurf_tot)
								ELSE
										PRINT *,' !!!! VDOP_surf -> sWGHTs:',swvsurf_tot,' !!!!' 
								ENDIF

								IF(swdvinsitu_tot.GT.1.)THEN
										PRINT *,' DVDOP_insitu -> sWGHTs:',swdvinsitu_tot
										PRINT *,'          rms_VAR(dTaft,dTfore):', &
													(sqrt(rms_var_vinsitu(i)/swadvinsitu_tot),i=1,2)
										PRINT *,'          -------(dRaft,dRfore):', &
													(sqrt(rms_var_vinsitu(i)/swadvinsitu_tot),i=3,4)
										PRINT *,'          -------(dPitch,dHdg):', &
													(sqrt(rms_var_vinsitu(i)/swadvinsitu_tot),i=5,6)
										PRINT *,'          -------(RDaft,RDfore):', &
													(sqrt(rms_var_vinsitu(i)/swadvinsitu_tot),i=7,8)
										PRINT *,'          -------(dXwe,dYsn,dZacft):', &
													(sqrt(rms_var_vinsitu(i)/swadvinsitu_tot),i=9,11)
										PRINT *,'          -------(dVHacft):', &
													sqrt(rms_var_vinsitu(12)/swadvinsitu_tot)
								ELSE
										PRINT *,' !!!! DVDOP_insitu -> sWGHTs:', &
														swdvinsitu_tot,' !!!!' 
								ENDIF

								!******************************************************************
								!**** NORMALIZED CORRELATION MATRIX BETWEEN THE NVAR VARIABLES
								!******************************************************************
								PRINT *,' '
								sp_zsvszi=swdzsurf_tot+swvsurf_tot+swdvinsitu_tot
								IF(sp_zsvszi.GT.1.)THEN
										PRINT *,' **********************************************'
										PRINT *,' ******* NORMALIZED CORRELATION MATRIX ********'
										PRINT *,' *******            (*1000)            ********'
										PRINT *,' *******   BETWEEN THE NVAR VARIABLES  ********'
										PRINT *,' **********************************************'
										PRINT *,' '
										PRINT *,'        dTa-dTf-dRa-dRf-dP-dH-RDa-RDf-dX-dY-dZ-dV '
										DO i=1,nvar
												IF(i.EQ.1)THEN
														PRINT *,' dTa  - ', &
														(int(1000.*corr_var(i,j)/amax1( 0.01,sqrt( corr_var(i,i)*corr_var(j,j)))),j=1,nvar)
												ELSEIF(i.EQ.2)THEN
														PRINT *,' dTf  - ', &
														(int(1000.*corr_var(i,j)/amax1( 0.01,sqrt( corr_var(i,i)*corr_var(j,j)))),j=1,nvar)
												ELSEIF(i.EQ.3)THEN
														PRINT *,' dRa  - ', &
														(int(1000.*corr_var(i,j)/amax1( 0.01,sqrt( corr_var(i,i)*corr_var(j,j)))),j=1,nvar)
												ELSEIF(i.EQ.4)THEN
														PRINT *,' dRf  - ', &
														(int(1000.*corr_var(i,j)/amax1( 0.01,sqrt( corr_var(i,i)*corr_var(j,j)))),j=1,nvar)
												ELSEIF(i.EQ.5)THEN
														PRINT *,' dP   - ', &
														(int(1000.*corr_var(i,j)/amax1( 0.01,sqrt( corr_var(i,i)*corr_var(j,j)))),j=1,nvar)
												ELSEIF(i.EQ.6)THEN
														PRINT *,' dH   - ', &
														(int(1000.*corr_var(i,j)/amax1( 0.01,sqrt( corr_var(i,i)*corr_var(j,j)))),j=1,nvar)
												ELSEIF(i.EQ.7)THEN
														PRINT *,' RDa - ', &
														(int(1000.*corr_var(i,j)/amax1( 0.01,sqrt( corr_var(i,i)*corr_var(j,j)))),j=1,nvar)
												ELSEIF(i.EQ.8)THEN
														PRINT *,' RDf - ', &
														(int(1000.*corr_var(i,j)/amax1( 0.01,sqrt( corr_var(i,i)*corr_var(j,j)))),j=1,nvar)
												ELSEIF(i.EQ.9)THEN
														PRINT *,' dX   - ', &
														(int(1000.*corr_var(i,j)/amax1( 0.01,sqrt( corr_var(i,i)*corr_var(j,j)))),j=1,nvar)
												ELSEIF(i.EQ.10)THEN
														PRINT *,' dY   - ', &
														(int(1000.*corr_var(i,j)/amax1( 0.01,sqrt( corr_var(i,i)*corr_var(j,j)))),j=1,nvar)
												ELSEIF(i.EQ.11)THEN
														PRINT *,' dZ   - ', &
														(int(1000.*corr_var(i,j)/amax1( 0.01,sqrt( corr_var(i,i)*corr_var(j,j)))),j=1,nvar)
												ELSEIF(i.EQ.12)THEN
														PRINT *,' dV   - ', &
														(int(1000.*corr_var(i,j)/amax1( 0.01,sqrt( corr_var(i,i)*corr_var(j,j)))),j=1,nvar)
												ENDIF
										ENDDO
								ELSE
										PRINT *,' !!!! SW_Zsurf+Vsurf+Vinsitu:',sp_zsvsvi,' !!!!'
								ENDIF

								PRINT *,' '

								!**************************************************************************************************
								!**** NORMALIZATION OF THE MATRICES AND VECTORS 
								!**** FOR DZ_surf, VDOP_surf et DVDOP_insitu
								!**** BY THE SUM OF THE POSITIVE VALUES OF THE OBSERVED ERRORS
								!**** THEN BUILD A UNIQUE MATRICE AND VECTOR BY SUMMING 
								!**************************************************************************************************
								PRINT *,' '
								PRINT *,' NORMALIZATION OF THE MATRICES AND VECTORS'
								PRINT *,' SumPosVal_DZsurf,VDOPsurf,DVDOPinsitu:', & 
													swadzsurf_tot,swavsurf_tot,swadvinsitu_tot
								itest_xmat=0
								DO i=1,nvar
										vect(i)=0.
										IF(kdzsurf.EQ.1.AND.swadzsurf_tot.GT.0.) &
										vect(i)=vect(i)+rw_dzsurf*vect_dzsurf(i)/swadzsurf_tot
										IF(kvsurf.EQ.1.AND.swavsurf_tot.GT.0.) &
										vect(i)=vect(i)+rw_vsurf*vect_vsurf(i)/swavsurf_tot
										IF(kdvinsitu.EQ.1.AND.swadvinsitu_tot.GT.0.) &
										vect(i)=vect(i)+rw_dvinsitu*vect_dvinsitu(i)/swadvinsitu_tot
										DO j=1,nvar
												xmat(i,j)=0.
												IF(kdzsurf.EQ.1.AND.swadzsurf_tot.GT.0.) &
												xmat(i,j)=xmat(i,j)+rw_dzsurf*xmat_dzsurf(i,j)/swadzsurf_tot
												IF(kvsurf.EQ.1.AND.swavsurf_tot.GT.0.) &
												xmat(i,j)=xmat(i,j)+rw_vsurf*xmat_vsurf(i,j)/swavsurf_tot
												IF(kdvinsitu.EQ.1.AND.swadvinsitu_tot.GT.0.) &
												xmat(i,j)=xmat(i,j)+rw_dvinsitu*xmat_dvinsitu(i,j)/swadvinsitu_tot
												IF(abs(xmat(i,j)).GT.0.)itest_xmat=1
										ENDDO
										!******************************************************************
										!**** CHECK THAT NO ELEMENT OF THE MATRIX' MAIN DIAGONAL IS NULL
										!******************************************************************
										IF(abs(xmat(i,i)).le.0.)THEN
												DO j=1,nvar
														xmat(i,j)=0.
														xmat(j,i)=0.
												ENDDO
												xmat(i,i)=1.
												vect(i)=0.
										ENDIF
								ENDDO
		
								PRINT *,' '

								!******************************************************************
								!**** INVERSION OF THE MATRIX
								!**** CALCULATION OF THE RESULTING VECTOR
								!******************************************************************
								IF(itest_xmat.EQ.1)THEN            
										CALL resoud(xmat,xinv,vect,res,nvar) ! <-- INGESTIGAR RESULTADO DE CFAC AQUI (ver res)
								ELSE
										PRINT *,' !!!! XMAT=0 !!!!'
								ENDIF              

								!******************************************************************
								!**** ASSIGNEMENT OF THE RESULTS
								!******************************************************************

								! CAI START -- add the function to WRITE out cfac files in SOLO format
								PRINT *,' '
								PRINT *,' '
								PRINT *,' '
								PRINT *,' /////////////////////////////////////////////'
								PRINT *,'     CORRECTIONS FOR NAVIGATIONAL ERRORS'
								PRINT *,' //////////// (add these values)  ////////////'
								PRINT *,' /////////////////////////////////////////////'
								PRINT *,' '
								PRINT *,' '

								WRITE(10,"(' /////////////////////////////////////////////',A)")
								WRITE(10,"('    CORRECTIONS FOR NAVIGATIONAL ERRORS',A)")
								WRITE(10,"(' //////////// (add these values)  ////////////',A)")
								WRITE(10,"(' /////////////////////////////////////////////',A)")

								IF(idtiltaft.EQ.1)THEN
										dtiltaft_res=res(1)
										PRINT *,' D_TILT_aft (deg) guess,residual,total : ', &
															dtiltaft_guess,dtiltaft_res, &
															dtiltaft_guess+dtiltaft_res
										WRITE(10,"(' D_TILT_aft (deg) guess,residual,total : ',3f7.3,/)") &
																dtiltaft_guess,dtiltaft_res, &
																dtiltaft_guess+dtiltaft_res
										! CAI
										tilt_corr_aft = dtiltaft_guess+dtiltaft_res
								ELSE
										dtiltaft_res=0.
										tilt_corr_aft = 0.0
								ENDIF

								IF(idtiltfore.EQ.1)THEN
										dtiltfore_res=res(2)
										PRINT *,' D_TILT_fore (deg) guess,residual,total : ', &
															dtiltfore_guess, dtiltfore_res, dtiltfore_guess+dtiltfore_res
										WRITE(10,"(' D_TILT_fore  (deg) guess,residual,total : ',3f7.3,/)") &
															dtiltfore_guess,dtiltfore_res, dtiltfore_guess+dtiltfore_res
										! CAI
										tilt_corr_fore = dtiltfore_guess+dtiltfore_res
								ELSE
										dtiltfore_res=0.
										tilt_corr_fore = 0.0
								ENDIF

								IF(idrotaaft.EQ.1)THEN
										drotaaft_res=res(3)
										PRINT *,' D_ROTA_aft (deg) guess,residual,total : ', &
															drotaaft_guess,drotaaft_res, &
															drotaaft_guess+drotaaft_res
										WRITE(10,"(' D_dROTA_aft (deg) guess,residual,total : ',3f7.3,/)") &
															drotaaft_guess,drotaaft_res, &
															drotaaft_guess+drotaaft_res
										! CAI
										rot_angle_corr_aft = drotaaft_guess+drotaaft_res
								ELSE
										drotaaft_res=0.
										rot_angle_corr_aft = 0.0
								ENDIF

								IF(idrotafore.EQ.1)THEN
										drotafore_res=res(4)
										PRINT *,' D_ROTA_fore (deg) guess,residual,total : ', &
															drotafore_guess,drotafore_res, &
															drotafore_guess+drotafore_res
										WRITE(10,"(' D_dROTA_fore (deg) guess,residual,total : ',3f7.3,/)") &
															drotafore_guess,drotafore_res, &
															drotafore_guess+drotafore_res
										! CAI
										rot_angle_corr_fore = drotafore_guess+drotafore_res
								ELSE
										drotafore_res=0.
										rot_angle_corr_fore = 0.0
								ENDIF

								IF(idpitch.EQ.1)THEN
										dpitch_res=res(5)
										PRINT *,' D_PITCH (deg) guess,residual,total : ', &
															dpitch_guess,dpitch_res,dpitch_guess+dpitch_res
										WRITE(10,"(' D_PITCH (deg) guess,residual,total : ',3f7.3,/)") &
															dpitch_guess,dpitch_res,dpitch_guess+dpitch_res
										! CAI
										pitch_corr_cfac = dpitch_guess+dpitch_res
								ELSE
										dpitch_res=0.
										pitch_corr_cfac = 0.0
								ENDIF

								IF(idhdg.EQ.1)THEN
										dhdg_res=res(6)
										PRINT *,' D_HEADING (deg) guess,residual,total : ', &
															dhdg_guess,dhdg_res,dhdg_guess+dhdg_res
										WRITE(10,"(' D_HEADING (deg) guess,residual,total : ',3f7.3,/)") &
															dhdg_guess,dhdg_res,dhdg_guess+dhdg_res
										! CAI
										drift_corr_cfac = dhdg_guess+dhdg_res
								ELSE
										dhdg_res=0.
										drift_corr_cfac = 0.0
								ENDIF

								IF(irdaft.EQ.1)THEN
										rdaft_res=100.*res(7)
										PRINT *,' RANGE_DELAY_AFT (m) guess,residual,total : ', &
															1000.*rdaft_guess,rdaft_res, &
															1000.*rdaft_guess+rdaft_res
										WRITE(10,"(' RANGE_DELAY_AFT (m) guess,residual,total : ',3f6.0,/)") &
															1000.*rdaft_guess,rdaft_res, &
															1000.*rdaft_guess+rdaft_res
										! CAI
										rangevalue_delay_corr_aft = 1000.*rdaft_guess+rdaft_res
								ELSE
										rdaft_res=0.
										rangevalue_delay_corr_aft = 0.0
								ENDIF

								IF(irdfore.EQ.1)THEN
										rdfore_res=100.*res(8)
										PRINT *,' RANGE_DELAY_FORE (m) guess,residual,total : ', &
														1000.*rdfore_guess,rdfore_res, &
														1000.*rdfore_guess+rdfore_res
										WRITE(10,"(' RANGE_DELAY_FORE (m) guess,residual,total : ',3f6.0,/)") &
														1000.*rdfore_guess,rdfore_res, &
														1000.*rdfore_guess+rdfore_res
										! CAI
										rangevalue_delay_corr_fore = 1000.*rdfore_guess+rdfore_res
								ELSE
										rdfore_res=0.
										rangevalue_delay_corr_fore = 0.0
								ENDIF

								IF(idxwe.EQ.1)THEN
										dxwe_res=100.*res(9)
										PRINT *,' D_XWE (m) guess,residual,total : ', &
														1000.*dxwe_guess,dxwe_res, &
														1000.*dxwe_guess+dxwe_res
										WRITE(10,"(' D_XWE (m) guess,residual,total : ',3f6.0,/)") &
														1000.*dxwe_guess,dxwe_res, &
														1000.*dxwe_guess+dxwe_res
								ELSE
										dxwe_res=0.
								ENDIF

								IF(idysn.EQ.1)THEN
										dysn_res=100.*res(10)
										PRINT *,' D_YSN (m) guess,residual,total : ', &
														1000.*dysn_guess,dysn_res, &
														1000.*dysn_guess+dysn_res
										WRITE(10,"(' D_YSN (m) guess,residual,total : ',3f6.0,/)") &
														1000.*dysn_guess,dysn_res, &
														1000.*dysn_guess+dysn_res
								ELSE
										dxwe_res=0.
								ENDIF

								IF(idzacft.EQ.1)THEN
										dzacft_res=100.*res(11)
										PRINT *,' D_ZACFT (m) guess,residual,total : ', &
														1000.*dzacft_guess,dzacft_res, &
														1000.*dzacft_guess+dzacft_res
										WRITE(10,"(' D_ZACFT (m) guess,residual,total : ',3f6.0,/)") &
														1000.*dzacft_guess,dzacft_res, &
														1000.*dzacft_guess+dzacft_res
										! CAI
										pressure_alt_corr = 1000.*dzacft_guess+dzacft_res
								ELSE
										dzacft_res=0.
										pressure_alt_corr = 0.0
								ENDIF

								IF(idvh.EQ.1)THEN
										dvh_res=res(12)
										PRINT *,' D_VHACFT (m/s) guess,residual,total : ', &
														dvh_guess,dvh_res,dvh_guess+dvh_res
										WRITE(10,"(' D_VHACFT (m/s) guess,residual,total : ',3f6.2,/)") &
														dvh_guess,dvh_res,dvh_guess+dvh_res
										! CAI
										ew_gndspd_corr = dvh_guess+dvh_res
								ELSE
										dvh_res=0.
										ew_gndspd_corr = 0.0
								ENDIF

								PRINT *,' '
								PRINT *,' '
								PRINT *,' '
								PRINT *,' '
								PRINT *,' //////////////////////////////////////////////////'
								PRINT *,' //////////////////////////////////////////////////'
								PRINT *,' //////////////////////////////////////////////////'

								WRITE(10,"(' /////////////////////////////////////////////',A)")

						ELSE    !!  of  !!  IF(ssurfins.GT.ssurfins_min)THEN  !!

								WRITE(10,"(' /////////////////////////////////////////////',A)")
								WRITE(10,"('    NO CORRECTIONS FOR NAVIGATIONAL ERRORS',A)")
								WRITE(10,"(' //////////// (not enough points) ////////////',A)")
								WRITE(10,"(' /////////////////////////////////////////////',A)")

								PRINT *,' '
								PRINT *,' '
								PRINT *,' '
								PRINT *,' /////////////////////////////////////////////'
								PRINT *,'    NO CORRECTIONS FOR NAVIGATIONAL ERRORS'
								PRINT *,' //////////// (not enough points) ////////////'
								PRINT *,' /////////////////////////////////////////////'
								PRINT *,' '

						ENDIF    !!  of  !!  IF(ssurfins.GT.ssurfins_min)THEN  !!

					  
						PRINT *,' '
						PRINT *,' END OF "CORNAV_EL_*" FILE #10 :',directory(1:ndir)//'/'//fich_cornav
		
						CLOSE(10)

						! CAI
						!******************************************************************
						!             Write the cfac files using SOLO format
						!******************************************************************

						! Write the aft cafc file					
						OPEN(11,file=directory(1:ndir)//'/'//'cfac.aft',form='formatted',status='unknown')
						WRITE(11,400) adjustl(cfac_text01),"=", 0.0
						WRITE(11,400) adjustl(cfac_text02),"=",0.0
						WRITE(11,400) adjustl(cfac_text03),"=",rangevalue_delay_corr_aft
						WRITE(11,400) adjustl(cfac_text04),"=",0.0
						WRITE(11,400) adjustl(cfac_text05),"=",0.0
						WRITE(11,400) adjustl(cfac_text06),"=",pressure_alt_corr
						WRITE(11,400) adjustl(cfac_text07),"=",0.0
						WRITE(11,400) adjustl(cfac_text08),"=",ew_gndspd_corr
						WRITE(11,400) adjustl(cfac_text09),"=",0.0
						WRITE(11,400) adjustl(cfac_text10),"=", 0.0
						WRITE(11,400) adjustl(cfac_text11),"=", 0.0
						WRITE(11,400) adjustl(cfac_text12),"=", 0.0
						WRITE(11,400) adjustl(cfac_text13),"=", pitch_corr_cfac
						WRITE(11,400) adjustl(cfac_text14),"=", drift_corr_cfac
						WRITE(11,400) adjustl(cfac_text15),"=", rot_angle_corr_aft
						WRITE(11,400) adjustl(cfac_text16),"=", tilt_corr_aft
						400		format (A24,A,F8.3)						
						CLOSE(11)

						! Write the fore cafc file
						OPEN(12,file=directory(1:ndir)//'/'//'cfac.fore',form='formatted',status='unknown')
						WRITE(12,"('azimuth_corr			=',f8.3)")0.0
						WRITE(12,"('elevation_corr		=',f8.3)")0.0
						WRITE(12,"('rangevalue_delay_corr	=',f8.3)")rangevalue_delay_corr_fore
						WRITE(12,"('longitude_corr		=',f8.3)")0.0
						WRITE(12,"('latitude_corr			=',f8.3)")0.0
						WRITE(12,"('pressure_alt_corr	=',f8.3)")pressure_alt_corr
						WRITE(12,"('radar_alt_corr		=',f8.3)")0.0
						WRITE(12,"('ew_gndspd_corr	=',f8.3)")ew_gndspd_corr
						WRITE(12,"('ns_gndspd_corr	=',f8.3)")0.0
						WRITE(12,"('vert_vel_corr			=',f8.3)")0.0
						WRITE(12,"('heading_corr			=',f8.3)")0.0
						WRITE(12,"('roll_corr					=',f8.3)")0.0
						WRITE(12,"('pitch_corr				=',f8.3)")pitch_corr_cfac
						WRITE(12,"('drift_corr					=',f8.3)")drift_corr_cfac
						WRITE(12,"('rot_angle_corr		=',f8.3)")rot_angle_corr_fore
						WRITE(12,"('tilt_corr					=',f8.3)")tilt_corr_fore
						CLOSE(12)
						! CAI ******  End of writing the cfac files  ******************

						!******************************************************************
						!**** WRITES THE "SURF_EL*" FILE #30 (IF IWRISURFILE=1)
						!******************************************************************
						IF(iwrisurfile.EQ.1)THEN
								PRINT *,' '
								PRINT *,' WRITES THE "SURF_EL_*" FILE #30 :', &
													directory(1:ndir)//'/'//wrisurfile(1:nsf)
								PRINT *,' INTERPOLATION OF THE RADAR-DERIVED SURFACE MAP'
								CALL inter(swdzsurf_wri,SW_or_altsurf_wri, &
													nx_wrisurf,ny_wrisurf,nxysurfmax)
								nwrisurf_ok=0
								DO j_wrisurf=1,ny_wrisurf
										DO i_wrisurf=1,nx_wrisurf
												IF(abs(SW_or_altsurf_wri(i_wrisurf,j_wrisurf)).lt.10.)THEN
														ialtsurf_wri(i_wrisurf)=(1000.*SW_or_altsurf_wri(i_wrisurf,j_wrisurf))
														nwrisurf_ok=nwrisurf_ok+1
												ELSE
														ialtsurf_wri(i_wrisurf)=-9999
												ENDIF
										ENDDO
										WRITE(30,222)(ialtsurf_wri(i_wrisurf),i_wrisurf=1,nx_wrisurf)
										222  format(500i6)
								ENDDO
								PRINT *,' -> NPTS WRITTEN ON THE "SURF_EL_*" FILE #30',nwrisurf_ok
								CLOSE(30)
						ENDIF

						!******************************************************************
						!**** END OF "SIS_EL_*" FILE #50 
						!******************************************************************
						PRINT *,' '
						PRINT *,' END OF "SIS_EL_*" FILE #50 :',directory(1:ndir)//'/'//fich_sis
						WRITE(50)999,999,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.
						WRITE(50)-999
						WRITE(50)-999
						WRITE(50)-999
						CLOSE(50)

				        !******************************************************************
				        !**** END OF "OUPUT_LOG_*" FILE #60 
				        !******************************************************************
				        IF (ioutputlog.EQ.1) THEN
				        		CLOSE(60)
				        ENDIF
				        
				        !******************************************************************
				        !**** END OF "OUPUT_SURF_*" FILE #70 
				        !******************************************************************
				        CLOSE(70)
				        				        
						!******************************************************************
						!**** END OF PROGRAMM
						!******************************************************************
						PRINT *,' '
						PRINT *,' **************************'
						PRINT *,' **** END OF PROGRAMM  ****'
						PRINT *,' **************************'

						PRINT *,'N1,N2,N3,N4,N5,N6,N7,N8= ',nb1,nb2,nb3,nb4,nb5,nb6,nb7,nb8
						PRINT *,'NSUP=', nsup
						PRINT *,'NTOTAL_OK= ',nbtotals
						PRINT *,'NBON, NMAUVAIS= ',nbon,nmauvais

						GOTO 3

				ENDIF    !!  of  !! IF(iend.EQ.2)THEN  !!


				!******************************************************************
				!**** INITIALIZATIONS AT THE BEGINNING OF A SWEEP (IF IEND=1)
				!******************************************************************
				istart_sweep(iradar_ray)=0
				xsweeps(iradar_ray)=xsweeps(iradar_ray)+1.
				nb_ray(iradar_ray)=0
				stilt(iradar_ray)=0.
				stilt2(iradar_ray)=0.
				rota_prev(iradar_ray)=-999.
				rota_start(iradar_ray)=-999.
				rota_end(iradar_ray)=-999.
				sxa(iradar_ray)=0.
				sya(iradar_ray)=0.
				sza(iradar_ray)=0.
				sacfthspd(iradar_ray)=0.
				stime(iradar_ray)=0.
				ssc(iradar_ray)=0.
				scc(iradar_ray)=0.
				xp_acft(iradar_ray)=0.
				su_acft(iradar_ray)=0.
				sv_acft(iradar_ray)=0.
				SW_acft(iradar_ray)=0.
				xp_wind(iradar_ray)=0.
				su_wind(iradar_ray)=0.
				sv_wind(iradar_ray)=0.
				SW_wind(iradar_ray)=0.
				n_dvinsitu(iradar_ray)=0
				n_dzsurf(iradar_ray)=0
				n_vsurf(iradar_ray)=0
				ndismiss_vhacft(iradar_ray)=0
				ndismiss_vdopcorr(iradar_ray)=0
				ndismiss_vdopsurf(iradar_ray)=0

				DO n=1,500
						zs_rot(iradar_ray,n)=0.
						zs_el(iradar_ray,n)=0.
						zs_az(iradar_ray,n)=0.
						zs_dsurf(iradar_ray,n)=0.
						zs_dhor(iradar_ray,n)=0.
						zs_zsurf(iradar_ray,n)=0.
						zs_hsurf(iradar_ray,n)=0.
						vs_dhor(iradar_ray,n)=0.
						vs_vdopsurf(iradar_ray,n)=0.
						vi_dhor(iradar_ray,n)=0.
						vi_vdop(iradar_ray,n)=0.
						vi_vinsitu(iradar_ray,n)=0.
				ENDDO

				swdzsurf_sweep(iradar_ray)=0.
				dzsurfsweep_mean(iradar_ray)=0.
				dzsurfsweep_rms(iradar_ray)=0.
				swvsurf_sweep(iradar_ray)=0.
				vsurfsweep_mean(iradar_ray)=0.
				vsurfsweep_rms(iradar_ray)=0.
				nsurf_wri(iradar_ray)=0
				swinsitu_sweep(iradar_ray)=0.
				dvinsitusweep_mean(iradar_ray)=0.
				dvinsitusweep_rms(iradar_ray)=0.
				  
				DO jgd=1,2
						s_vpv(iradar_ray,jgd)=0.
						sv_vpv(iradar_ray,jgd)=0.
						svv_vpv(iradar_ray,jgd)=0.
				ENDDO

      ENDIF  !!  of  !! IF(iend.ge.1)THEN

		!************************************************************************
		!**** NEW RAY
		!************************************************************************
		nb_ray(iradar_ray)=nb_ray(iradar_ray)+1

		!******************************************************************
		!**** FRENCH->ENGLISH TRANSLATIONS 
		!**** CONSTANT CORRECTIONS READ ON THE TAPE
		!******************************************************************
		azeast_ray=azimuth+corr_azest(iradar_ray)	! Mod Oliv
		elhor_ray=elevation+corr_elhor(iradar_ray)		!
		xlon_acft=longitude+corr_lon(iradar_ray)			!	radar coord in RADD descriptor
		xlat_acft=latitude+corr_lat(iradar_ray)				!	radar coord in RADD descriptor
		p_alt_acft=altitude+corr_p_alt(iradar_ray)			!
		r_alt_acft=altitude_agl+corr_r_alt(iradar_ray)	!
		roll_acft=roll+corr_roul(iradar_ray)						!
		pitch_acft=pitch+corr_tang(iradar_ray)				!
		hdg_acft=heading+corr_cap(iradar_ray)				!
		drift_acft=drift+corr_derv(iradar_ray)					!
		rota_ray=rotation+corr_rota(iradar_ray)				!
		tilt_ray=tilt+corr_incl(iradar_ray)							!

		

		!******************************************************************
		!**** EARTH-RELATIVE ANGLES AND
		!**** PARAMETERS FOR THE ANALYSIS
		!******************************************************************
		IF(iaftfore.EQ.-1)THEN
				dtilt_guess=dtiltaft_guess
				drota_guess=drotaaft_guess
		ELSE
				dtilt_guess=dtiltfore_guess
				drota_guess=drotafore_guess
		ENDIF
		
		!------------------------------------------------------------------
		!---- ( IF ISIM=1 ) -> SIMULATED TRUE NAVIGATION WITHOUT dXXX_GUESS
		!------------------------------------------------------------------
		IF(isim.EQ.1)THEN
				CALL azel(rota_ray+roll_acft, 							&
										tilt_ray, 										&
										hdg_acft,drift_acft, 					&
										pitch_acft, 								&
										azeast_ray_true,elhor_ray_true, 	&
										cxa_true,cya_true,cza_true, 		&
										cwe_true,csn_true,cnz_true)
				caze_true=caze
				saze_true=saze
				celh_true=celh
				selh_true=selh
				dcwe_dt_true=+crr*sti*spit*shdg-srr*sti*chdg+cti*cpit*shdg
				dcwe_dr_true=+srr*cti*spit*shdg+crr*cti*chdg
				dcwe_dp_true=-crr*cti*cpit*shdg-sti*spit*shdg
				dcwe_dh_true=-crr*cti*spit*chdg-srr*cti*shdg+sti*cpit*chdg
				dcsn_dt_true=+crr*sti*spit*chdg+srr*sti*shdg+cti*cpit*chdg
				dcsn_dr_true=+srr*cti*spit*chdg-crr*cti*shdg
				dcsn_dp_true=-crr*cti*cpit*chdg-sti*spit*chdg
				dcsn_dh_true=+crr*cti*spit*shdg-srr*cti*chdg-sti*cpit*shdg
				dcnz_dt_true=-crr*sti*cpit+cti*spit
				dcnz_dr_true=-srr*cti*cpit
				dcnz_dp_true=-crr*cti*spit+sti*cpit
				dcnz_dh_true=0.
				duacft_dv_true=+shdg*cdri+chdg*sdri
				dvacft_dv_true=+chdg*cdri-shdg*sdri
		ENDIF
		!------------------------------------------------------------------
		
		CALL azel( rota_ray+drota_guess+roll_acft, 		&
								tilt_ray+dtilt_guess, 					&
								hdg_acft+dhdg_guess,drift_acft, 	&
								pitch_acft+dpitch_guess, 			&
								azeast_ray,elhor_ray, 				&
								cxa,cya,cza,cwe,csn,cnz)
		
		IF(sin(conv*(rota_ray+drota_guess+roll_acft)).lt.0.)THEN
				side=-1.
				ilr=1
		ELSE
				side=+1.
				ilr=2
		ENDIF
		
		! Derivatives of radar (Earth-relative) cartesian coordinates (G00, Apendix)
		dcwe_dt=+crr*sti*spit*shdg-srr*sti*chdg+cti*cpit*shdg		! dXrad/dTilt
		dcwe_dr=+srr*cti*spit*shdg+crr*cti*chdg							! dXrad/dR 
		dcwe_dp=-crr*cti*cpit*shdg-sti*spit*shdg							! dXrad/dP
		dcwe_dh=-crr*cti*spit*chdg-srr*cti*shdg+sti*cpit*chdg		! dXrad/dH
		dcsn_dt=+crr*sti*spit*chdg+srr*sti*shdg+cti*cpit*chdg		! dYrad/dTilt
		dcsn_dr=+srr*cti*spit*chdg-crr*cti*shdg							! dYrad/dR
		dcsn_dp=-crr*cti*cpit*chdg-sti*spit*chdg							! dYrad/dP
		dcsn_dh=+crr*cti*spit*shdg-srr*cti*chdg-sti*cpit*shdg		! dYrad/dH
		dcnz_dt=-crr*sti*cpit+cti*spit											! dZrad/dTilt
		dcnz_dr=-srr*cti*cpit														! dZrad/dR
		dcnz_dp=-crr*cti*spit+sti*cpit											! dZrad/dP
		dcnz_dh=0.																	! dZrad/dH
		
		! Earth-relarive Doppler velocity derivatives  (G00, Apendix)
		duacft_dv=+shdg*cdri+chdg*sdri	! duac/dV u-com
		dvacft_dv=+chdg*cdri-shdg*sdri	! dvac/dV v-com

		!******************************************************************
		!**** DISTANCE OF THE RANGE GATES
		!******************************************************************
		IF(iaftfore.EQ.-1)THEN
				d_dgate_guess=rdaft_guess
		ELSE
				d_dgate_guess=rdfore_guess
		ENDIF
		ngates=nb_portes
		
		!------------------------------------------------------------------
		!---- ( IF ISIM=1 ) -> SIMULATED TRUE RANGE GATES WITHOUT dXXX_GUESS
		!------------------------------------------------------------------
		IF(isim.EQ.1)THEN
				DO ig=1,ngates
						! CAI-START
						!          dgate_true(ig)=d_porte(iradar*maxport+ig)
						!    &                        +corr_dist(iradar+1)
						! It seems that the above code IF wrong, since iradar has not been assigned values yet,
						! therefore, following are the new code:
						dgate_true(ig)=d_porte((iradar_ray-1)*maxport+ig)+corr_dist(iradar_ray)
				ENDDO
		ENDIF
		!------------------------------------------------------------------
		DO ig=1,ngates
				dgate_corr(ig)=d_porte((iradar_ray-1)*maxport+ig) & ! Mod Oliv
											+corr_dist(iradar_ray)				&
											+d_dgate_guess 
		ENDDO
		ddg=dgate_corr(2)-dgate_corr(1)

		!******************************************************************
		!**** AIRCRAFT POSITION, (PRESSURE OR RADAR) ALTITUDE AND HEADING
		!******************************************************************
		ylat=(xlat_acft+orig_lat)/2.
		deg_lon=deg_lon0*cos(conv*ylat) !	length (km) of 1 deg of longitude at latitude ylat
		
		!------------------------------------------------------------------
		!---- ( IF ISIM=1 ) -> SIMULATED TRUE AIRCRAFT POSITION WITHOUT dXXX_GUESS
		!------------------------------------------------------------------
		IF(isim.EQ.1)THEN
				x_acft_true=(xlon_acft-orig_lon)*deg_lon
				y_acft_true=(xlat_acft-orig_lat)*deg_lat
				IF(ipr_alt.EQ.1)THEN
						z_acft_true=p_alt_acft
				ELSE
						z_acft_true=r_alt_acft
				ENDIF
		ENDIF
		!------------------------------------------------------------------
		
		x_acft=(xlon_acft-orig_lon)*deg_lon+dxwe_guess ! [km], x position relative to origin; positive eastward
		y_acft=(orig_lat-xlat_acft)*deg_lat+dysn_guess	! [km], y position relative to origin; positive southward
		IF(ipr_alt.EQ.1)THEN
				z_acft=p_alt_acft+dzacft_guess	! [km]
		ELSE
				z_acft=r_alt_acft+dzacft_guess	! [km]
		ENDIF

		!	PRINT *,'Z_ACFT,P_ALT,D_GUESS= ',z_acft,p_alt_acft,dzacft_guess

		!******************************************************************
		!**** ADD TO THE MEAN PARAMETERS FOR THE CURRENT SWEEP
		!******************************************************************
		stilt(iradar_ray)=stilt(iradar_ray)+tilt_ray											! sum of tilt
		stilt2(iradar_ray)=stilt2(iradar_ray)+tilt_ray*tilt_ray						! sum of tilt squared
		IF(nb_ray(iradar_ray).EQ.1)rota_start(iradar_ray)=rota_ray		
		sxa(iradar_ray)=sxa(iradar_ray)+x_acft											! sum of x aircraft position
		sya(iradar_ray)=sya(iradar_ray)+y_acft											! sum of y aircraft position
		sza(iradar_ray)=sza(iradar_ray)+z_acft											! sum of z aircraft position
		stime(iradar_ray)=stime(iradar_ray)+time_ks 								! sum of time
		ssc(iradar_ray)=ssc(iradar_ray)+shdg											! sum of sin(heading)
		scc(iradar_ray)=scc(iradar_ray)+chdg											! sum of cos(heading)
		dmax=amin1(dmax0,dgate_corr(ngates))

		!******************************************************************
		!**** AIRCRAFT SPEED
		!******************************************************************
		IF( (abs(acftspd_we).LT.10..AND.abs(acftspd_sn).LT.10.) .OR. &
				(abs(acftspd_we).GT.200..OR.abs(acftspd_sn).GT.200.) )THEN
				PRINT *,' !!!! NO_RAY:',nb_ray, &
								' -> U,V,W_acft:',acftspd_we,acftspd_sn,acftspd_nz,' !!!!'
				!WRITE(60,'(A)') '    line 2131		NO RAY!!'
				GOTO 1
		ENDIF
		!----------------------------------------------------------------------
		!---- ( IF ISIM=1 ) -> SIMULATED TRUE AIRCRAFT SPEED WITHOUT dXXX_GUESS
		!----------------------------------------------------------------------
		IF(isim.EQ.1)THEN
				acftspd_we_true=acftspd_we
				acftspd_sn_true=acftspd_sn
		ENDIF
		!----------------------------------------------------------------------
		acftspd_we=acftspd_we+duacft_dv*dvh_guess											! [m/s] (dvh_guess=0)
		acftspd_sn=acftspd_sn+dvacft_dv*dvh_guess											! [m/s] (dvh_guess=0)
		acftspd_hor=sqrt(acftspd_we*acftspd_we+acftspd_sn*acftspd_sn)	! [m/s]
		sacfthspd(iradar_ray)=sacfthspd(iradar_ray)+acftspd_hor						! sum of horizontal acft speed
		xp_acft(iradar_ray)=xp_acft(iradar_ray)+1.													! counter for aircraft
		su_acft(iradar_ray)=su_acft(iradar_ray)+acftspd_we									! sum of acft u-comp speed
		sv_acft(iradar_ray)=sv_acft(iradar_ray)+acftspd_sn									! sum of acft v-comp speed
		SW_acft(iradar_ray)=SW_acft(iradar_ray)+acftspd_nz								! sum of acft vertical speed
		proj_acftspd=acftspd_we*cwe+acftspd_sn*csn+acftspd_nz*cnz			! 3D projected acft speed

		!******************************************************************
		!**** FLIGHT-LEVEL WIND
		!******************************************************************

		!----------------------------------------------------------------------
		!---- ( IF ISIM=1 ) -> SIMULATED TRUE AIRCRAFT SPEED WITHOUT dXXX_GUESS
		!----------------------------------------------------------------------
		IF(isim.EQ.1)THEN
				proj_wind=wind_we*cwe_true+wind_sn*csn_true+wind_nz*cnz_true
				wa_we_true=wind_we-acftspd_we_true
				wa_sn_true=wind_sn-acftspd_sn_true
		ENDIF
		!----------------------------------------------------------------------
		IF( (abs(wind_we).LE.0..AND.abs(wind_sn).LE.0.) .OR. &
				(abs(wind_we).GT.100..OR.abs(wind_sn).GT.100.))THEN
				PRINT *,' !!!! NO_RAY:',nb_ray(iradar_ray),' -> Uwe,Vsn_wind:', &
									wind_we,wind_sn,' !!!!'
				!WRITE(60,'(A)') '    line 2169		NO RAY!!'									
				GOTO 1
		ENDIF
      
		IF(abs(wind_nz).LE. 0. .OR.abs(wind_nz).GT.50.)wind_nz=0.
		xp_wind(iradar_ray)=xp_wind(iradar_ray)+1.						! counter of wind speed
		su_wind(iradar_ray)=su_wind(iradar_ray)+wind_we			! sum of u-wind spd
		sv_wind(iradar_ray)=sv_wind(iradar_ray)+wind_sn			! sum of v-wind spd
		SW_wind(iradar_ray)=SW_wind(iradar_ray)+wind_nz			! sum of w-wind spd
		proj_wind=wind_we*cwe+wind_sn*csn+wind_nz*cnz		! 3D projected wind speed
		wa_we=wind_we-acftspd_we													! difference btw u-wind spd and u-acft speed
		wa_sn=wind_sn-acftspd_sn														! difference btw v-wind spd and v-acft speed
		wa_nz=wind_nz-acftspd_nz														! difference btw w-wind spd and w-acft speed

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!!!! ELIMINATION OF GATES CONTAMINATED WITH GROUND-CLUTTER
		!!!! ONLY FOR TOGA-COARE DATA !!!!
		!!!!  -> aft_SWEEP : dROTA=+6 deg
		!!!!  -> fore_SWEEP : dROTA=+3 deg)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!!!!
		!!!!     IF(tilt_ray.lt.-15.)THEN
		!!!!       rota_sidelobe=rota_ray+roll_acft+6.
		!!!!     ELSEIF(tilt_ray.GT.15.)THEN
		!!!!       rota_sidelobe=rota_ray+roll_acft+3.
		!!!!     ENDIF
		!!!!      IF(a_don.le.1993.AND.cos(conv*rota_sidelobe).lt.0.)THEN
		!!!!       dmax_sidelobe=(z_min-z_acft)/cos(conv*rota_sidelobe)
		!!!!       dmax=amin1(dmax0,dmax_sidelobe)
		!!!!     ELSE
		!!!!       dmax=dmax0
		!!!!     ENDIF
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!******************************************************************
		!**** DISMISS THE SPECIFIED RANGE GATES
		!******************************************************************
		DO iig=1,15
				IF(ig_dismiss(iig).GT.0)THEN
						ig=ig_dismiss(iig)
						ZE(ig)=-999.
						VR(ig)=-999.
						vs(ig)=-999.	!Olivier
						vl(ig)=-999.	!Olivier
						vg(ig)=-999.
						vu(ig)=-999.
				ENDIF
		ENDDO

		!******************************************************************
		!**** RANGE GATES FOR COMPARISONS WITH FLIGHT-LEVEL (IN SITU) DATA
		!******************************************************************
		ngates_insitu_max=-999
		IF(abs(selh).LT.selhinsitu_max)THEN
				ig=1
				DO WHILE (ig.LT.maxport.AND.dgate_corr(ig).LT.dmax_insitu)
						ngates_insitu_max=ig
						ig=ig+1
				ENDDO
		ENDIF

		!******************************************************************
		!**** CHECK THE NCP, SW AND REFLECTIVITY VALUES
		!******************************************************************
		ngates_max=1
		DO ig=1,ngates
				! ref_min=ref_min0+20.*(log10(dgate_corr(ig))-1.) ! ?? (RV)
				IF(dgate_corr(ig).LT.dmin.OR.	dgate_corr(ig).GT.dmax.OR. &
								NCP(ig).LT.xncp_min.OR. &
								SW(ig).GT.SW_max .OR. &
								ZE(ig).LT.ref_min0.OR.	ZE(ig).GT.ref_max)THEN
						ZE(ig)=-999.
						VR(ig)=-999.
						vs(ig)=-999.	!Olivier
						vl(ig)=-999.	!Olivier
						vg(ig)=-999.
						vu(ig)=-999.
				ELSE
						ngates_max=ig														! farthest gate with non NaN value 
						nref_ok(iradar_ray)=nref_ok(iradar_ray)+1
				ENDIF
		ENDDO
		
		dum=1 ! dummy for debugging
		
10	CONTINUE

		!******************************************************************
		!**** CHOOSE WHICH DOPPLER VELOCITY WILL BE USED (FOLLOWING ICHOICE_VDOP)
		!****   -> 1:RAW(VR), 2:CORRECTED FOR VACFT(VG), 3:UNFOLDED(VU)
		!******************************************************************
		DO ig=1,ngates_max
				vdop_READ=-999.
				vdop_corr(ig)=-999.
				IF(ZE(ig).GT.-900.)THEN
						! CAI-START: get rid of all Oliver's modification, since it is for P3
						!          IF(ichoice_vdop.EQ.1)THEN		! Olivier
						!            IF(     abs(VR(ig)).GT.0.		
						!    &          .AND.abs(VR(ig)).lt.vdop_max	
						!    &          .AND.abs(vs(ig)).GT.0.		
						!    &          .AND.abs(vs(ig)).lt.vdop_max	
						!    &          .AND.abs(vl(ig)).GT.0.		
						!    &          .AND.abs(vl(ig)).lt.vdop_max	
						!    &          .AND.proj_acftspd.GT.-900.   )THEN
						!
						!               d_vs_vl=vs(ig)-vl(ig)		! Olivier
						!               kvsl=ifix((d_vs_vl/vnyq_el)*0.5)+5	! Olivier
						!
						!               IF(kvsl.ge.1.AND.kvsl.le.9)THEN		! Olivier
						!                 vs_depl=vs(ig)+xms(kvsl)*vnyq_el	
						!                 vl_depl=vl(ig)+xml(kvsl)*vnyq_el	
						!                 vsl_depl=(vs_depl+vl_depl)/2.		! Olivier
						!
						!               IF(    abs(vs_depl-vl_depl).lt.vnyq_el/2.
						!    &		    .AND.abs(VR(ig)-vsl_depl).lt.vnyq_el/2.     )THEN 
						!          	    vdop_READ=VR(ig)+proj_acftspd
						!                   nbon=nbon+1
						!               ELSE
						!                   nmauvais=nmauvais+1
						!               ENDIF
						!
						!                  nbtotals=nbtotals+1
						!
						!              ENDIF
						!
						!             ENDIF					! Olivier
						!            ENDIF					! Olivier
						!
						!	      IF(     abs(VR(ig)).GT.0.			
						!    &		 .AND.abs(VR(ig)).lt.vdop_max		
						!    &	         .AND.proj_acftspd.GT.-900.     )THEN	! Olivier
						!		 vdop_READ=VR(ig)+proj_acftspd		! Olivier
						!	      ENDIF					! Olivier
						!	    ENDIF					! Olivier

						IF(ichoice_vdop.EQ.1.AND.abs(VR(ig)).GT.0..AND.abs(VR(ig)).LT.vdop_max &
						     .AND.proj_acftspd.GT.-900.)  THEN 
						     vdop_READ=VR(ig)+proj_acftspd
								! CAI-STOP
						ELSEIF(ichoice_vdop.EQ.2.AND.abs(vg(ig)).GT.0..AND.abs(vg(ig)).LT.vdop_max) THEN
						     vdop_READ=vg(ig)
						ELSEIF(ichoice_vdop.EQ.3.AND.abs(vu(ig)).GT.0..AND.abs(vu(ig)).LT.vdop_max) THEN
							 vdop_READ=vu(ig)
						ENDIF
						
						IF(vdop_READ.GT.-900.)THEN
								ndop_ok(iradar_ray)=ndop_ok(iradar_ray)+1
								vdop_corr(ig)=vdop_READ
						ENDIF
				ENDIF
		ENDDO

		dum=1 ! dummy for debugging

		!******************************************************************
		!**** ( IF     ( KZsurf=1  or  KVsurf=1 )
		!****      AND  Z_ACFT > Z_ACFTmin
		!****      AND SIN(ELEV_HOR) < SELH_SURF 
		!****      AND  VFF_AV>0 )
		!****  -> DETERMINE ALTITUDE (THEN DOPPLER VELOCITY)
		!****     OF THE SURFACE FOR THIS RAY
		!******************************************************************
		IF((kdzsurf+kvsurf).GE.1.AND.selh.LT.selh_surf .AND.z_acft.GT.zacftmin_surf)THEN
					  
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				IF(iprintscreen.GT.0) THEN
						IF (nb_ray(iradar_ray).EQ.iprintscreen*(nb_ray(iradar_ray)/iprintscreen))THEN
								PRINT *,' '
								PRINT *,' ',1000*ihhmmss+ims_ray, &
												' IRADAR:',iradar_ray, &
												' NO_RAY:',nb_ray(iradar_ray)
								PRINT *,'    ROTA,TILT_RAY:',rota_ray,tilt_ray
								PRINT *,'    ROLL,PITCH,HDG,DRIFT_ACFT:',roll_acft, &
												pitch_acft,hdg_acft,drift_acft
								PRINT *,'    AZ_EAST:',azeast_ray,' EL_HOR:',elhor_ray
								PRINT *,'    CWE,CSN,CNZ:',cwe,csn,cnz
								PRINT *,'    U,V,W_ACFT:',acftspd_we,acftspd_sn,acftspd_nz, &
												' PROJ_VACFT:',proj_acftspd
								PRINT *,'HOLA_2341'
						ENDIF
				ENDIF

				!=====WRITE TO OUTPUT LOG FILE====
				IF (ioutputlog.EQ.1) THEN
						WRITE(60,600) 'line_2341'
						WRITE(60,601) 1000*ihhmmss+ims_ray, &
													' IRADAR:',iradar_ray, &
													' NO_RAY:',nb_ray(iradar_ray)
						WRITE(60,602) '    ROTA,TILT_RAY:',rota_ray,tilt_ray
						WRITE(60,603) '    ROLL,PITCH,HDG,DRIFT_ACFT:',roll_acft, &
															pitch_acft,hdg_acft,drift_acft
						WRITE(60,604) '    AZ_EAST:',azeast_ray,' EL_HOR:',elhor_ray
						WRITE(60,605) '    CWE,CSN,CNZ:',cwe,csn,cnz
						WRITE(60,606) '    U,V,W_ACFT:',acftspd_we,acftspd_sn,acftspd_nz, &
												' PROJ_VACFT:',proj_acftspd						
						605   format(A,3F10.2)						
						606   format(A,3F10.2,A,F10.2)			
				ENDIF					
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				dsurf_ray=-999.
				dhorsurf_ray=-999.
				hsurf_ray=-999.
				refsurf_ray=0.
				gradrefsurf_ray=0.
				refsurf_min=refsurf_min0*((abs(selh))**0.7)						! Georgis et al (2000)				
				gradrefsurf_min=gradrefsurf_min0*((abs(selh))**0.7)		! Georgis et al (2000)
				wghtsurf_ray=0.
				refmax_ray=-999.
				ig_refmax=999		! Olivier (real->entier)
				d_refmax=-999.
				h_refmax=-999.
				z_refmax=-999.
				gradrefmax_ray=-999.
				ig_gradrefmax=999	! Olivier (real->entier)
				d_gradrefmax=-999.
				h_gradrefmax=-999.
				z_gradrefmax=-999.

				!******************************************************************
				!**** DETERMINE REFmax AND dREF/dD|max
				!******************************************************************
				distmax=(z_acft-altdtm_min+1.)/abs(selh)	! distance btw radar and ground echo
																				! T95 Rg Eq. 4 (why +1?)
				dhor_prevgate=0.	! Mod Oliv
				z_prevgate=0.		! Mod Oliv
				DO ig=igstart_surf,ngates_max
						IF(dgate_corr(ig).LE.distmax)THEN		! if range gate is less than distance to ground
								d_ig=dgate_corr(ig)
								dver_ig=d_ig*selh									! vertical location of beam volume relative to acft
								dhor_ig=d_ig*celh									! horizontal location of beam volume relative to acft
								frac1=2.*(z_acft+dver_ig)/rayter
								frac2=(z_acft*z_acft+d_ig*d_ig+2.*z_acft*dver_ig)/(rayter*rayter)
								z_ig=rayter*(sqrt(1.+frac1+frac2)-1.)	! [km] altitude of the beam volume 
								theta=atan(dhor_ig/(rayter+z_acft+dver_ig))
								dhor_ig=rayter*theta							! [km] arc distance (L=theta*radius)
								IF(ZE(ig).GT.-900.)THEN
										IF(ZE(ig).GT.refmax_ray)THEN
												refmax_ray=ZE(ig)				! max dBZ
												ig_refmax=ig						! gate number of max dBZ
												d_refmax=d_ig					! [km] gate range of max dBZ
												dhor_refmax=dhor_ig			! horizontal distance of max dBZ
												z_refmax=z_ig					! altitude of max dBZ
										ENDIF
										IF(ig.GT.1.AND.ZE(ig-1).GT.-900.)THEN
												gradref=(ZE(ig)-ZE(ig-1))/(d_ig-dgate_corr(ig-1))
												IF(gradref.GT.gradrefmax_ray)THEN
														gradrefmax_ray=gradref
														ig_gradrefmax=ig
														d_gradrefmax=(d_ig+dgate_corr(ig-1))/2.
														dhor_gradrefmax=(dhor_prevgate+dhor_ig)/2.
														z_gradrefmax=(z_prevgate+z_ig)/2.
												ENDIF
										ENDIF
								ENDIF
								z_prevgate=z_ig
								dhor_prevgate=dhor_ig
						ENDIF  
				ENDDO

				!*************************************************************************************
				!**** WEIGHT ASSOCIATED WITH THE OBTAINED SURFACE POINT
				!*************************************************************************************
				! If the ray contains a surface gate
				IF(refmax_ray.GT.refsurf_min.AND.gradrefmax_ray.GT.gradrefsurf_min)THEN
						IF((d_refmax.GT.d_gradrefmax).AND.abs(z_refmax-z_gradrefmax).LT.1.)THEN
								wght_ref=1.+(refmax_ray-refsurf_min)/refsurf_min
								wght_grad=1.+(gradrefmax_ray-gradrefsurf_min)/(gradrefsurf_min)
								wghtsurf_ray=sqrt(wght_ref*wght_grad)
								dsurf_ray=d_refmax
								hsurf_ray=z_refmax
								dhorsurf_ray=dhor_refmax
								xsurf_ray=x_acft+dhorsurf_ray*caze
								ysurf_ray=y_acft+dhorsurf_ray*saze
								!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
								!!!!      IF( nb_ray(iradar_ray).EQ.10*(nb_ray(iradar_ray)/10) )THEN
								!!!!        PRINT *,'     -> SURF: d,dhor,z:'
								!!!!     &         ,dsurf_ray,dhorsurf_ray,hsurf_ray
								!!!!        PRINT *,'     -> X,Y,H_SURF:',xsurf_ray,ysurf_ray,hsurf_ray
								!!!!        PRINT *,'        WGHTSURF_ray :',wghtsurf_ray
								!!!!      ENDIF
								!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
								!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
								!!!!          ELSE
								!!!!      IF( nb_ray(iradar_ray).EQ.10*(nb_ray(iradar_ray)/10) )
								!!!!     &  PRINT *,'    !!!! VALUES OK, BUT PB ON d_REF AND/OR d_GRAD !!!'
								!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						ENDIF
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						!!!!        ELSE
						!!!!      IF( nb_ray(iradar_ray).EQ.10*(nb_ray(iradar_ray)/10) )
						!!!!     &  PRINT *,'    !!!! PB ON REF AND/OR GRAD VALUES !!!'
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				ENDIF

				!--------------------------------------------------------------
				!---- ( IF ISIM=1 ) -> (XYH)SURF_RAY FROM dXXX_GUESS & DTM
				!--------------------------------------------------------------
				IF(isim.EQ.1)THEN
						!---- *_app -> with dXXX_GUESS ("apparent" navigation)
						!---- *_true -> without ("true" navigation)
						hsurf_ray=-999.
						wghtsurf_ray=0.
						igsurf_ray=-999
						DO ig=igstart_surf,ngates_max
								IF(dgate_corr(ig).LE.distmax)THEN
										d_app=dgate_corr(ig)
										dver_app=d_app*selh
										dhor_app=d_app*celh
										frac1=2.*(z_acft+dver_app)/rayter
										frac2=(z_acft*z_acft+d_app*d_app+2.*z_acft*dver_app)/(rayter*rayter)
										z_app=rayter*(sqrt(1.+frac1+frac2)-1.)
										theta=atan(dhor_app/(rayter+z_acft+dver_app))
										dhor_app=rayter*theta
										x_app=x_acft+dhor_app*caze
										y_app=y_acft+dhor_app*saze
										d_true=dgate_true(ig)
										dver_true=d_true*selh_true
										dhor_true=d_true*celh_true
										frac1=2.*(z_acft_true+dver_true)/rayter
										frac2=(z_acft_true*z_acft_true+d_true*d_true+2.*z_acft_true*dver_true)/(rayter*rayter)
										z_true=rayter*(sqrt(1.+frac1+frac2)-1.)
										theta=atan(dhor_true/(rayter+z_acft_true+dver_true))
										dhor_true=rayter*theta
										x_true=x_acft_true+dhor_true*caze_true
										y_true=y_acft_true+dhor_true*saze_true
										IF(igsurf_ray.EQ.-999.AND.x_true.GT.xmin_dtm &
													.AND.x_true.lt.xmax_dtm.AND.y_true.GT.ymin_dtm &
													.AND.y_true.lt.ymax_dtm)THEN
												isurf_true=(x_true-xmin_dtm)/hx_dtm+1
												jsurf_true=(y_true-ymin_dtm)/hy_dtm+1
												aa=alt_dtm(isurf_true,jsurf_true)
												bb=(-alt_dtm(isurf_true,jsurf_true)+alt_dtm(isurf_true+1,jsurf_true))/hx_dtm
												cc=(-alt_dtm(isurf_true,jsurf_true)+alt_dtm(isurf_true,jsurf_true+1))/hy_dtm
												dd=(+alt_dtm(isurf_true,jsurf_true) 				&
														-alt_dtm(isurf_true+1,jsurf_true) 			&
														-alt_dtm(isurf_true,jsurf_true+1) 			&
														+alt_dtm(isurf_true+1,jsurf_true+1)) 	&
														/(hx_dtm*hy_dtm)
												x_dtm=xmin_dtm+float(isurf_true-1)*hx_dtm
												dx=x_true-x_dtm
												y_dtm=ymin_dtm+float(jsurf_true-1)*hy_dtm
												dy=y_true-y_dtm
												hsurf_dtm=aa+bb*dx+cc*dy+dd*dx*dy
												IF(hsurf_dtm.GE.z_true)THEN
														xsurf_true=x_true
														ysurf_true=y_true
														hsurf_true=z_true
														dxh_dtm=bb+dd*dy
														dyh_dtm=cc+dd*dx
														igsurf_ray=ig
														dsurf_ray=d_app
														xsurf_ray=x_app
														ysurf_ray=y_app
														hsurf_ray=z_app
														wghtsurf_ray=1.
												ENDIF
										ENDIF
								ENDIF
						ENDDO
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						!IF( nb_ray(iradar_ray).EQ.10*(nb_ray(iradar_ray)/10) )THEN
								PRINT *,'     -> X,Y,H_SURF-TRUE:', xsurf_true,ysurf_true,hsurf_true
								PRINT *,'     ->I,J_SURF_true :',isurf_true,jsurf_true, ' dxH,dyH :',dxh_dtm,dyh_dtm
								PRINT *,'     -> X,Y,H_SURF-RAY:', xsurf_ray,ysurf_ray,hsurf_ray
								PRINT *,'HOLA_2529'																
						!ENDIF
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				ENDIF
				!--------------------------------------------------------------

				!******************************************************************
				!**** IF THIS SURFACE POINT IS CORRECT (IF WGHTSURF_ray > 0)
				!**** THEN COMPARE WITH THE SURFACE POINT DERIVED FROM THE DTM
				!******************************************************************
				IF(hsurf_ray.GT.-900..AND.wghtsurf_ray.GT.0.)THEN

						IF(iprintscreen.GT.0) THEN
								IF (nb_ray(iradar_ray).EQ.iprintscreen*(nb_ray(iradar_ray)/iprintscreen))THEN
										PRINT *,'    -> X,Y,H_SURF-RAY:',xsurf_ray,ysurf_ray,hsurf_ray
										PRINT *,'HOLA_2542'
								ENDIF
						ENDIF
		
						!=====WRITE TO OUTPUT LOG FILE====
						IF (ioutputlog.EQ.1) THEN
								WRITE(60, '(A,2F10.7,X,F10.5)') 'line_2542  -> X,Y,H_SURF-RAY:',xsurf_ray,ysurf_ray,hsurf_ray
						ENDIF		

						!******************************************************************************************************
						!**** INTERPOLATION OF ALT_DTM(x,y) [READ ON SURF_DTM_* OR CONSTANT]
						!******************************************************************************************************
						! If surface ray position is within the DTM domain (RV)
						IF(xsurf_ray.GT.xmin_dtm.AND.xsurf_ray.LT.xmax_dtm.AND. &
							ysurf_ray.GT.ymin_dtm.AND.ysurf_ray.LT.ymax_dtm)THEN
								isurf_ray=(xsurf_ray-xmin_dtm)/hx_dtm+1
								jsurf_ray=(ysurf_ray-ymin_dtm)/hy_dtm+1						
								aa=alt_dtm(isurf_ray,jsurf_ray)
								bb=(-alt_dtm(isurf_ray,jsurf_ray)+alt_dtm(isurf_ray+1,jsurf_ray))/hx_dtm
								cc=(-alt_dtm(isurf_ray,jsurf_ray)+alt_dtm(isurf_ray,jsurf_ray+1))/hy_dtm
								dd=(+alt_dtm(isurf_ray,jsurf_ray) &
										-alt_dtm(isurf_ray+1,jsurf_ray) &
										-alt_dtm(isurf_ray,jsurf_ray+1) &
										+alt_dtm(isurf_ray+1,jsurf_ray+1)) &
										/(hx_dtm*hy_dtm)
								x_dtm=xmin_dtm+float(isurf_ray-1)*hx_dtm
								dx=xsurf_ray-x_dtm
								y_dtm=ymin_dtm+float(jsurf_ray-1)*hy_dtm
								dy=ysurf_ray-y_dtm
								hsurf_dtm=aa+bb*dx+cc*dy+dd*dx*dy
								d_hsurf=hsurf_ray-hsurf_dtm
								dxh_dtm=bb+dd*dy
								dyh_dtm=cc+dd*dx
								!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!							    
								IF(iprintscreen.GT.0) THEN
										IF (nb_ray(iradar_ray).EQ.iprintscreen*(nb_ray(iradar_ray)/iprintscreen))THEN
												PRINT *,'        I,J_SURF_ray :',isurf_ray,jsurf_ray, ' dxH,dyH :',dxh_dtm,dyh_dtm
												PRINT *,'     -> H_SURF-DTM:',hsurf_dtm, '  =>> D_HSURF :',d_hsurf
												PRINT *,'HOLA_2575'
										ENDIF
								ENDIF
														
								!=====WRITE TO OUTPUT LOG FILE====
								IF (ioutputlog.EQ.1) THEN
										WRITE(60, '(A,2I4,A,2F5.2)') 'line_2575  -> I,J_SURF_ray :',isurf_ray,jsurf_ray,' dxH,dyH :',dxh_dtm,dyh_dtm
										WRITE(60, '(A,F8.3,A,F8.3)') 'line_2575  -> H_SURF-DTM:',hsurf_dtm,'  =>> D_HSURF :',d_hsurf
										WRITE(60, '(A,4F8.3)') 'line_2575  -> aa,bb,cc,dd:',aa,bb,cc,dd 
								ENDIF																					    
								!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

								!******************************************************************
								!**** IF ( ABS(HSURF_RADAR-HSURF_DTM) < DHSURF_MAX ) THEN
								!**** !!!! DHSURF_MAX=999. !!!! -> NOT IN USE !!!!
								!******************************************************************
								IF(abs(d_hsurf).LT.dhsurf_max)THEN
								
										ssurfins=ssurfins+wghtsurf_ray
										
										
										!*********** SAVE SURF POINTS TO TEXT FILE HERE ************************************************
										WRITE(70,700) counter, 										&
																isurf_ray,jsurf_ray, 						&
																alt_dtm(isurf_ray,jsurf_ray)*1000, 	&
																xsurf_ray,ysurf_ray, 						&
																x_acft,y_acft, 								&
																dhorsurf_ray, azeast_ray, 				&
																dhor_ig,dver_ig,z_acft, 					&
																roll,d_ig
										700   format(I5,2I4,12F7.2)
										
										
										!******************************************************************
										!**** CASE "DZ_surf"
										!******************************************************************
										IF(kdzsurf.EQ.1)THEN
												!----------------------------------------------------------------------
												!---- ( IF ISIM=1 ) -> SIMULATED DZ_surf FROM dXXX_GUESS
												!----------------------------------------------------------------------
												IF(isim.EQ.1)THEN
												d_hsurf_dxxx=-dsurf_ray*(-dcnz_dt+dxh_dtm*dcwe_dt+dyh_dtm*dcsn_dt) 						&
																								  *dtilt_guess*conv-dsurf_ray 											&
																								 *(-dcnz_dr+dxh_dtm*dcwe_dr+dyh_dtm*dcsn_dr) 			&
																								  *drota_guess*conv-dsurf_ray 										&
																								 *(-dcnz_dp+dxh_dtm*dcwe_dp+dyh_dtm*dcsn_dp) 			&
																								  *dpitch_guess*conv-dsurf_ray 										&
																								 *(-dcnz_dh+dxh_dtm*dcwe_dh+dyh_dtm*dcsn_dh) 			&
																								  *dhdg_guess*conv-(-cnz+dxh_dtm*cwe+dyh_dtm*csn)	&
																								 *d_dgate_guess-dxh_dtm*dxwe_guess 							&
																								-dyh_dtm*dysn_guess+dzacft_guess
												!!!!            d_hsurf=d_hsurf_dxxx
												ENDIF
												!----------------------------------------------------------------------
												!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
												IF(iprintscreen.GT.0) THEN
														IF (nb_ray(iradar_ray).EQ.iprintscreen*(nb_ray(iradar_ray)/iprintscreen))THEN
																PRINT *,'     -> D_HSURF_dXXX :',d_hsurf_dxxx
																PRINT *,'HOLA_2624'
														ENDIF
												ENDIF
																														
														!=====WRITE TO OUTPUT LOG FILE====
												IF (ioutputlog.EQ.1) THEN
														WRITE(60, '(A,F9.5)') 'line_2624  -> D_HSURF_dXXX :',d_hsurf_dxx
												ENDIF												
												!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

												!******************************************************************
												!**** ADD WEIGHTS AND DZ_surf
												!****************************************************************** 
												n_dzsurf(iradar_ray)=n_dzsurf(iradar_ray)+1
												swdzsurf_sweep(iradar_ray)=swdzsurf_sweep(iradar_ray)+wghtsurf_ray
												dzsurfsweep_mean(iradar_ray)=dzsurfsweep_mean(iradar_ray) &
																											  +wghtsurf_ray*d_hsurf
												dzsurfsweep_rms(iradar_ray)=dzsurfsweep_rms(iradar_ray) &
																											+wghtsurf_ray*d_hsurf*d_hsurf
												swdzsurf_tot=swdzsurf_tot+wghtsurf_ray
												swdzmsurf_tot=swdzmsurf_tot+wghtsurf_ray*d_hsurf
												swdz2surf_tot=swdz2surf_tot+wghtsurf_ray*d_hsurf*d_hsurf	   
												swadzsurf_tot=swadzsurf_tot+wghtsurf_ray*abs(d_hsurf)

												!******************************************************************
												!**** VALUES OF VAR(1->NVAR) FOR FIELD "DZ_surf"
												!****  - VAR(1->6) -> [dT_aft,dT_fore,dR_aft,dR_fore,dP,dH] in DEGREES
												!****  - VAR(7->11) -> [RD_aft,RD_fore,dXwe,dYsn,dZ] in HECTOMETERS
												!****  - VAR(12) -> [dVH] in METER/SECOND
												!******************************************************************
												IF(iaftfore.EQ.-1)THEN
														IF(idtiltaft.EQ.1)THEN
																var(1)=dsurf_ray*(-dcnz_dt+dxh_dtm*dcwe_dt+dyh_dtm*dcsn_dt)*conv
														ELSE
																var(1)=0.
																xmat_dzsurf(1,1)=xmat_dzsurf(1,1)+wghtsurf_ray
														ENDIF
														var(2)=0.
												ELSE
														var(1)=0.
														IF(idtiltfore.EQ.1)THEN
																var(2)=dsurf_ray*(-dcnz_dt+dxh_dtm*dcwe_dt+dyh_dtm*dcsn_dt)*conv
														ELSE
																var(2)=0.
																xmat_dzsurf(2,2)=xmat_dzsurf(2,2)+wghtsurf_ray
														ENDIF
												ENDIF

												IF(iaftfore.EQ.-1)THEN
														IF(idrotaaft.EQ.1)THEN
																var(3)=dsurf_ray*(-dcnz_dr+dxh_dtm*dcwe_dr+dyh_dtm*dcsn_dr)*conv
														ELSE
																var(3)=0.
																xmat_dzsurf(3,3)=xmat_dzsurf(3,3)+wghtsurf_ray
														ENDIF
														var(4)=0.
												ELSE
														var(3)=0.
														IF(idrotafore.EQ.1)THEN
																var(4)=dsurf_ray*(-dcnz_dr+dxh_dtm*dcwe_dr+dyh_dtm*dcsn_dr)*conv
														ELSE
																var(4)=0.
																xmat_dzsurf(4,4)=xmat_dzsurf(2,2)+wghtsurf_ray
														ENDIF
												ENDIF

												IF(idpitch.EQ.1)THEN
														var(5)=dsurf_ray*(-dcnz_dp+dxh_dtm*dcwe_dp+dyh_dtm*dcsn_dp)*conv
												ELSE
														var(5)=0.
														xmat_dzsurf(5,5)=xmat_dzsurf(4,4)+wghtsurf_ray
												ENDIF

												IF(idhdg.EQ.1)THEN
														var(6)=dsurf_ray*(+dxh_dtm*dcwe_dh+dyh_dtm*dcsn_dh)*conv
												ELSE
														var(6)=0.
														xmat_dzsurf(6,6)=xmat_dzsurf(5,5)+wghtsurf_ray
												ENDIF

												IF(iaftfore.EQ.-1)THEN
														IF(irdaft.EQ.1)THEN
																var(7)=(-cnz+dxh_dtm*cwe+dyh_dtm*csn)*0.1
														ELSE
																var(7)=0.
																xmat_dzsurf(7,7)=xmat_dzsurf(6,6)+wghtsurf_ray
														ENDIF
														var(8)=0.
												ELSE
														var(7)=0.
														IF(irdfore.EQ.1)THEN
																var(8)=(-cnz+dxh_dtm*cwe+dyh_dtm*csn)*0.1
														ELSE
																var(8)=0.
																xmat_dzsurf(8,8)=xmat_dzsurf(8,8)+wghtsurf_ray
														ENDIF
												ENDIF

												IF(idxwe.EQ.1)THEN
														var(9)=dxh_dtm*0.1
												ELSE
														var(9)=0.
														xmat_dzsurf(9,9)=xmat_dzsurf(9,9)+wghtsurf_ray
												ENDIF

												IF(idysn.EQ.1)THEN
														var(10)=dyh_dtm*0.1
												ELSE
														var(10)=0.
														xmat_dzsurf(10,10)=xmat_dzsurf(10,10)+wghtsurf_ray
												ENDIF

												IF(idzacft.EQ.1)THEN
														var(11)=-0.1
												ELSE
														var(11)=0.
														xmat_dzsurf(11,11)=xmat_dzsurf(11,11)+wghtsurf_ray
												ENDIF

												 var(12)=0.
	

												!******************************************************************
												!**** ADD TO XMAT_dzsurf(1->NVAR,1->NVAR) AND VECT_dzsurf(1->NVAR)
												!******************************************************************
												DO i=1,nvar
														DO j=1,nvar
																xmat_dzsurf(i,j)=xmat_dzsurf(i,j)+wghtsurf_ray*var(i)*var(j)
														ENDDO
														vect_dzsurf(i)=vect_dzsurf(i)+wghtsurf_ray*var(i)*d_hsurf
												ENDDO

												!******************************************************************
												!**** ADD TO COVARIANCE MATRIX FOR FIELD "DZ_surf"
												!******************************************************************
												DO i=1,nvar
														rms_var_zsurf(i)=rms_var_zsurf(i)+wghtsurf_ray*var(i)*var(i)
														DO j=1,nvar
																 corr_var(i,j)=corr_var(i,j)+wghtsurf_ray*var(i)*var(j)
														ENDDO
												ENDDO

												!******************************************************************
												!**** CASE "DZ_surf" ONLY -> D_VH CANNOT BE CALCULATED
												!******************************************************************
												IF(rw_vsurf+rw_dvinsitu.le.0.)THEN
														xmat_vsurf(12,12)=xmat_vsurf(12,12)+wghtsurf_ray

												!******************************************************************
												!**** CASE "FLAT SURFACE" -> D_HEADING,D_XWE,D_YSN CANNOT BE OBTAINED
												!******************************************************************
														IF(altdtm_min.ge.altdtm_max)THEN
																xmat_vsurf(6,6)=xmat_vsurf(6,6)+wghtsurf_ray
																xmat_vsurf(9,9)=xmat_vsurf(9,9)+wghtsurf_ray
																xmat_vsurf(10,10)=xmat_vsurf(10,10)+wghtsurf_ray
														ENDIF
												ENDIF

												!******************************************************************
												!**** ARRAYS FOR "SIS_EL_*" FILE #50
												!******************************************************************
												zs_rot(iradar_ray,n_dzsurf(iradar_ray))=rota_ray
												zs_el(iradar_ray,n_dzsurf(iradar_ray))=elhor_ray
												zs_az(iradar_ray,n_dzsurf(iradar_ray))=azeast_ray
												zs_dsurf(iradar_ray,n_dzsurf(iradar_ray))=dsurf_ray
												zs_dhor(iradar_ray,n_dzsurf(iradar_ray))=side*dsurf_ray*celh
												zs_zsurf(iradar_ray,n_dzsurf(iradar_ray))=hsurf_ray
												zs_hsurf(iradar_ray,n_dzsurf(iradar_ray))=hsurf_dtm
										ENDIF     !!  of !! IF(kdzsurf.EQ.1)THEN

										!******************************************************************
										!**** (IF IWRISURFILE=1)
										!**** WEIGHTED SUM FOR ALT_SURF(x,y) 
										!**** TO BE WRITTEN ON "SURF_EL_*" FILE #30
										!******************************************************************
										IF(iwrisurfile.EQ.1)THEN
												IF(xsurf_ray.GT.xmin_wrisurf-hxy_wrisurf &
																		.AND.xsurf_ray.LT.xmax_wrisurf+hxy_wrisurf &
																		.AND.ysurf_ray.GT.ymin_wrisurf-hxy_wrisurf &
																		.AND.ysurf_ray.LT.ymax_wrisurf+hxy_wrisurf &
																		.AND.hsurf_ray.GT.zsurfrad_min &
																		.AND.hsurf_ray.LT.zsurfrad_max)THEN
														nsurf_wri(iradar_ray)=nsurf_wri(iradar_ray)+1
														i_wrisurf=(xsurf_ray-xmin_wrisurf)/hxy_wrisurf+1
														IF(xsurf_ray.lt.xmin_wrisurf)i_wrisurf=i_wrisurf-1
														j_wrisurf=(ysurf_ray-ymin_wrisurf)/hxy_wrisurf+1
														IF(ysurf_ray.lt.ymin_wrisurf)j_wrisurf=j_wrisurf-1
														DO ii=max0(i_wrisurf,1),min0(i_wrisurf+1,nx_wrisurf)
																xi=xmin_wrisurf+float(ii-1)*hxy_wrisurf
																dx=(xsurf_ray-xi)/hxy_wrisurf
																DO jj=max0(j_wrisurf,1),min0(j_wrisurf+1,ny_wrisurf)
																		yj=ymin_wrisurf+float(jj-1)*hxy_wrisurf
																		dy=(ysurf_ray-yj)/hxy_wrisurf
																		d2=dx*dx+dy*dy
																		wghtsurf_wri=wghtsurf_ray*((4.-d2)/(4.+d2))
																		swdzsurf_wri(ii,jj)=swdzsurf_wri(ii,jj)+wghtsurf_wri
																		SW_or_altsurf_wri(ii,jj)=SW_or_altsurf_wri(ii,jj)+wghtsurf_wri*hsurf_ray
																ENDDO
														ENDDO
												ENDIF
										ENDIF

										!******************************************************************
										!**** CASE "VDOP_surf"
										!******************************************************************
										IF(kvsurf.EQ.1.AND.acftspd_hor.GT.0.)THEN

												IF(vdop_corr(ig_refmax).GT.-900..or.isim.EQ.1)THEN
												 
														vdopsurf_ray=vdop_corr(ig_refmax)
														IF(abs(vdopsurf_ray).le.1.)nb1=nb1+1
														IF(abs(vdopsurf_ray).le.2..AND.abs(vdopsurf_ray).GT.1.)nb2=nb2+1
														IF(abs(vdopsurf_ray).le.3..AND.abs(vdopsurf_ray).GT.2.)nb3=nb3+1
														IF(abs(vdopsurf_ray).le.4..AND.abs(vdopsurf_ray).GT.3.)nb4=nb4+1
														IF(abs(vdopsurf_ray).le.5..AND.abs(vdopsurf_ray).GT.4.)nb5=nb5+1
														IF(abs(vdopsurf_ray).le.6..AND.abs(vdopsurf_ray).GT.5.)nb6=nb6+1
														IF(abs(vdopsurf_ray).le.7..AND.abs(vdopsurf_ray).GT.6.)nb7=nb7+1
														IF(abs(vdopsurf_ray).le.8..AND.abs(vdopsurf_ray).GT.7.)nb8=nb8+1
														IF(abs(vdopsurf_ray).GT.8.)nsup=nsup+1

														!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
														!  IF( nb_ray(iradar_ray).EQ.10*(nb_ray(iradar_ray)/10) )THEN
														!	IF(ig_refmax.NE.-999)THEN	! Olivier
														!
														!	v_ctrl=-999.				! Olivier
														!
														!	d_vs_vl=vs(ig_refmax)-vl(ig_refmax)	! Olivier
														!	kvsl=ifix((d_vs_vl/vnyq_el)*0.5)+5	! Olivier
														!        PRINT *,'    d_VS_VL :',d_vs_vl,' -> KVSL :',kvsl	
														! IF(kvsl.ge.1.AND.kvsl.le.9)THEN		! Olivier
														!	  vs_depl=vs(ig_refmax)+xms(kvsl)*vnyq_el	! Olivier
														!	  vl_depl=vl(ig_refmax)+xml(kvsl)*vnyq_el	! Olivier
														!	  vsl_depl=(vs_depl+vl_depl)/2.		! Olivier
														!
														!	  IF(    abs(vs_depl-vl_depl).lt.vnyq_el/2.
														!     &	    .AND.abs(VR(ig_refmax)-vsl_depl).lt.vnyq_el/2.)THEN ! Oliv
														!		PRINT *,'IG_REFMAX= ',ig_refmax
														!		PRINT *,'VR= ',VR(ig_refmax)
														!		PRINT *,'VS, VL= ', vs(ig_refmax),vl(ig_refmax)
														!		PRINT *,'VS_depl,VL_depl= ',vs_depl,vl_depl
														!		PRINT *,'VSL_depl= ',vsl_depl
														!	          v_ctrl=VR(ig_refmax)			! Olivier
														!		PRINT *,'VDOP_CTRL= ',v_ctrl
														!
														!	      IF(proj_acftspd.GT.-900.)THEN
														!		v_corr=v_ctrl+proj_acftspd
														!		PRINT *,' VDOP_CORR= ',v_corr
														!	      ENDIF
														!
														!	  ENDIF
														!
														!	ENDIF
														!
														!	ENDIF
														!  ENDIF

														!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
														!--------------------------------------------------------------
														!---- ( IF ISIM=1 ) -> SIMULATED VDOPSURF_RAY FROM dXXX_GUESS
														!--------------------------------------------------------------
														IF(isim.EQ.1)THEN
																vdopsurf_ray=-(-acftspd_we_true*dcwe_dt_true		&
																				 -acftspd_sn_true*dcsn_dt_true				&
																				 -acftspd_nz*dcnz_dt_true)						&
																				*dtilt_guess*conv									&
																			  -(-acftspd_we_true*dcwe_dr_true				&
																				 -acftspd_sn_true*dcsn_dr_true				&
																				 -acftspd_nz*dcnz_dr_true)						&
																				*drota_guess*conv									&
																			  -(-acftspd_we_true*dcwe_dp_true				&
																				 -acftspd_sn_true*dcsn_dp_true				&
																				 -acftspd_nz*dcnz_dp_true)						&
																				*dpitch_guess*conv								&
																			  -(-acftspd_we_true*dcwe_dh_true				&
																				 -acftspd_sn_true*dcsn_dh_true				&
																				 -acftspd_nz*dcnz_dh_true)						&
																				*dhdg_guess*conv									&
																			  -(-cwe_true*duacft_dv_true						&
																				 -csn_true*dvacft_dv_true)*dvh_guess
														ENDIF
														!--------------------------------------------------------------

														!IF(abs(vdopsurf_ray).lt.vdopsurf_max)THEN !!!! <--- IT IS COMMENTED IN THE ORIGINAL (RV)

														IF(iprintscreen.GT.0) THEN
																IF (nb_ray(iradar_ray).EQ.iprintscreen*(nb_ray(iradar_ray)/iprintscreen))THEN
																		PRINT *,'vdopsurf_ray =',vdopsurf_ray
																		PRINT *,'HOLA_2904'
																ENDIF
														ENDIF
																																
																!=====WRITE TO OUTPUT LOG FILE====
														IF (ioutputlog.EQ.1) THEN
																WRITE(60, '(A,F9.5)') 'line_2904  -> VDOPSURF_RAY :',vdopsurf_ray
														ENDIF
																		
														IF(abs(vdopsurf_ray).LT.6.)THEN

																!******************************************************************
																!**** ADD WEIGHTS AND VDOP_surf
																!******************************************************************
																n_vsurf(iradar_ray)=n_vsurf(iradar_ray)+1
																swvsurf_sweep(iradar_ray)=swvsurf_sweep(iradar_ray)+wghtsurf_ray
																!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
																
																IF(iprintscreen.GT.0) THEN
																		IF (nb_ray(iradar_ray).EQ.iprintscreen*(nb_ray(iradar_ray)/iprintscreen))THEN
																				PRINT *,'     -> VDOPSURF_RAY :',vdopsurf_ray
																				!!!!        PRINT *,'        SWVSURF_SWEEP(',iradar_ray,') :'
																				!!!!     &         ,swvsurf_sweep(iradar_ray)
																				PRINT *,'HOLA_2919'
																		ENDIF
																ENDIF
																																		
																		!=====WRITE TO OUTPUT LOG FILE====
																IF (ioutputlog.EQ.1) THEN
																		WRITE(60, '(A,F8.6)') 'line_2919  -> VDOPSURF_RAY :',vdopsurf_ray
																ENDIF
			
																
																!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
																vsurfsweep_mean(iradar_ray)=vsurfsweep_mean(iradar_ray) &
																															+wghtsurf_ray*vdopsurf_ray
																vsurfsweep_rms(iradar_ray)=vsurfsweep_rms(iradar_ray) &
																														+wghtsurf_ray*vdopsurf_ray*vdopsurf_ray
																swvsurf_tot=swvsurf_tot+wghtsurf_ray
																swvmsurf_tot=swvmsurf_tot+wghtsurf_ray*vdopsurf_ray
																swv2surf_tot=swv2surf_tot+wghtsurf_ray*vdopsurf_ray*vdopsurf_ray
																swavsurf_tot=swavsurf_tot +wghtsurf_ray*abs(vdopsurf_ray)
																!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
																!!!!      PRINT *,'    -> WGHTSURF_RAY:',wghtsurf_ray,' VSURF:',vdopsurf_ray
																!!!!      PRINT *,'        N_VSURF:',n_vsurf(iradar_ray)
																!!!!     &       ,' SWV,SV,SV2:',swvsurf_sweep(iradar_ray)
																!!!!     &       ,vsurfsweep_mean(iradar_ray),vsurfsweep_rms(iradar_ray)
																!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

																!******************************************************************
																!**** VALUES OF VAR(1->NVAR) FOR FIELD "V_surf"
																!****  - VAR(1->6) -> [dT_aft,dT_fore,dR_aft,dR_fore,dP,dH] in DEGREES
																!****  - VAR(7->11) -> [RD_aft,RD_fore,dXwe,dYsn,dZ] in HECTOMETERS
																!****  - VAR(12) -> [dVH] in METER/SECOND
																!******************************************************************
																IF(iaftfore.EQ.-1)THEN
																		IF(idtiltaft.EQ.1)THEN
																				var(1)=(-acftspd_we*dcwe_dt-acftspd_sn*dcsn_dt-acftspd_nz*dcnz_dt)*conv
																		ELSE
																				var(1)=0.
																				xmat_vsurf(1,1)=xmat_vsurf(1,1)+wghtsurf_ray
																		ENDIF
																		var(2)=0.
																ELSE
																		var(1)=0.
																		IF(idtiltfore.EQ.1)THEN
																				var(2)=(-acftspd_we*dcwe_dt-acftspd_sn*dcsn_dt-acftspd_nz*dcnz_dt)*conv
																		ELSE
																				var(2)=0.
																				xmat_vsurf(2,2)=xmat_vsurf(2,2)+wghtsurf_ray
																		ENDIF
																ENDIF

																IF(iaftfore.EQ.-1)THEN
																		IF(idrotaaft.EQ.1)THEN
																				var(3)=(-acftspd_we*dcwe_dr-acftspd_sn*dcsn_dr-acftspd_nz*dcnz_dr)*conv
																		ELSE
																				var(3)=0.
																				xmat_vsurf(3,3)=xmat_vsurf(3,3)+wghtsurf_ray
																		ENDIF
																		var(4)=0.
																ELSE
																		var(3)=0.
																		IF(idrotafore.EQ.1)THEN
																				var(4)=(-acftspd_we*dcwe_dr-acftspd_sn*dcsn_dr-acftspd_nz*dcnz_dr)*conv
																		ELSE
																				var(4)=0.
																				xmat_vsurf(4,4)=xmat_vsurf(4,4)+wghtsurf_ray
																		ENDIF
																ENDIF

																IF(idpitch.EQ.1)THEN
																		var(5)=(-acftspd_we*dcwe_dp-acftspd_sn*dcsn_dp-acftspd_nz*dcnz_dp)*conv
																ELSE
																		var(5)=0.
																		xmat_vsurf(5,5)=xmat_vsurf(5,5)+wghtsurf_ray
																ENDIF

																IF(idhdg.EQ.1)THEN
																		var(6)=(-acftspd_we*dcwe_dh-acftspd_sn*dcsn_dh)*conv
																ELSE
																		var(6)=0.
																		xmat_vsurf(6,6)=xmat_vsurf(6,6)+wghtsurf_ray
																ENDIF

																var(7)=0.
																var(8)=0.
																var(9)=0.
																var(10)=0.
																var(11)=0.

																IF(idvh.EQ.1)THEN
																		var(12)=-duacft_dv*cwe-dvacft_dv*csn
																ELSE
																		var(12)=0.
																		xmat_vsurf(12,12)=xmat_vsurf(12,12)+wghtsurf_ray
																ENDIF
																!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
																!!!!      PRINT *,'    VAR_VSURF(1->12):',(var(i),i=1,12)
																!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

																!******************************************************************
																!**** ADD TO XMAT_vsurf(1->NVAR,1->NVAR) AND VECT_vsurf(1->NVAR)
																!******************************************************************
																DO i=1,nvar
																		DO j=1,nvar
																				xmat_vsurf(i,j)=xmat_vsurf(i,j)+wghtsurf_ray*var(i)*var(j)
																		ENDDO
																		vect_vsurf(i)=vect_vsurf(i)+wghtsurf_ray*var(i)*vdopsurf_ray
																ENDDO

																!******************************************************************
																!**** ADD TO COVARIANCE MATRIX FOR FIELD "VDOP_surf"
																!******************************************************************
																DO i=1,nvar
																		rms_var_vsurf(i)=rms_var_vsurf(i)+wghtsurf_ray*var(i)*var(i)
																		DO j=1,nvar
																				corr_var(i,j)=corr_var(i,j)+wghtsurf_ray*var(i)*var(j)
																		ENDDO
																ENDDO

																!******************************************************************
																!**** CASE "VDOP_surf" AND/or "DVDOP_insitu" ONLY :
																!**** -> RGE-DLY_aft,RGE-DLY_aft,D_XWE,D_YSN,D_ZACFT CANNOT BE CALCULATED
																!******************************************************************
																IF(rw_dzsurf.le.0.)THEN
																		DO ij=7,11
																				 xmat_vsurf(ij,ij)=xmat_vsurf(ij,ij)+wghtsurf_ray
																		ENDDO
																ENDIF

																!******************************************************************
																!**** ARRAYS FOR "SIS_EL_*" FILE #50
																!******************************************************************
																vs_dhor(iradar_ray,n_vsurf(iradar_ray))=side*dsurf_ray*celh
																vs_vdopsurf(iradar_ray,n_vsurf(iradar_ray))=vdopsurf_ray

														ELSE  !!  of  !! IF(abs(vdopsurf_ray).lt.vdopsurf_max) !! <--  someone replaced it
																																								 !!       for IF(abs(vdopsurf_ray).lt.6.) (RV)
														
																ndismiss_vdopsurf(iradar_ray)=ndismiss_vdopsurf(iradar_ray)+1
														
														ENDIF  !! of !! IF(abs(vdopsurf_ray).lt.vdopsurf_max) !!

												ELSE  !!  of  !!  IF(vdop_corr(ig_refmax).GT.-900.) !!
												
														ndismiss_vdopcorr(iradar_ray)=ndismiss_vdopcorr(iradar_ray)+1

												ENDIF  !!  of  !!  IF(vdop_corr(ig_refmax).GT.-900.) !!

										ELSE  !!  of  !! IF(kvsurf.EQ.1.AND.acftspd_hor.GT.0.)  !!
												IF(acftspd_hor.le.0.)ndismiss_vhacft(iradar_ray)=ndismiss_vhacft(iradar_ray)+1
										ENDIF  !!  of  !! IF(kvsurf.EQ.1.AND.acftspd_hor.GT.0.)  !!
								ENDIF  !!  of  !! IF(abs(d_hsurf).lt.dhsurf_max)  !!
						ENDIF  !!  of  !!  IF(xsurf_ray.GT.xmin_dtm ... )  !!
				ENDIF  !!  of  !!  IF(hsurf_ray.GT.-999..AND.wghtsurf_ray.GT.0.)  !!
		ENDIF  !!  of  !!  IF(kdzsurf+kvsurf.ge.1 ... )  !!


      
		!******************************************************************
		!**** CASE "DVDOP_insitu"
		!**** (IF D<DMAX_insitu AND ||sin(ELEV_HOR)||<0.1)
		!******************************************************************
		IF(kdvinsitu.EQ.1.AND.ngates_insitu_max.GT.1)THEN
		
				!******************************************************************
				!**** CONTROL CONTINUITY ALONG THE RAY ( IF ICTRL_CONTRAY=1 )
				!**** DISMISS VDOP IF |VDOP-VDOP_PREV|>dVDOP_MAX AFTER UNFOLDING
				!******************************************************************
				IF(ictrl_contray.EQ.1)THEN
						init=0
						DO ig=1,ngates_insitu_max
								d_ig=dgate_corr(ig)
								IF(ZE(ig).GT.-900..AND.vdop_corr(ig).GT.-900.)THEN
										xis=0.
										svis=0.
										xrad=0.
										svrad=0.
										IF(init.EQ.0)THEN
												init=1
												xis=xpmin_contray+1.
												svis=xis*proj_wind
										ELSE
												init=2
												IF(d_ig.lt.dmax_insitu)THEN
														xis=(dmax_insitu-d_ig)/ddg
														svis=xis*proj_wind
														igmin=1
												ELSE
														igmin=((d_ig-dmax_insitu)/ddg)
												ENDIF
												DO jg=igmin,max0(1,ig-1)
														IF(abs(vdop_corr(jg)).lt.vdop_max)THEN
																xrad=xrad+1.
																svrad=svrad+vdop_corr(jg)
														ENDIF
												ENDDO
										ENDIF
										xctrl=xis+xrad
										IF(xctrl.ge.xpmin_contray)THEN
												vctrl=(svis+svrad)/xctrl
												dv=vdop_corr(ig)-vctrl
												idepl=0
												IF(ichoice_vdop.EQ.1.or.ichoice_vdop.EQ.2)THEN
														IF(abs(dv).GT.vnyq)THEN
																idepl=1
																DO WHILE (dv.GT.+vnyq)
																		vdop_corr(ig)=vdop_corr(ig)-2.*vnyq
																		dv=vdop_corr(ig)-vctrl
																ENDDO
																DO WHILE (dv.lt.-vnyq)
																		vdop_corr(ig)=vdop_corr(ig)+2.*vnyq
																		dv=vdop_corr(ig)-vctrl
																ENDDO
														ENDIF
												ENDIF
												IF(abs(dv).GT.dvdop_max)THEN
														vdop_corr(ig)=-999.
														IF(init.EQ.1)init=0
												ENDIF
										ENDIF
								ENDIF
						ENDDO
				ENDIF    !!!!  OF IF(ictrl_contray.EQ.1)

				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				itest=1
				DO ig=1,ngates_insitu_max
						IF(ZE(ig).GT.-900..AND.vdop_corr(ig).GT.-900.)itest=1
				ENDDO
				
				! print to screen
				IF(iprintscreen.GT.0) THEN
						IF (nb_ray(iradar_ray).EQ.iprintscreen*(nb_ray(iradar_ray)/iprintscreen))THEN
							PRINT *,' '
							PRINT *,' ',1000*ihhmmss+ims_ray, &
												' IRADAR:',iradar_ray, &
												' NO_RAY:',nb_ray(iradar_ray)
							PRINT *,'    ROTA,TILT_RAY:',rota_ray,tilt_ray
							PRINT *,'    ROLL,PITCH,HDG,DRIFT_ACFT:',roll_acft, &
													pitch_acft,hdg_acft,drift_acft
							PRINT *,'    AZ_EAST:',azeast_ray,' EL_HOR:',elhor_ray
							PRINT *,'HOLA_3150'
						ENDIF
				ENDIF
				!=====WRITE TO OUTPUT LOG FILE====
				IF(ioutputlog .EQ.1) THEN
						WRITE(60,600) ' line_3150'
						WRITE(60,601) 1000*ihhmmss+ims_ray, &
													' IRADAR:',iradar_ray, &
													' NO_RAY:',nb_ray(iradar_ray)
						WRITE(60,602) '    ROTA,TILT_RAY:',rota_ray,tilt_ray
						WRITE(60,603) '    ROLL,PITCH,HDG,DRIFT_ACFT:',roll_acft, &
															pitch_acft,hdg_acft,drift_acft
						WRITE(60,604) '    AZ_EAST:',azeast_ray,' EL_HOR:',elhor_ray
						600   format(A)
						601   format(I10,X,A,I10,A,I10)
						602   format(A,2F10.2)
						603   format(A,4F10.2)
						604   format(A,F10.2,A,F10.2)
				ENDIF
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				DO ig=1,ngates_insitu_max
						!	IF(     ZE(ig).GT.10.
						!     &     .AND.abs(VR(ig)).GT.0.
						!     &     .AND.abs(VR(ig)).lt.vdop_max
						!     &     .AND.abs(vs(ig)).GT.0.
						!     &     .AND.abs(vs(ig)).lt.vdop_max
						!     &     .AND.abs(vl(ig)).GT.0.
						!     &     .AND.abs(vl(ig)).lt.vdop_max
						!     &     .AND.proj_acftspd.GT.-900.   
						!     &     .AND.dgate_corr(ig).lt.10.
						!     &     .AND.elhor_ray.lt.5.
						!     &     .AND.elhor_ray.GT.-5.     )THEN	! 111111
						!
						!       d_vs_vl=vs(ig)-vl(ig)     ! Olivier
						!       kvsl=ifix((d_vs_vl/vnyq_el)*0.5)+5      ! Olivier
						!       IF(kvsl.ge.1.AND.kvsl.le.9)THEN         ! Olivier
						!         vs_depl=vs(ig)+xms(kvsl)*vnyq_el       ! Olivier
						!         vl_depl=vl(ig)+xml(kvsl)*vnyq_el       ! Olivier
						!         vsl_depl=(vs_depl+vl_depl)/2.         ! Olivier
						!
						!         IF(     abs(vs_depl-vl_depl).lt.vnyq_el/2. 
						! c    &      .AND.abs(VR(ig)-vsl_depl).lt.vnyq_el/2.     )THEN ! 112 Oliv
						!		PRINT *,'IG= ',ig
						!		PRINT *,'VR= ',VR(ig)
						!		PRINT *,'VS,VL= ',vs(ig),vl(ig)
						!		PRINT *,'VS_depl,VL_depl= ',vs_depl,vl_depl
						!		PRINT *,'VSL_depl= ',vsl_depl!!
						!
						!	    IF(proj_acftspd.GT.-900.)THEN
						!	      v_corr=VR(ig)+proj_acftspd                ! Olivier
						!	      vdop_corr(ig)=v_corr
						!              PRINT *,'    -> VDOP_CORR :',v_corr       ! Olivier
						!	    ENDIF
						!
						!         ENDIF
						!
						!       ENDIF
						!	ENDIF                               ! Olivier

						!------------------------------------------------------------------
						!---- ( IF ISIM=1 ) -> SIMULATED dV_dopinsitu WITH dXXX_GUESS
						!------------------------------------------------------------------
						IF(ig.EQ.1.AND.isim.EQ.1)THEN
								ZE(1)=999.
								dv_dopinsitu=-( wa_we_true*dcwe_dt_true						&
															+wa_sn_true*dcsn_dt_true					&
															+wa_nz*dcnz_dt_true)						&
															*dtilt_guess*conv								&
															-( wa_we_true*dcwe_dr_true				&
															+wa_sn_true*dcsn_dr_true					&
															+wa_nz*dcnz_dr_true)						&
															*drota_guess*conv								&
															-( wa_we_true*dcwe_dp_true				&
															+wa_sn_true*dcsn_dp_true					&
															+wa_nz*dcnz_dp_true)						&
															*dpitch_guess*conv							&
															-( wa_we_true*dcwe_dh_true				&
															+wa_sn_true*dcsn_dh_true					&
															+wa_nz*dcnz_dh_true)						&
															*dhdg_guess*conv								&
															-(-cwe_true*duacft_dv_true					&
															-csn_true*dvacft_dv_true)*dvh_guess
								vdop_corr(1)=dv_dopinsitu+proj_wind_true
								DO iig=2,ngates_insitu_max
										ZE(iig)=-999.
										vdop_corr(iig)=-999.
								ENDDO
						ENDIF
				
						!------------------------------------------------------------------
						IF(ZE(ig).GT.-900. .AND.vdop_corr(ig).GT.-900.)THEN

								wghtinsitu_ig=1.-0.5*dgate_corr(ig)/dmax_insitu
								dv_dopinsitu=vdop_corr(ig)-proj_wind

								IF(abs(dv_dopinsitu).LT.dvdopinsitu_max)THEN

										!******************************************************************
										!**** ADD WEIGHTS AND DVDOP_insitu
										!******************************************************************
										n_dvinsitu(iradar_ray)=n_dvinsitu(iradar_ray)+1
										ssurfins=ssurfins+wghtinsitu_ig
										swinsitu_sweep(iradar_ray)=swinsitu_sweep(iradar_ray)+wghtinsitu_ig
										dvinsitusweep_mean(iradar_ray)=dvinsitusweep_mean(iradar_ray) &
																										+wghtinsitu_ig*dv_dopinsitu
										dvinsitusweep_rms(iradar_ray)=dvinsitusweep_rms(iradar_ray) &
																									+wghtinsitu_ig*dv_dopinsitu*dv_dopinsitu
										swdvinsitu_tot=swdvinsitu_tot+wghtinsitu_ig
										swdvminsitu_tot=swdvminsitu_tot+wghtinsitu_ig*dv_dopinsitu
										swdv2insitu_tot=swdv2insitu_tot+wghtinsitu_ig*dv_dopinsitu*dv_dopinsitu   
										swadvinsitu_tot=swadvinsitu_tot+wghtinsitu_ig*abs(dv_dopinsitu)

										s_vpv(iradar_ray,ilr)=s_vpv(iradar_ray,ilr)+wghtinsitu_ig
										sv_vpv(iradar_ray,ilr)=sv_vpv(iradar_ray,ilr)+wghtinsitu_ig*dv_dopinsitu
										svv_vpv(iradar_ray,ilr)=svv_vpv(iradar_ray,ilr)+wghtinsitu_ig*dv_dopinsitu*dv_dopinsitu
										x_vpv(iradar_ray,ilr)=x_vpv(iradar_ray,ilr)+wghtinsitu_ig
										xv_vpv(iradar_ray,ilr)=xv_vpv(iradar_ray,ilr)+wghtinsitu_ig*dv_dopinsitu
										xvv_vpv(iradar_ray,ilr)=xvv_vpv(iradar_ray,ilr)+wghtinsitu_ig*dv_dopinsitu*dv_dopinsitu

										!******************************************************************
										!**** VALUES OF VAR(1->NVAR) FOR FIELD "DV_insitu"
										!****  - VAR(1->6) -> [dT_aft,dT_fore,dR_aft,dR_fore,dP,dH] in DEGREES
										!****  - VAR(7->11) -> [RD_aft,RD_fore,dXwe,dYsn,dZ] in HECTOMETERS
										!****  - VAR(12) -> [dVH] in METER/SECOND
										!******************************************************************
										IF(iaftfore.EQ.-1)THEN
												IF(idtiltaft.EQ.1)THEN
														var(1)=( wa_we*dcwe_dt+wa_sn*dcsn_dt+wa_nz*dcnz_dt)*conv	
												ELSE
														var(1)=0.
														xmat_dvinsitu(1,1)=xmat_dvinsitu(1,1)+wghtinsitu_ig
												ENDIF
												var(2)=0.
										ELSE
												var(1)=0.
												IF(idtiltfore.EQ.1)THEN
														var(2)=( wa_we*dcwe_dt+wa_sn*dcsn_dt+wa_nz*dcnz_dt)*conv
												ELSE
														var(2)=0.
														xmat_dvinsitu(2,2)=xmat_dvinsitu(2,2)+wghtinsitu_ig
												ENDIF
										ENDIF

										IF(iaftfore.EQ.-1)THEN
												IF(idrotaaft.EQ.1)THEN
														var(3)=( wa_we*dcwe_dr+wa_sn*dcsn_dr+wa_nz*dcnz_dr)*conv
												ELSE
														var(3)=0.
														xmat_dvinsitu(3,3)=xmat_dvinsitu(3,3)+wghtinsitu_ig
												ENDIF
												var(4)=0.
										ELSE
												var(3)=0.
												IF(idrotafore.EQ.1)THEN
														var(4)=( wa_we*dcwe_dr+wa_sn*dcsn_dr+wa_nz*dcnz_dr)*conv
												ELSE
														var(4)=0.
														xmat_dvinsitu(4,4)=xmat_dvinsitu(4,4)+wghtinsitu_ig
												ENDIF
										ENDIF

										IF(idpitch.EQ.1)THEN
												var(5)=( wa_we*dcwe_dp+wa_sn*dcsn_dp+wa_nz*dcnz_dp)*conv
										ELSE
												var(5)=0.
												xmat_dvinsitu(5,5)=xmat_dvinsitu(5,5)+wghtinsitu_ig
										ENDIF

										IF(idhdg.EQ.1)THEN
												var(6)=(wa_we*dcwe_dh+wa_sn*dcsn_dh)*conv
										ELSE
												var(6)=0.
												xmat_dvinsitu(6,6)=xmat_dvinsitu(6,6)+wghtinsitu_ig
										ENDIF

										var(7)=0.
										var(8)=0.
										var(9)=0.
										var(10)=0.
										var(11)=0.

										IF(idvh.EQ.1)THEN
												var(12)=-duacft_dv*cwe-dvacft_dv*csn
										ELSE
												var(12)=0.
												xmat_dvinsitu(12,12)=xmat_dvinsitu(12,12)+wghtinsitu_ig
										ENDIF
								

										!******************************************************************
										!**** ADD TO XMAT_dvinsitu(1->NVAR,1->NVAR) AND VECT_dvinsitu(1->NVAR)
										!******************************************************************
										DO i=1,nvar
												DO j=1,nvar
														xmat_dvinsitu(i,j)=xmat_dvinsitu(i,j)+wghtinsitu_ig*var(i)*var(j)
												ENDDO
												vect_dvinsitu(i)=vect_dvinsitu(i)+wghtinsitu_ig*var(i)*dv_dopinsitu
										ENDDO

										!******************************************************************
										!**** ADD TO COVARIANCE MATRIX FOR FIELD "DVDOP_insitu"
										!******************************************************************
										DO i=1,nvar
												rms_var_vinsitu(i)=rms_var_vinsitu(i)+wghtinsitu_ig*var(i)*var(i)
												DO j=1,nvar
														corr_var(i,j)=corr_var(i,j)+wghtinsitu_ig*var(i)*var(j)
												ENDDO
										ENDDO

										!******************************************************************
										!**** CASE "VDOP_surf" AND/or "DVDOP_insitu" ONLY :
										!**** -> RGE-DLY_aft,RGE-DLY_aft,D_XWE,D_YSN,D_ZACFT CANNOT BE CALCULATED
										!******************************************************************
										IF(rw_dzsurf.le.0.)THEN
												DO ij=7,11
														xmat_vsurf(ij,ij)=xmat_vsurf(ij,ij)+wghtinsitu_ig
												ENDDO
										ENDIF

										!******************************************************************
										!**** ARRAYS FOR "SIS_EL_*" FILE
										!******************************************************************
										vi_dhor(iradar_ray,n_dvinsitu(iradar_ray))=side*dgate_corr(ig)*celh
										vi_vdop(iradar_ray,n_dvinsitu(iradar_ray))=vdop
										vi_vinsitu(iradar_ray,n_dvinsitu(iradar_ray))=proj_wind

								ENDIF  !!  of  !!  IF(abs(dv_dopinsitu).lt.dvdopinsitu_max)  !!
						ENDIF  !!  of  !!  IF(ZE(ig).GT.-900. ... )  !!
				ENDDO !!  of  !!  DO ig=1,ngates_insitu_max  !!
		ENDIF  !!  of  !!  IF(kdvinsitu.EQ.1.AND.ngates_insitu_max.GT.1)  !!

		!******************************************************************
		!**** STORE FOR NEXT RAY
		!******************************************************************
		istart_sweep(iradar_ray)=1
		swp_prev(iradar_ray)=swp(iradar_ray)
		vnyq_prev=vnyq
		rota_prev(iradar_ray)=rota_ray
		tilt_prev=tilt_ray

		GOTO 1 ! if iopen=1 read next CfRadial line

3		STOP
	
		END ! finish program
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!******************************************************************
!
!**** CALCULATE DIRECTOR COSINES, AZIM_EST ET DE ELEV_HOR
!**** FROM LEE ET AL. (JTech, 1994, 11, 572-578) 
!
!******************************************************************
      SUBROUTINE azel(rotaroll,tilt_ray, hdg_acft,drift_acft,pitch_acft, &
					                     azeast_ray,elhor_ray, cxa,cya,cza,cwe,csn,cnz)
					                     
      COMMON	/cosinang/		crr,srr,cti,sti, &
															chdg,shdg,cdri,sdri,cpit,spit, &
															caze,saze,celh,selh

      conv=3.14159/180.

      crr=cos(conv*rotaroll)
      srr=sin(conv*rotaroll)
      cti=cos(conv*tilt_ray)
      sti=sin(conv*tilt_ray)
      IF(srr.GT.0.)THEN
        side=+1.
      ELSE
        side=-1.
      ENDIF

      chdg=cos(conv*hdg_acft)
      shdg=sin(conv*hdg_acft)
      cdri=cos(conv*drift_acft)
      sdri=sin(conv*drift_acft)
      cpit=cos(conv*pitch_acft)
      spit=sin(conv*pitch_acft)

      !cxa=+srr*cti
      !cya=-cti*crr*spit+sti*cpit
      !cza=+cti*crr*cpit+sti*spit
      
      !cwe=+chdg*cxa+shdg*cya
      !csn=-shdg*cxa+chdg*cya
      !cnz=+cza
      
		! Earth-relative coordinates (Eq. 15)
      x=-(crr*shdg*cti*spit)+(chdg*srr*cti)+(shdg*cpit*sti) 
      y=-(crr*chdg*cti*spit)-(shdg*srr*cti)+(cpit*chdg*sti)	
      z= cpit*cti*crr+spit*sti															
      
      cwe=x
      csn=y
      cnz=z
      
      !azeast_ray=atan2(csn,cwe)/conv
      azeast_ray=atan2(x,y)/conv			! azimuth angle (earth-relative, Eq 16)
      caze=cos(conv*azeast_ray)
      saze=sin(conv*azeast_ray)
      DO WHILE (azeast_ray.LE.0.)
         azeast_ray=azeast_ray+360.
      ENDDO
      
      !chor=amax1(0.1,sqrt(cwe*cwe+csn*csn))
      chor=amax1(0.1,sqrt(x*x+y*y))	!	range
      elhor_ray=atan(z/chor)/conv			!	elevation angle over horizontal radar plane (Eq 17)
      																!	original Eq 17 uses asin but isn't correct
      celh=cos(conv*elhor_ray)
      selh=sin(conv*elhor_ray)

      RETURN
      END

!****************************************************************************
!
!**** INVERSION OF MATRIX (NVAR,NVAR)
!
!****************************************************************************
      SUBROUTINE resoud(xmat,xinv,vect,res,nvar)

      dimension xmat(nvar,nvar),vect(nvar),res(nvar)

      call chol_inv(xmat,res,vect,nvar,ierr)

      return
      end

!****************************************************************************
!
!**** INTERPOLATION TO FILL THE HOLES IN THE RADAR-DERIVED SURFACE MAP
!
!****************************************************************************
      SUBROUTINE inter(sp,sz,nx,ny,nxysurfmax)
      dimension sp(nxysurfmax,nxysurfmax),sz(nxysurfmax,nxysurfmax), &
     							x(1000),y(1000),s(1000),d(1000)
      nsaut=5
      nmin=5
      spmin=1.
      nin=0
      nintx=0
      ninty=0
      nout=0

      DO j=1,ny
         DO i=1,nx
            IF(sp(i,j).GT.spmin)THEN
              sz(i,j)=sz(i,j)/sp(i,j)
              nin=nin+1
	    ELSE       
	      sz(i,j)=-999.
            ENDIF
         ENDDO
      ENDDO

      PRINT *,'     -> ALONG X'
      DO j=1,ny
         imax=1
  1      imin=imax
	 iant=0
	 n=0
	 DO i=imin,nx
	    imax=i
	    IF(sz(i,j).GT.-900.)THEN
	      IF(iant.NE.0.AND.(i-iant).GT.nsaut+1)go to 2
	      iant=i
	      n=n+1
	      x(n)=float(i)
	      y(n)=sz(i,j)
            ENDIF
         ENDDO 
  2      IF(n.ge.nmin)THEN
	   q1=(y(2)-y(1))/(x(2)-x(1))
	   qn=(y(n)-y(n-1))/(x(n)-x(n-1))
	   call spline(x,y,s,d,q1,qn,n)
	   DO i=imin,imax
              IF(sz(i,j).lt.-900.)THEN
	        xi=float(i)
	        val=splin(xi,x,y,s,d,q1,qn,n)
                IF(val.GT.0.)THEN
                  sz(i,j)=val
                  nintx=nintx+1
                ENDIF
              ENDIF
           ENDDO
         ENDIF
	 IF(imax.le.(nx-nsaut+1))go to 1
      ENDDO

      PRINT *,'     -> ALONG Y'
      DO i=1,nx
         jmax=1
  3      jmin=jmax
         jant=0
 	 n=0
	 DO j=jmin,ny
	    jmax=j
	    IF(sz(i,j).GT.-900.)THEN
	      IF(jant.NE.0.AND.(j-jant).GT.nsaut+1)go to 4
	      jant=j
	      n=n+1
	      x(n)=float(j)
	      y(n)=sz(i,j)
            ENDIF
         ENDDO 
  4      IF(n.ge.nmin)THEN
	   q1=(y(2)-y(1))/(x(2)-x(1))
	   qn=(y(n)-y(n-1))/(x(n)-x(n-1))
	   call spline(x,y,s,d,q1,qn,n)
	   DO j=jmin,jmax
              IF(sz(i,j).lt.-900.)THEN
	        yj=float(j)
	        val=splin(yj,x,y,s,d,q1,qn,n)
                IF(val.GT.0.)THEN
                  sz(i,j)=val
                  ninty=ninty+1
                  nput=nout+1
                ENDIF
              ELSE
                nout=nout+1
              ENDIF
           ENDDO
         ENDIF
	 IF(jmax.le.(ny-nmin+1))go to 3
      ENDDO

      PRINT *,'     -> N_in,int_X,int_Y,out :',nin,nintx,ninty,nout
      
      RETURN
      END

!****************************************************************************
!
!**** SUBROUTINE SPLINE
!
!****************************************************************************
      SUBROUTINE spline(x,u,s,del,q1,qn,n)
      dimension x(1000),u(1000),s(1000),del(1000)
      dimension a(1000),v(1000)

      del(2)=x(2)-x(1)
      v(1)=6.*(((u(2)-u(1))/del(2))-q1)
      n1=n-1
      DO i=2,n1
         del(i+1)=x(i+1)-x(i)
         v(i)=((u(i-1)/del(i))-u(i)*((1./del(i))+(1./del(i+1)))+(u(i+1)/del(i+1)))*6.
      ENDDO
      v(n)=(qn+(u(n1)-u(n))/del(n))*6.

      a(1)=2.*del(2)
      a(2)=1.5*del(2)+2.*del(3)
      v(2)=v(2)-.5*v(1)
      DO i=3,n1
         c=del(i)/a(i-1)
         a(i)=2.*(del(i)+del(i+1))-c*del(i)
         v(i)=v(i)-c*v(i-1)
      ENDDO
      c=del(n)/a(n1)
      a(n)=2.*del(n)-c*del(n)
      v(n)=v(n)-c*v(n1)
 
      s(n)=v(n)/a(n)
      DO j=1,n1
         i=n-j
         s(i)=(v(i)-del(i+1)*s(i+1))/a(i)
      ENDDO
 
       RETURN
      END

!****************************************************************************
!
!**** FUNCTION SPLIN
!
!****************************************************************************

      FUNCTION splin(v,x,u,s,del,q1,qn,n)
      dimension x(1000),u(1000),s(1000),del(1000)

      IF(v-x(1))50,10,20
  10  splin=u(1)
      return
  20  DO k=2,n
         IF(v-x(k))30,30,40
  30     k1=k-1
         ff1=s(k1)*(x(k)-v)**3.
         ff2=s(k)*(v-x(k1))**3.
         ff3=1./(6.*del(k))
         f1=(ff1+ff2)*ff3
         f2=(v-x(k1))*(u(k)/del(k)-s(k)*del(k)/6.)
         f3=(x(k)-v)*(u(k1)/del(k)-s(k1)*del(k)/6.)
         splin=f1+f2+f3
         return
  40     continue
      ENDDO
  50  splin=0.

      RETURN
      END
