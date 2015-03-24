
	PROGRAM readnetcdf 
	!
	!   Program to read netcdf aircraft sweepfiles.
	!   Huaqing Cai April 2010
	!   Niles Oien  March 2010
	!	 Raul Valenzuela October 2014

	! //////////////////////////////////////////////////////////////////////////////////////////
	!   This program creates a list with netcdf (CfRadial) files
	!   contained in the current folder and output a text file.
	!
	!   The executable needs to be located in the same CfRadial Folder.
	!
	!   Syntax:
	!	 within the folder containing CfRadial files run
	!   $ ./netcdf2text DZ VG
	!
	!   RV 2014
	!//////////////////////////////////////////////////////////////////////////////////////////

	IMPLICIT none 
	INCLUDE '/usr/include/netcdf.inc'

	!  Variable Declaration 
	INTEGER NBEAMS    ! Number of times rays were collected --- Number of rays
	INTEGER nranges   ! Number of range gates per time  --- Number of gates for each ray
	INTEGER nsweeps   ! Number of sweep in each netcdf file; usually it is one

	INTEGER*4 counter,T_RAYS  ! Ray number and total rays after each sweep
	INTEGER RETVAL          ! Return value for netcdf calls.
	INTEGER nid             ! NETCDF ID
	INTEGER i,j,k,ii,jj,kk  ! General purpose
	INTEGER NUM_MISSING     ! Number of gates that had missing data

	LOGICAL FIRST     ! Used to get min, max

	!  Scaler variable for coccrection factors

	REAL azimuth_correction   
	REAL elevation_correction
	REAL range_correction
	REAL longitude_correction
	REAL latitude_correction
	REAL pressure_altitude_correction
	! Change radar_altitude_correction to altitude_correction
	!   REAL radar_altitude_correction
	REAL altitude_correction
	! Change ew_gound_speed_correction to eastward_velocity_correction
	!   REAL ew_gound_speed_correction
	REAL eastward_velocity_correction
	! Change ns_ground_speed_correction to northward_velocity_correction
	!   REAL ns_ground_speed_correction
	REAL northward_velocity_correction
	REAL vertical_velocity_correction
	REAL heading_correction
	REAL roll_correction
	REAL pitch_correction
	REAL drift_correction
	REAL rotation_correction
	REAL tilt_correction

	! Miscellaneous variables

	REAL MAX,MIN      ! max, min values

	REAL BAD_DATA_VALUE
	PARAMETER(BAD_DATA_VALUE=-999.0) ! Value we use for bad data

	! One dimensional dynamical array of azimuths, one per time

	INTEGER,ALLOCATABLE,DIMENSION(:) :: sweep_number

	REAL*8,ALLOCATABLE,DIMENSION(:) :: time   
	REAL,ALLOCATABLE,DIMENSION(:) :: range 
	REAL,ALLOCATABLE,DIMENSION(:) :: azimuth      
	REAL,ALLOCATABLE,DIMENSION(:) :: elevation   
	REAL*8,ALLOCATABLE,DIMENSION(:) :: latitude
	REAL*8,ALLOCATABLE,DIMENSION(:) :: longitude
	REAL*8,ALLOCATABLE,DIMENSION(:) :: altitude
	REAL,ALLOCATABLE,DIMENSION(:) :: altitude_agl
	REAL,ALLOCATABLE,DIMENSION(:) :: heading 
	REAL,ALLOCATABLE,DIMENSION(:) :: roll
	REAL,ALLOCATABLE,DIMENSION(:) :: pitch 
	REAL,ALLOCATABLE,DIMENSION(:) :: drift
	REAL,ALLOCATABLE,DIMENSION(:) :: rotation 
	REAL,ALLOCATABLE,DIMENSION(:) :: tilt 
	! Change ew_velocity to eastward_velocity
	! Change ns_velocity to northward_velocity
	!   REAL,ALLOCATABLE,DIMENSION(:) :: ew_velocity 
	!   REAL,ALLOCATABLE,DIMENSION(:) :: ns_velocity
	REAL,ALLOCATABLE,DIMENSION(:) :: eastward_velocity 
	REAL,ALLOCATABLE,DIMENSION(:) :: northward_velocity 

	REAL,ALLOCATABLE,DIMENSION(:) :: vertical_velocity 
	! Change ew_wind to eastward_wind, ns_wind to northward_wind
	!   REAL,ALLOCATABLE,DIMENSION(:) :: ew_wind
	!   REAL,ALLOCATABLE,DIMENSION(:) :: ns_wind 
	REAL,ALLOCATABLE,DIMENSION(:) :: eastward_wind 
	REAL,ALLOCATABLE,DIMENSION(:) :: northward_wind 

	REAL,ALLOCATABLE,DIMENSION(:) :: vertical_wind 

	! Two dimensional dynamical array of DBZ, VR, SW, NCP, etc 

	REAL,ALLOCATABLE,DIMENSION(:,:) :: FIELD1 
	REAL,ALLOCATABLE,DIMENSION(:,:) :: FIELD2
	REAL,ALLOCATABLE,DIMENSION(:,:) :: FIELDNAN
	! REAL,ALLOCATABLE,DIMENSION(:,:) :: NCP 
	! REAL,ALLOCATABLE,DIMENSION(:,:) :: SW 

	! Two dimensional dynamical array of 2 byte ints as stored in file
	INTEGER*2,ALLOCATABLE,DIMENSION(:,:) :: I2DATA    
	REAL SCALE  ! Scale used to convert ints to reals
	REAL OFFSET ! Offset used to convert ints to reals
	INTEGER*2 MISSING ! Value used to indicate missing data in the netcdf file

	! Variables for start and end time for one netcdf file

	! Change start_time to time_coverage_start
	! Change end_time to time_coverage_end

	!   CHARACTER(len=64) start_time,end_time
	CHARACTER(len=64) time_coverage_start,time_coverage_end 

	! Variables for input file list
	CHARACTER(len=80) filename,outfilename , argu
	INTEGER  nfile     ! number of netcdf files

	! field names
	CHARACTER(len=2) FIELD1_name,FIELD2_name

	! system command
	CHARACTER(len=40) command

	!   Initialization
	counter = 0
	T_RAYS = 0

	! get field names from command line
	CALL GETARG(1, FIELD1_name)
	CALL GETARG(2, FIELD2_name)

	! create a temporal text file with the number 
	! of CfRadial files within the working folder (RV)
	command='ls -l *.nc | wc -l > nfiles.txt'
	CALL system(command)
	OPEN(10, file='nfiles.txt', status='old')
	READ(10,*) nfile
	CLOSE(10)
	command='rm nfiles.txt'
	CALL system(command)
	
	! create a temporal text file with a list
	! of CfRadial files within the folder (RV)
	command='ls *.nc > filenames.txt'
	CALL system(command)
	OPEN(20, file='filenames.txt', status='old')

	DO ii=1,nfile
			READ(20,*) filename 
			PRINT *, 'Total Files =',nfile,'file=',ii,'  ',filename

			! Open the output file for storing ray information
			WRITE(outfilename,'(I0.3)') ii ! file name of form 001, 002, etc (RV)
			outfilename =  TRIM(ADJUSTL(outfilename)) // '.txt'

			OPEN(30, file=outfilename, status='new')
			!
			! Open the netcdf file.
			!
			! In the following functions, RETVAL contains the status of
			! the CfRadial file. The variable itself is saved in the last argument
			! of the function (RV)
			!
			RETVAL = NF_OPEN(TRIM(filename), NF_NOWRITE, nid)
			CALL CHECK_NETCDF(RETVAL, 'Failed to open input file')
			!
			! Read the number of sweeps, times and ranges. These are
			! dimensions rather than variables and so are stored differently.
			!
			RETVAL = NF_INQ_DIMID(nid, 'sweep', I)
			CALL CHECK_NETCDF(RETVAL, 'Failed to find a dimension named sweep')
			RETVAL = NF_INQ_DIMLEN(NID, I, NSWEEPS)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read number of sweeps')

			RETVAL = NF_INQ_DIMID(nid, 'time', I)
			CALL CHECK_NETCDF(RETVAL, 'Failed to find a dimension named time')
			RETVAL = NF_INQ_DIMLEN(NID, I, NBEAMS)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read number of times')

			RETVAL = NF_INQ_DIMID(NID, 'range', I)
			CALL CHECK_NETCDF(RETVAL, 'Failed to find a dimension named range')
			RETVAL = NF_INQ_DIMLEN(NID, I, NRANGES)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read number of ranges')

			PRINT *, NSWEEPS,' sweeps ', NBEAMS, ' beams found, each with ', NRANGES, ' gates.'
			!
			! Read scalar variables - all the correction factors.
			!
			RETVAL = NF_INQ_VARID(NID, 'azimuth_correction', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named elevation_correction found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, azimuth_correction)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read azimuth correction.')

			RETVAL = NF_INQ_VARID(NID, 'elevation_correction', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named elevation_correction found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, elevation_correction)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read elevation correction.')

			RETVAL = NF_INQ_VARID(NID, 'range_correction', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named range_correction found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, range_correction)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read range correction.')

			! Change range_correction unit from meters to km:
			range_correction = range_correction/1000.0

			RETVAL = NF_INQ_VARID(NID, 'longitude_correction', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named longitude_correction found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, longitude_correction)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read longitude correction.')

			RETVAL = NF_INQ_VARID(NID, 'latitude_correction', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named latitude_correction found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, latitude_correction)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read  latitude correction.')

			RETVAL = NF_INQ_VARID(NID, 'pressure_altitude_correction', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named pressure_altitude_correction found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, pressure_altitude_correction)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read pressure altitude correction.')

			! change  pressure_altitude_correction unit from meters to km:
			pressure_altitude_correction = pressure_altitude_correction/1000.0

			RETVAL = NF_INQ_VARID(NID, 'altitude_correction', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named altitude_correction found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, altitude_correction)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read altitude correction.')

			! change  altitude_correction unit from meters to km:
			altitude_correction = altitude_correction/1000.0

			RETVAL = NF_INQ_VARID(NID, 'eastward_velocity_correction', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named eastward_velocity_correction found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, eastward_velocity_correction)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read eastward_velocity_correction.')

			RETVAL = NF_INQ_VARID(NID, 'northward_velocity_correction', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named northward_velocity_correction found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, northward_velocity_correction)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read northward_velocity_correction.')

			RETVAL = NF_INQ_VARID(NID, 'vertical_velocity_correction', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named vertical_velocity_correction found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, vertical_velocity_correction)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read vertical_velocity correction.')

			RETVAL = NF_INQ_VARID(NID, 'heading_correction', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named heading_correction found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, heading_correction)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read heading correction.')

			RETVAL = NF_INQ_VARID(NID, 'roll_correction', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named roll_correction found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, roll_correction)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read roll correction.')

			RETVAL = NF_INQ_VARID(NID, 'pitch_correction', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named pitch_correction found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, pitch_correction)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read pitch correction.')

			RETVAL = NF_INQ_VARID(NID, 'drift_correction', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named drift_correction found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, drift_correction)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read drift correction.')

			RETVAL = NF_INQ_VARID(NID, 'rotation_correction', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named rotation_correction found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, rotation_correction)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read rotation correction.')

			RETVAL = NF_INQ_VARID(NID, 'tilt_correction', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named tilt_correction found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, tilt_correction)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read tilt correction.')

			!  Get the start_time  of the sweep
			!
			RETVAL = NF_INQ_VARID(NID, 'time_coverage_start', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named time_coverage_start found.')
			RETVAL = NF_GET_VAR_TEXT(NID, I, time_coverage_start)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read time_coverage_start.')

			PRINT *,'start_time of the sweep: ', TRIM(time_coverage_start(1:20)) 

			PRINT *, 'Correction Factors : '
			PRINT *, 'azimuth_correction =',azimuth_correction
			PRINT *, 'elevation_correction =',elevation_correction
			PRINT *, 'range_correction =',range_correction
			PRINT *, 'longitude_correction =',longitude_correction
			PRINT *, 'latitude_correction =',latitude_correction
			PRINT *, 'pressure_altitude_correction =',pressure_altitude_correction
			PRINT *, 'altitude_correction =',altitude_correction
			PRINT *, 'eastward_velocity_correction =',eastward_velocity_correction
			PRINT *, 'northward_velocity_correction =',northward_velocity_correction
			PRINT *, 'vertical_velocity_correction =',vertical_velocity_correction
			PRINT *, 'heading_correction =',heading_correction
			PRINT *, 'roll_correction =',roll_correction
			PRINT *, 'pitch_correction =',pitch_correction
			PRINT *, 'drift_correction =',drift_correction
			PRINT *, 'rotation_correction =',rotation_correction
			PRINT *, 'tilt_correction =',tilt_correction


			!
			! Read one dimensional array - the time, range, azimuths, elevation, etc
			!
			!  ========== Sweep ============
			ALLOCATE(sweep_number(nsweeps)) ! Dynamic memory allocation
			RETVAL = NF_INQ_VARID(NID, 'sweep_number', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named sweep_number found.')
			RETVAL = NF_GET_VAR_INT(NID, I, sweep_number)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read sweep_number.')

			!  ========== Time ============
			ALLOCATE(time(NBEAMS)) ! Dynamic memory allocation
			RETVAL = NF_INQ_VARID(NID, 'time', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named time found.')
			RETVAL = NF_GET_VAR_DOUBLE(NID, I, time)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read time.')

			!  ========== range ============
			ALLOCATE(range(nranges)) ! Dynamic memory allocation
			RETVAL = NF_INQ_VARID(NID, 'range', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named range found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, range)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read range.')


			! Change range from meters to km, since in the netcdf file, range is in 
			! meters, not km anymore
			DO I=1, NRANGES
					range(I) = range(I)/1000.0
			END DO

			!  ========== azimuth ============
			ALLOCATE(azimuth(NBEAMS)) ! Dynamic memory allocation
			RETVAL = NF_INQ_VARID(NID, 'azimuth', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named azimuth found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, azimuth)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read azimuth.')
			
			!  ========== elevation ============
			ALLOCATE(elevation(NBEAMS)) ! Dynamic memory allocation
			RETVAL = NF_INQ_VARID(NID, 'elevation', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named elevation found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, elevation)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read elevation.')

			!  ========== latitude ============
			ALLOCATE(latitude(NBEAMS)) ! Dynamic memory allocation
			RETVAL = NF_INQ_VARID(NID, 'latitude', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named latitude found.')
			RETVAL = NF_GET_VAR_DOUBLE(NID, I, latitude)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read latitude.')

			!  ========== longitude ============
			ALLOCATE(longitude(NBEAMS)) ! Dynamic memory allocation
			RETVAL = NF_INQ_VARID(NID, 'longitude', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named longitude found.')
			RETVAL = NF_GET_VAR_DOUBLE(NID, I, longitude)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read longitude.')

			!  ========== altitude ============
			ALLOCATE(altitude(NBEAMS)) ! Dynamic memory allocation
			RETVAL = NF_INQ_VARID(NID, 'altitude', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named altitude found.')
			RETVAL = NF_GET_VAR_DOUBLE(NID, I, altitude)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read altitude.')

			! change unit of altitude from meter to km:
			DO I=1, NBEAMS
					altitude(I) = altitude(I)/1000.0
			END DO

			!  ========== altitude_agl ============
			ALLOCATE(altitude_agl(NBEAMS)) ! Dynamic memory allocation
			RETVAL = NF_INQ_VARID(NID, 'altitude_agl', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named altitude_agl found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, altitude_agl)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read altitude_agl.')


			! change unit of altitude_agl from meter to km:
			DO I=1, NBEAMS
					altitude_agl(I) = altitude_agl(I)/1000.0
			END DO


			!  ========== heading ============
			ALLOCATE(heading(NBEAMS)) ! Dynamic memory allocation
			RETVAL = NF_INQ_VARID(NID, 'heading', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named heading found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, heading)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read heading.')

			!  ========== roll ============
			ALLOCATE(roll(NBEAMS)) ! Dynamic memory allocation
			RETVAL = NF_INQ_VARID(NID, 'roll', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named roll found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, roll)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read roll.')

			!  ========== pitch ============
			ALLOCATE(pitch(NBEAMS)) ! Dynamic memory allocation
			RETVAL = NF_INQ_VARID(NID, 'pitch', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named pitch found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, pitch)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read pitch.')

			!  ========== drift ============
			ALLOCATE(drift(NBEAMS)) ! Dynamic memory allocation
			RETVAL = NF_INQ_VARID(NID, 'drift', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named drift found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, drift)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read drift.')

			!  ========== rotation ============
			ALLOCATE(rotation(NBEAMS)) ! Dynamic memory allocation
			RETVAL = NF_INQ_VARID(NID, 'rotation', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named rotation found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, rotation)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read rotation.')

			!  ========== tilt ============
			ALLOCATE(tilt(NBEAMS)) ! Dynamic memory allocation
			RETVAL = NF_INQ_VARID(NID, 'tilt', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named tilt found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, tilt)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read tilt.')

			!  ========== eastward_velocity ============
			ALLOCATE(eastward_velocity(NBEAMS)) ! Dynamic memory allocation
			RETVAL = NF_INQ_VARID(NID, 'eastward_velocity', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named eastward_velocity found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, eastward_velocity)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read eastward_velocity.')

			!  ========== northward_velocity ============
			ALLOCATE(northward_velocity(NBEAMS)) ! Dynamic memory allocation
			RETVAL = NF_INQ_VARID(NID, 'northward_velocity', I)
			CALL CHECK_NETCDF(RETVAL,'No variable named northward_velocity found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, northward_velocity)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read northward_velocity.')
			
			!  ========== vertical_velocity ============
			ALLOCATE(vertical_velocity(NBEAMS)) ! Dynamic memory allocation
			RETVAL = NF_INQ_VARID(NID, 'vertical_velocity', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named vertical_velocity found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, vertical_velocity)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read vertical_velocity.')

			!  ========== eastward_wind ============
			ALLOCATE(eastward_wind(NBEAMS)) ! Dynamic memory allocation
			RETVAL = NF_INQ_VARID(NID, 'eastward_wind', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named eastward_wind found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, eastward_wind)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read eastward_wind.')

			!  ========== northward_wind ============
			ALLOCATE(northward_wind(NBEAMS)) ! Dynamic memory allocation
			RETVAL = NF_INQ_VARID(NID, 'northward_wind', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named northward_wind found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, northward_wind)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read northward_wind.')

			!  ========== vertical_wind ============
			ALLOCATE(vertical_wind(NBEAMS)) ! Dynamic memory allocation
			RETVAL = NF_INQ_VARID(NID, 'vertical_wind', I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named vertical_wind found.')
			RETVAL = NF_GET_VAR_REAL(NID, I, vertical_wind)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read vertical_wind.')
			
			!
			! Read two-dimensional arrays such as DBZ. Read DBZ integers and scale them into real dbz values.
			!
			!  ================================= FIELD1  =================================
			ALLOCATE(FIELD1(NRANGES,NBEAMS))
			ALLOCATE(I2DATA(NRANGES,NBEAMS))

			RETVAL = NF_INQ_VARID(NID, FIELD1_name, I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named ' // FIELD1_name // ' found.')
			RETVAL = NF_GET_VAR_INT2(NID, I, I2DATA)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read ' //  FIELD1_name)
			!
			! Also need to get scale, offset which are attributes for this variable.
			!

			! Change attibutes: from missing_value to: _FillValue

			RETVAL = NF_GET_ATT_REAL(NID,I,'scale_factor', SCALE)
			CALL CHECK_NETCDF(RETVAL, 'Failed to get scale attribute for ' // FIELD1_name)
			RETVAL = NF_GET_ATT_REAL(NID,I,'add_offset', OFFSET)
			CALL CHECK_NETCDF(RETVAL, 'Failed to get offset attribute for ' // FIELD1_name)
			!       RETVAL = NF_GET_ATT_INT2(NID,I,'missing_value', MISSING)
			RETVAL = NF_GET_ATT_INT2(NID,I,'_FillValue', MISSING)
			CALL CHECK_NETCDF(RETVAL, 'Failed to get missing attribute for ' // FIELD1_name)
			!
			! Apply scale, offset to get dbz values. Count how many were missing as well.
			!
			NUM_MISSING = 0
			FIRST = .TRUE.
			DO J=1,NBEAMS
					DO K=1,NRANGES
							IF (I2DATA(K,J) .EQ. MISSING) THEN
									FIELD1(K,J) = BAD_DATA_VALUE
									NUM_MISSING = NUM_MISSING + 1
							ELSE
									FIELD1(K,J) = SCALE*I2DATA(K,J) + OFFSET
									IF (FIRST) THEN
											MIN = FIELD1(K,J)
											MAX = MIN
											FIRST = .FALSE.
									ELSE 
											IF (FIELD1(K,J) .LT. MIN) MIN = FIELD1(K,J)
											IF (FIELD1(K,J) .GT. MAX) MAX = FIELD1(K,J)
									END IF
							END IF
					END DO
			END DO

			DEALLOCATE(I2DATA) ! We are done with the interger array

			PRINT *, NUM_MISSING, ' of ', NBEAMS*NRANGES, ' gates had missing dbz values.'
			PRINT *, FIELD1_name // ' Min, Max: ', MIN, ' to ', MAX
			  

			!  ================================= FIELD2 =================================
			ALLOCATE(FIELD2(NRANGES,NBEAMS))
			ALLOCATE(I2DATA(NRANGES,NBEAMS))

			RETVAL = NF_INQ_VARID(NID, FIELD2_name, I)
			CALL CHECK_NETCDF(RETVAL, 'No variable named ' // FIELD2_name // ' found.')
			RETVAL = NF_GET_VAR_INT2(NID, I, I2DATA)
			CALL CHECK_NETCDF(RETVAL, 'Failed to read ' // FIELD2_name )
			!
			! Also need to get scale, offset which are attributes for this variable.
			!
			RETVAL = NF_GET_ATT_REAL(NID,I,'scale_factor', SCALE)
			CALL CHECK_NETCDF(RETVAL, 'Failed to get scale attribute for ' // FIELD2_name)
			RETVAL = NF_GET_ATT_REAL(NID,I,'add_offset', OFFSET)
			CALL CHECK_NETCDF(RETVAL, 'Failed to get offset attribute for ' // FIELD2_name)
			RETVAL = NF_GET_ATT_INT2(NID,I,'_FillValue', MISSING)
			CALL CHECK_NETCDF(RETVAL, 'Failed to get missing attribute for ' // FIELD2_name)
			!
			! Apply scale, offset to get dbz values. Count how many were missing as well.
			!
			NUM_MISSING = 0
			FIRST = .TRUE.
			DO J=1,NBEAMS
					DO K=1,NRANGES
							IF (I2DATA(K,J) .EQ. MISSING) THEN
									FIELD2(K,J) = BAD_DATA_VALUE
									NUM_MISSING = NUM_MISSING + 1
							ELSE
									FIELD2(K,J) = SCALE*I2DATA(K,J) + OFFSET
									IF (FIRST) THEN
											MIN = FIELD2(K,J)
											MAX = MIN
											FIRST = .FALSE.
									ELSE
											IF (FIELD2(K,J) .LT. MIN) MIN = FIELD2(K,J)
											IF (FIELD2(K,J) .GT. MAX) MAX = FIELD2(K,J)
									END IF
							END IF
					END DO
			END DO

			DEALLOCATE(I2DATA) ! We are done with the interger array

			PRINT *, NUM_MISSING, ' of ', NBEAMS*NRANGES, ' gates had missing vr values.'
			PRINT *, FIELD2_name // ' Min, Max: ', MIN, ' to ', MAX


			!  ==================FAKE NAN FIELD FOR NOAA P3===================
			ALLOCATE(FIELDNAN(NRANGES,NBEAMS))		
			NUM_MISSING = 0
			FIRST = .TRUE.
			DO J=1,NBEAMS
					DO K=1,NRANGES
									FIELDNAN(K,J) =  1.0
					END DO
			END DO

			! ====================================
			! Write the ray informatin to the output file
			! =====================================        
			DO I=1,NBEAMS
					counter = I+T_RAYS
					WRITE(30,101)counter,filename,sweep_number(nsweeps) &
					,NBEAMS,NRANGES &
					,time_coverage_start(1:4),time_coverage_start(6:7) &
					,time_coverage_start(9:10) &
					,time_coverage_start(12:13),time_coverage_start(15:16) &
					,time_coverage_start(18:19),time(I) & 
					,azimuth(I),elevation(I),latitude(I),longitude(I),altitude(I) & 
					,altitude_agl(I),heading(I),roll(I),pitch(I),drift(I) & 
					,rotation(I),tilt(I),eastward_velocity(I) & 
					,northward_velocity(I),vertical_velocity(I),eastward_wind(I) &
					,northward_wind(I) & 
					,vertical_wind(I),azimuth_correction,elevation_correction & 
					,range_correction,longitude_correction,latitude_correction & 
					,pressure_altitude_correction,altitude_correction & 
					,eastward_velocity_correction,northward_velocity_correction & 
					,vertical_velocity_correction,heading_correction & 
					,roll_correction,pitch_correction,drift_correction & 
					,rotation_correction,tilt_correction 
					WRITE(30,102)counter,(range(J), J=1,nranges)     
					WRITE(30,102)counter,(FIELD1(J,I),J=1,nranges)
					WRITE(30,102)counter,(FIELDNAN(J,I),J=1,nranges)
					WRITE(30,102)counter,(FIELD2(J,I),J=1,nranges)
					WRITE(30,102)counter,(FIELDNAN(J,I),J=1,nranges)
			ENDDO

			101   format(I10,2x,A50,3I10,A5,5A3,D20.8,2F10.4,3D20.8,29F10.3)
			102   format(I10,800F10.4) 


			! Increase the ray counter by NBEAMS; this is done after each sweep
			T_RAYS = T_RAYS + NBEAMS          

			! Close the file just read
			RETVAL = NF_CLOSE(NID)
			IF (RETVAL .NE. 0) THEN
					PRINT *, 'INFO : Unable to close ', TRIM(FILENAME) ! Really just info, not abnormal termination
			END IF


			! ====== DEALLOCATE all dynamical arrays  ============
			!  One dimensional arrays
			DEALLOCATE(sweep_number)

			DEALLOCATE(time)
			DEALLOCATE(range)
			DEALLOCATE(azimuth)
			DEALLOCATE(elevation)
			DEALLOCATE(latitude)
			DEALLOCATE(longitude)
			DEALLOCATE(altitude)
			DEALLOCATE(altitude_agl)
			DEALLOCATE(heading)
			DEALLOCATE(roll)
			DEALLOCATE(pitch)
			DEALLOCATE(drift)
			DEALLOCATE(rotation)
			DEALLOCATE(tilt)
			DEALLOCATE(eastward_velocity)
			DEALLOCATE(northward_velocity)
			DEALLOCATE(vertical_velocity)
			DEALLOCATE(eastward_wind)
			DEALLOCATE(northward_wind)
			DEALLOCATE(vertical_wind)

			! two dimensional arrays
			DEALLOCATE(FIELD1)
			DEALLOCATE(FIELD2)
			DEALLOCATE(FIELDNAN)
			  !DEALLOCATE(VR)
			  !DEALLOCATE(SW)
			  
			! close the output file
			CLOSE(30)

	ENDDO       !   end of netcdf file loop, go to next file

	command='rm filenames.txt'
	CALL system(command)

	CLOSE(20)
	
	STOP
	END

	!---------------------------------------------------
	!
	! Routine to check if a netcdf call has gone well.
	!

	SUBROUTINE CHECK_NETCDF(RETVAL, ERRSTR)

	INTEGER RETVAL
	CHARACTER(*) ERRSTR

	IF (RETVAL .NE. 0) THEN
	PRINT *, TRIM(ERRSTR)
	PRINT *, 'Netcdf return value was ', RETVAL
	PRINT *, 'Abnormal termination.'
	STOP
	ENDIF

	RETURN
	END

