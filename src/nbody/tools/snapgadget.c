/*
 * SNAPGADGET.C: convert binary N-body file to GADGET binary format.
 *         KYLE OMAN
 */

#include "stdinc.h"
#include "getparam.h"
#include "vectmath.h"
#include "filestruct.h"
#include "bodytags.h"

#  define PURPOSE		";Convert binary N-body file to GADGET binary"

#include "proto.h"
#include "allvars.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

string defv[] = {
  PURPOSE,
  "in=???",                   ";Input file name",
  "out=???",                  ";Output file name",
  "rmax=???",                 ";Rescale r coordinate",
  "vmax=???",                 ";Rescale v coordinate",
  "options=" MassTag "," PosTag "," VelTag,
  //"gadgetformat=" 1,
  ";Others are " PhiTag "," AuxTag "," KeyTag ", use KeyTag only (for particle IDs)",
  "VERSION=2.3",		";Josh Barnes  11 June 1998, mod by Kyle Oman 2013",
  ";mod is certainly work in progress, works for a specific use-case, not general!",
  NULL,
};

local string options;

local stream instr;

local float rmax, vmax;

static long int n_type[6];

//TODO review reqd functions
local void snapgadget(void);
local void out_header(int, real);
local void out_idata(int *, int);
local void out_rdata(real *, int);
local void out_vdata(real *, int);
local void out_phase(real *, int);
local void out_int(int);
local void out_3int(int, int, int);
local void out_real(real);
local void out_vect(vector);

local void fill_write_buffer_zeno(enum iofields, long int*, long int, int, int*, int*, real*, real*, real*, real*, real*, real*);

int main(int argc, string argv[])
{
  initparam(argv, defv);
  instr = stropen(getparam("in"), "r");
  rmax = atof(getparam("rmax"));
  vmax = atof(getparam("vmax"));
  get_history(instr);
  options = getparam("options"); //TODO mod to take only m, xyz, vxyz?
  //All.SnapFormat = getparam("gadgetformat"); TODO make this work
  All.SnapFormat = 1;

  snapgadget();
  return (0);
}

local void snapgadget(void)
{
  long int nbody;
  int *key_buf, *type_buf;
  real tsnap, *mass_buf, *pos_buf, *vel_buf, *phase_buf, *phi_buf, *aux_buf;
  int allocated_key = 0, allocated_type = 0, allocated_mass = 0, allocated_pos = 0, allocated_vel = 0, allocated_phase = 0, allocated_phi = 0, allocated_aux = 0;

  while (get_tag_ok(instr, SnapShotTag)) {	/* loop reading data frames */
    get_set(instr, SnapShotTag);
    get_set(instr, ParametersTag);
    if (get_tag_ok(instr, TimeTag))
      get_data(instr, TimeTag, RealType, &tsnap, 0);
    else
      tsnap = 0.0;
    get_data(instr, NBodyTag, IntType, &nbody, 0);
    get_tes(instr, ParametersTag);
    get_set(instr, ParticlesTag);
    if (scanopt(options, MassTag)) {	/* mass data to convert?    */
      mass_buf = (real *) allocate_long((unsigned long int) nbody * (unsigned long int) sizeof(real));
      allocated_mass = 1;
      get_data(instr, MassTag, RealType, mass_buf, nbody, 0);
    }
    if (scanopt(options, PosTag)) {		/* position data?           */
      pos_buf = (real *) allocate_long((unsigned long int) nbody * (unsigned long int) NDIM * (unsigned long int) sizeof(real));
      allocated_pos = 1;
      get_data(instr, PosTag, RealType, pos_buf, nbody, NDIM, 0);
    }
    if (scanopt(options, VelTag)) {		/* velocity data?           */
      vel_buf = (real *) allocate_long((unsigned long int) nbody * (unsigned long int) NDIM * (unsigned long int) sizeof(real));
      allocated_vel = 1;
      get_data(instr, VelTag, RealType, vel_buf, nbody, NDIM, 0);
    }
    if (scanopt(options, PhaseTag)) {	/* phase-space data?        */
      phase_buf = (real *) allocate_long((unsigned long int) nbody * (unsigned long int) 2 * (unsigned long int) NDIM * (unsigned long int) sizeof(real));
      allocated_phase = 1;
      get_data(instr, PhaseTag, RealType, phase_buf, nbody, 2, NDIM, 0);
    }
    if (scanopt(options, PhiTag)) {		/* potential data?          */
      phi_buf = (real *) allocate_long((unsigned long int) nbody * (unsigned long int) sizeof(real));
      allocated_phi = 1;
      get_data(instr, PhiTag, RealType, phi_buf, nbody, 0);
    }
    if (scanopt(options, AuxTag)) {		/* auxiliary data?          */
      aux_buf = (real *) allocate_long((unsigned long int) nbody * (unsigned long int) sizeof(real));
      allocated_aux = 1;
      get_data(instr, AuxTag, RealType, aux_buf, nbody, 0);
    }
    if (scanopt(options, KeyTag)) {		/* key data?                */
      key_buf = (int *) allocate_long((unsigned long int) nbody * (unsigned long int) sizeof(int));
      allocated_key = 1;
      get_data(instr, KeyTag, IntType, key_buf, nbody, 0);
    }
    get_tes(instr, ParticlesTag);
    get_tes(instr, SnapShotTag);
  }
  
  //add a buffer for gadget type (here DM only)
  
  type_buf = (int *) allocate_long((unsigned long int) nbody * (unsigned long int) sizeof(int));
  allocated_type = 1;
  long int loop;
  for(loop = 0; loop < nbody; loop++)
    type_buf[loop] = 1;
  
  //check for IDs, if none, allocate

  if(allocated_key == 0)
    {
      key_buf = (int *) allocate_long((unsigned long int) nbody * (unsigned long int) sizeof(int));
      allocated_key = 1;
      for(loop = 0; loop < nbody; loop++)
	key_buf[loop] = loop;
    }

  //now have all data read into bufs, proceed to write a gadget binary

  All.BufferSize = 100000;
  
  //prepare the gadget header struct
  
  long int npart, bytes_per_blockelement;
  int type, nextblock;
  long int typelist[6];
  long int n_for_this_task, ntask, n, p, pc, offset = 0, task;
  int blockmaxlen, ntot_type[6], nn[6];
  enum iofields blocknr;
  int blksize;
  //MPI_Status status;
  FILE *fd = 0;
  
#define SKIP  {my_fwrite(&blksize,sizeof(int),1,fd);}
  
  /* determine particle numbers of each tiype in file */ //TODO for now just 1 particle type, type 1 for halo DM particles
  
  for(n = 0; n < 6; n++)
    {
      ntot_type[n] = 0;
    }
  ntot_type[1] = nbody;
  
  for(n = 0; n < 6; n++)
    {
      n_type[n] = 0;
    }
  n_type[1] = nbody;

  /* fill file header */
  
  for(n = 0; n < 6; n++)
    {
      header.npart[n] = ntot_type[n];
      //header.npartTotal[n] = (unsigned int) ntot_type_all[n]; //TODO support for multiple files?
      header.npartTotal[n] = (unsigned int) ntot_type[n]; //TODO this supports single file only
      //header.npartTotalHighWord[n] = (unsigned int) (ntot_type_all[n] >> 32);
      header.npartTotalHighWord[n] = 0;
    }

  printf("NPARTICLES: %i\n",nbody);

  for(n = 0; n < 6; n++)
    {
      All.MassTable[n] = 0.0;
      if(n == 1)
	{
	  int write_masses = 0;
	  long int m;
	  for(m = 1; m < nbody; m++)
	    {
	      if(mass_buf[m] != mass_buf[m-1])
		write_masses = 1;
	    }
	  if(write_masses == 0)
	    All.MassTable[n] = mass_buf[0] / (vmax * vmax * rmax);
	  else
	    All.MassTable[n] = 0;
	}
    }

  for(n = 0; n < 6; n++)
    header.mass[n] = All.MassTable[n];
  
  All.Time = tsnap;
  header.time = All.Time;
  
  All.ComovingIntegrationOn = 0; //ZENO doesn't do cosmological ICs, I think - support for comoving integration has been REMOVED!
  header.redshift = 0;
  
  header.flag_sfr = 0;
  header.flag_feedback = 0;
  header.flag_cooling = 0;
  header.flag_stellarage = 0;
  header.flag_metals = 0;
  
  /* COOLING, SFR, STELLARAGE, METALS not defined in ZENO, typically... esp. not for simple DM only case*/
  
#ifdef COOLING
  header.flag_cooling = 1;
#endif
#ifdef SFR
  header.flag_sfr = 1;
  header.flag_feedback = 1;
#ifdef STELLARAGE
  header.flag_stellarage = 1;
#endif
#ifdef METALS
  header.flag_metals = 1;
#endif
#endif

  //TODO recheck these parameters
  All.NumFilesPerSnapshot = 1;
  All.BoxSize = 0.0; //not needed for non-periodic BCs, right?
  All.Omega0 = 0;
  All.OmegaLambda = 0; //comoving is off
  All.HubbleParam = 1.0;
  
  header.num_files = All.NumFilesPerSnapshot;
  header.BoxSize = All.BoxSize;
  header.Omega0 = All.Omega0;
  header.OmegaLambda = All.OmegaLambda;
  header.HubbleParam = All.HubbleParam;
  
    /* Write header data to file */
  
  if(!(fd = fopen( (char*) getparam("out"), "w")))
    {
      printf("can't open file `%s' for writing snapshot.\n", (char*) getparam("out"));
      exit(123);
    }
  
  if(All.SnapFormat == 2)
    {
      blksize = sizeof(int) + 4 * sizeof(char);
      SKIP;
      my_fwrite("HEAD", sizeof(char), 4, fd);
      nextblock = sizeof(header) + 2 * sizeof(int);
      my_fwrite(&nextblock, sizeof(int), 1, fd);
      SKIP;
    }
  
  blksize = sizeof(header);
  SKIP;
  my_fwrite(&header, sizeof(header), 1, fd);
  SKIP;
  printf("wrote header\n");

  /* Write particle data to file. */
  
  for(blocknr = 0; blocknr < IO_NBLOCKS; blocknr++)
    {
      if(blockpresent(blocknr))
	{
	  bytes_per_blockelement = get_bytes_per_blockelement(blocknr);
	  blockmaxlen = ((int) (All.BufferSize * 1024 * 1024)) / bytes_per_blockelement;
	  npart = get_particles_in_block(blocknr, &typelist[0]);
	  
	  if(npart > 0)
	    {
	      if(All.SnapFormat == 1 || All.SnapFormat == 2)
		{
		  if(All.SnapFormat == 2)
		    {
		      blksize = sizeof(int) + 4 * sizeof(char);
		      SKIP;
		      my_fwrite(Tab_IO_Labels[blocknr], sizeof(char), 4, fd);
		      nextblock = npart * bytes_per_blockelement + 2 * sizeof(int);
		      my_fwrite(&nextblock, sizeof(int), 1, fd);
		      SKIP;
		    }
		  
		  blksize = npart * bytes_per_blockelement;
		  SKIP;
		}
	      
	      for(type = 0; type < 6; type++)
		{
		  if(typelist[type])
		    {
		      n_for_this_task = n_type[type];
		      while(n_for_this_task > 0)
			{
			  pc = n_for_this_task;
			  
			  if(pc > blockmaxlen)
			    pc = blockmaxlen;
			  
			  printf("allocate comm: %u\n", blockmaxlen);
			  //THIS REGION SOMEWHAT BROKEN? ALWAYS BEEN A PAIN IN THE ASS. LOOKS LIKE IS CAPPED AT 100 000 000 PARTICLES FOR NOW.
			  CommBuffer = (void*) allocate_long((unsigned long int) nbody * (unsigned long int) NDIM * (unsigned long int) sizeof(real));
			  //CommBuffer = (void*) allocate_long((unsigned long int) blockmaxlen);
			  printf("all allocated!\n");
			  
			  printf("buffering\n");
			  fill_write_buffer_zeno(blocknr, &offset, pc, type, key_buf, type_buf, mass_buf, pos_buf, vel_buf, phase_buf, phi_buf, aux_buf);
			  printf("done buffering\n");
			  
			  my_fwrite(CommBuffer, bytes_per_blockelement, pc, fd);
			  
			  free(CommBuffer);
			  
			  n_for_this_task -= pc;
			}
		    }
		}
	      
	      if(All.SnapFormat == 1 || All.SnapFormat == 2)
		SKIP;
	    }
	}
    }
  
  fclose(fd);
  
  //TODO re-check that these free's make sense with allocate's above
  if(allocated_mass)
    free(mass_buf);
  if(allocated_pos)
    free(pos_buf);
  if(allocated_vel)
    free(vel_buf);
  if(allocated_phase)
    free(phase_buf);
  if(allocated_phi)
    free(phi_buf);
  if(allocated_aux)
    free(aux_buf);
  if(allocated_key)
    free(key_buf);
  if(allocated_type)
    free(type_buf);
}

/*! This function tells whether or not a given block in the output file is
 *  present, depending on the type of simulation run and the compile-time
 *  options. If one wants to add a new output-block, this function should be
 *  augmented accordingly.
 */
int blockpresent(enum iofields blocknr)
{

      //TODO review makefile for these flags

#ifndef OUTPUTPOTENTIAL
  if(blocknr == IO_POT)
    return 0;
#endif

#ifndef OUTPUTACCELERATION
  if(blocknr == IO_ACCEL)
    return 0;
#endif

#ifndef OUTPUTCHANGEOFENTROPY
  if(blocknr == IO_DTENTR)
    return 0;
#endif

#ifndef OUTPUTTIMESTEP
  if(blocknr == IO_TSTP)
    return 0;
#endif

  return 1;			/* default: present */
}

/*! This function tells the size of one data entry in each of the blocks
 *  defined for the output file. If one wants to add a new output-block, this
 *  function should be augmented accordingly.
 */
int get_bytes_per_blockelement(enum iofields blocknr)
{
  int bytes_per_blockelement = 0;

  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ACCEL:
      bytes_per_blockelement = 3 * sizeof(float);
      break;

    case IO_ID:
#ifdef LONGIDS 
      bytes_per_blockelement = sizeof(long long);
#else
      bytes_per_blockelement = sizeof(int);
#endif
      break;

    case IO_MASS:
    case IO_U:
    case IO_RHO:
    case IO_HSML:
    case IO_POT:
    case IO_DTENTR:
    case IO_TSTP:
      bytes_per_blockelement = sizeof(float);
      break;
    }

  return bytes_per_blockelement;
}

/*! This function determines how many particles there are in a given block,
 *  based on the information in the header-structure.  It also flags particle
 *  types that are present in the block in the typelist array. If one wants to
 *  add a new output-block, this function should be augmented accordingly.
 */
long int get_particles_in_block(enum iofields blocknr, long int *typelist)
{
  long int i, nall, ntot_withmasses, ngas, nstars;

  nall = 0;
  ntot_withmasses = 0;

  for(i = 0; i < 6; i++)
    {
      typelist[i] = 0;

      if(header.npart[i] > 0)
	{
	  nall += header.npart[i];
	  typelist[i] = 1;
	}

      if(All.MassTable[i] == 0)
	ntot_withmasses += header.npart[i];
    }

  ngas = header.npart[0];
  nstars = header.npart[4];


  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ACCEL:
    case IO_TSTP:
    case IO_ID:
    case IO_POT:
      return nall;
      break;

    case IO_MASS:
      for(i = 0; i < 6; i++)
	{
	  typelist[i] = 0;
	  if(All.MassTable[i] == 0 && header.npart[i] > 0)
	    typelist[i] = 1;
	}
      return ntot_withmasses;
      break;

    case IO_U:
    case IO_RHO:
    case IO_HSML:
    case IO_DTENTR:
      for(i = 1; i < 6; i++)
	typelist[i] = 0;
      return ngas;
      break;
    }

  exit(212);
  return 0;
}

/*! This catches I/O errors occuring for my_fwrite(). In this case we
 *  better stop.
 */
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nwritten;

  if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
    {
      printf("I/O error (fwrite) on task=%d has occured: %s\n", ThisTask, strerror(errno));
      fflush(stdout);
      exit(777);
    }
  return nwritten;
}

/*! This function fills the write buffer with particle data. New output blocks
 *  can in principle be added here.
 */
 void fill_write_buffer_zeno(enum iofields blocknr, long int *startindex, long int pc, int type, int* key_buf, int* type_buf, real* mass_buf, real* pos_buf, real* vel_buf, real* phase_buf, real* phi_buf, real* aux_buf)
{
  long int n, pindex;
  int k;
  float *fp;

#ifdef LONGIDS
  long long *ip;
#else
  int *ip;
#endif

#ifdef PERIODIC
  FLOAT boxSize;
#endif
  double dt_gravkick, dt_hydrokick, a3inv = 1, fac1, fac2;

  a3inv = fac1 = fac2 = 1;

  fp = CommBuffer;
  ip = CommBuffer;

  pindex = 0;

  switch (blocknr)
    {
    case IO_POS:		/* positions */
      for(n = 0; n < pc; pindex++)
	if(type_buf[pindex] == type)
	  {
	    for(k = 0; k < 3; k++)
	      {
		fp[k] = pos_buf[pindex * NDIM + k] / rmax;	
/*#ifdef PERIODIC
		boxSize = All.BoxSize;
#ifdef LONG_X
		if(k == 0)
		  boxSize = All.BoxSize * LONG_X;
#endif
#ifdef LONG_Y
		if(k == 1)
		  boxSize = All.BoxSize * LONG_Y;
#endif
#ifdef LONG_Z
		if(k == 2)
		  boxSize = All.BoxSize * LONG_Z;
#endif
		while(fp[k] < 0)
		  fp[k] += boxSize;
		while(fp[k] >= boxSize)
		  fp[k] -= boxSize;
#endif*/
	      }
	    n++;
	    fp += 3;
	  }
      break;

    case IO_VEL:		/* velocities */
      for(n = 0; n < pc; pindex++)
	if(type_buf[pindex] == type)
	  {
	    //removed gravkick & hydrokick stuff, only relevant for comoving integration
	    for(k = 0; k < 3; k++)
	      {
		fp[k] = vel_buf[pindex * NDIM + k] / vmax;
	      }
	    for(k = 0; k < 3; k++)
	      fp[k] *= sqrt(a3inv); //no effect

	    n++;
	    fp += 3;
	  }
      break;

    case IO_ID:		/* particle ID */
      for(n = 0; n < pc; pindex++)
	if(type_buf[pindex] == type)
	  {
	    *ip++ = key_buf[pindex];
	    n++;
	  }
      break;

    case IO_MASS:		/* particle mass */
      for(n = 0; n < pc; pindex++)
	if(type_buf[pindex] == type)
	  {
	    *fp++ = mass_buf[pindex] / (vmax * vmax * rmax);
	    n++;
	  }
      break;
      /*      
    case IO_U:			// internal energy
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
#ifdef ISOTHERM_EQS
	    *fp++ = SphP[pindex].Entropy;
#else
	    *fp++ =
	      dmax(All.MinEgySpec,
		   SphP[pindex].Entropy / GAMMA_MINUS1 * pow(SphP[pindex].Density * a3inv, GAMMA_MINUS1));
#endif
	    n++;
	  }
      break;
      */
      /*
    case IO_RHO:		// density
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Density;
	    n++;
	  }
      break;
      */
      /*
    case IO_HSML:		// SPH smoothing length
      for(n = 0; n < pc; pindex++)
	if(type_buf[pindex] == type)
	  {
	    *fp++ = SphP[pindex].Hsml; //won't execute for DM only, TODO if hydro supported
	    n++;
	  }
      break;
      */

    case IO_POT:		/* gravitational potential */
#ifdef OUTPUTPOTENTIAL
      for(n = 0; n < pc; pindex++)
	if(type_buf[pindex] == type)
	  {
	    *fp++ = phi_buf[pindex];
	    n++;
	  }
#endif
      break;
      /*
    case IO_ACCEL:		// acceleration
#ifdef OUTPUTACCELERATION //not supported by zeno, don't define!
      for(n = 0; n < pc; pindex++)
	if(type_buf[pindex] == type)
	  {
	    for(k = 0; k < 3; k++)
	      fp[k] = fac1 * P[pindex].GravAccel[k];

	    if(type_buf[pindex] == 0)
	      for(k = 0; k < 3; k++)
		fp[k] += fac2 * SphP[pindex].HydroAccel[k];
	    fp += 3;
	    n++;
	  }
#endif
      break;
      */
      /*
    case IO_DTENTR:		// rate of change of entropy
#ifdef OUTPUTCHANGEOFENTROPY //not supported by zeno, don't define!
      for(n = 0; n < pc; pindex++)
	if(type_buf[pindex] == type)
	  {
	    *fp++ = SphP[pindex].DtEntropy;
	    n++;
	  }
#endif
      break;
      */
      /*
    case IO_TSTP:		// timestep 
#ifdef OUTPUTTIMESTEP //not supported by zeno, don't define!

      for(n = 0; n < pc; pindex++)
	if(type_buf[pindex] == type)
	  {
	    *fp++ = (P[pindex].Ti_endstep - P[pindex].Ti_begstep) * All.Timebase_interval;
	    n++;
	  }
#endif
      break;
      */
    }

  *startindex = pindex;
}
