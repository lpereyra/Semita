#ifndef LEESNAP_H
#define LEESNAP_H

/* Header for the standard 
 * file GADGET format.
 */
struct io_header
{
  int npart[N_part_types];             /*!< number of particles of each type in this file */
  double mass[N_part_types];           /*!< mass of particles of each type. If 0, then the masses are explicitly
                                            stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;                         /*!< time of snapshot file */
  double redshift;                     /*!< redshift of snapshot file */
  int flag_sfr;                        /*!< flags whether the simulation was including star formation */
  int flag_feedback;                   /*!< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[N_part_types];          /*!< total number of particles of each type in this snapshot. This can be
                                            different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;                    /*!< flags whether cooling was included  */
  int num_files;                       /*!< number of files in multi-file snapshot */
  double BoxSize;                      /*!< box-size of simulation in case periodic boundaries were used */
  double Omega0;                       /*!< matter density in units of critical density */
  double OmegaLambda;                  /*!< cosmological constant parameter */
  double HubbleParam;                  /*!< Hubble parameter in units of 100 km/sec/Mpc */
  int flag_stellarage;                 /*!< flags whether the file contains formation times of star particles */
  int flag_metals;                     /*!< flags whether the file contains metallicity values for gas and star particles */
  unsigned int npartTotalHighWord[N_part_types];  /*!< High word of the total number of particles of each type */
  int  flag_entropy_instead_u;         /*!< flags that IC-file contains entropy instead of u */
  char fill[60];	                     /*!< fills to 256 Bytes */
};                                     /*!< holds header for snapshot files */


/* Generic read function */
extern void read_gadget(void);
#ifdef CHANGE_POSITION
  extern void change_positions(TYPE_INT);
#endif

#endif
