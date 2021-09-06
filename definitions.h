#define  PHYSICS                        HD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     VECTOR
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            4

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  THETA                          0
#define  MPOLY                          1
#define  AMP                            2
#define  SEED                           3

/* [Beg] user-defined constants (do not change this line) */

// Adimensioned system :
// L=1; V=sqrt(L/T0)=1; rho=rho_0=1
#define UNIT_DENSITY  1.0
#define UNIT_LENGTH   1.0
#define UNIT_VELOCITY 1.0

/* [End] user-defined constants (do not change this line) */
