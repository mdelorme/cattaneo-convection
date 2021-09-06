/* ///////////////////////////////////////////////////////////////////// */
/*! 
  Slab cartésien pour la convection inspiré de Catteneo et al 1991.
  
  Le domaine est initialisé sur un modèle polytropique : 
  T(z)   = (1.0 + theta*z)
  rho(z) = (1.0 + theta*z)**m
  P(z)   = (1.0 + theta*z)**(m+1)

  supplémenté d'une équation d'état de gaz parfaits : P = rho * T

  avec :
   . z la profondeur
   . m l'indice polytropique
   . theta le gradient de température.

  dans ce setup, la gravité s'exprime comme g = theta*(m+1). Elle est imposée
  de manière constante sur le domaine pour préserver l'équilibre hydrostatique
  à l'initialisation.

  L'instabilité est déclenchée en imposant une perturbation aléatoire sur la pression

  Conditions au bord : 
   . En haut du domaine (z=0) -> Température Ttop imposée.
   . En bas du domaine  (z=1) -> Gradient de température theta imposé.
   . Dans les deux cas : 
     * Gradient de vitesses horizontales nulles
     * Vitesses verticales nulles
     * Pression obtenue par l'équilibre hydrostatique
     * Densité obtenue par l'équation d'état.
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include <assert.h>

int first = 1;

const double Ttop = 1.0;

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3) {

  // Paramètres du run
  int    seed  = g_inputParam[SEED];
  double mpoly = g_inputParam[MPOLY];
  double theta = g_inputParam[THETA];
  double amp   = g_inputParam[AMP];

  // One-time init par processus
  if (first) {
    g_gamma = 5.0/3.0;
    srand(seed*prank);
    first = 0;
  }

  // Choix de la dimension représentant la verticale
  double z = (DIMENSIONS == 2 ? x2 : x3);

  // Perturbation aléatoire
  double pert = (((float)rand() / RAND_MAX) - 0.5) * amp;

  // ICs
  v[RHO] = pow(1.0 + theta*z, mpoly);
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
  v[PRS] = pow(1.0 + theta*z, mpoly+1.0) * (1.0 + pert);
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{
}

void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) {
  /**
   * Notes sur les tests effectués sur les conditions au bord.
   * 
   * Chaque test de setup est validé par :
   *  . Un plot de l'évolution de la masse et des énergies dans le domaine. 
   *  . Un plot de l'évolution du profil de température tous les 50 pas de temps
   *  . Un film de l'évolution du domaine présentant rho, rho-rho_mean, T, T-T_mean, uz, uh
   *
   * Jusqu'à présent TOUS les setups montrent :
   *  . La préservation du gradient de température en bas du domaine fonctionne correctement
   *  . La dérive et l'augmentation de la température en haut du domaine
   *  . Une variation de la masse totale dans le domaine.
   *
   * Il existe une ambiguité sur la définition des conditions en haut du domaine. En 
   * théorie ce que l'on veut c'est que la température à l'interface soit Ttop. 
   * Puisque l'on passe par les cellules fantômes, la définition de la température 
   * dans ces zones n'est pas claire. Il existe donc plusieurs setups possibles, et le 
   * résultat va beaucoup varier en fonction de ce que l'on choisit. L'implémentation 
   * actuelle prend comme valeur dans les fantômes la valeur de Ttop.
   */
  int   i, j, k, nv;
  double  *x1, *x2, *x3;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  // Paramètres du run
  double theta = g_inputParam[THETA];
  double m     = g_inputParam[MPOLY];
  double gval  = theta * (m+1);
  
  if (side == X2_BEG && DIMENSIONS == 2) {  // z = 0 -> Haut du domaine
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){
	// État de référence à l'intérieur du domaine
	double rho_ref = d->Vc[RHO][k][JBEG][i];
	double P_ref   = d->Vc[PRS][k][JBEG][i];
	double T_ref   = P_ref / rho_ref;

	// Échelle de pression et altitude
	double Hp  = T_ref / gval;
	double dz  = x2[JBEG] - x2[j];
	
	// Calcul des quantités locales
	double P   = P_ref * exp(-dz / Hp);
	double T   = Ttop;
	double rho = P / T;

	// Précautions : On vérifie que tout tient la route
	if (P > P_ref) {
	  printf("Top : P=%lf; Pref=%lf; dz=%lf\n", P, P_ref, dz);
	  MPI_Abort(MPI_COMM_WORLD, 1);
	}
	else if (T > T_ref || T < 0.0) {
	  printf("Top : rho=%lf; rho_ref=%lf; P=%lf; Pref=%lf; T=%lf; T_ref=%lf; dz=%lf\n",
		 rho, rho_ref, P, P_ref, T, T_ref, dz);
	  MPI_Abort(MPI_COMM_WORLD, 2);
	}

	// Cellule de ref pour les vitesses : cellule miroir à la limite.
	int jref = 2*JBEG - (j+1);

	// Et remplissage de la condition au bord
	d->Vc[PRS][k][j][i] = P;
	d->Vc[RHO][k][j][i] = rho;
	d->Vc[VX1][k][j][i] =  d->Vc[VX1][k][jref][i];
	d->Vc[VX2][k][j][i] = -d->Vc[VX2][k][jref][i];
      }
    }
  }

  if (side == X2_END && DIMENSIONS == 2){  // z = 1 -> Bas du domaine
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){
	// État de référence à l'intérieur du domaine
	double rho_ref = d->Vc[RHO][k][JEND][i];
	double P_ref   = d->Vc[PRS][k][JEND][i];
	double T_ref   = P_ref / rho_ref;

	// Échelle de pression et altitude.
	double Hp  = T_ref / gval; 
	double dz  = x2[j] - x2[JEND];

	// Calcul des quantités locales
	double P   = P_ref * exp(dz / Hp);
	double T   = T_ref + dz * g_inputParam[THETA];
	double rho = P / T;

	// Précautions
	if (P < P_ref) {
	  printf("Bottom : P=%lf; Pref=%lf; dz=%lf\n", P, P_ref, dz);
	  MPI_Abort(MPI_COMM_WORLD, 1);
	}
	else if (T < T_ref) {
	  printf("Top : rho=%lf; rho_ref=%lf; P=%lf; Pref=%lf; T=%lf; T_ref=%lf; dz=%lf\n",
		 rho, rho_ref, P, P_ref, T, T_ref, dz);
	  MPI_Abort(MPI_COMM_WORLD, 2);
	}

	// Référence pour les vitesses
	int jref = 2*JEND - (j-1);

	// Et remplissage de la condition au bord
	d->Vc[PRS][k][j][i] = P;
	d->Vc[RHO][k][j][i] = rho;
	d->Vc[VX1][k][j][i] =  d->Vc[VX1][k][jref][i];
	d->Vc[VX2][k][j][i] = -d->Vc[VX2][k][jref][i];
      }
    }
  }

  if (side == X3_BEG && DIMENSIONS == 3){  // z = 0 -> Haut du domaine
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){

      }
    }
  }

  if (side == X3_END){  // z = 1 -> Bas du domaine
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){

      }
    }
  }
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3) {

  // Calcul de l'accélération gravitationnelle g = theta*(m+1)
  // On inverse le signe si on est en dehors du domaine pour aider avec les conditions
  // de barrière impénétrable.
  // Note : Cela ne fait que très peu de différences sur le résultat.
  double sgn = 1.0;
  double z = (DIMENSIONS == 2 ? x2 : x3);
  if (z > 1.0 || z < 0.0)
    sgn = -1.0;
  double gval = g_inputParam[THETA] * (g_inputParam[MPOLY] + 1.0) * sgn;
  
#if DIMENSIONS==2
  g[IDIR] = 0.0;
  g[JDIR] = gval;
  g[KDIR] = 0.0;
#else
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = gval;
#endif
  
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  return 0.0;
}
#endif
