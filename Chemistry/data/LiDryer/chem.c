
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#if defined(BL_FORT_USE_UPPERCASE)
#define CKINDX CKINDX
#define CKINIT CKINIT
#define CKFINALIZE CKFINALIZE
#define CKXNUM CKXNUM
#define CKSYME CKSYME
#define CKSYMS CKSYMS
#define CKRP CKRP
#define CKPX CKPX
#define CKPY CKPY
#define CKPC CKPC
#define CKRHOX CKRHOX
#define CKRHOY CKRHOY
#define CKRHOC CKRHOC
#define CKWT CKWT
#define CKAWT CKAWT
#define CKMMWY CKMMWY
#define CKMMWX CKMMWX
#define CKMMWC CKMMWC
#define CKYTX CKYTX
#define CKYTCP CKYTCP
#define CKYTCR CKYTCR
#define CKXTY CKXTY
#define CKXTCP CKXTCP
#define CKXTCR CKXTCR
#define CKCTX CKCTX
#define CKCTY CKCTY
#define CKCPOR CKCPOR
#define CKHORT CKHORT
#define CKSOR CKSOR
#define CKCVML CKCVML
#define CKCPML CKCPML
#define CKUML CKUML
#define CKHML CKHML
#define CKGML CKGML
#define CKAML CKAML
#define CKSML CKSML
#define CKCVMS CKCVMS
#define CKCPMS CKCPMS
#define CKUMS CKUMS
#define CKHMS CKHMS
#define CKGMS CKGMS
#define CKAMS CKAMS
#define CKSMS CKSMS
#define CKCPBL CKCPBL
#define CKCPBS CKCPBS
#define CKCVBL CKCVBL
#define CKCVBS CKCVBS
#define CKHBML CKHBML
#define CKHBMS CKHBMS
#define CKUBML CKUBML
#define CKUBMS CKUBMS
#define CKSBML CKSBML
#define CKSBMS CKSBMS
#define CKGBML CKGBML
#define CKGBMS CKGBMS
#define CKABML CKABML
#define CKABMS CKABMS
#define CKWC CKWC
#define CKWYP CKWYP
#define CKWXP CKWXP
#define CKWYR CKWYR
#define CKWXR CKWXR
#define CKQC CKQC
#define CKKFKR CKKFKR
#define CKQYP CKQYP
#define CKQXP CKQXP
#define CKQYR CKQYR
#define CKQXR CKQXR
#define CKNU CKNU
#define CKNCF CKNCF
#define CKABE CKABE
#define CKEQC CKEQC
#define CKEQYP CKEQYP
#define CKEQXP CKEQXP
#define CKEQYR CKEQYR
#define CKEQXR CKEQXR
#define DWDOT DWDOT
#define VCKHMS VCKHMS
#define VCKPY VCKPY
#define VCKWYR VCKWYR
#define VCKYTX VCKYTX
#define GET_T_GIVEN_EY GET_T_GIVEN_EY
#define GET_T_GIVEN_HY GET_T_GIVEN_HY
#define GET_REACTION_MAP GET_REACTION_MAP
#elif defined(BL_FORT_USE_LOWERCASE)
#define CKINDX ckindx
#define CKINIT ckinit
#define CKFINALIZE ckfinalize
#define CKXNUM ckxnum
#define CKSYME cksyme
#define CKSYMS cksyms
#define CKRP ckrp
#define CKPX ckpx
#define CKPY ckpy
#define CKPC ckpc
#define CKRHOX ckrhox
#define CKRHOY ckrhoy
#define CKRHOC ckrhoc
#define CKWT ckwt
#define CKAWT ckawt
#define CKMMWY ckmmwy
#define CKMMWX ckmmwx
#define CKMMWC ckmmwc
#define CKYTX ckytx
#define CKYTCP ckytcp
#define CKYTCR ckytcr
#define CKXTY ckxty
#define CKXTCP ckxtcp
#define CKXTCR ckxtcr
#define CKCTX ckctx
#define CKCTY ckcty
#define CKCPOR ckcpor
#define CKHORT ckhort
#define CKSOR cksor
#define CKCVML ckcvml
#define CKCPML ckcpml
#define CKUML ckuml
#define CKHML ckhml
#define CKGML ckgml
#define CKAML ckaml
#define CKSML cksml
#define CKCVMS ckcvms
#define CKCPMS ckcpms
#define CKUMS ckums
#define CKHMS ckhms
#define CKGMS ckgms
#define CKAMS ckams
#define CKSMS cksms
#define CKCPBL ckcpbl
#define CKCPBS ckcpbs
#define CKCVBL ckcvbl
#define CKCVBS ckcvbs
#define CKHBML ckhbml
#define CKHBMS ckhbms
#define CKUBML ckubml
#define CKUBMS ckubms
#define CKSBML cksbml
#define CKSBMS cksbms
#define CKGBML ckgbml
#define CKGBMS ckgbms
#define CKABML ckabml
#define CKABMS ckabms
#define CKWC ckwc
#define CKWYP ckwyp
#define CKWXP ckwxp
#define CKWYR ckwyr
#define CKWXR ckwxr
#define CKQC ckqc
#define CKKFKR ckkfkr
#define CKQYP ckqyp
#define CKQXP ckqxp
#define CKQYR ckqyr
#define CKQXR ckqxr
#define CKNU cknu
#define CKNCF ckncf
#define CKABE ckabe
#define CKEQC ckeqc
#define CKEQYP ckeqyp
#define CKEQXP ckeqxp
#define CKEQYR ckeqyr
#define CKEQXR ckeqxr
#define DWDOT dwdot
#define VCKHMS vckhms
#define VCKPY vckpy
#define VCKWYR vckwyr
#define VCKYTX vckytx
#define GET_T_GIVEN_EY get_t_given_ey
#define GET_T_GIVEN_HY get_t_given_hy
#define GET_REACTION_MAP get_reaction_map
#elif defined(BL_FORT_USE_UNDERSCORE)
#define CKINDX ckindx_
#define CKINIT ckinit_
#define CKFINALIZE ckfinalize_
#define CKXNUM ckxnum_
#define CKSYME cksyme_
#define CKSYMS cksyms_
#define CKRP ckrp_
#define CKPX ckpx_
#define CKPY ckpy_
#define CKPC ckpc_
#define CKRHOX ckrhox_
#define CKRHOY ckrhoy_
#define CKRHOC ckrhoc_
#define CKWT ckwt_
#define CKAWT ckawt_
#define CKMMWY ckmmwy_
#define CKMMWX ckmmwx_
#define CKMMWC ckmmwc_
#define CKYTX ckytx_
#define CKYTCP ckytcp_
#define CKYTCR ckytcr_
#define CKXTY ckxty_
#define CKXTCP ckxtcp_
#define CKXTCR ckxtcr_
#define CKCTX ckctx_
#define CKCTY ckcty_
#define CKCPOR ckcpor_
#define CKHORT ckhort_
#define CKSOR cksor_
#define CKCVML ckcvml_
#define CKCPML ckcpml_
#define CKUML ckuml_
#define CKHML ckhml_
#define CKGML ckgml_
#define CKAML ckaml_
#define CKSML cksml_
#define CKCVMS ckcvms_
#define CKCPMS ckcpms_
#define CKUMS ckums_
#define CKHMS ckhms_
#define CKGMS ckgms_
#define CKAMS ckams_
#define CKSMS cksms_
#define CKCPBL ckcpbl_
#define CKCPBS ckcpbs_
#define CKCVBL ckcvbl_
#define CKCVBS ckcvbs_
#define CKHBML ckhbml_
#define CKHBMS ckhbms_
#define CKUBML ckubml_
#define CKUBMS ckubms_
#define CKSBML cksbml_
#define CKSBMS cksbms_
#define CKGBML ckgbml_
#define CKGBMS ckgbms_
#define CKABML ckabml_
#define CKABMS ckabms_
#define CKWC ckwc_
#define CKWYP ckwyp_
#define CKWXP ckwxp_
#define CKWYR ckwyr_
#define CKWXR ckwxr_
#define CKQC ckqc_
#define CKKFKR ckkfkr_
#define CKQYP ckqyp_
#define CKQXP ckqxp_
#define CKQYR ckqyr_
#define CKQXR ckqxr_
#define CKNU cknu_
#define CKNCF ckncf_
#define CKABE ckabe_
#define CKEQC ckeqc_
#define CKEQYP ckeqyp_
#define CKEQXP ckeqxp_
#define CKEQYR ckeqyr_
#define CKEQXR ckeqxr_
#define DWDOT dwdot_
#define VCKHMS vckhms_
#define VCKPY vckpy_
#define VCKWYR vckwyr_
#define VCKYTX vckytx_
#define GET_T_GIVEN_EY get_t_given_ey_
#define GET_T_GIVEN_HY get_t_given_hy_
#define GET_REACTION_MAP get_reaction_map_
#endif

/*function declarations */
void atomicWeight(double * restrict awt);
void molecularWeight(double * restrict wt);
void gibbs(double * restrict species, double * restrict tc);
void helmholtz(double * restrict species, double * restrict tc);
void speciesInternalEnergy(double * restrict species, double * restrict tc);
void speciesEnthalpy(double * restrict species, double * restrict tc);
void speciesEntropy(double * restrict species, double * restrict tc);
void cp_R(double * restrict species, double * restrict tc);
void cv_R(double * restrict species, double * restrict tc);
void equilibriumConstants(double * restrict kc, double * restrict g_RT, double T);
void productionRate(double * restrict wdot, double * restrict sc, double T);
void comp_k_f(double * restrict tc, double invT, double * restrict k_f);
void comp_Kc(double * restrict tc, double invT, double * restrict Kc);
void comp_qfqr(double * restrict q_f, double * restrict q_r, double * restrict sc, double * restrict tc, double invT);
void progressRate(double * restrict qdot, double * restrict speciesConc, double T);
void progressRateFR(double * restrict q_f, double * restrict q_r, double * restrict speciesConc, double T);
void CKINIT();
void CKFINALIZE();
void CKINDX(int * iwrk, double * restrict rwrk, int * mm, int * kk, int * ii, int * nfit );
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double * restrict rval, int * kerr, int lenline);
void CKSNUM(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double * restrict rval, int * kerr, int lenline, int lenkray);
void CKSYME(int * kname, int * lenkname);
void CKSYMS(int * kname, int * lenkname);
void CKRP(int * ickwrk, double * restrict rckwrk, double * restrict ru, double * restrict ruc, double * restrict pa);
void CKPX(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict P);
void CKPY(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P);
void CKPC(double * restrict rho, double * restrict T, double * restrict c, int * iwrk, double * restrict rwrk, double * restrict P);
void CKRHOX(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict rho);
void CKRHOY(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict rho);
void CKRHOC(double * restrict P, double * restrict T, double * restrict c, int * iwrk, double * restrict rwrk, double * restrict rho);
void CKWT(int * iwrk, double * restrict rwrk, double * restrict wt);
void CKAWT(int * iwrk, double * restrict rwrk, double * restrict awt);
void CKMMWY(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wtm);
void CKMMWX(double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wtm);
void CKMMWC(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict wtm);
void CKYTX(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict x);
void CKYTCP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c);
void CKYTCR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c);
void CKXTY(double * restrict x, int * iwrk, double * restrict rwrk, double * restrict y);
void CKXTCP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c);
void CKXTCR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c);
void CKCTX(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict x);
void CKCTY(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict y);
void CKCPOR(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cpor);
void CKHORT(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hort);
void CKSOR(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sor);
void CKCVML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvml);
void CKCPML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvml);
void CKUML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict uml);
void CKHML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict uml);
void CKGML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict gml);
void CKAML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict aml);
void CKSML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sml);
void CKCVMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvms);
void CKCPMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvms);
void CKUMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums);
void CKHMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums);
void CKGMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict gms);
void CKAMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ams);
void CKSMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sms);
void CKCPBL(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cpbl);
void CKCPBS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cpbs);
void CKCVBL(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cpbl);
void CKCVBS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cpbs);
void CKHBML(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict hbml);
void CKHBMS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict hbms);
void CKUBML(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict ubml);
void CKUBMS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict ubms);
void CKSBML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict sbml);
void CKSBMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict sbms);
void CKGBML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict gbml);
void CKGBMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict gbms);
void CKABML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict abml);
void CKABMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict abms);
void CKWC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict wdot);
void CKWYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot);
void CKWXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot);
void CKWYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot);
void CKWXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot);
void CKQC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict qdot);
void CKKFKR(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict q_f, double * restrict q_r);
void CKQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot);
void CKQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot);
void CKQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot);
void CKQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot);
void CKNU(int * kdim, int * iwrk, double * restrict rwrk, int * nuki);
void CKNCF(int * mdim, int * iwrk, double * restrict rwrk, int * ncf);
void CKABE(int * iwrk, double * restrict rwrk, double * restrict a, double * restrict b, double * restrict e );
void CKEQC(double * restrict T, double * restrict C , int * iwrk, double * restrict rwrk, double * restrict eqcon );
void CKEQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon);
void CKEQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon);
void CKEQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon);
void CKEQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon);
void DWDOT(double * restrict J, double * restrict sc, double * restrict T, int * consP);
void aJacobian(double * restrict J, double * restrict sc, double T, int consP);
void dcvpRdT(double * restrict species, double * restrict tc);
void GET_T_GIVEN_EY(double * restrict e, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict t, int *ierr);
void GET_T_GIVEN_HY(double * restrict h, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict t, int *ierr);
void GET_REACTION_MAP(int * restrict rmap);
/*vector version */
void vproductionRate(int npt, double * restrict wdot, double * restrict c, double * restrict T);
void VCKHMS(int * restrict np, double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums);
void VCKPY(int * restrict np, double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P);
void VCKWYR(int * restrict np, double * restrict rho, double * restrict T,
            double * restrict y, int * restrict iwrk, double * restrict rwrk,
            double * restrict wdot);
void VCKYTX(int * restrict np, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict x);
void vcomp_k_f(int npt, double * restrict k_f_s, double * restrict tc, double * restrict invT);
void vcomp_gibbs(int npt, double * restrict g_RT, double * restrict tc);
void vcomp_Kc(int npt, double * restrict Kc_s, double * restrict g_RT, double * restrict invT);
void vcomp_wdot(int npt, double * restrict wdot, double * restrict mixture, double * restrict sc,
                double * restrict k_f_s, double * restrict Kc_s,
                double * restrict tc, double * restrict invT, double * restrict T);

/* Inverse molecular weights */
static const double imw[9] = {
    1.0 / 2.015940,  /*H2 */
    1.0 / 31.998800,  /*O2 */
    1.0 / 18.015340,  /*H2O */
    1.0 / 1.007970,  /*H */
    1.0 / 15.999400,  /*O */
    1.0 / 17.007370,  /*OH */
    1.0 / 33.006770,  /*HO2 */
    1.0 / 34.014740,  /*H2O2 */
    1.0 / 28.013400};  /*N2 */



static double fwd_A[21], fwd_beta[21], fwd_Ea[21];
static double low_A[21], low_beta[21], low_Ea[21];
static double rev_A[21], rev_beta[21], rev_Ea[21];
static double troe_a[21],troe_Ts[21], troe_Tss[21], troe_Tsss[21];
static double sri_a[21], sri_b[21], sri_c[21], sri_d[21], sri_e[21];
static double activation_units[21], prefactor_units[21], phase_units[21];
static int is_PD[21], troe_len[21], sri_len[21], nTB[21], *TBid[21];
static double *TB[21];

static double fwd_A_DEF[21], fwd_beta_DEF[21], fwd_Ea_DEF[21];
static double low_A_DEF[21], low_beta_DEF[21], low_Ea_DEF[21];
static double rev_A_DEF[21], rev_beta_DEF[21], rev_Ea_DEF[21];
static double troe_a_DEF[21],troe_Ts_DEF[21], troe_Tss_DEF[21], troe_Tsss_DEF[21];
static double sri_a_DEF[21], sri_b_DEF[21], sri_c_DEF[21], sri_d_DEF[21], sri_e_DEF[21];
static double activation_units_DEF[21], prefactor_units_DEF[21], phase_units_DEF[21];
static int is_PD_DEF[21], troe_len_DEF[21], sri_len_DEF[21], nTB_DEF[21], *TBid_DEF[21];
static double *TB_DEF[21];
static int rxn_map[21] = {6,7,8,9,2,3,4,5,0,10,11,12,13,14,15,1,16,17,18,19,20};

void GET_REACTION_MAP(int *rmap)
{
    for (int i=0; i<21; ++i) {
        rmap[i] = rxn_map[i];
    }
}


#include <ReactionData.H>
double* GetParamPtr(int                reaction_id,
                    REACTION_PARAMETER param_id,
                    int                species_id,
                    int                get_default)
{
  double* ret = 0;
  if (reaction_id<0 || reaction_id>=21) {
    printf("Bad reaction id = %d",reaction_id);
    abort();
  };
  int mrid = rxn_map[reaction_id];

  if (param_id == THIRD_BODY) {
    if (species_id<0 || species_id>=9) {
      printf("GetParamPtr: Bad species id = %d",species_id);
      abort();
    }
    if (get_default) {
      for (int i=0; i<nTB_DEF[mrid]; ++i) {
        if (species_id == TBid_DEF[mrid][i]) {
          ret = &(TB_DEF[mrid][i]);
        }
      }
    }
    else {
      for (int i=0; i<nTB[mrid]; ++i) {
        if (species_id == TBid[mrid][i]) {
          ret = &(TB[mrid][i]);
        }
      }
    }
    if (ret == 0) {
      printf("GetParamPtr: No TB for reaction id = %d",reaction_id);
      abort();
    }
  }
  else {
    if (     param_id == FWD_A)     {ret = (get_default ? &(fwd_A_DEF[mrid]) : &(fwd_A[mrid]));}
      else if (param_id == FWD_BETA)  {ret = (get_default ? &(fwd_beta_DEF[mrid]) : &(fwd_beta[mrid]));}
      else if (param_id == FWD_EA)    {ret = (get_default ? &(fwd_Ea_DEF[mrid]) : &(fwd_Ea[mrid]));}
      else if (param_id == LOW_A)     {ret = (get_default ? &(low_A_DEF[mrid]) : &(low_A[mrid]));}
      else if (param_id == LOW_BETA)  {ret = (get_default ? &(low_beta_DEF[mrid]) : &(low_beta[mrid]));}
      else if (param_id == LOW_EA)    {ret = (get_default ? &(low_Ea_DEF[mrid]) : &(low_Ea[mrid]));}
      else if (param_id == REV_A)     {ret = (get_default ? &(rev_A_DEF[mrid]) : &(rev_A[mrid]));}
      else if (param_id == REV_BETA)  {ret = (get_default ? &(rev_beta_DEF[mrid]) : &(rev_beta[mrid]));}
      else if (param_id == REV_EA)    {ret = (get_default ? &(rev_Ea_DEF[mrid]) : &(rev_Ea[mrid]));}
      else if (param_id == TROE_A)    {ret = (get_default ? &(troe_a_DEF[mrid]) : &(troe_a[mrid]));}
      else if (param_id == TROE_TS)   {ret = (get_default ? &(troe_Ts_DEF[mrid]) : &(troe_Ts[mrid]));}
      else if (param_id == TROE_TSS)  {ret = (get_default ? &(troe_Tss_DEF[mrid]) : &(troe_Tss[mrid]));}
      else if (param_id == TROE_TSSS) {ret = (get_default ? &(troe_Tsss_DEF[mrid]) : &(troe_Tsss[mrid]));}
      else if (param_id == SRI_A)     {ret = (get_default ? &(sri_a_DEF[mrid]) : &(sri_a[mrid]));}
      else if (param_id == SRI_B)     {ret = (get_default ? &(sri_b_DEF[mrid]) : &(sri_b[mrid]));}
      else if (param_id == SRI_C)     {ret = (get_default ? &(sri_c_DEF[mrid]) : &(sri_c[mrid]));}
      else if (param_id == SRI_D)     {ret = (get_default ? &(sri_d_DEF[mrid]) : &(sri_d[mrid]));}
      else if (param_id == SRI_E)     {ret = (get_default ? &(sri_e_DEF[mrid]) : &(sri_e[mrid]));}
    else {
      printf("GetParamPtr: Unknown parameter id");
      abort();
    }
  }
  return ret;
}

void ResetAllParametersToDefault()
{
    for (int i=0; i<21; i++) {
        if (nTB[i] != 0) {
            nTB[i] = 0;
            free(TB[i]);
            free(TBid[i]);
        }

        fwd_A[i]    = fwd_A_DEF[i];
        fwd_beta[i] = fwd_beta_DEF[i];
        fwd_Ea[i]   = fwd_Ea_DEF[i];

        low_A[i]    = low_A_DEF[i];
        low_beta[i] = low_beta_DEF[i];
        low_Ea[i]   = low_Ea_DEF[i];

        rev_A[i]    = rev_A_DEF[i];
        rev_beta[i] = rev_beta_DEF[i];
        rev_Ea[i]   = rev_Ea_DEF[i];

        troe_a[i]    = troe_a_DEF[i];
        troe_Ts[i]   = troe_Ts_DEF[i];
        troe_Tss[i]  = troe_Tss_DEF[i];
        troe_Tsss[i] = troe_Tsss_DEF[i];

        sri_a[i] = sri_a_DEF[i];
        sri_b[i] = sri_b_DEF[i];
        sri_c[i] = sri_c_DEF[i];
        sri_d[i] = sri_d_DEF[i];
        sri_e[i] = sri_e_DEF[i];

        is_PD[i]    = is_PD_DEF[i];
        troe_len[i] = troe_len_DEF[i];
        sri_len[i]  = sri_len_DEF[i];

        activation_units[i] = activation_units_DEF[i];
        prefactor_units[i]  = prefactor_units_DEF[i];
        phase_units[i]      = phase_units_DEF[i];

        nTB[i]  = nTB_DEF[i];
        if (nTB[i] != 0) {
           TB[i] = (double *) malloc(sizeof(double) * nTB[i]);
           TBid[i] = (int *) malloc(sizeof(int) * nTB[i]);
           for (int j=0; j<nTB[i]; j++) {
             TB[i][j] = TB_DEF[i][j];
             TBid[i][j] = TBid_DEF[i][j];
           }
        }
    }
}

void SetAllDefaults()
{
    for (int i=0; i<21; i++) {
        if (nTB_DEF[i] != 0) {
            nTB_DEF[i] = 0;
            free(TB_DEF[i]);
            free(TBid_DEF[i]);
        }

        fwd_A_DEF[i]    = fwd_A[i];
        fwd_beta_DEF[i] = fwd_beta[i];
        fwd_Ea_DEF[i]   = fwd_Ea[i];

        low_A_DEF[i]    = low_A[i];
        low_beta_DEF[i] = low_beta[i];
        low_Ea_DEF[i]   = low_Ea[i];

        rev_A_DEF[i]    = rev_A[i];
        rev_beta_DEF[i] = rev_beta[i];
        rev_Ea_DEF[i]   = rev_Ea[i];

        troe_a_DEF[i]    = troe_a[i];
        troe_Ts_DEF[i]   = troe_Ts[i];
        troe_Tss_DEF[i]  = troe_Tss[i];
        troe_Tsss_DEF[i] = troe_Tsss[i];

        sri_a_DEF[i] = sri_a[i];
        sri_b_DEF[i] = sri_b[i];
        sri_c_DEF[i] = sri_c[i];
        sri_d_DEF[i] = sri_d[i];
        sri_e_DEF[i] = sri_e[i];

        is_PD_DEF[i]    = is_PD[i];
        troe_len_DEF[i] = troe_len[i];
        sri_len_DEF[i]  = sri_len[i];

        activation_units_DEF[i] = activation_units[i];
        prefactor_units_DEF[i]  = prefactor_units[i];
        phase_units_DEF[i]      = phase_units[i];

        nTB_DEF[i]  = nTB[i];
        if (nTB_DEF[i] != 0) {
           TB_DEF[i] = (double *) malloc(sizeof(double) * nTB_DEF[i]);
           TBid_DEF[i] = (int *) malloc(sizeof(int) * nTB_DEF[i]);
           for (int j=0; j<nTB_DEF[i]; j++) {
             TB_DEF[i][j] = TB[i][j];
             TBid_DEF[i][j] = TBid[i][j];
           }
        }
    }
}

/* Finalizes parameter database */
void CKFINALIZE()
{
  for (int i=0; i<21; ++i) {
    free(TB[i]); TB[i] = 0; 
    free(TBid[i]); TBid[i] = 0;
    nTB[i] = 0;

    free(TB_DEF[i]); TB_DEF[i] = 0; 
    free(TBid_DEF[i]); TBid_DEF[i] = 0;
    nTB_DEF[i] = 0;
  }
}

/* Initializes parameter database */
void CKINIT()
{
    // (0):  H + O2 <=> O + OH
    fwd_A[6]     = 3547000000000000;
    fwd_beta[6]  = -0.40600000000000003;
    fwd_Ea[6]    = 16599;
    prefactor_units[6]  = 1.0000000000000002e-06;
    activation_units[6] = 0.50321666580471969;
    phase_units[6]      = 1e-12;
    is_PD[6] = 0;
    nTB[6] = 0;

    // (1):  O + H2 <=> H + OH
    fwd_A[7]     = 50800;
    fwd_beta[7]  = 2.6699999999999999;
    fwd_Ea[7]    = 6290;
    prefactor_units[7]  = 1.0000000000000002e-06;
    activation_units[7] = 0.50321666580471969;
    phase_units[7]      = 1e-12;
    is_PD[7] = 0;
    nTB[7] = 0;

    // (2):  H2 + OH <=> H2O + H
    fwd_A[8]     = 216000000;
    fwd_beta[8]  = 1.51;
    fwd_Ea[8]    = 3430;
    prefactor_units[8]  = 1.0000000000000002e-06;
    activation_units[8] = 0.50321666580471969;
    phase_units[8]      = 1e-12;
    is_PD[8] = 0;
    nTB[8] = 0;

    // (3):  O + H2O <=> OH + OH
    fwd_A[9]     = 2970000;
    fwd_beta[9]  = 2.02;
    fwd_Ea[9]    = 13400;
    prefactor_units[9]  = 1.0000000000000002e-06;
    activation_units[9] = 0.50321666580471969;
    phase_units[9]      = 1e-12;
    is_PD[9] = 0;
    nTB[9] = 0;

    // (4):  H2 + M <=> H + H + M
    fwd_A[2]     = 4.577e+19;
    fwd_beta[2]  = -1.3999999999999999;
    fwd_Ea[2]    = 104380;
    prefactor_units[2]  = 1.0000000000000002e-06;
    activation_units[2] = 0.50321666580471969;
    phase_units[2]      = 1e-6;
    is_PD[2] = 0;
    nTB[2] = 2;
    TB[2] = (double *) malloc(2 * sizeof(double));
    TBid[2] = (int *) malloc(2 * sizeof(int));
    TBid[2][0] = 0; TB[2][0] = 2.5; // H2
    TBid[2][1] = 2; TB[2][1] = 12; // H2O

    // (5):  O + O + M <=> O2 + M
    fwd_A[3]     = 6165000000000000;
    fwd_beta[3]  = -0.5;
    fwd_Ea[3]    = 0;
    prefactor_units[3]  = 1.0000000000000002e-12;
    activation_units[3] = 0.50321666580471969;
    phase_units[3]      = 1e-12;
    is_PD[3] = 0;
    nTB[3] = 2;
    TB[3] = (double *) malloc(2 * sizeof(double));
    TBid[3] = (int *) malloc(2 * sizeof(int));
    TBid[3][0] = 0; TB[3][0] = 2.5; // H2
    TBid[3][1] = 2; TB[3][1] = 12; // H2O

    // (6):  O + H + M <=> OH + M
    fwd_A[4]     = 4.714e+18;
    fwd_beta[4]  = -1;
    fwd_Ea[4]    = 0;
    prefactor_units[4]  = 1.0000000000000002e-12;
    activation_units[4] = 0.50321666580471969;
    phase_units[4]      = 1e-12;
    is_PD[4] = 0;
    nTB[4] = 2;
    TB[4] = (double *) malloc(2 * sizeof(double));
    TBid[4] = (int *) malloc(2 * sizeof(int));
    TBid[4][0] = 0; TB[4][0] = 2.5; // H2
    TBid[4][1] = 2; TB[4][1] = 12; // H2O

    // (7):  H + OH + M <=> H2O + M
    fwd_A[5]     = 3.8000000000000004e+22;
    fwd_beta[5]  = -2;
    fwd_Ea[5]    = 0;
    prefactor_units[5]  = 1.0000000000000002e-12;
    activation_units[5] = 0.50321666580471969;
    phase_units[5]      = 1e-12;
    is_PD[5] = 0;
    nTB[5] = 2;
    TB[5] = (double *) malloc(2 * sizeof(double));
    TBid[5] = (int *) malloc(2 * sizeof(int));
    TBid[5][0] = 0; TB[5][0] = 2.5; // H2
    TBid[5][1] = 2; TB[5][1] = 12; // H2O

    // (8):  H + O2 (+M) <=> HO2 (+M)
    fwd_A[0]     = 1475000000000;
    fwd_beta[0]  = 0.59999999999999998;
    fwd_Ea[0]    = 0;
    low_A[0]     = 6.366e+20;
    low_beta[0]  = -1.72;
    low_Ea[0]    = 524.79999999999995;
    troe_a[0]    = 0.80000000000000004;
    troe_Tsss[0] = 1.0000000000000001e-30;
    troe_Ts[0]   = 1e+30;
    troe_len[0]  = 3;
    prefactor_units[0]  = 1.0000000000000002e-06;
    activation_units[0] = 0.50321666580471969;
    phase_units[0]      = 1e-12;
    is_PD[0] = 1;
    nTB[0] = 3;
    TB[0] = (double *) malloc(3 * sizeof(double));
    TBid[0] = (int *) malloc(3 * sizeof(int));
    TBid[0][0] = 0; TB[0][0] = 2; // H2
    TBid[0][1] = 2; TB[0][1] = 11; // H2O
    TBid[0][2] = 1; TB[0][2] = 0.78000000000000003; // O2

    // (9):  HO2 + H <=> H2 + O2
    fwd_A[10]     = 16600000000000;
    fwd_beta[10]  = 0;
    fwd_Ea[10]    = 823;
    prefactor_units[10]  = 1.0000000000000002e-06;
    activation_units[10] = 0.50321666580471969;
    phase_units[10]      = 1e-12;
    is_PD[10] = 0;
    nTB[10] = 0;

    // (10):  HO2 + H <=> OH + OH
    fwd_A[11]     = 70790000000000;
    fwd_beta[11]  = 0;
    fwd_Ea[11]    = 295;
    prefactor_units[11]  = 1.0000000000000002e-06;
    activation_units[11] = 0.50321666580471969;
    phase_units[11]      = 1e-12;
    is_PD[11] = 0;
    nTB[11] = 0;

    // (11):  HO2 + O <=> O2 + OH
    fwd_A[12]     = 32500000000000;
    fwd_beta[12]  = 0;
    fwd_Ea[12]    = 0;
    prefactor_units[12]  = 1.0000000000000002e-06;
    activation_units[12] = 0.50321666580471969;
    phase_units[12]      = 1e-12;
    is_PD[12] = 0;
    nTB[12] = 0;

    // (12):  HO2 + OH <=> H2O + O2
    fwd_A[13]     = 28900000000000;
    fwd_beta[13]  = 0;
    fwd_Ea[13]    = -497;
    prefactor_units[13]  = 1.0000000000000002e-06;
    activation_units[13] = 0.50321666580471969;
    phase_units[13]      = 1e-12;
    is_PD[13] = 0;
    nTB[13] = 0;

    // (13):  HO2 + HO2 <=> H2O2 + O2
    fwd_A[14]     = 420000000000000;
    fwd_beta[14]  = 0;
    fwd_Ea[14]    = 11982;
    prefactor_units[14]  = 1.0000000000000002e-06;
    activation_units[14] = 0.50321666580471969;
    phase_units[14]      = 1e-12;
    is_PD[14] = 0;
    nTB[14] = 0;

    // (14):  HO2 + HO2 <=> H2O2 + O2
    fwd_A[15]     = 130000000000;
    fwd_beta[15]  = 0;
    fwd_Ea[15]    = -1629.3;
    prefactor_units[15]  = 1.0000000000000002e-06;
    activation_units[15] = 0.50321666580471969;
    phase_units[15]      = 1e-12;
    is_PD[15] = 0;
    nTB[15] = 0;

    // (15):  H2O2 (+M) <=> OH + OH (+M)
    fwd_A[1]     = 295100000000000;
    fwd_beta[1]  = 0;
    fwd_Ea[1]    = 48430;
    low_A[1]     = 1.202e+17;
    low_beta[1]  = 0;
    low_Ea[1]    = 45500;
    troe_a[1]    = 0.5;
    troe_Tsss[1] = 1.0000000000000001e-30;
    troe_Ts[1]   = 1e+30;
    troe_len[1]  = 3;
    prefactor_units[1]  = 1;
    activation_units[1] = 0.50321666580471969;
    phase_units[1]      = 1e-6;
    is_PD[1] = 1;
    nTB[1] = 2;
    TB[1] = (double *) malloc(2 * sizeof(double));
    TBid[1] = (int *) malloc(2 * sizeof(int));
    TBid[1][0] = 0; TB[1][0] = 2.5; // H2
    TBid[1][1] = 2; TB[1][1] = 12; // H2O

    // (16):  H2O2 + H <=> H2O + OH
    fwd_A[16]     = 24100000000000;
    fwd_beta[16]  = 0;
    fwd_Ea[16]    = 3970;
    prefactor_units[16]  = 1.0000000000000002e-06;
    activation_units[16] = 0.50321666580471969;
    phase_units[16]      = 1e-12;
    is_PD[16] = 0;
    nTB[16] = 0;

    // (17):  H2O2 + H <=> HO2 + H2
    fwd_A[17]     = 48200000000000;
    fwd_beta[17]  = 0;
    fwd_Ea[17]    = 7950;
    prefactor_units[17]  = 1.0000000000000002e-06;
    activation_units[17] = 0.50321666580471969;
    phase_units[17]      = 1e-12;
    is_PD[17] = 0;
    nTB[17] = 0;

    // (18):  H2O2 + O <=> OH + HO2
    fwd_A[18]     = 9550000;
    fwd_beta[18]  = 2;
    fwd_Ea[18]    = 3970;
    prefactor_units[18]  = 1.0000000000000002e-06;
    activation_units[18] = 0.50321666580471969;
    phase_units[18]      = 1e-12;
    is_PD[18] = 0;
    nTB[18] = 0;

    // (19):  H2O2 + OH <=> HO2 + H2O
    fwd_A[19]     = 1000000000000;
    fwd_beta[19]  = 0;
    fwd_Ea[19]    = 0;
    prefactor_units[19]  = 1.0000000000000002e-06;
    activation_units[19] = 0.50321666580471969;
    phase_units[19]      = 1e-12;
    is_PD[19] = 0;
    nTB[19] = 0;

    // (20):  H2O2 + OH <=> HO2 + H2O
    fwd_A[20]     = 580000000000000;
    fwd_beta[20]  = 0;
    fwd_Ea[20]    = 9557;
    prefactor_units[20]  = 1.0000000000000002e-06;
    activation_units[20] = 0.50321666580471969;
    phase_units[20]      = 1e-12;
    is_PD[20] = 0;
    nTB[20] = 0;

    SetAllDefaults();
}



/*A few mechanism parameters */
void CKINDX(int * iwrk, double * restrict rwrk, int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 3;
    *kk = 9;
    *ii = 21;
    *nfit = -1; /*Why do you need this anyway ?  */
}



/* ckxnum... for parsing strings  */
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double * restrict rval, int * kerr, int lenline )
{
    int n,i; /*Loop Counters */
    char cstr[1000];
    char *saveptr;
    char *p; /*String Tokens */
    /* Strip Comments  */
    for (i=0; i<lenline; ++i) {
        if (line[i]=='!') {
            break;
        }
        cstr[i] = line[i];
    }
    cstr[i] = '\0';

    p = strtok_r(cstr," ", &saveptr);
    if (!p) {
        *nval = 0;
        *kerr = 1;
        return;
    }
    for (n=0; n<*nexp; ++n) {
        rval[n] = atof(p);
        p = strtok_r(NULL, " ", &saveptr);
        if (!p) break;
    }
    *nval = n+1;
    if (*nval < *nexp) *kerr = 1;
    return;
}


/* cksnum... for parsing strings  */
void CKSNUM(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double * restrict rval, int * kerr, int lenline, int lenkray)
{
    /*Not done yet ... */
}


/* Returns the char strings of element names */
void CKSYME(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*3; i++) {
        kname[i] = ' ';
    }

    /* H  */
    kname[ 0*lenkname + 0 ] = 'H';
    kname[ 0*lenkname + 1 ] = ' ';

    /* O  */
    kname[ 1*lenkname + 0 ] = 'O';
    kname[ 1*lenkname + 1 ] = ' ';

    /* N  */
    kname[ 2*lenkname + 0 ] = 'N';
    kname[ 2*lenkname + 1 ] = ' ';

}


/* Returns the char strings of species names */
void CKSYMS(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*9; i++) {
        kname[i] = ' ';
    }

    /* H2  */
    kname[ 0*lenkname + 0 ] = 'H';
    kname[ 0*lenkname + 1 ] = '2';
    kname[ 0*lenkname + 2 ] = ' ';

    /* O2  */
    kname[ 1*lenkname + 0 ] = 'O';
    kname[ 1*lenkname + 1 ] = '2';
    kname[ 1*lenkname + 2 ] = ' ';

    /* H2O  */
    kname[ 2*lenkname + 0 ] = 'H';
    kname[ 2*lenkname + 1 ] = '2';
    kname[ 2*lenkname + 2 ] = 'O';
    kname[ 2*lenkname + 3 ] = ' ';

    /* H  */
    kname[ 3*lenkname + 0 ] = 'H';
    kname[ 3*lenkname + 1 ] = ' ';

    /* O  */
    kname[ 4*lenkname + 0 ] = 'O';
    kname[ 4*lenkname + 1 ] = ' ';

    /* OH  */
    kname[ 5*lenkname + 0 ] = 'O';
    kname[ 5*lenkname + 1 ] = 'H';
    kname[ 5*lenkname + 2 ] = ' ';

    /* HO2  */
    kname[ 6*lenkname + 0 ] = 'H';
    kname[ 6*lenkname + 1 ] = 'O';
    kname[ 6*lenkname + 2 ] = '2';
    kname[ 6*lenkname + 3 ] = ' ';

    /* H2O2  */
    kname[ 7*lenkname + 0 ] = 'H';
    kname[ 7*lenkname + 1 ] = '2';
    kname[ 7*lenkname + 2 ] = 'O';
    kname[ 7*lenkname + 3 ] = '2';
    kname[ 7*lenkname + 4 ] = ' ';

    /* N2  */
    kname[ 8*lenkname + 0 ] = 'N';
    kname[ 8*lenkname + 1 ] = '2';
    kname[ 8*lenkname + 2 ] = ' ';

}


/* Returns R, Rc, Patm */
void CKRP(int * ickwrk, double * restrict rckwrk, double * restrict ru, double * restrict ruc, double * restrict pa)
{
     *ru  = 8.31451e+07; 
     *ruc = 1.98721558317399615845; 
     *pa  = 1.01325e+06; 
}


/*Compute P = rhoRT/W(x) */
void CKPX(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict P)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    *P = *rho * 8.31451e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void CKPY(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*H */
    YOW += y[4]*imw[4]; /*O */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    *P = *rho * 8.31451e+07 * (*T) * YOW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void VCKPY(int * restrict np, double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P)
{
    double YOW[*np];
    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<9; n++) {
        for (int i=0; i<(*np); i++) {
            YOW[i] += y[n*(*np)+i] * imw[n];
        }
    }

    for (int i=0; i<(*np); i++) {
        P[i] = rho[i] * 8.31451e+07 * T[i] * YOW[i]; /*P = rho*R*T/W */
    }

    return;
}


/*Compute P = rhoRT/W(c) */
void CKPC(double * restrict rho, double * restrict T, double * restrict c, int * iwrk, double * restrict rwrk, double * restrict P)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*2.015940; /*H2 */
    W += c[1]*31.998800; /*O2 */
    W += c[2]*18.015340; /*H2O */
    W += c[3]*1.007970; /*H */
    W += c[4]*15.999400; /*O */
    W += c[5]*17.007370; /*OH */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.013400; /*N2 */

    for (id = 0; id < 9; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.31451e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
void CKRHOX(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict rho)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    *rho = *P * XW / (8.31451e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
void CKRHOY(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict rho)
{
    double YOW = 0;
    double tmp[9];

    for (int i = 0; i < 9; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 9; i++)
    {
        YOW += tmp[i];
    }

    *rho = *P / (8.31451e+07 * (*T) * YOW);/*rho = P*W/(R*T) */
    return;
}


/*Compute rho = P*W(c)/(R*T) */
void CKRHOC(double * restrict P, double * restrict T, double * restrict c, int * iwrk, double * restrict rwrk, double * restrict rho)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*2.015940; /*H2 */
    W += c[1]*31.998800; /*O2 */
    W += c[2]*18.015340; /*H2O */
    W += c[3]*1.007970; /*H */
    W += c[4]*15.999400; /*O */
    W += c[5]*17.007370; /*OH */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.013400; /*N2 */

    for (id = 0; id < 9; ++id) {
        sumC += c[id];
    }
    *rho = *P * W / (sumC * (*T) * 8.31451e+07); /*rho = PW/(R*T) */

    return;
}


/*get molecular weight for all species */
void CKWT(int * iwrk, double * restrict rwrk, double * restrict wt)
{
    molecularWeight(wt);
}


/*get atomic weight for all elements */
void CKAWT(int * iwrk, double * restrict rwrk, double * restrict awt)
{
    atomicWeight(awt);
}


/*given y[species]: mass fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWY(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wtm)
{
    double YOW = 0;
    double tmp[9];

    for (int i = 0; i < 9; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 9; i++)
    {
        YOW += tmp[i];
    }

    *wtm = 1.0 / YOW;
    return;
}


/*given x[species]: mole fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWX(double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wtm)
{
    double XW = 0;/* see Eq 4 in CK Manual */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    *wtm = XW;

    return;
}


/*given c[species]: molar concentration */
/*returns mean molecular weight (gm/mole) */
void CKMMWC(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict wtm)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*2.015940; /*H2 */
    W += c[1]*31.998800; /*O2 */
    W += c[2]*18.015340; /*H2O */
    W += c[3]*1.007970; /*H */
    W += c[4]*15.999400; /*O */
    W += c[5]*17.007370; /*OH */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.013400; /*N2 */

    for (id = 0; id < 9; ++id) {
        sumC += c[id];
    }
    /* CK provides no guard against divison by zero */
    *wtm = W/sumC;

    return;
}


/*convert y[species] (mass fracs) to x[species] (mole fracs) */
void CKYTX(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict x)
{
    double YOW = 0;
    double tmp[9];

    for (int i = 0; i < 9; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 9; i++)
    {
        YOW += tmp[i];
    }

    double YOWINV = 1.0/YOW;

    for (int i = 0; i < 9; i++)
    {
        x[i] = y[i]*imw[i]*YOWINV;
    }
    return;
}


/*convert y[npoints*species] (mass fracs) to x[npoints*species] (mole fracs) */
void VCKYTX(int * restrict np, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict x)
{
    double YOW[*np];
    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<9; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] = y[n*(*np)+i] * imw[n];
            YOW[i] += x[n*(*np)+i];
        }
    }

    for (int i=0; i<(*np); i++) {
        YOW[i] = 1.0/YOW[i];
    }

    for (int n=0; n<9; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] *=  YOW[i];
        }
    }
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c)
{
    double YOW = 0;
    double PWORT;

    /*Compute inverse of mean molecular wt first */
    for (int i = 0; i < 9; i++)
    {
        c[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 9; i++)
    {
        YOW += c[i];
    }

    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*Now compute conversion */

    for (int i = 0; i < 9; i++)
    {
        c[i] = PWORT * y[i] * imw[i];
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c)
{
    for (int i = 0; i < 9; i++)
    {
        c[i] = (*rho)  * y[i] * imw[i];
    }
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
void CKXTY(double * restrict x, int * iwrk, double * restrict rwrk, double * restrict y)
{
    double XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    /*Now compute conversion */
    double XWinv = 1.0/XW;
    y[0] = x[0]*2.015940*XWinv; 
    y[1] = x[1]*31.998800*XWinv; 
    y[2] = x[2]*18.015340*XWinv; 
    y[3] = x[3]*1.007970*XWinv; 
    y[4] = x[4]*15.999400*XWinv; 
    y[5] = x[5]*17.007370*XWinv; 
    y[6] = x[6]*33.006770*XWinv; 
    y[7] = x[7]*34.014740*XWinv; 
    y[8] = x[8]*28.013400*XWinv; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.31451e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 9; ++id) {
        c[id] = x[id]*PORT;
    }

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c)
{
    int id; /*loop counter */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 9; ++id) {
        c[id] = x[id]*ROW;
    }

    return;
}


/*convert c[species] (molar conc) to x[species] (mole fracs) */
void CKCTX(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict x)
{
    int id; /*loop counter */
    double sumC = 0; 

    /*compute sum of c  */
    for (id = 0; id < 9; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    double sumCinv = 1.0/sumC;
    for (id = 0; id < 9; ++id) {
        x[id] = c[id]*sumCinv;
    }

    return;
}


/*convert c[species] (molar conc) to y[species] (mass fracs) */
void CKCTY(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict y)
{
    double CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*2.015940; /*H2 */
    CW += c[1]*31.998800; /*O2 */
    CW += c[2]*18.015340; /*H2O */
    CW += c[3]*1.007970; /*H */
    CW += c[4]*15.999400; /*O */
    CW += c[5]*17.007370; /*OH */
    CW += c[6]*33.006770; /*HO2 */
    CW += c[7]*34.014740; /*H2O2 */
    CW += c[8]*28.013400; /*N2 */
    /*Now compute conversion */
    double CWinv = 1.0/CW;
    y[0] = c[0]*2.015940*CWinv; 
    y[1] = c[1]*31.998800*CWinv; 
    y[2] = c[2]*18.015340*CWinv; 
    y[3] = c[3]*1.007970*CWinv; 
    y[4] = c[4]*15.999400*CWinv; 
    y[5] = c[5]*17.007370*CWinv; 
    y[6] = c[6]*33.006770*CWinv; 
    y[7] = c[7]*34.014740*CWinv; 
    y[8] = c[8]*28.013400*CWinv; 

    return;
}


/*get Cp/R as a function of T  */
/*for all species (Eq 19) */
void CKCPOR(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cpor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpor, tc);
}


/*get H/RT as a function of T  */
/*for all species (Eq 20) */
void CKHORT(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hort)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEnthalpy(hort, tc);
}


/*get S/R as a function of T  */
/*for all species (Eq 21) */
void CKSOR(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sor, tc);
}


/*get specific heat at constant volume as a function  */
/*of T for all species (molar units) */
void CKCVML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        cvml[id] *= 8.31451e+07;
    }
}


/*get specific heat at constant pressure as a  */
/*function of T for all species (molar units) */
void CKCPML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cpml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        cpml[id] *= 8.31451e+07;
    }
}


/*get internal energy as a function  */
/*of T for all species (molar units) */
void CKUML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict uml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        uml[id] *= RT;
    }
}


/*get enthalpy as a function  */
/*of T for all species (molar units) */
void CKHML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        hml[id] *= RT;
    }
}


/*get standard-state Gibbs energy as a function  */
/*of T for all species (molar units) */
void CKGML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict gml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    gibbs(gml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        gml[id] *= RT;
    }
}


/*get standard-state Helmholtz free energy as a  */
/*function of T for all species (molar units) */
void CKAML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict aml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    helmholtz(aml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        aml[id] *= RT;
    }
}


/*Returns the standard-state entropies in molar units */
void CKSML(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        sml[id] *= 8.31451e+07;
    }
}


/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
void CKCVMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cvms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvms, tc);
    /*multiply by R/molecularweight */
    cvms[0] *= 4.124383662212169e+07; /*H2 */
    cvms[1] *= 2.598381814318037e+06; /*O2 */
    cvms[2] *= 4.615239012974499e+06; /*H2O */
    cvms[3] *= 8.248767324424338e+07; /*H */
    cvms[4] *= 5.196763628636074e+06; /*O */
    cvms[5] *= 4.888768810227566e+06; /*OH */
    cvms[6] *= 2.519031701678171e+06; /*HO2 */
    cvms[7] *= 2.444384405113783e+06; /*H2O2 */
    cvms[8] *= 2.968047434442088e+06; /*N2 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
void CKCPMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cpms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 4.124383662212169e+07; /*H2 */
    cpms[1] *= 2.598381814318037e+06; /*O2 */
    cpms[2] *= 4.615239012974499e+06; /*H2O */
    cpms[3] *= 8.248767324424338e+07; /*H */
    cpms[4] *= 5.196763628636074e+06; /*O */
    cpms[5] *= 4.888768810227566e+06; /*OH */
    cpms[6] *= 2.519031701678171e+06; /*HO2 */
    cpms[7] *= 2.444384405113783e+06; /*H2O2 */
    cpms[8] *= 2.968047434442088e+06; /*N2 */
}


/*Returns internal energy in mass units (Eq 30.) */
void CKUMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    for (int i = 0; i < 9; i++)
    {
        ums[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
void CKHMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hms, tc);
    for (int i = 0; i < 9; i++)
    {
        hms[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int * restrict np, double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hms)
{
    double tc[5], h[9];

    for (int i=0; i<(*np); i++) {
        tc[0] = 0.0;
        tc[1] = T[i];
        tc[2] = T[i]*T[i];
        tc[3] = T[i]*T[i]*T[i];
        tc[4] = T[i]*T[i]*T[i]*T[i];

        speciesEnthalpy(h, tc);

        hms[0*(*np)+i] = h[0];
        hms[1*(*np)+i] = h[1];
        hms[2*(*np)+i] = h[2];
        hms[3*(*np)+i] = h[3];
        hms[4*(*np)+i] = h[4];
        hms[5*(*np)+i] = h[5];
        hms[6*(*np)+i] = h[6];
        hms[7*(*np)+i] = h[7];
        hms[8*(*np)+i] = h[8];
    }

    for (int n=0; n<9; n++) {
        for (int i=0; i<(*np); i++) {
            hms[n*(*np)+i] *= 8.31451e+07 * T[i] * imw[n];
        }
    }
}


/*Returns gibbs in mass units (Eq 31.) */
void CKGMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict gms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    gibbs(gms, tc);
    for (int i = 0; i < 9; i++)
    {
        gms[i] *= RT*imw[i];
    }
}


/*Returns helmholtz in mass units (Eq 32.) */
void CKAMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ams)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    helmholtz(ams, tc);
    for (int i = 0; i < 9; i++)
    {
        ams[i] *= RT*imw[i];
    }
}


/*Returns the entropies in mass units (Eq 28.) */
void CKSMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict sms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sms, tc);
    /*multiply by R/molecularweight */
    sms[0] *= 4.124383662212169e+07; /*H2 */
    sms[1] *= 2.598381814318037e+06; /*O2 */
    sms[2] *= 4.615239012974499e+06; /*H2O */
    sms[3] *= 8.248767324424338e+07; /*H */
    sms[4] *= 5.196763628636074e+06; /*O */
    sms[5] *= 4.888768810227566e+06; /*OH */
    sms[6] *= 2.519031701678171e+06; /*HO2 */
    sms[7] *= 2.444384405113783e+06; /*H2O2 */
    sms[8] *= 2.968047434442088e+06; /*N2 */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[9]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 9; ++id) {
        result += x[id]*cpor[id];
    }

    *cpbl = result * 8.31451e+07;
}


/*Returns the mean specific heat at CP (Eq. 34) */
void CKCPBS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cpbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[9], tresult[9]; /* temporary storage */
    cp_R(cpor, tc);
    for (int i = 0; i < 9; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 9; i++)
    {
        result += tresult[i];
    }

    *cpbs = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
void CKCVBL(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cvbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[9]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 9; ++id) {
        result += x[id]*cvor[id];
    }

    *cvbl = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 36) */
void CKCVBS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cvbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[9]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*H2 */
    result += cvor[1]*y[1]*imw[1]; /*O2 */
    result += cvor[2]*y[2]*imw[2]; /*H2O */
    result += cvor[3]*y[3]*imw[3]; /*H */
    result += cvor[4]*y[4]*imw[4]; /*O */
    result += cvor[5]*y[5]*imw[5]; /*OH */
    result += cvor[6]*y[6]*imw[6]; /*HO2 */
    result += cvor[7]*y[7]*imw[7]; /*H2O2 */
    result += cvor[8]*y[8]*imw[8]; /*N2 */

    *cvbs = result * 8.31451e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[9]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 9; ++id) {
        result += x[id]*hml[id];
    }

    *hbml = result * RT;
}


/*Returns mean enthalpy of mixture in mass units */
void CKHBMS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict hbms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[9], tmp[9]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);
    int id;
    for (id = 0; id < 9; ++id) {
        tmp[id] = y[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 9; ++id) {
        result += tmp[id];
    }

    *hbms = result * RT;
}


/*get mean internal energy in molar units */
void CKUBML(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict ubml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double uml[9]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 9; ++id) {
        result += x[id]*uml[id];
    }

    *ubml = result * RT;
}


/*get mean internal energy in mass units */
void CKUBMS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict ubms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double ums[9]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*H2 */
    result += y[1]*ums[1]*imw[1]; /*O2 */
    result += y[2]*ums[2]*imw[2]; /*H2O */
    result += y[3]*ums[3]*imw[3]; /*H */
    result += y[4]*ums[4]*imw[4]; /*O */
    result += y[5]*ums[5]*imw[5]; /*OH */
    result += y[6]*ums[6]*imw[6]; /*HO2 */
    result += y[7]*ums[7]*imw[7]; /*H2O2 */
    result += y[8]*ums[8]*imw[8]; /*N2 */

    *ubms = result * RT;
}


/*get mixture entropy in molar units */
void CKSBML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict sbml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[9]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 9; ++id) {
        result += x[id]*(sor[id]-log((x[id]+1e-100))-logPratio);
    }

    *sbml = result * 8.31451e+07;
}


/*get mixture entropy in mass units */
void CKSBMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict sbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[9]; /* temporary storage */
    double x[9]; /* need a ytx conversion */
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*H */
    YOW += y[4]*imw[4]; /*O */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(31.998800*YOW); 
    x[2] = y[2]/(18.015340*YOW); 
    x[3] = y[3]/(1.007970*YOW); 
    x[4] = y[4]/(15.999400*YOW); 
    x[5] = y[5]/(17.007370*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(28.013400*YOW); 
    speciesEntropy(sor, tc);
    /*Perform computation in Eq 42 and 43 */
    result += x[0]*(sor[0]-log((x[0]+1e-100))-logPratio);
    result += x[1]*(sor[1]-log((x[1]+1e-100))-logPratio);
    result += x[2]*(sor[2]-log((x[2]+1e-100))-logPratio);
    result += x[3]*(sor[3]-log((x[3]+1e-100))-logPratio);
    result += x[4]*(sor[4]-log((x[4]+1e-100))-logPratio);
    result += x[5]*(sor[5]-log((x[5]+1e-100))-logPratio);
    result += x[6]*(sor[6]-log((x[6]+1e-100))-logPratio);
    result += x[7]*(sor[7]-log((x[7]+1e-100))-logPratio);
    result += x[8]*(sor[8]-log((x[8]+1e-100))-logPratio);
    /*Scale by R/W */
    *sbms = result * 8.31451e+07 * YOW;
}


/*Returns mean gibbs free energy in molar units */
void CKGBML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict gbml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double gort[9]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 9; ++id) {
        result += x[id]*(gort[id]+log((x[id]+1e-100))+logPratio);
    }

    *gbml = result * RT;
}


/*Returns mixture gibbs free energy in mass units */
void CKGBMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict gbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double gort[9]; /* temporary storage */
    double x[9]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*H */
    YOW += y[4]*imw[4]; /*O */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(31.998800*YOW); 
    x[2] = y[2]/(18.015340*YOW); 
    x[3] = y[3]/(1.007970*YOW); 
    x[4] = y[4]/(15.999400*YOW); 
    x[5] = y[5]/(17.007370*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(28.013400*YOW); 
    gibbs(gort, tc);
    /*Perform computation in Eq 44 */
    result += x[0]*(gort[0]+log((x[0]+1e-100))+logPratio);
    result += x[1]*(gort[1]+log((x[1]+1e-100))+logPratio);
    result += x[2]*(gort[2]+log((x[2]+1e-100))+logPratio);
    result += x[3]*(gort[3]+log((x[3]+1e-100))+logPratio);
    result += x[4]*(gort[4]+log((x[4]+1e-100))+logPratio);
    result += x[5]*(gort[5]+log((x[5]+1e-100))+logPratio);
    result += x[6]*(gort[6]+log((x[6]+1e-100))+logPratio);
    result += x[7]*(gort[7]+log((x[7]+1e-100))+logPratio);
    result += x[8]*(gort[8]+log((x[8]+1e-100))+logPratio);
    /*Scale by RT/W */
    *gbms = result * RT * YOW;
}


/*Returns mean helmholtz free energy in molar units */
void CKABML(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict abml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double aort[9]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 9; ++id) {
        result += x[id]*(aort[id]+log((x[id]+1e-100))+logPratio);
    }

    *abml = result * RT;
}


/*Returns mixture helmholtz free energy in mass units */
void CKABMS(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict abms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double aort[9]; /* temporary storage */
    double x[9]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*H */
    YOW += y[4]*imw[4]; /*O */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(31.998800*YOW); 
    x[2] = y[2]/(18.015340*YOW); 
    x[3] = y[3]/(1.007970*YOW); 
    x[4] = y[4]/(15.999400*YOW); 
    x[5] = y[5]/(17.007370*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(28.013400*YOW); 
    helmholtz(aort, tc);
    /*Perform computation in Eq 44 */
    result += x[0]*(aort[0]+log((x[0]+1e-100))+logPratio);
    result += x[1]*(aort[1]+log((x[1]+1e-100))+logPratio);
    result += x[2]*(aort[2]+log((x[2]+1e-100))+logPratio);
    result += x[3]*(aort[3]+log((x[3]+1e-100))+logPratio);
    result += x[4]*(aort[4]+log((x[4]+1e-100))+logPratio);
    result += x[5]*(aort[5]+log((x[5]+1e-100))+logPratio);
    result += x[6]*(aort[6]+log((x[6]+1e-100))+logPratio);
    result += x[7]*(aort[7]+log((x[7]+1e-100))+logPratio);
    result += x[8]*(aort[8]+log((x[8]+1e-100))+logPratio);
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
void CKWC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 9; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*H */
    YOW += y[4]*imw[4]; /*O */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]*imw[0]; 
    c[1] = PWORT * y[1]*imw[1]; 
    c[2] = PWORT * y[2]*imw[2]; 
    c[3] = PWORT * y[3]*imw[3]; 
    c[4] = PWORT * y[4]*imw[4]; 
    c[5] = PWORT * y[5]*imw[5]; 
    c[6] = PWORT * y[6]*imw[6]; 
    c[7] = PWORT * y[7]*imw[7]; 
    c[8] = PWORT * y[8]*imw[8]; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 9; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void CKWYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]*imw[0]; 
    c[1] = 1e6 * (*rho) * y[1]*imw[1]; 
    c[2] = 1e6 * (*rho) * y[2]*imw[2]; 
    c[3] = 1e6 * (*rho) * y[3]*imw[3]; 
    c[4] = 1e6 * (*rho) * y[4]*imw[4]; 
    c[5] = 1e6 * (*rho) * y[5]*imw[5]; 
    c[6] = 1e6 * (*rho) * y[6]*imw[6]; 
    c[7] = 1e6 * (*rho) * y[7]*imw[7]; 
    c[8] = 1e6 * (*rho) * y[8]*imw[8]; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void VCKWYR(int * restrict np, double * restrict rho, double * restrict T,
	    double * restrict y, int * restrict iwrk, double * restrict rwrk,
	    double * restrict wdot)
{
    double c[9*(*np)]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    for (int n=0; n<9; n++) {
        for (int i=0; i<(*np); i++) {
            c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];
        }
    }

    /*call productionRate */
    vproductionRate(*np, wdot, c, T);

    /*convert to chemkin units */
    for (int i=0; i<9*(*np); i++) {
        wdot[i] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 9; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 9; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 21; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKKFKR(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict q_f, double * restrict q_r)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 9; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*H */
    YOW += y[4]*imw[4]; /*O */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]*imw[0]; 
    c[1] = PWORT * y[1]*imw[1]; 
    c[2] = PWORT * y[2]*imw[2]; 
    c[3] = PWORT * y[3]*imw[3]; 
    c[4] = PWORT * y[4]*imw[4]; 
    c[5] = PWORT * y[5]*imw[5]; 
    c[6] = PWORT * y[6]*imw[6]; 
    c[7] = PWORT * y[7]*imw[7]; 
    c[8] = PWORT * y[8]*imw[8]; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 9; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]*imw[0]; 
    c[1] = 1e6 * (*rho) * y[1]*imw[1]; 
    c[2] = 1e6 * (*rho) * y[2]*imw[2]; 
    c[3] = 1e6 * (*rho) * y[3]*imw[3]; 
    c[4] = 1e6 * (*rho) * y[4]*imw[4]; 
    c[5] = 1e6 * (*rho) * y[5]*imw[5]; 
    c[6] = 1e6 * (*rho) * y[6]*imw[6]; 
    c[7] = 1e6 * (*rho) * y[7]*imw[7]; 
    c[8] = 1e6 * (*rho) * y[8]*imw[8]; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 9; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the stoichiometric coefficients */
/*of the reaction mechanism. (Eq 50) */
void CKNU(int * kdim, int * iwrk, double * restrict rwrk, int * nuki)
{
    int id; /*loop counter */
    int kd = (*kdim); 
    /*Zero nuki */
    for (id = 0; id < 9 * kd; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    nuki[ 3 * kd + 0 ] += -1 ;
    nuki[ 1 * kd + 0 ] += -1 ;
    nuki[ 6 * kd + 0 ] += +1 ;

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    nuki[ 7 * kd + 1 ] += -1 ;
    nuki[ 5 * kd + 1 ] += +1 ;
    nuki[ 5 * kd + 1 ] += +1 ;

    /*reaction 3: H2 + M <=> H + H + M */
    nuki[ 0 * kd + 2 ] += -1 ;
    nuki[ 3 * kd + 2 ] += +1 ;
    nuki[ 3 * kd + 2 ] += +1 ;

    /*reaction 4: O + O + M <=> O2 + M */
    nuki[ 4 * kd + 3 ] += -1 ;
    nuki[ 4 * kd + 3 ] += -1 ;
    nuki[ 1 * kd + 3 ] += +1 ;

    /*reaction 5: O + H + M <=> OH + M */
    nuki[ 4 * kd + 4 ] += -1 ;
    nuki[ 3 * kd + 4 ] += -1 ;
    nuki[ 5 * kd + 4 ] += +1 ;

    /*reaction 6: H + OH + M <=> H2O + M */
    nuki[ 3 * kd + 5 ] += -1 ;
    nuki[ 5 * kd + 5 ] += -1 ;
    nuki[ 2 * kd + 5 ] += +1 ;

    /*reaction 7: H + O2 <=> O + OH */
    nuki[ 3 * kd + 6 ] += -1 ;
    nuki[ 1 * kd + 6 ] += -1 ;
    nuki[ 4 * kd + 6 ] += +1 ;
    nuki[ 5 * kd + 6 ] += +1 ;

    /*reaction 8: O + H2 <=> H + OH */
    nuki[ 4 * kd + 7 ] += -1 ;
    nuki[ 0 * kd + 7 ] += -1 ;
    nuki[ 3 * kd + 7 ] += +1 ;
    nuki[ 5 * kd + 7 ] += +1 ;

    /*reaction 9: H2 + OH <=> H2O + H */
    nuki[ 0 * kd + 8 ] += -1 ;
    nuki[ 5 * kd + 8 ] += -1 ;
    nuki[ 2 * kd + 8 ] += +1 ;
    nuki[ 3 * kd + 8 ] += +1 ;

    /*reaction 10: O + H2O <=> OH + OH */
    nuki[ 4 * kd + 9 ] += -1 ;
    nuki[ 2 * kd + 9 ] += -1 ;
    nuki[ 5 * kd + 9 ] += +1 ;
    nuki[ 5 * kd + 9 ] += +1 ;

    /*reaction 11: HO2 + H <=> H2 + O2 */
    nuki[ 6 * kd + 10 ] += -1 ;
    nuki[ 3 * kd + 10 ] += -1 ;
    nuki[ 0 * kd + 10 ] += +1 ;
    nuki[ 1 * kd + 10 ] += +1 ;

    /*reaction 12: HO2 + H <=> OH + OH */
    nuki[ 6 * kd + 11 ] += -1 ;
    nuki[ 3 * kd + 11 ] += -1 ;
    nuki[ 5 * kd + 11 ] += +1 ;
    nuki[ 5 * kd + 11 ] += +1 ;

    /*reaction 13: HO2 + O <=> O2 + OH */
    nuki[ 6 * kd + 12 ] += -1 ;
    nuki[ 4 * kd + 12 ] += -1 ;
    nuki[ 1 * kd + 12 ] += +1 ;
    nuki[ 5 * kd + 12 ] += +1 ;

    /*reaction 14: HO2 + OH <=> H2O + O2 */
    nuki[ 6 * kd + 13 ] += -1 ;
    nuki[ 5 * kd + 13 ] += -1 ;
    nuki[ 2 * kd + 13 ] += +1 ;
    nuki[ 1 * kd + 13 ] += +1 ;

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    nuki[ 6 * kd + 14 ] += -1 ;
    nuki[ 6 * kd + 14 ] += -1 ;
    nuki[ 7 * kd + 14 ] += +1 ;
    nuki[ 1 * kd + 14 ] += +1 ;

    /*reaction 16: HO2 + HO2 <=> H2O2 + O2 */
    nuki[ 6 * kd + 15 ] += -1 ;
    nuki[ 6 * kd + 15 ] += -1 ;
    nuki[ 7 * kd + 15 ] += +1 ;
    nuki[ 1 * kd + 15 ] += +1 ;

    /*reaction 17: H2O2 + H <=> H2O + OH */
    nuki[ 7 * kd + 16 ] += -1 ;
    nuki[ 3 * kd + 16 ] += -1 ;
    nuki[ 2 * kd + 16 ] += +1 ;
    nuki[ 5 * kd + 16 ] += +1 ;

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    nuki[ 7 * kd + 17 ] += -1 ;
    nuki[ 3 * kd + 17 ] += -1 ;
    nuki[ 6 * kd + 17 ] += +1 ;
    nuki[ 0 * kd + 17 ] += +1 ;

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    nuki[ 7 * kd + 18 ] += -1 ;
    nuki[ 4 * kd + 18 ] += -1 ;
    nuki[ 5 * kd + 18 ] += +1 ;
    nuki[ 6 * kd + 18 ] += +1 ;

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    nuki[ 7 * kd + 19 ] += -1 ;
    nuki[ 5 * kd + 19 ] += -1 ;
    nuki[ 6 * kd + 19 ] += +1 ;
    nuki[ 2 * kd + 19 ] += +1 ;

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    nuki[ 7 * kd + 20 ] += -1 ;
    nuki[ 5 * kd + 20 ] += -1 ;
    nuki[ 6 * kd + 20 ] += +1 ;
    nuki[ 2 * kd + 20 ] += +1 ;
}


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim, int * iwrk, double * restrict rwrk, int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < kd * 9; ++ id) {
         ncf[id] = 0; 
    }

    /*H2 */
    ncf[ 0 * kd + 0 ] = 2; /*H */

    /*O2 */
    ncf[ 1 * kd + 1 ] = 2; /*O */

    /*H2O */
    ncf[ 2 * kd + 0 ] = 2; /*H */
    ncf[ 2 * kd + 1 ] = 1; /*O */

    /*H */
    ncf[ 3 * kd + 0 ] = 1; /*H */

    /*O */
    ncf[ 4 * kd + 1 ] = 1; /*O */

    /*OH */
    ncf[ 5 * kd + 1 ] = 1; /*O */
    ncf[ 5 * kd + 0 ] = 1; /*H */

    /*HO2 */
    ncf[ 6 * kd + 0 ] = 1; /*H */
    ncf[ 6 * kd + 1 ] = 2; /*O */

    /*H2O2 */
    ncf[ 7 * kd + 0 ] = 2; /*H */
    ncf[ 7 * kd + 1 ] = 2; /*O */

    /*N2 */
    ncf[ 8 * kd + 2 ] = 2; /*N */

}


/*Returns the arrehenius coefficients  */
/*for all reactions */
void CKABE(int * iwrk, double * restrict rwrk, double * restrict a, double * restrict b, double * restrict e)
{
    for (int i=0; i<21; ++i) {
        a[i] = fwd_A[i];
        b[i] = fwd_beta[i];
        e[i] = fwd_Ea[i];
    }

    return;
}


/*Returns the equil constants for each reaction */
void CKEQC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[9]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[1] *= 1e-06; 

    /*reaction 3: H2 + M <=> H + H + M */
    eqcon[2] *= 1e-06; 

    /*reaction 4: O + O + M <=> O2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: O + H + M <=> OH + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H + OH + M <=> H2O + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + O2 <=> O + OH */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + H2 <=> H + OH */
    /*eqcon[7] *= 1;  */

    /*reaction 9: H2 + OH <=> H2O + H */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + H2O <=> OH + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: HO2 + H <=> H2 + O2 */
    /*eqcon[10] *= 1;  */

    /*reaction 12: HO2 + H <=> OH + OH */
    /*eqcon[11] *= 1;  */

    /*reaction 13: HO2 + O <=> O2 + OH */
    /*eqcon[12] *= 1;  */

    /*reaction 14: HO2 + OH <=> H2O + O2 */
    /*eqcon[13] *= 1;  */

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[14] *= 1;  */

    /*reaction 16: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: H2O2 + H <=> H2O + OH */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[20] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[9]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[1] *= 1e-06; 

    /*reaction 3: H2 + M <=> H + H + M */
    eqcon[2] *= 1e-06; 

    /*reaction 4: O + O + M <=> O2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: O + H + M <=> OH + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H + OH + M <=> H2O + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + O2 <=> O + OH */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + H2 <=> H + OH */
    /*eqcon[7] *= 1;  */

    /*reaction 9: H2 + OH <=> H2O + H */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + H2O <=> OH + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: HO2 + H <=> H2 + O2 */
    /*eqcon[10] *= 1;  */

    /*reaction 12: HO2 + H <=> OH + OH */
    /*eqcon[11] *= 1;  */

    /*reaction 13: HO2 + O <=> O2 + OH */
    /*eqcon[12] *= 1;  */

    /*reaction 14: HO2 + OH <=> H2O + O2 */
    /*eqcon[13] *= 1;  */

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[14] *= 1;  */

    /*reaction 16: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: H2O2 + H <=> H2O + OH */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[20] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[9]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[1] *= 1e-06; 

    /*reaction 3: H2 + M <=> H + H + M */
    eqcon[2] *= 1e-06; 

    /*reaction 4: O + O + M <=> O2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: O + H + M <=> OH + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H + OH + M <=> H2O + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + O2 <=> O + OH */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + H2 <=> H + OH */
    /*eqcon[7] *= 1;  */

    /*reaction 9: H2 + OH <=> H2O + H */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + H2O <=> OH + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: HO2 + H <=> H2 + O2 */
    /*eqcon[10] *= 1;  */

    /*reaction 12: HO2 + H <=> OH + OH */
    /*eqcon[11] *= 1;  */

    /*reaction 13: HO2 + O <=> O2 + OH */
    /*eqcon[12] *= 1;  */

    /*reaction 14: HO2 + OH <=> H2O + O2 */
    /*eqcon[13] *= 1;  */

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[14] *= 1;  */

    /*reaction 16: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: H2O2 + H <=> H2O + OH */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[20] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[9]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[1] *= 1e-06; 

    /*reaction 3: H2 + M <=> H + H + M */
    eqcon[2] *= 1e-06; 

    /*reaction 4: O + O + M <=> O2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: O + H + M <=> OH + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H + OH + M <=> H2O + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + O2 <=> O + OH */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + H2 <=> H + OH */
    /*eqcon[7] *= 1;  */

    /*reaction 9: H2 + OH <=> H2O + H */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + H2O <=> OH + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: HO2 + H <=> H2 + O2 */
    /*eqcon[10] *= 1;  */

    /*reaction 12: HO2 + H <=> OH + OH */
    /*eqcon[11] *= 1;  */

    /*reaction 13: HO2 + O <=> O2 + OH */
    /*eqcon[12] *= 1;  */

    /*reaction 14: HO2 + OH <=> H2O + O2 */
    /*eqcon[13] *= 1;  */

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[14] *= 1;  */

    /*reaction 16: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: H2O2 + H <=> H2O + OH */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[20] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[9]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[1] *= 1e-06; 

    /*reaction 3: H2 + M <=> H + H + M */
    eqcon[2] *= 1e-06; 

    /*reaction 4: O + O + M <=> O2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: O + H + M <=> OH + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H + OH + M <=> H2O + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: H + O2 <=> O + OH */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + H2 <=> H + OH */
    /*eqcon[7] *= 1;  */

    /*reaction 9: H2 + OH <=> H2O + H */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + H2O <=> OH + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: HO2 + H <=> H2 + O2 */
    /*eqcon[10] *= 1;  */

    /*reaction 12: HO2 + H <=> OH + OH */
    /*eqcon[11] *= 1;  */

    /*reaction 13: HO2 + O <=> O2 + OH */
    /*eqcon[12] *= 1;  */

    /*reaction 14: HO2 + OH <=> H2O + O2 */
    /*eqcon[13] *= 1;  */

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[14] *= 1;  */

    /*reaction 16: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: H2O2 + H <=> H2O + OH */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[20] *= 1;  */
}

static double T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static double k_f_save[21];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static double Kc_save[21];
#ifdef _OPENMP
#pragma omp threadprivate(Kc_save)
#endif


/*compute the production rate for each species */
void productionRate(double * restrict wdot, double * restrict sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }

    double qdot, q_f[21], q_r[21];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 9; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[1]-q_r[1];
    wdot[5] += qdot;
    wdot[5] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[2]-q_r[2];
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[3] += qdot;

    qdot = q_f[3]-q_r[3];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[4] -= qdot;

    qdot = q_f[4]-q_r[4];
    wdot[3] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[5]-q_r[5];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[5] -= qdot;

    qdot = q_f[6]-q_r[6];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[5] += qdot;

    qdot = q_f[7]-q_r[7];
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[8]-q_r[8];
    wdot[0] -= qdot;
    wdot[2] += qdot;
    wdot[3] += qdot;
    wdot[5] -= qdot;

    qdot = q_f[9]-q_r[9];
    wdot[2] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[5] += qdot;

    qdot = q_f[10]-q_r[10];
    wdot[0] += qdot;
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[6] -= qdot;

    qdot = q_f[11]-q_r[11];
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[12]-q_r[12];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[13]-q_r[13];
    wdot[1] += qdot;
    wdot[2] += qdot;
    wdot[5] -= qdot;
    wdot[6] -= qdot;

    qdot = q_f[14]-q_r[14];
    wdot[1] += qdot;
    wdot[6] -= qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[15]-q_r[15];
    wdot[1] += qdot;
    wdot[6] -= qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[16]-q_r[16];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[17]-q_r[17];
    wdot[0] += qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[18]-q_r[18];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[19]-q_r[19];
    wdot[2] += qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[20]-q_r[20];
    wdot[2] += qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    return;
}

void comp_k_f(double * restrict tc, double invT, double * restrict k_f)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<21; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(double * restrict tc, double invT, double * restrict Kc)
{
    /*compute the Gibbs free energy */
    double g_RT[9];
    gibbs(g_RT, tc);

    Kc[0] = g_RT[1] + g_RT[3] - g_RT[6];
    Kc[1] = -g_RT[5] - g_RT[5] + g_RT[7];
    Kc[2] = g_RT[0] - g_RT[3] - g_RT[3];
    Kc[3] = -g_RT[1] + g_RT[4] + g_RT[4];
    Kc[4] = g_RT[3] + g_RT[4] - g_RT[5];
    Kc[5] = -g_RT[2] + g_RT[3] + g_RT[5];
    Kc[6] = g_RT[1] + g_RT[3] - g_RT[4] - g_RT[5];
    Kc[7] = g_RT[0] - g_RT[3] + g_RT[4] - g_RT[5];
    Kc[8] = g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5];
    Kc[9] = g_RT[2] + g_RT[4] - g_RT[5] - g_RT[5];
    Kc[10] = -g_RT[0] - g_RT[1] + g_RT[3] + g_RT[6];
    Kc[11] = g_RT[3] - g_RT[5] - g_RT[5] + g_RT[6];
    Kc[12] = -g_RT[1] + g_RT[4] - g_RT[5] + g_RT[6];
    Kc[13] = -g_RT[1] - g_RT[2] + g_RT[5] + g_RT[6];
    Kc[14] = -g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7];
    Kc[15] = -g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7];
    Kc[16] = -g_RT[2] + g_RT[3] - g_RT[5] + g_RT[7];
    Kc[17] = -g_RT[0] + g_RT[3] - g_RT[6] + g_RT[7];
    Kc[18] = g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7];
    Kc[19] = -g_RT[2] + g_RT[5] - g_RT[6] + g_RT[7];
    Kc[20] = -g_RT[2] + g_RT[5] - g_RT[6] + g_RT[7];

#ifdef __INTEL_COMPILER
     #pragma simd
#endif
    for (int i=0; i<21; ++i) {
        Kc[i] = exp(Kc[i]);
    };

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 * invT;
    double refCinv = 1 / refC;

    Kc[0] *= refCinv;
    Kc[1] *= refC;
    Kc[2] *= refC;
    Kc[3] *= refCinv;
    Kc[4] *= refCinv;
    Kc[5] *= refCinv;

    return;
}

void comp_qfqr(double * restrict qf, double * restrict qr, double * restrict sc, double * restrict tc, double invT)
{

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    qf[0] = sc[1]*sc[3];
    qr[0] = sc[6];

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    qf[1] = sc[7];
    qr[1] = sc[5]*sc[5];

    /*reaction 3: H2 + M <=> H + H + M */
    qf[2] = sc[0];
    qr[2] = sc[3]*sc[3];

    /*reaction 4: O + O + M <=> O2 + M */
    qf[3] = sc[4]*sc[4];
    qr[3] = sc[1];

    /*reaction 5: O + H + M <=> OH + M */
    qf[4] = sc[3]*sc[4];
    qr[4] = sc[5];

    /*reaction 6: H + OH + M <=> H2O + M */
    qf[5] = sc[3]*sc[5];
    qr[5] = sc[2];

    /*reaction 7: H + O2 <=> O + OH */
    qf[6] = sc[1]*sc[3];
    qr[6] = sc[4]*sc[5];

    /*reaction 8: O + H2 <=> H + OH */
    qf[7] = sc[0]*sc[4];
    qr[7] = sc[3]*sc[5];

    /*reaction 9: H2 + OH <=> H2O + H */
    qf[8] = sc[0]*sc[5];
    qr[8] = sc[2]*sc[3];

    /*reaction 10: O + H2O <=> OH + OH */
    qf[9] = sc[2]*sc[4];
    qr[9] = sc[5]*sc[5];

    /*reaction 11: HO2 + H <=> H2 + O2 */
    qf[10] = sc[3]*sc[6];
    qr[10] = sc[0]*sc[1];

    /*reaction 12: HO2 + H <=> OH + OH */
    qf[11] = sc[3]*sc[6];
    qr[11] = sc[5]*sc[5];

    /*reaction 13: HO2 + O <=> O2 + OH */
    qf[12] = sc[4]*sc[6];
    qr[12] = sc[1]*sc[5];

    /*reaction 14: HO2 + OH <=> H2O + O2 */
    qf[13] = sc[5]*sc[6];
    qr[13] = sc[1]*sc[2];

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    qf[14] = sc[6]*sc[6];
    qr[14] = sc[1]*sc[7];

    /*reaction 16: HO2 + HO2 <=> H2O2 + O2 */
    qf[15] = sc[6]*sc[6];
    qr[15] = sc[1]*sc[7];

    /*reaction 17: H2O2 + H <=> H2O + OH */
    qf[16] = sc[3]*sc[7];
    qr[16] = sc[2]*sc[5];

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    qf[17] = sc[3]*sc[7];
    qr[17] = sc[0]*sc[6];

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    qf[18] = sc[4]*sc[7];
    qr[18] = sc[5]*sc[6];

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    qf[19] = sc[5]*sc[7];
    qr[19] = sc[2]*sc[6];

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    qf[20] = sc[5]*sc[7];
    qr[20] = sc[2]*sc[6];

    double T = tc[1];

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int i = 0; i < 9; ++i) {
        mixture += sc[i];
    }

    double Corr[21];
    for (int i = 0; i < 21; ++i) {
        Corr[i] = 1.0;
    }

    /* troe */
    {
        double alpha[2];
        alpha[0] = mixture + (TB[0][0] - 1)*sc[0] + (TB[0][1] - 1)*sc[2] + (TB[0][2] - 1)*sc[1];
        alpha[1] = mixture + (TB[1][0] - 1)*sc[0] + (TB[1][1] - 1)*sc[2];
        for (int i=0; i<2; i++)
        {
            double redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;
            redP = alpha[i-0] / k_f_save[i] * phase_units[i] * low_A[i] * exp(low_beta[i] * tc[0] - activation_units[i] * low_Ea[i] *invT);
            F = redP / (1.0 + redP);
            logPred = log10(redP);
            logFcent = log10(
                (fabs(troe_Tsss[i]) > 1.e-100 ? (1.-troe_a[i])*exp(-T/troe_Tsss[i]) : 0.) 
                + (fabs(troe_Ts[i]) > 1.e-100 ? troe_a[i] * exp(-T/troe_Ts[i]) : 0.) 
                + (troe_len[i] == 4 ? exp(-troe_Tss[i] * invT) : 0.) );
            troe_c = -.4 - .67 * logFcent;
            troe_n = .75 - 1.27 * logFcent;
            troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
            F_troe = pow(10., logFcent / (1.0 + troe*troe));
            Corr[i] = F * F_troe;
        }
    }

    /* simple three-body correction */
    {
        double alpha;
        alpha = mixture + (TB[2][0] - 1)*sc[0] + (TB[2][1] - 1)*sc[2];
        Corr[2] = alpha;
        alpha = mixture + (TB[3][0] - 1)*sc[0] + (TB[3][1] - 1)*sc[2];
        Corr[3] = alpha;
        alpha = mixture + (TB[4][0] - 1)*sc[0] + (TB[4][1] - 1)*sc[2];
        Corr[4] = alpha;
        alpha = mixture + (TB[5][0] - 1)*sc[0] + (TB[5][1] - 1)*sc[2];
        Corr[5] = alpha;
    }

    for (int i=0; i<21; i++)
    {
        qf[i] *= Corr[i] * k_f_save[i];
        qr[i] *= Corr[i] * k_f_save[i] / Kc_save[i];
    }

    return;
}


/*compute the production rate for each species */
void vproductionRate(int npt, double * restrict wdot, double * restrict sc, double * restrict T)
{
    double k_f_s[21*npt], Kc_s[21*npt], mixture[npt], g_RT[9*npt];
    double tc[5*npt], invT[npt];

#ifdef __INTEL_COMPILER
     #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        tc[0*npt+i] = log(T[i]);
        tc[1*npt+i] = T[i];
        tc[2*npt+i] = T[i]*T[i];
        tc[3*npt+i] = T[i]*T[i]*T[i];
        tc[4*npt+i] = T[i]*T[i]*T[i]*T[i];
        invT[i] = 1.0 / T[i];
    }

    for (int i=0; i<npt; i++) {
        mixture[i] = 0.0;
    }

    for (int n=0; n<9; n++) {
        for (int i=0; i<npt; i++) {
            mixture[i] += sc[n*npt+i];
            wdot[n*npt+i] = 0.0;
        }
    }

    vcomp_k_f(npt, k_f_s, tc, invT);

    vcomp_gibbs(npt, g_RT, tc);

    vcomp_Kc(npt, Kc_s, g_RT, invT);

    vcomp_wdot(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);
}

void vcomp_k_f(int npt, double * restrict k_f_s, double * restrict tc, double * restrict invT)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        k_f_s[0*npt+i] = prefactor_units[0] * fwd_A[0] * exp(fwd_beta[0] * tc[i] - activation_units[0] * fwd_Ea[0] * invT[i]);
        k_f_s[1*npt+i] = prefactor_units[1] * fwd_A[1] * exp(fwd_beta[1] * tc[i] - activation_units[1] * fwd_Ea[1] * invT[i]);
        k_f_s[2*npt+i] = prefactor_units[2] * fwd_A[2] * exp(fwd_beta[2] * tc[i] - activation_units[2] * fwd_Ea[2] * invT[i]);
        k_f_s[3*npt+i] = prefactor_units[3] * fwd_A[3] * exp(fwd_beta[3] * tc[i] - activation_units[3] * fwd_Ea[3] * invT[i]);
        k_f_s[4*npt+i] = prefactor_units[4] * fwd_A[4] * exp(fwd_beta[4] * tc[i] - activation_units[4] * fwd_Ea[4] * invT[i]);
        k_f_s[5*npt+i] = prefactor_units[5] * fwd_A[5] * exp(fwd_beta[5] * tc[i] - activation_units[5] * fwd_Ea[5] * invT[i]);
        k_f_s[6*npt+i] = prefactor_units[6] * fwd_A[6] * exp(fwd_beta[6] * tc[i] - activation_units[6] * fwd_Ea[6] * invT[i]);
        k_f_s[7*npt+i] = prefactor_units[7] * fwd_A[7] * exp(fwd_beta[7] * tc[i] - activation_units[7] * fwd_Ea[7] * invT[i]);
        k_f_s[8*npt+i] = prefactor_units[8] * fwd_A[8] * exp(fwd_beta[8] * tc[i] - activation_units[8] * fwd_Ea[8] * invT[i]);
        k_f_s[9*npt+i] = prefactor_units[9] * fwd_A[9] * exp(fwd_beta[9] * tc[i] - activation_units[9] * fwd_Ea[9] * invT[i]);
        k_f_s[10*npt+i] = prefactor_units[10] * fwd_A[10] * exp(fwd_beta[10] * tc[i] - activation_units[10] * fwd_Ea[10] * invT[i]);
        k_f_s[11*npt+i] = prefactor_units[11] * fwd_A[11] * exp(fwd_beta[11] * tc[i] - activation_units[11] * fwd_Ea[11] * invT[i]);
        k_f_s[12*npt+i] = prefactor_units[12] * fwd_A[12] * exp(fwd_beta[12] * tc[i] - activation_units[12] * fwd_Ea[12] * invT[i]);
        k_f_s[13*npt+i] = prefactor_units[13] * fwd_A[13] * exp(fwd_beta[13] * tc[i] - activation_units[13] * fwd_Ea[13] * invT[i]);
        k_f_s[14*npt+i] = prefactor_units[14] * fwd_A[14] * exp(fwd_beta[14] * tc[i] - activation_units[14] * fwd_Ea[14] * invT[i]);
        k_f_s[15*npt+i] = prefactor_units[15] * fwd_A[15] * exp(fwd_beta[15] * tc[i] - activation_units[15] * fwd_Ea[15] * invT[i]);
        k_f_s[16*npt+i] = prefactor_units[16] * fwd_A[16] * exp(fwd_beta[16] * tc[i] - activation_units[16] * fwd_Ea[16] * invT[i]);
        k_f_s[17*npt+i] = prefactor_units[17] * fwd_A[17] * exp(fwd_beta[17] * tc[i] - activation_units[17] * fwd_Ea[17] * invT[i]);
        k_f_s[18*npt+i] = prefactor_units[18] * fwd_A[18] * exp(fwd_beta[18] * tc[i] - activation_units[18] * fwd_Ea[18] * invT[i]);
        k_f_s[19*npt+i] = prefactor_units[19] * fwd_A[19] * exp(fwd_beta[19] * tc[i] - activation_units[19] * fwd_Ea[19] * invT[i]);
        k_f_s[20*npt+i] = prefactor_units[20] * fwd_A[20] * exp(fwd_beta[20] * tc[i] - activation_units[20] * fwd_Ea[20] * invT[i]);
    }
}

void vcomp_gibbs(int npt, double * restrict g_RT, double * restrict tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        double tg[5], g[9];
        tg[0] = tc[0*npt+i];
        tg[1] = tc[1*npt+i];
        tg[2] = tc[2*npt+i];
        tg[3] = tc[3*npt+i];
        tg[4] = tc[4*npt+i];

        gibbs(g, tg);

        g_RT[0*npt+i] = g[0];
        g_RT[1*npt+i] = g[1];
        g_RT[2*npt+i] = g[2];
        g_RT[3*npt+i] = g[3];
        g_RT[4*npt+i] = g[4];
        g_RT[5*npt+i] = g[5];
        g_RT[6*npt+i] = g[6];
        g_RT[7*npt+i] = g[7];
        g_RT[8*npt+i] = g[8];
    }
}

void vcomp_Kc(int npt, double * restrict Kc_s, double * restrict g_RT, double * restrict invT)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
        double refC = (101325. / 8.31451) * invT[i];
        double refCinv = 1.0 / refC;

        Kc_s[0*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[3*npt+i]) - (g_RT[6*npt+i]));
        Kc_s[1*npt+i] = refC * exp((g_RT[7*npt+i]) - (g_RT[5*npt+i] + g_RT[5*npt+i]));
        Kc_s[2*npt+i] = refC * exp((g_RT[0*npt+i]) - (g_RT[3*npt+i] + g_RT[3*npt+i]));
        Kc_s[3*npt+i] = refCinv * exp((g_RT[4*npt+i] + g_RT[4*npt+i]) - (g_RT[1*npt+i]));
        Kc_s[4*npt+i] = refCinv * exp((g_RT[3*npt+i] + g_RT[4*npt+i]) - (g_RT[5*npt+i]));
        Kc_s[5*npt+i] = refCinv * exp((g_RT[3*npt+i] + g_RT[5*npt+i]) - (g_RT[2*npt+i]));
        Kc_s[6*npt+i] = exp((g_RT[1*npt+i] + g_RT[3*npt+i]) - (g_RT[4*npt+i] + g_RT[5*npt+i]));
        Kc_s[7*npt+i] = exp((g_RT[0*npt+i] + g_RT[4*npt+i]) - (g_RT[3*npt+i] + g_RT[5*npt+i]));
        Kc_s[8*npt+i] = exp((g_RT[0*npt+i] + g_RT[5*npt+i]) - (g_RT[2*npt+i] + g_RT[3*npt+i]));
        Kc_s[9*npt+i] = exp((g_RT[2*npt+i] + g_RT[4*npt+i]) - (g_RT[5*npt+i] + g_RT[5*npt+i]));
        Kc_s[10*npt+i] = exp((g_RT[3*npt+i] + g_RT[6*npt+i]) - (g_RT[0*npt+i] + g_RT[1*npt+i]));
        Kc_s[11*npt+i] = exp((g_RT[3*npt+i] + g_RT[6*npt+i]) - (g_RT[5*npt+i] + g_RT[5*npt+i]));
        Kc_s[12*npt+i] = exp((g_RT[4*npt+i] + g_RT[6*npt+i]) - (g_RT[1*npt+i] + g_RT[5*npt+i]));
        Kc_s[13*npt+i] = exp((g_RT[5*npt+i] + g_RT[6*npt+i]) - (g_RT[1*npt+i] + g_RT[2*npt+i]));
        Kc_s[14*npt+i] = exp((g_RT[6*npt+i] + g_RT[6*npt+i]) - (g_RT[1*npt+i] + g_RT[7*npt+i]));
        Kc_s[15*npt+i] = exp((g_RT[6*npt+i] + g_RT[6*npt+i]) - (g_RT[1*npt+i] + g_RT[7*npt+i]));
        Kc_s[16*npt+i] = exp((g_RT[3*npt+i] + g_RT[7*npt+i]) - (g_RT[2*npt+i] + g_RT[5*npt+i]));
        Kc_s[17*npt+i] = exp((g_RT[3*npt+i] + g_RT[7*npt+i]) - (g_RT[0*npt+i] + g_RT[6*npt+i]));
        Kc_s[18*npt+i] = exp((g_RT[4*npt+i] + g_RT[7*npt+i]) - (g_RT[5*npt+i] + g_RT[6*npt+i]));
        Kc_s[19*npt+i] = exp((g_RT[5*npt+i] + g_RT[7*npt+i]) - (g_RT[2*npt+i] + g_RT[6*npt+i]));
        Kc_s[20*npt+i] = exp((g_RT[5*npt+i] + g_RT[7*npt+i]) - (g_RT[2*npt+i] + g_RT[6*npt+i]));
    }
}

void vcomp_wdot(int npt, double * restrict wdot, double * restrict mixture, double * restrict sc,
		double * restrict k_f_s, double * restrict Kc_s,
		double * restrict tc, double * restrict invT, double * restrict T)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        double qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        double alpha;
        double redP, F;
        double logPred;
        double logFcent, troe_c, troe_n, troe, F_troe;

        /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
        phi_f = sc[1*npt+i]*sc[3*npt+i];
        alpha = mixture[i] + (TB[0][0] - 1)*sc[0*npt+i] + (TB[0][1] - 1)*sc[2*npt+i] + (TB[0][2] - 1)*sc[1*npt+i];
        k_f = k_f_s[0*npt+i];
        redP = alpha / k_f * phase_units[0] * low_A[0] * exp(low_beta[0] * tc[i] - activation_units[0] * low_Ea[0] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[0]) > 1.e-100 ? (1.-troe_a[0])*exp(-T[i]/troe_Tsss[0]) : 0.) 
            + (fabs(troe_Ts[0]) > 1.e-100 ? troe_a[0] * exp(-T[i]/troe_Ts[0]) : 0.) 
            + (troe_len[0] == 4 ? exp(-troe_Tss[0] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i];
        Kc = Kc_s[0*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
        phi_f = sc[7*npt+i];
        alpha = mixture[i] + (TB[1][0] - 1)*sc[0*npt+i] + (TB[1][1] - 1)*sc[2*npt+i];
        k_f = k_f_s[1*npt+i];
        redP = alpha / k_f * phase_units[1] * low_A[1] * exp(low_beta[1] * tc[i] - activation_units[1] * low_Ea[1] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[1]) > 1.e-100 ? (1.-troe_a[1])*exp(-T[i]/troe_Tsss[1]) : 0.) 
            + (fabs(troe_Ts[1]) > 1.e-100 ? troe_a[1] * exp(-T[i]/troe_Ts[1]) : 0.) 
            + (troe_len[1] == 4 ? exp(-troe_Tss[1] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[5*npt+i];
        Kc = Kc_s[1*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[5*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 3: H2 + M <=> H + H + M */
        phi_f = sc[0*npt+i];
        alpha = mixture[i] + (TB[2][0] - 1)*sc[0*npt+i] + (TB[2][1] - 1)*sc[2*npt+i];
        k_f = alpha * k_f_s[2*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[3*npt+i];
        Kc = Kc_s[2*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[3*npt+i] += qdot;

        /*reaction 4: O + O + M <=> O2 + M */
        phi_f = sc[4*npt+i]*sc[4*npt+i];
        alpha = mixture[i] + (TB[3][0] - 1)*sc[0*npt+i] + (TB[3][1] - 1)*sc[2*npt+i];
        k_f = alpha * k_f_s[3*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i];
        Kc = Kc_s[3*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;

        /*reaction 5: O + H + M <=> OH + M */
        phi_f = sc[3*npt+i]*sc[4*npt+i];
        alpha = mixture[i] + (TB[4][0] - 1)*sc[0*npt+i] + (TB[4][1] - 1)*sc[2*npt+i];
        k_f = alpha * k_f_s[4*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i];
        Kc = Kc_s[4*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 6: H + OH + M <=> H2O + M */
        phi_f = sc[3*npt+i]*sc[5*npt+i];
        alpha = mixture[i] + (TB[5][0] - 1)*sc[0*npt+i] + (TB[5][1] - 1)*sc[2*npt+i];
        k_f = alpha * k_f_s[5*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i];
        Kc = Kc_s[5*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;

        /*reaction 7: H + O2 <=> O + OH */
        phi_f = sc[1*npt+i]*sc[3*npt+i];
        k_f = k_f_s[6*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[5*npt+i];
        Kc = Kc_s[6*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 8: O + H2 <=> H + OH */
        phi_f = sc[0*npt+i]*sc[4*npt+i];
        k_f = k_f_s[7*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[5*npt+i];
        Kc = Kc_s[7*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 9: H2 + OH <=> H2O + H */
        phi_f = sc[0*npt+i]*sc[5*npt+i];
        k_f = k_f_s[8*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[3*npt+i];
        Kc = Kc_s[8*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;

        /*reaction 10: O + H2O <=> OH + OH */
        phi_f = sc[2*npt+i]*sc[4*npt+i];
        k_f = k_f_s[9*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[5*npt+i];
        Kc = Kc_s[9*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 11: HO2 + H <=> H2 + O2 */
        phi_f = sc[3*npt+i]*sc[6*npt+i];
        k_f = k_f_s[10*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[1*npt+i];
        Kc = Kc_s[10*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 12: HO2 + H <=> OH + OH */
        phi_f = sc[3*npt+i]*sc[6*npt+i];
        k_f = k_f_s[11*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[5*npt+i];
        Kc = Kc_s[11*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 13: HO2 + O <=> O2 + OH */
        phi_f = sc[4*npt+i]*sc[6*npt+i];
        k_f = k_f_s[12*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[5*npt+i];
        Kc = Kc_s[12*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 14: HO2 + OH <=> H2O + O2 */
        phi_f = sc[5*npt+i]*sc[6*npt+i];
        k_f = k_f_s[13*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[2*npt+i];
        Kc = Kc_s[13*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
        phi_f = sc[6*npt+i]*sc[6*npt+i];
        k_f = k_f_s[14*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[7*npt+i];
        Kc = Kc_s[14*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 16: HO2 + HO2 <=> H2O2 + O2 */
        phi_f = sc[6*npt+i]*sc[6*npt+i];
        k_f = k_f_s[15*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[7*npt+i];
        Kc = Kc_s[15*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 17: H2O2 + H <=> H2O + OH */
        phi_f = sc[3*npt+i]*sc[7*npt+i];
        k_f = k_f_s[16*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[5*npt+i];
        Kc = Kc_s[16*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 18: H2O2 + H <=> HO2 + H2 */
        phi_f = sc[3*npt+i]*sc[7*npt+i];
        k_f = k_f_s[17*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[6*npt+i];
        Kc = Kc_s[17*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 19: H2O2 + O <=> OH + HO2 */
        phi_f = sc[4*npt+i]*sc[7*npt+i];
        k_f = k_f_s[18*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[6*npt+i];
        Kc = Kc_s[18*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 20: H2O2 + OH <=> HO2 + H2O */
        phi_f = sc[5*npt+i]*sc[7*npt+i];
        k_f = k_f_s[19*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[6*npt+i];
        Kc = Kc_s[19*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 21: H2O2 + OH <=> HO2 + H2O */
        phi_f = sc[5*npt+i]*sc[7*npt+i];
        k_f = k_f_s[20*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[6*npt+i];
        Kc = Kc_s[20*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
    }
}

/*compute the reaction Jacobian */
void DWDOT(double * restrict J, double * restrict sc, double * restrict Tp, int * consP)
{
    double c[9];

    for (int k=0; k<9; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian(J, c, *Tp, *consP);

    /* dwdot[k]/dT */
    for (int k=0; k<9; k++) {
        J[90+k] *= 1.e-6;
    }

    /* dTdot/d[X] */
    for (int k=0; k<9; k++) {
        J[k*10+9] *= 1.e6;
    }

    return;
}

/*compute the reaction Jacobian */
void aJacobian(double * restrict J, double * restrict sc, double T, int consP)
{
    for (int i=0; i<100; i++) {
        J[i] = 0.0;
    }

    double wdot[9];
    for (int k=0; k<9; k++) {
        wdot[k] = 0.0;
    }

    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];
    double invT2 = invT * invT;

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;
    double refCinv = 1.0 / refC;

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int k = 0; k < 9; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[9];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[9];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[9];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[0][0] - 1)*sc[0] + (TB[0][1] - 1)*sc[2] + (TB[0][2] - 1)*sc[1];
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = prefactor_units[0] * fwd_A[0]
                * exp(fwd_beta[0] * tc[0] - activation_units[0] * fwd_Ea[0] * invT);
    dlnkfdT = fwd_beta[0] * invT + activation_units[0] * fwd_Ea[0] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[0] * exp(low_beta[0] * tc[0] - activation_units[0] * low_Ea[0] * invT);
    Pr = phase_units[0] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[0] * invT + activation_units[0] * low_Ea[0] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[0]) > 1.e-100 ? (1.-troe_a[0])*exp(-T/troe_Tsss[0]) : 0.);
    Fcent2 = (fabs(troe_Ts[0]) > 1.e-100 ? troe_a[0] * exp(-T/troe_Ts[0]) : 0.);
    Fcent3 = (troe_len[0] == 4 ? exp(-troe_Tss[0] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[0]) > 1.e-100 ? -Fcent1/troe_Tsss[0] : 0.)
      + (fabs(troe_Ts[0]) > 1.e-100 ? -Fcent2/troe_Ts[0] : 0.)
      + (troe_len[0] == 4 ? Fcent3*troe_Tss[0]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[6];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[6]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[3] -= q; /* H */
    wdot[6] += q; /* HO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[0][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[O2]/d[H2] */
        J[3] -= dqdci;                /* dwdot[H]/d[H2] */
        J[6] += dqdci;                /* dwdot[HO2]/d[H2] */
        /* d()/d[O2] */
        dqdci = (TB[0][2] - 1)*dcdc_fac + k_f*sc[3];
        J[11] -= dqdci;               /* dwdot[O2]/d[O2] */
        J[13] -= dqdci;               /* dwdot[H]/d[O2] */
        J[16] += dqdci;               /* dwdot[HO2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (TB[0][1] - 1)*dcdc_fac;
        J[21] -= dqdci;               /* dwdot[O2]/d[H2O] */
        J[23] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[26] += dqdci;               /* dwdot[HO2]/d[H2O] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[1];
        J[31] -= dqdci;               /* dwdot[O2]/d[H] */
        J[33] -= dqdci;               /* dwdot[H]/d[H] */
        J[36] += dqdci;               /* dwdot[HO2]/d[H] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[61] -= dqdci;               /* dwdot[O2]/d[HO2] */
        J[63] -= dqdci;               /* dwdot[H]/d[HO2] */
        J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    }
    else {
        dqdc[0] = TB[0][0]*dcdc_fac;
        dqdc[1] = TB[0][2]*dcdc_fac + k_f*sc[3];
        dqdc[2] = TB[0][1]*dcdc_fac;
        dqdc[3] = dcdc_fac + k_f*sc[1];
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac - k_r;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        for (int k=0; k<9; k++) {
            J[10*k+1] -= dqdc[k];
            J[10*k+3] -= dqdc[k];
            J[10*k+6] += dqdc[k];
        }
    }
    J[91] -= dqdT; /* dwdot[O2]/dT */
    J[93] -= dqdT; /* dwdot[H]/dT */
    J[96] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[1][0] - 1)*sc[0] + (TB[1][1] - 1)*sc[2];
    /* forward */
    phi_f = sc[7];
    k_f = prefactor_units[1] * fwd_A[1]
                * exp(fwd_beta[1] * tc[0] - activation_units[1] * fwd_Ea[1] * invT);
    dlnkfdT = fwd_beta[1] * invT + activation_units[1] * fwd_Ea[1] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[1] * exp(low_beta[1] * tc[0] - activation_units[1] * low_Ea[1] * invT);
    Pr = phase_units[1] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[1] * invT + activation_units[1] * low_Ea[1] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[1]) > 1.e-100 ? (1.-troe_a[1])*exp(-T/troe_Tsss[1]) : 0.);
    Fcent2 = (fabs(troe_Ts[1]) > 1.e-100 ? troe_a[1] * exp(-T/troe_Ts[1]) : 0.);
    Fcent3 = (troe_len[1] == 4 ? exp(-troe_Tss[1] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[1]) > 1.e-100 ? -Fcent1/troe_Tsss[1] : 0.)
      + (fabs(troe_Ts[1]) > 1.e-100 ? -Fcent2/troe_Ts[1] : 0.)
      + (troe_len[1] == 4 ? Fcent3*troe_Tss[1]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[5]*sc[5];
    Kc = refC * exp(-g_RT[5] - g_RT[5] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7]) + (2*h_RT[5]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[5] += 2 * q; /* OH */
    wdot[7] -= q; /* H2O2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[1][0] - 1)*dcdc_fac;
        J[5] += 2 * dqdci;            /* dwdot[OH]/d[H2] */
        J[7] -= dqdci;                /* dwdot[H2O2]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[1][1] - 1)*dcdc_fac;
        J[25] += 2 * dqdci;           /* dwdot[OH]/d[H2O] */
        J[27] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
        /* d()/d[OH] */
        dqdci =  - k_r*2*sc[5];
        J[55] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
        J[57] -= dqdci;               /* dwdot[H2O2]/d[OH] */
        /* d()/d[H2O2] */
        dqdci =  + k_f;
        J[75] += 2 * dqdci;           /* dwdot[OH]/d[H2O2] */
        J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    }
    else {
        dqdc[0] = TB[1][0]*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = TB[1][1]*dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac - k_r*2*sc[5];
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac + k_f;
        dqdc[8] = dcdc_fac;
        for (int k=0; k<9; k++) {
            J[10*k+5] += 2 * dqdc[k];
            J[10*k+7] -= dqdc[k];
        }
    }
    J[95] += 2 * dqdT; /* dwdot[OH]/dT */
    J[97] -= dqdT; /* dwdot[H2O2]/dT */

    /*reaction 3: H2 + M <=> H + H + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[2][0] - 1)*sc[0] + (TB[2][1] - 1)*sc[2];
    /* forward */
    phi_f = sc[0];
    k_f = prefactor_units[2] * fwd_A[2]
                * exp(fwd_beta[2] * tc[0] - activation_units[2] * fwd_Ea[2] * invT);
    dlnkfdT = fwd_beta[2] * invT + activation_units[2] * fwd_Ea[2] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[3];
    Kc = refC * exp(g_RT[0] - g_RT[3] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0]) + (2*h_RT[3]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[3] += 2 * q; /* H */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[2][0] - 1)*q_nocor + k_f;
        J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
        J[3] += 2 * dqdci;            /* dwdot[H]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[2][1] - 1)*q_nocor;
        J[20] -= dqdci;               /* dwdot[H2]/d[H2O] */
        J[23] += 2 * dqdci;           /* dwdot[H]/d[H2O] */
        /* d()/d[H] */
        dqdci =  - k_r*2*sc[3];
        J[30] -= dqdci;               /* dwdot[H2]/d[H] */
        J[33] += 2 * dqdci;           /* dwdot[H]/d[H] */
    }
    else {
        dqdc[0] = TB[2][0] + k_f;
        dqdc[1] = dcdc_fac;
        dqdc[2] = TB[2][1];
        dqdc[3] = dcdc_fac - k_r*2*sc[3];
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        for (int k=0; k<9; k++) {
            J[10*k+0] -= dqdc[k];
            J[10*k+3] += 2 * dqdc[k];
        }
    }
    J[90] -= dqdT; /* dwdot[H2]/dT */
    J[93] += 2 * dqdT; /* dwdot[H]/dT */

    /*reaction 4: O + O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[3][0] - 1)*sc[0] + (TB[3][1] - 1)*sc[2];
    /* forward */
    phi_f = sc[4]*sc[4];
    k_f = prefactor_units[3] * fwd_A[3]
                * exp(fwd_beta[3] * tc[0] - activation_units[3] * fwd_Ea[3] * invT);
    dlnkfdT = fwd_beta[3] * invT + activation_units[3] * fwd_Ea[3] * invT2;
    /* reverse */
    phi_r = sc[1];
    Kc = refCinv * exp(-g_RT[1] + g_RT[4] + g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[4]) + (h_RT[1]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[4] -= 2 * q; /* O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[3][0] - 1)*q_nocor;
        J[1] += dqdci;                /* dwdot[O2]/d[H2] */
        J[4] += -2 * dqdci;           /* dwdot[O]/d[H2] */
        /* d()/d[O2] */
        dqdci =  - k_r;
        J[11] += dqdci;               /* dwdot[O2]/d[O2] */
        J[14] += -2 * dqdci;          /* dwdot[O]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (TB[3][1] - 1)*q_nocor;
        J[21] += dqdci;               /* dwdot[O2]/d[H2O] */
        J[24] += -2 * dqdci;          /* dwdot[O]/d[H2O] */
        /* d()/d[O] */
        dqdci =  + k_f*2*sc[4];
        J[41] += dqdci;               /* dwdot[O2]/d[O] */
        J[44] += -2 * dqdci;          /* dwdot[O]/d[O] */
    }
    else {
        dqdc[0] = TB[3][0];
        dqdc[1] = dcdc_fac - k_r;
        dqdc[2] = TB[3][1];
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac + k_f*2*sc[4];
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        for (int k=0; k<9; k++) {
            J[10*k+1] += dqdc[k];
            J[10*k+4] += -2 * dqdc[k];
        }
    }
    J[91] += dqdT; /* dwdot[O2]/dT */
    J[94] += -2 * dqdT; /* dwdot[O]/dT */

    /*reaction 5: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[4][0] - 1)*sc[0] + (TB[4][1] - 1)*sc[2];
    /* forward */
    phi_f = sc[3]*sc[4];
    k_f = prefactor_units[4] * fwd_A[4]
                * exp(fwd_beta[4] * tc[0] - activation_units[4] * fwd_Ea[4] * invT);
    dlnkfdT = fwd_beta[4] * invT + activation_units[4] * fwd_Ea[4] * invT2;
    /* reverse */
    phi_r = sc[5];
    Kc = refCinv * exp(g_RT[3] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[4]) + (h_RT[5]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* H */
    wdot[4] -= q; /* O */
    wdot[5] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[4][0] - 1)*q_nocor;
        J[3] -= dqdci;                /* dwdot[H]/d[H2] */
        J[4] -= dqdci;                /* dwdot[O]/d[H2] */
        J[5] += dqdci;                /* dwdot[OH]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[4][1] - 1)*q_nocor;
        J[23] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[24] -= dqdci;               /* dwdot[O]/d[H2O] */
        J[25] += dqdci;               /* dwdot[OH]/d[H2O] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[4];
        J[33] -= dqdci;               /* dwdot[H]/d[H] */
        J[34] -= dqdci;               /* dwdot[O]/d[H] */
        J[35] += dqdci;               /* dwdot[OH]/d[H] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[3];
        J[43] -= dqdci;               /* dwdot[H]/d[O] */
        J[44] -= dqdci;               /* dwdot[O]/d[O] */
        J[45] += dqdci;               /* dwdot[OH]/d[O] */
        /* d()/d[OH] */
        dqdci =  - k_r;
        J[53] -= dqdci;               /* dwdot[H]/d[OH] */
        J[54] -= dqdci;               /* dwdot[O]/d[OH] */
        J[55] += dqdci;               /* dwdot[OH]/d[OH] */
    }
    else {
        dqdc[0] = TB[4][0];
        dqdc[1] = dcdc_fac;
        dqdc[2] = TB[4][1];
        dqdc[3] = dcdc_fac + k_f*sc[4];
        dqdc[4] = dcdc_fac + k_f*sc[3];
        dqdc[5] = dcdc_fac - k_r;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        for (int k=0; k<9; k++) {
            J[10*k+3] -= dqdc[k];
            J[10*k+4] -= dqdc[k];
            J[10*k+5] += dqdc[k];
        }
    }
    J[93] -= dqdT; /* dwdot[H]/dT */
    J[94] -= dqdT; /* dwdot[O]/dT */
    J[95] += dqdT; /* dwdot[OH]/dT */

    /*reaction 6: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[5][0] - 1)*sc[0] + (TB[5][1] - 1)*sc[2];
    /* forward */
    phi_f = sc[3]*sc[5];
    k_f = prefactor_units[5] * fwd_A[5]
                * exp(fwd_beta[5] * tc[0] - activation_units[5] * fwd_Ea[5] * invT);
    dlnkfdT = fwd_beta[5] * invT + activation_units[5] * fwd_Ea[5] * invT2;
    /* reverse */
    phi_r = sc[2];
    Kc = refCinv * exp(-g_RT[2] + g_RT[3] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[5]) + (h_RT[2]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H2O */
    wdot[3] -= q; /* H */
    wdot[5] -= q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[5][0] - 1)*q_nocor;
        J[2] += dqdci;                /* dwdot[H2O]/d[H2] */
        J[3] -= dqdci;                /* dwdot[H]/d[H2] */
        J[5] -= dqdci;                /* dwdot[OH]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[5][1] - 1)*q_nocor - k_r;
        J[22] += dqdci;               /* dwdot[H2O]/d[H2O] */
        J[23] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[25] -= dqdci;               /* dwdot[OH]/d[H2O] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[5];
        J[32] += dqdci;               /* dwdot[H2O]/d[H] */
        J[33] -= dqdci;               /* dwdot[H]/d[H] */
        J[35] -= dqdci;               /* dwdot[OH]/d[H] */
        /* d()/d[OH] */
        dqdci =  + k_f*sc[3];
        J[52] += dqdci;               /* dwdot[H2O]/d[OH] */
        J[53] -= dqdci;               /* dwdot[H]/d[OH] */
        J[55] -= dqdci;               /* dwdot[OH]/d[OH] */
    }
    else {
        dqdc[0] = TB[5][0];
        dqdc[1] = dcdc_fac;
        dqdc[2] = TB[5][1] - k_r;
        dqdc[3] = dcdc_fac + k_f*sc[5];
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac + k_f*sc[3];
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        for (int k=0; k<9; k++) {
            J[10*k+2] += dqdc[k];
            J[10*k+3] -= dqdc[k];
            J[10*k+5] -= dqdc[k];
        }
    }
    J[92] += dqdT; /* dwdot[H2O]/dT */
    J[93] -= dqdT; /* dwdot[H]/dT */
    J[95] -= dqdT; /* dwdot[OH]/dT */

    /*reaction 7: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = prefactor_units[6] * fwd_A[6]
                * exp(fwd_beta[6] * tc[0] - activation_units[6] * fwd_Ea[6] * invT);
    dlnkfdT = fwd_beta[6] * invT + activation_units[6] * fwd_Ea[6] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[5];
    Kc = exp(g_RT[1] + g_RT[3] - g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[4] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[3] -= q; /* H */
    wdot[4] += q; /* O */
    wdot[5] += q; /* OH */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[3];
    J[11] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[13] -= dqdci;               /* dwdot[H]/d[O2] */
    J[14] += dqdci;               /* dwdot[O]/d[O2] */
    J[15] += dqdci;               /* dwdot[OH]/d[O2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[1];
    J[31] -= dqdci;               /* dwdot[O2]/d[H] */
    J[33] -= dqdci;               /* dwdot[H]/d[H] */
    J[34] += dqdci;               /* dwdot[O]/d[H] */
    J[35] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[5];
    J[41] -= dqdci;               /* dwdot[O2]/d[O] */
    J[43] -= dqdci;               /* dwdot[H]/d[O] */
    J[44] += dqdci;               /* dwdot[O]/d[O] */
    J[45] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[4];
    J[51] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[53] -= dqdci;               /* dwdot[H]/d[OH] */
    J[54] += dqdci;               /* dwdot[O]/d[OH] */
    J[55] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[O2]/dT */
    J[93] -= dqdT;                /* dwdot[H]/dT */
    J[94] += dqdT;                /* dwdot[O]/dT */
    J[95] += dqdT;                /* dwdot[OH]/dT */

    /*reaction 8: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[4];
    k_f = prefactor_units[7] * fwd_A[7]
                * exp(fwd_beta[7] * tc[0] - activation_units[7] * fwd_Ea[7] * invT);
    dlnkfdT = fwd_beta[7] * invT + activation_units[7] * fwd_Ea[7] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[5];
    Kc = exp(g_RT[0] - g_RT[3] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[4]) + (h_RT[3] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[3] += q; /* H */
    wdot[4] -= q; /* O */
    wdot[5] += q; /* OH */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[4];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[3] += dqdci;                /* dwdot[H]/d[H2] */
    J[4] -= dqdci;                /* dwdot[O]/d[H2] */
    J[5] += dqdci;                /* dwdot[OH]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[5];
    J[30] -= dqdci;               /* dwdot[H2]/d[H] */
    J[33] += dqdci;               /* dwdot[H]/d[H] */
    J[34] -= dqdci;               /* dwdot[O]/d[H] */
    J[35] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[0];
    J[40] -= dqdci;               /* dwdot[H2]/d[O] */
    J[43] += dqdci;               /* dwdot[H]/d[O] */
    J[44] -= dqdci;               /* dwdot[O]/d[O] */
    J[45] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[3];
    J[50] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[53] += dqdci;               /* dwdot[H]/d[OH] */
    J[54] -= dqdci;               /* dwdot[O]/d[OH] */
    J[55] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[90] -= dqdT;                /* dwdot[H2]/dT */
    J[93] += dqdT;                /* dwdot[H]/dT */
    J[94] -= dqdT;                /* dwdot[O]/dT */
    J[95] += dqdT;                /* dwdot[OH]/dT */

    /*reaction 9: H2 + OH <=> H2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[5];
    k_f = prefactor_units[8] * fwd_A[8]
                * exp(fwd_beta[8] * tc[0] - activation_units[8] * fwd_Ea[8] * invT);
    dlnkfdT = fwd_beta[8] * invT + activation_units[8] * fwd_Ea[8] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[3];
    Kc = exp(g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[2] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[2] += q; /* H2O */
    wdot[3] += q; /* H */
    wdot[5] -= q; /* OH */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[5];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[2] += dqdci;                /* dwdot[H2O]/d[H2] */
    J[3] += dqdci;                /* dwdot[H]/d[H2] */
    J[5] -= dqdci;                /* dwdot[OH]/d[H2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[20] -= dqdci;               /* dwdot[H2]/d[H2O] */
    J[22] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[23] += dqdci;               /* dwdot[H]/d[H2O] */
    J[25] -= dqdci;               /* dwdot[OH]/d[H2O] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[2];
    J[30] -= dqdci;               /* dwdot[H2]/d[H] */
    J[32] += dqdci;               /* dwdot[H2O]/d[H] */
    J[33] += dqdci;               /* dwdot[H]/d[H] */
    J[35] -= dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[0];
    J[50] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[52] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[53] += dqdci;               /* dwdot[H]/d[OH] */
    J[55] -= dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[90] -= dqdT;                /* dwdot[H2]/dT */
    J[92] += dqdT;                /* dwdot[H2O]/dT */
    J[93] += dqdT;                /* dwdot[H]/dT */
    J[95] -= dqdT;                /* dwdot[OH]/dT */

    /*reaction 10: O + H2O <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[4];
    k_f = prefactor_units[9] * fwd_A[9]
                * exp(fwd_beta[9] * tc[0] - activation_units[9] * fwd_Ea[9] * invT);
    dlnkfdT = fwd_beta[9] * invT + activation_units[9] * fwd_Ea[9] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[5];
    Kc = exp(g_RT[2] + g_RT[4] - g_RT[5] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[4]) + (2*h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* H2O */
    wdot[4] -= q; /* O */
    wdot[5] += 2 * q; /* OH */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[4];
    J[22] -= dqdci;               /* dwdot[H2O]/d[H2O] */
    J[24] -= dqdci;               /* dwdot[O]/d[H2O] */
    J[25] += 2 * dqdci;           /* dwdot[OH]/d[H2O] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[2];
    J[42] -= dqdci;               /* dwdot[H2O]/d[O] */
    J[44] -= dqdci;               /* dwdot[O]/d[O] */
    J[45] += 2 * dqdci;           /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*2*sc[5];
    J[52] -= dqdci;               /* dwdot[H2O]/d[OH] */
    J[54] -= dqdci;               /* dwdot[O]/d[OH] */
    J[55] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[92] -= dqdT;                /* dwdot[H2O]/dT */
    J[94] -= dqdT;                /* dwdot[O]/dT */
    J[95] += 2 * dqdT;            /* dwdot[OH]/dT */

    /*reaction 11: HO2 + H <=> H2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[6];
    k_f = prefactor_units[10] * fwd_A[10]
                * exp(fwd_beta[10] * tc[0] - activation_units[10] * fwd_Ea[10] * invT);
    dlnkfdT = fwd_beta[10] * invT + activation_units[10] * fwd_Ea[10] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[1];
    Kc = exp(-g_RT[0] - g_RT[1] + g_RT[3] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (h_RT[0] + h_RT[1]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] += q; /* O2 */
    wdot[3] -= q; /* H */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[1];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[O2]/d[H2] */
    J[3] -= dqdci;                /* dwdot[H]/d[H2] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H2] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[0];
    J[10] += dqdci;               /* dwdot[H2]/d[O2] */
    J[11] += dqdci;               /* dwdot[O2]/d[O2] */
    J[13] -= dqdci;               /* dwdot[H]/d[O2] */
    J[16] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[30] += dqdci;               /* dwdot[H2]/d[H] */
    J[31] += dqdci;               /* dwdot[O2]/d[H] */
    J[33] -= dqdci;               /* dwdot[H]/d[H] */
    J[36] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[3];
    J[60] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[61] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[63] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[90] += dqdT;                /* dwdot[H2]/dT */
    J[91] += dqdT;                /* dwdot[O2]/dT */
    J[93] -= dqdT;                /* dwdot[H]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 12: HO2 + H <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[6];
    k_f = prefactor_units[11] * fwd_A[11]
                * exp(fwd_beta[11] * tc[0] - activation_units[11] * fwd_Ea[11] * invT);
    dlnkfdT = fwd_beta[11] * invT + activation_units[11] * fwd_Ea[11] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[5];
    Kc = exp(g_RT[3] - g_RT[5] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (2*h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* H */
    wdot[5] += 2 * q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[33] -= dqdci;               /* dwdot[H]/d[H] */
    J[35] += 2 * dqdci;           /* dwdot[OH]/d[H] */
    J[36] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2*sc[5];
    J[53] -= dqdci;               /* dwdot[H]/d[OH] */
    J[55] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[56] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[3];
    J[63] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[65] += 2 * dqdci;           /* dwdot[OH]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[93] -= dqdT;                /* dwdot[H]/dT */
    J[95] += 2 * dqdT;            /* dwdot[OH]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 13: HO2 + O <=> O2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[6];
    k_f = prefactor_units[12] * fwd_A[12]
                * exp(fwd_beta[12] * tc[0] - activation_units[12] * fwd_Ea[12] * invT);
    dlnkfdT = fwd_beta[12] * invT + activation_units[12] * fwd_Ea[12] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[5];
    Kc = exp(-g_RT[1] + g_RT[4] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[6]) + (h_RT[1] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[4] -= q; /* O */
    wdot[5] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[5];
    J[11] += dqdci;               /* dwdot[O2]/d[O2] */
    J[14] -= dqdci;               /* dwdot[O]/d[O2] */
    J[15] += dqdci;               /* dwdot[OH]/d[O2] */
    J[16] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[6];
    J[41] += dqdci;               /* dwdot[O2]/d[O] */
    J[44] -= dqdci;               /* dwdot[O]/d[O] */
    J[45] += dqdci;               /* dwdot[OH]/d[O] */
    J[46] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1];
    J[51] += dqdci;               /* dwdot[O2]/d[OH] */
    J[54] -= dqdci;               /* dwdot[O]/d[OH] */
    J[55] += dqdci;               /* dwdot[OH]/d[OH] */
    J[56] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[4];
    J[61] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[64] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[65] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] += dqdT;                /* dwdot[O2]/dT */
    J[94] -= dqdT;                /* dwdot[O]/dT */
    J[95] += dqdT;                /* dwdot[OH]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 14: HO2 + OH <=> H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[6];
    k_f = prefactor_units[13] * fwd_A[13]
                * exp(fwd_beta[13] * tc[0] - activation_units[13] * fwd_Ea[13] * invT);
    dlnkfdT = fwd_beta[13] * invT + activation_units[13] * fwd_Ea[13] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[2];
    Kc = exp(-g_RT[1] - g_RT[2] + g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[6]) + (h_RT[1] + h_RT[2]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[2] += q; /* H2O */
    wdot[5] -= q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[2];
    J[11] += dqdci;               /* dwdot[O2]/d[O2] */
    J[12] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[15] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[16] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[1];
    J[21] += dqdci;               /* dwdot[O2]/d[H2O] */
    J[22] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[25] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[26] -= dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[6];
    J[51] += dqdci;               /* dwdot[O2]/d[OH] */
    J[52] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[55] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[56] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[5];
    J[61] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[62] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[65] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] += dqdT;                /* dwdot[O2]/dT */
    J[92] += dqdT;                /* dwdot[H2O]/dT */
    J[95] -= dqdT;                /* dwdot[OH]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[6];
    k_f = prefactor_units[14] * fwd_A[14]
                * exp(fwd_beta[14] * tc[0] - activation_units[14] * fwd_Ea[14] * invT);
    dlnkfdT = fwd_beta[14] * invT + activation_units[14] * fwd_Ea[14] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[7];
    Kc = exp(-g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[6]) + (h_RT[1] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[11] += dqdci;               /* dwdot[O2]/d[O2] */
    J[16] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[17] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2*sc[6];
    J[61] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[67] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[1];
    J[71] += dqdci;               /* dwdot[O2]/d[H2O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[H2O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[91] += dqdT;                /* dwdot[O2]/dT */
    J[96] += -2 * dqdT;           /* dwdot[HO2]/dT */
    J[97] += dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 16: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[6];
    k_f = prefactor_units[15] * fwd_A[15]
                * exp(fwd_beta[15] * tc[0] - activation_units[15] * fwd_Ea[15] * invT);
    dlnkfdT = fwd_beta[15] * invT + activation_units[15] * fwd_Ea[15] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[7];
    Kc = exp(-g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[6]) + (h_RT[1] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[11] += dqdci;               /* dwdot[O2]/d[O2] */
    J[16] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[17] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2*sc[6];
    J[61] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[67] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[1];
    J[71] += dqdci;               /* dwdot[O2]/d[H2O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[H2O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[91] += dqdT;                /* dwdot[O2]/dT */
    J[96] += -2 * dqdT;           /* dwdot[HO2]/dT */
    J[97] += dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 17: H2O2 + H <=> H2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = prefactor_units[16] * fwd_A[16]
                * exp(fwd_beta[16] * tc[0] - activation_units[16] * fwd_Ea[16] * invT);
    dlnkfdT = fwd_beta[16] * invT + activation_units[16] * fwd_Ea[16] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[5];
    Kc = exp(-g_RT[2] + g_RT[3] - g_RT[5] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[2] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H2O */
    wdot[3] -= q; /* H */
    wdot[5] += q; /* OH */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[5];
    J[22] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[23] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[25] += dqdci;               /* dwdot[OH]/d[H2O] */
    J[27] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[32] += dqdci;               /* dwdot[H2O]/d[H] */
    J[33] -= dqdci;               /* dwdot[H]/d[H] */
    J[35] += dqdci;               /* dwdot[OH]/d[H] */
    J[37] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[2];
    J[52] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[53] -= dqdci;               /* dwdot[H]/d[OH] */
    J[55] += dqdci;               /* dwdot[OH]/d[OH] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[72] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[73] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[75] += dqdci;               /* dwdot[OH]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[92] += dqdT;                /* dwdot[H2O]/dT */
    J[93] -= dqdT;                /* dwdot[H]/dT */
    J[95] += dqdT;                /* dwdot[OH]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = prefactor_units[17] * fwd_A[17]
                * exp(fwd_beta[17] * tc[0] - activation_units[17] * fwd_Ea[17] * invT);
    dlnkfdT = fwd_beta[17] * invT + activation_units[17] * fwd_Ea[17] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[6];
    Kc = exp(-g_RT[0] + g_RT[3] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[0] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[3] -= q; /* H */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[6];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[3] -= dqdci;                /* dwdot[H]/d[H2] */
    J[6] += dqdci;                /* dwdot[HO2]/d[H2] */
    J[7] -= dqdci;                /* dwdot[H2O2]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[30] += dqdci;               /* dwdot[H2]/d[H] */
    J[33] -= dqdci;               /* dwdot[H]/d[H] */
    J[36] += dqdci;               /* dwdot[HO2]/d[H] */
    J[37] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[0];
    J[60] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[63] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[70] += dqdci;               /* dwdot[H2]/d[H2O2] */
    J[73] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[90] += dqdT;                /* dwdot[H2]/dT */
    J[93] -= dqdT;                /* dwdot[H]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[7];
    k_f = prefactor_units[18] * fwd_A[18]
                * exp(fwd_beta[18] * tc[0] - activation_units[18] * fwd_Ea[18] * invT);
    dlnkfdT = fwd_beta[18] * invT + activation_units[18] * fwd_Ea[18] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[6];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[7]) + (h_RT[5] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* O */
    wdot[5] += q; /* OH */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[7];
    J[44] -= dqdci;               /* dwdot[O]/d[O] */
    J[45] += dqdci;               /* dwdot[OH]/d[O] */
    J[46] += dqdci;               /* dwdot[HO2]/d[O] */
    J[47] -= dqdci;               /* dwdot[H2O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[6];
    J[54] -= dqdci;               /* dwdot[O]/d[OH] */
    J[55] += dqdci;               /* dwdot[OH]/d[OH] */
    J[56] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[5];
    J[64] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[65] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[4];
    J[74] -= dqdci;               /* dwdot[O]/d[H2O2] */
    J[75] += dqdci;               /* dwdot[OH]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[94] -= dqdT;                /* dwdot[O]/dT */
    J[95] += dqdT;                /* dwdot[OH]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[7];
    k_f = prefactor_units[19] * fwd_A[19]
                * exp(fwd_beta[19] * tc[0] - activation_units[19] * fwd_Ea[19] * invT);
    dlnkfdT = fwd_beta[19] * invT + activation_units[19] * fwd_Ea[19] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[6];
    Kc = exp(-g_RT[2] + g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[2] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H2O */
    wdot[5] -= q; /* OH */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[22] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[25] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[26] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[27] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[52] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[55] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[56] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[2];
    J[62] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[65] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[5];
    J[72] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[75] -= dqdci;               /* dwdot[OH]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[92] += dqdT;                /* dwdot[H2O]/dT */
    J[95] -= dqdT;                /* dwdot[OH]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[7];
    k_f = prefactor_units[20] * fwd_A[20]
                * exp(fwd_beta[20] * tc[0] - activation_units[20] * fwd_Ea[20] * invT);
    dlnkfdT = fwd_beta[20] * invT + activation_units[20] * fwd_Ea[20] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[6];
    Kc = exp(-g_RT[2] + g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[2] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H2O */
    wdot[5] -= q; /* OH */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[22] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[25] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[26] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[27] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[52] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[55] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[56] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[2];
    J[62] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[65] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[5];
    J[72] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[75] -= dqdci;               /* dwdot[OH]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[92] += dqdT;                /* dwdot[H2O]/dT */
    J[95] -= dqdT;                /* dwdot[OH]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    double c_R[9], dcRdT[9], e_RT[9];
    double * eh_RT;
    if (consP) {
        cp_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        eh_RT = &h_RT[0];
    }
    else {
        cv_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        speciesInternalEnergy(e_RT, tc);
        eh_RT = &e_RT[0];
    }

    double cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;
    for (int k = 0; k < 9; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[90+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 9; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 9; ++m) {
            dehmixdc += eh_RT[m]*J[k*10+m];
        }
        J[k*10+9] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[99] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}


/*compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void dcvpRdT(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=300 kelvin */
    if (T < 300) {
        /*species 0: H2 */
        species[0] =
            +8.24944174e-04
            -1.62860306e-06 * tc[1]
            -2.84263030e-10 * tc[2]
            +1.65394890e-12 * tc[3];
        /*species 1: O2 */
        species[1] =
            +1.12748635e-03
            -1.15123009e-06 * tc[1]
            +3.94163169e-09 * tc[2]
            -3.50742157e-12 * tc[3];
        /*species 2: H2O */
        species[2] =
            +3.47498246e-03
            -1.27093927e-05 * tc[1]
            +2.09057438e-08 * tc[2]
            -1.00263539e-11 * tc[3];
        /*species 3: H */
        species[3] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 4: O */
        species[4] =
            -1.63816649e-03
            +4.84206340e-06 * tc[1]
            -4.80852957e-09 * tc[2]
            +1.55627854e-12 * tc[3];
        /*species 5: OH */
        species[5] =
            -3.22544939e-03
            +1.30552938e-05 * tc[1]
            -1.73956093e-08 * tc[2]
            +8.24949516e-12 * tc[3];
        /*species 6: HO2 */
        species[6] =
            -4.74912051e-03
            +4.23165782e-05 * tc[1]
            -7.28291682e-08 * tc[2]
            +3.71690050e-11 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            +6.56922581e-03
            -2.97002516e-07 * tc[1]
            -1.38774166e-08 * tc[2]
            +9.88605900e-12 * tc[3];
        /*species 8: N2 */
        species[8] =
            +1.40824000e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77942000e-12 * tc[3];
    } else {
        /*species 0: H2 */
        species[0] =
            +7.00064411e-04
            -1.12676574e-07 * tc[1]
            -2.76947345e-11 * tc[2]
            +6.33100716e-15 * tc[3];
        /*species 1: O2 */
        species[1] =
            +6.13519689e-04
            -2.51768398e-07 * tc[1]
            +5.32584444e-11 * tc[2]
            -4.54574124e-15 * tc[3];
        /*species 2: H2O */
        species[2] =
            +3.05629289e-03
            -1.74605202e-06 * tc[1]
            +3.60298917e-10 * tc[2]
            -2.55664715e-14 * tc[3];
        /*species 3: H */
        species[3] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 4: O */
        species[4] =
            -2.75506191e-05
            -6.20560670e-09 * tc[1]
            +1.36532023e-11 * tc[2]
            -1.74722060e-15 * tc[3];
        /*species 5: OH */
        species[5] =
            +1.05650448e-03
            -5.18165516e-07 * tc[1]
            +9.15656022e-11 * tc[2]
            -5.32783504e-15 * tc[3];
        /*species 6: HO2 */
        species[6] =
            +2.23982013e-03
            -1.26731630e-06 * tc[1]
            +3.42739110e-10 * tc[2]
            -4.31634140e-14 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            +4.33613639e-03
            -2.94937764e-06 * tc[1]
            +7.04671071e-10 * tc[2]
            -5.72661424e-14 * tc[3];
        /*species 8: N2 */
        species[8] =
            +1.48797700e-03
            -1.13695220e-06 * tc[1]
            +3.02911200e-10 * tc[2]
            -2.70134040e-14 * tc[3];
    }
    return;
}


/*compute the progress rate for each reaction */
void progressRate(double * restrict qdot, double * restrict sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }

    double q_f[21], q_r[21];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 21; ++i) {
        qdot[i] = q_f[i] - q_r[i];
    }

    return;
}


/*compute the progress rate for each reaction */
void progressRateFR(double * restrict q_f, double * restrict q_r, double * restrict sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }

    comp_qfqr(q_f, q_r, sc, tc, invT);

    return;
}


/*compute the equilibrium constants for each reaction */
void equilibriumConstants(double * restrict kc, double * restrict g_RT, double T)
{
    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    kc[0] = 1.0 / (refC) * exp((g_RT[3] + g_RT[1]) - (g_RT[6]));

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    kc[1] = refC * exp((g_RT[7]) - (g_RT[5] + g_RT[5]));

    /*reaction 3: H2 + M <=> H + H + M */
    kc[2] = refC * exp((g_RT[0]) - (g_RT[3] + g_RT[3]));

    /*reaction 4: O + O + M <=> O2 + M */
    kc[3] = 1.0 / (refC) * exp((g_RT[4] + g_RT[4]) - (g_RT[1]));

    /*reaction 5: O + H + M <=> OH + M */
    kc[4] = 1.0 / (refC) * exp((g_RT[4] + g_RT[3]) - (g_RT[5]));

    /*reaction 6: H + OH + M <=> H2O + M */
    kc[5] = 1.0 / (refC) * exp((g_RT[3] + g_RT[5]) - (g_RT[2]));

    /*reaction 7: H + O2 <=> O + OH */
    kc[6] = exp((g_RT[3] + g_RT[1]) - (g_RT[4] + g_RT[5]));

    /*reaction 8: O + H2 <=> H + OH */
    kc[7] = exp((g_RT[4] + g_RT[0]) - (g_RT[3] + g_RT[5]));

    /*reaction 9: H2 + OH <=> H2O + H */
    kc[8] = exp((g_RT[0] + g_RT[5]) - (g_RT[2] + g_RT[3]));

    /*reaction 10: O + H2O <=> OH + OH */
    kc[9] = exp((g_RT[4] + g_RT[2]) - (g_RT[5] + g_RT[5]));

    /*reaction 11: HO2 + H <=> H2 + O2 */
    kc[10] = exp((g_RT[6] + g_RT[3]) - (g_RT[0] + g_RT[1]));

    /*reaction 12: HO2 + H <=> OH + OH */
    kc[11] = exp((g_RT[6] + g_RT[3]) - (g_RT[5] + g_RT[5]));

    /*reaction 13: HO2 + O <=> O2 + OH */
    kc[12] = exp((g_RT[6] + g_RT[4]) - (g_RT[1] + g_RT[5]));

    /*reaction 14: HO2 + OH <=> H2O + O2 */
    kc[13] = exp((g_RT[6] + g_RT[5]) - (g_RT[2] + g_RT[1]));

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    kc[14] = exp((g_RT[6] + g_RT[6]) - (g_RT[7] + g_RT[1]));

    /*reaction 16: HO2 + HO2 <=> H2O2 + O2 */
    kc[15] = exp((g_RT[6] + g_RT[6]) - (g_RT[7] + g_RT[1]));

    /*reaction 17: H2O2 + H <=> H2O + OH */
    kc[16] = exp((g_RT[7] + g_RT[3]) - (g_RT[2] + g_RT[5]));

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    kc[17] = exp((g_RT[7] + g_RT[3]) - (g_RT[6] + g_RT[0]));

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    kc[18] = exp((g_RT[7] + g_RT[4]) - (g_RT[5] + g_RT[6]));

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    kc[19] = exp((g_RT[7] + g_RT[5]) - (g_RT[6] + g_RT[2]));

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    kc[20] = exp((g_RT[7] + g_RT[5]) - (g_RT[6] + g_RT[2]));

    return;
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void gibbs(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=300 kelvin */
    if (T < 300) {
        /*species 0: H2 */
        species[0] =
            -1.012520870000000e+03 * invT
            +6.592218400000000e+00
            -3.298124310000000e+00 * tc[0]
            -4.124720870000000e-04 * tc[1]
            +1.357169215000000e-07 * tc[2]
            +7.896195275000000e-12 * tc[3]
            -2.067436120000000e-14 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.005249020000000e+03 * invT
            -2.821801190000000e+00
            -3.212936400000000e+00 * tc[0]
            -5.637431750000000e-04 * tc[1]
            +9.593584116666666e-08 * tc[2]
            -1.094897691666667e-10 * tc[3]
            +4.384276960000000e-14 * tc[4];
        /*species 2: H2O */
        species[2] =
            -3.020811330000000e+04 * invT
            +7.966096399999998e-01
            -3.386842490000000e+00 * tc[0]
            -1.737491230000000e-03 * tc[1]
            +1.059116055000000e-06 * tc[2]
            -5.807151058333333e-10 * tc[3]
            +1.253294235000000e-13 * tc[4];
        /*species 3: H */
        species[3] =
            +2.547162700000000e+04 * invT
            +2.960117608000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.914764450000000e+04 * invT
            -1.756619999999964e-02
            -2.946428780000000e+00 * tc[0]
            +8.190832450000000e-04 * tc[1]
            -4.035052833333333e-07 * tc[2]
            +1.335702658333333e-10 * tc[3]
            -1.945348180000000e-14 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.346309130000000e+03 * invT
            +4.815738570000000e+00
            -4.125305610000000e+00 * tc[0]
            +1.612724695000000e-03 * tc[1]
            -1.087941151666667e-06 * tc[2]
            +4.832113691666666e-10 * tc[3]
            -1.031186895000000e-13 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +2.948080400000000e+02 * invT
            +5.851355599999999e-01
            -4.301798010000000e+00 * tc[0]
            +2.374560255000000e-03 * tc[1]
            -3.526381516666666e-06 * tc[2]
            +2.023032450000000e-09 * tc[3]
            -4.646125620000001e-13 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.766314650000000e+04 * invT
            -3.396609550000000e+00
            -3.388753650000000e+00 * tc[0]
            -3.284612905000000e-03 * tc[1]
            +2.475020966666667e-08 * tc[2]
            +3.854837933333333e-10 * tc[3]
            -1.235757375000000e-13 * tc[4];
        /*species 8: N2 */
        species[8] =
            -1.020900000000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041200000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427500000000e-13 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            -8.350339970000000e+02 * invT
            +4.346533540000000e+00
            -2.991423370000000e+00 * tc[0]
            -3.500322055000000e-04 * tc[1]
            +9.389714483333333e-09 * tc[2]
            +7.692981816666667e-13 * tc[3]
            -7.913758950000000e-17 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.233930180000000e+03 * invT
            +5.084126000000002e-01
            -3.697578190000000e+00 * tc[0]
            -3.067598445000000e-04 * tc[1]
            +2.098069983333333e-08 * tc[2]
            -1.479401233333333e-12 * tc[3]
            +5.682176550000000e-17 * tc[4];
        /*species 2: H2O */
        species[2] =
            -2.989920900000000e+04 * invT
            -4.190671200000001e+00
            -2.672145610000000e+00 * tc[0]
            -1.528146445000000e-03 * tc[1]
            +1.455043351666667e-07 * tc[2]
            -1.000830325000000e-11 * tc[3]
            +3.195808935000000e-16 * tc[4];
        /*species 3: H */
        species[3] =
            +2.547162700000000e+04 * invT
            +2.960117638000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.923080270000000e+04 * invT
            -2.378248450000000e+00
            -2.542059660000000e+00 * tc[0]
            +1.377530955000000e-05 * tc[1]
            +5.171338916666667e-10 * tc[2]
            -3.792556183333333e-13 * tc[3]
            +2.184025750000000e-17 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.683628750000000e+03 * invT
            -2.836911870000000e+00
            -2.864728860000000e+00 * tc[0]
            -5.282522400000000e-04 * tc[1]
            +4.318045966666667e-08 * tc[2]
            -2.543488950000000e-12 * tc[3]
            +6.659793800000000e-17 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +1.118567130000000e+02 * invT
            +2.321087500000001e-01
            -4.017210900000000e+00 * tc[0]
            -1.119910065000000e-03 * tc[1]
            +1.056096916666667e-07 * tc[2]
            -9.520530833333334e-12 * tc[3]
            +5.395426750000000e-16 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.800696090000000e+04 * invT
            +4.072029891000000e+00
            -4.573166850000000e+00 * tc[0]
            -2.168068195000000e-03 * tc[1]
            +2.457814700000000e-07 * tc[2]
            -1.957419641666667e-11 * tc[3]
            +7.158267800000000e-16 * tc[4];
        /*species 8: N2 */
        species[8] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439885000000000e-04 * tc[1]
            +9.474601666666666e-08 * tc[2]
            -8.414199999999999e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
    }
    return;
}


/*compute the a/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void helmholtz(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=300 kelvin */
    if (T < 300) {
        /*species 0: H2 */
        species[0] =
            -1.01252087e+03 * invT
            +5.59221840e+00
            -3.29812431e+00 * tc[0]
            -4.12472087e-04 * tc[1]
            +1.35716922e-07 * tc[2]
            +7.89619527e-12 * tc[3]
            -2.06743612e-14 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.00524902e+03 * invT
            -3.82180119e+00
            -3.21293640e+00 * tc[0]
            -5.63743175e-04 * tc[1]
            +9.59358412e-08 * tc[2]
            -1.09489769e-10 * tc[3]
            +4.38427696e-14 * tc[4];
        /*species 2: H2O */
        species[2] =
            -3.02081133e+04 * invT
            -2.03390360e-01
            -3.38684249e+00 * tc[0]
            -1.73749123e-03 * tc[1]
            +1.05911606e-06 * tc[2]
            -5.80715106e-10 * tc[3]
            +1.25329424e-13 * tc[4];
        /*species 3: H */
        species[3] =
            +2.54716270e+04 * invT
            +1.96011761e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.91476445e+04 * invT
            -1.01756620e+00
            -2.94642878e+00 * tc[0]
            +8.19083245e-04 * tc[1]
            -4.03505283e-07 * tc[2]
            +1.33570266e-10 * tc[3]
            -1.94534818e-14 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.34630913e+03 * invT
            +3.81573857e+00
            -4.12530561e+00 * tc[0]
            +1.61272470e-03 * tc[1]
            -1.08794115e-06 * tc[2]
            +4.83211369e-10 * tc[3]
            -1.03118689e-13 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +2.94808040e+02 * invT
            -4.14864440e-01
            -4.30179801e+00 * tc[0]
            +2.37456025e-03 * tc[1]
            -3.52638152e-06 * tc[2]
            +2.02303245e-09 * tc[3]
            -4.64612562e-13 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.76631465e+04 * invT
            -4.39660955e+00
            -3.38875365e+00 * tc[0]
            -3.28461290e-03 * tc[1]
            +2.47502097e-08 * tc[2]
            +3.85483793e-10 * tc[3]
            -1.23575738e-13 * tc[4];
        /*species 8: N2 */
        species[8] =
            -1.02090000e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120000e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242750e-13 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            -8.35033997e+02 * invT
            +3.34653354e+00
            -2.99142337e+00 * tc[0]
            -3.50032206e-04 * tc[1]
            +9.38971448e-09 * tc[2]
            +7.69298182e-13 * tc[3]
            -7.91375895e-17 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.23393018e+03 * invT
            -4.91587400e-01
            -3.69757819e+00 * tc[0]
            -3.06759845e-04 * tc[1]
            +2.09806998e-08 * tc[2]
            -1.47940123e-12 * tc[3]
            +5.68217655e-17 * tc[4];
        /*species 2: H2O */
        species[2] =
            -2.98992090e+04 * invT
            -5.19067120e+00
            -2.67214561e+00 * tc[0]
            -1.52814644e-03 * tc[1]
            +1.45504335e-07 * tc[2]
            -1.00083033e-11 * tc[3]
            +3.19580894e-16 * tc[4];
        /*species 3: H */
        species[3] =
            +2.54716270e+04 * invT
            +1.96011764e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.92308027e+04 * invT
            -3.37824845e+00
            -2.54205966e+00 * tc[0]
            +1.37753096e-05 * tc[1]
            +5.17133892e-10 * tc[2]
            -3.79255618e-13 * tc[3]
            +2.18402575e-17 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.68362875e+03 * invT
            -3.83691187e+00
            -2.86472886e+00 * tc[0]
            -5.28252240e-04 * tc[1]
            +4.31804597e-08 * tc[2]
            -2.54348895e-12 * tc[3]
            +6.65979380e-17 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +1.11856713e+02 * invT
            -7.67891250e-01
            -4.01721090e+00 * tc[0]
            -1.11991006e-03 * tc[1]
            +1.05609692e-07 * tc[2]
            -9.52053083e-12 * tc[3]
            +5.39542675e-16 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.80069609e+04 * invT
            +3.07202989e+00
            -4.57316685e+00 * tc[0]
            -2.16806820e-03 * tc[1]
            +2.45781470e-07 * tc[2]
            -1.95741964e-11 * tc[3]
            +7.15826780e-16 * tc[4];
        /*species 8: N2 */
        species[8] =
            -9.22797700e+02 * invT
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988500e-04 * tc[1]
            +9.47460167e-08 * tc[2]
            -8.41420000e-12 * tc[3]
            +3.37667550e-16 * tc[4];
    }
    return;
}


/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void cv_R(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=300 kelvin */
    if (T < 300) {
        /*species 0: H2 */
        species[0] =
            +2.29812431e+00
            +8.24944174e-04 * tc[1]
            -8.14301529e-07 * tc[2]
            -9.47543433e-11 * tc[3]
            +4.13487224e-13 * tc[4];
        /*species 1: O2 */
        species[1] =
            +2.21293640e+00
            +1.12748635e-03 * tc[1]
            -5.75615047e-07 * tc[2]
            +1.31387723e-09 * tc[3]
            -8.76855392e-13 * tc[4];
        /*species 2: H2O */
        species[2] =
            +2.38684249e+00
            +3.47498246e-03 * tc[1]
            -6.35469633e-06 * tc[2]
            +6.96858127e-09 * tc[3]
            -2.50658847e-12 * tc[4];
        /*species 3: H */
        species[3] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +1.94642878e+00
            -1.63816649e-03 * tc[1]
            +2.42103170e-06 * tc[2]
            -1.60284319e-09 * tc[3]
            +3.89069636e-13 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +2.38875365e+00
            +6.56922581e-03 * tc[1]
            -1.48501258e-07 * tc[2]
            -4.62580552e-09 * tc[3]
            +2.47151475e-12 * tc[4];
        /*species 8: N2 */
        species[8] =
            +2.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            +1.99142337e+00
            +7.00064411e-04 * tc[1]
            -5.63382869e-08 * tc[2]
            -9.23157818e-12 * tc[3]
            +1.58275179e-15 * tc[4];
        /*species 1: O2 */
        species[1] =
            +2.69757819e+00
            +6.13519689e-04 * tc[1]
            -1.25884199e-07 * tc[2]
            +1.77528148e-11 * tc[3]
            -1.13643531e-15 * tc[4];
        /*species 2: H2O */
        species[2] =
            +1.67214561e+00
            +3.05629289e-03 * tc[1]
            -8.73026011e-07 * tc[2]
            +1.20099639e-10 * tc[3]
            -6.39161787e-15 * tc[4];
        /*species 3: H */
        species[3] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +1.54205966e+00
            -2.75506191e-05 * tc[1]
            -3.10280335e-09 * tc[2]
            +4.55106742e-12 * tc[3]
            -4.36805150e-16 * tc[4];
        /*species 5: OH */
        species[5] =
            +1.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.57316685e+00
            +4.33613639e-03 * tc[1]
            -1.47468882e-06 * tc[2]
            +2.34890357e-10 * tc[3]
            -1.43165356e-14 * tc[4];
        /*species 8: N2 */
        species[8] =
            +1.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }
    return;
}


/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void cp_R(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=300 kelvin */
    if (T < 300) {
        /*species 0: H2 */
        species[0] =
            +3.29812431e+00
            +8.24944174e-04 * tc[1]
            -8.14301529e-07 * tc[2]
            -9.47543433e-11 * tc[3]
            +4.13487224e-13 * tc[4];
        /*species 1: O2 */
        species[1] =
            +3.21293640e+00
            +1.12748635e-03 * tc[1]
            -5.75615047e-07 * tc[2]
            +1.31387723e-09 * tc[3]
            -8.76855392e-13 * tc[4];
        /*species 2: H2O */
        species[2] =
            +3.38684249e+00
            +3.47498246e-03 * tc[1]
            -6.35469633e-06 * tc[2]
            +6.96858127e-09 * tc[3]
            -2.50658847e-12 * tc[4];
        /*species 3: H */
        species[3] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.94642878e+00
            -1.63816649e-03 * tc[1]
            +2.42103170e-06 * tc[2]
            -1.60284319e-09 * tc[3]
            +3.89069636e-13 * tc[4];
        /*species 5: OH */
        species[5] =
            +4.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.38875365e+00
            +6.56922581e-03 * tc[1]
            -1.48501258e-07 * tc[2]
            -4.62580552e-09 * tc[3]
            +2.47151475e-12 * tc[4];
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            +2.99142337e+00
            +7.00064411e-04 * tc[1]
            -5.63382869e-08 * tc[2]
            -9.23157818e-12 * tc[3]
            +1.58275179e-15 * tc[4];
        /*species 1: O2 */
        species[1] =
            +3.69757819e+00
            +6.13519689e-04 * tc[1]
            -1.25884199e-07 * tc[2]
            +1.77528148e-11 * tc[3]
            -1.13643531e-15 * tc[4];
        /*species 2: H2O */
        species[2] =
            +2.67214561e+00
            +3.05629289e-03 * tc[1]
            -8.73026011e-07 * tc[2]
            +1.20099639e-10 * tc[3]
            -6.39161787e-15 * tc[4];
        /*species 3: H */
        species[3] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.54205966e+00
            -2.75506191e-05 * tc[1]
            -3.10280335e-09 * tc[2]
            +4.55106742e-12 * tc[3]
            -4.36805150e-16 * tc[4];
        /*species 5: OH */
        species[5] =
            +2.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +4.57316685e+00
            +4.33613639e-03 * tc[1]
            -1.47468882e-06 * tc[2]
            +2.34890357e-10 * tc[3]
            -1.43165356e-14 * tc[4];
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }
    return;
}


/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesInternalEnergy(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=300 kelvin */
    if (T < 300) {
        /*species 0: H2 */
        species[0] =
            +2.29812431e+00
            +4.12472087e-04 * tc[1]
            -2.71433843e-07 * tc[2]
            -2.36885858e-11 * tc[3]
            +8.26974448e-14 * tc[4]
            -1.01252087e+03 * invT;
        /*species 1: O2 */
        species[1] =
            +2.21293640e+00
            +5.63743175e-04 * tc[1]
            -1.91871682e-07 * tc[2]
            +3.28469308e-10 * tc[3]
            -1.75371078e-13 * tc[4]
            -1.00524902e+03 * invT;
        /*species 2: H2O */
        species[2] =
            +2.38684249e+00
            +1.73749123e-03 * tc[1]
            -2.11823211e-06 * tc[2]
            +1.74214532e-09 * tc[3]
            -5.01317694e-13 * tc[4]
            -3.02081133e+04 * invT;
        /*species 3: H */
        species[3] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716270e+04 * invT;
        /*species 4: O */
        species[4] =
            +1.94642878e+00
            -8.19083245e-04 * tc[1]
            +8.07010567e-07 * tc[2]
            -4.00710797e-10 * tc[3]
            +7.78139272e-14 * tc[4]
            +2.91476445e+04 * invT;
        /*species 5: OH */
        species[5] =
            +3.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.34630913e+03 * invT;
        /*species 6: HO2 */
        species[6] =
            +3.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 7: H2O2 */
        species[7] =
            +2.38875365e+00
            +3.28461290e-03 * tc[1]
            -4.95004193e-08 * tc[2]
            -1.15645138e-09 * tc[3]
            +4.94302950e-13 * tc[4]
            -1.76631465e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +2.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
    } else {
        /*species 0: H2 */
        species[0] =
            +1.99142337e+00
            +3.50032206e-04 * tc[1]
            -1.87794290e-08 * tc[2]
            -2.30789455e-12 * tc[3]
            +3.16550358e-16 * tc[4]
            -8.35033997e+02 * invT;
        /*species 1: O2 */
        species[1] =
            +2.69757819e+00
            +3.06759845e-04 * tc[1]
            -4.19613997e-08 * tc[2]
            +4.43820370e-12 * tc[3]
            -2.27287062e-16 * tc[4]
            -1.23393018e+03 * invT;
        /*species 2: H2O */
        species[2] =
            +1.67214561e+00
            +1.52814644e-03 * tc[1]
            -2.91008670e-07 * tc[2]
            +3.00249098e-11 * tc[3]
            -1.27832357e-15 * tc[4]
            -2.98992090e+04 * invT;
        /*species 3: H */
        species[3] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716270e+04 * invT;
        /*species 4: O */
        species[4] =
            +1.54205966e+00
            -1.37753096e-05 * tc[1]
            -1.03426778e-09 * tc[2]
            +1.13776685e-12 * tc[3]
            -8.73610300e-17 * tc[4]
            +2.92308027e+04 * invT;
        /*species 5: OH */
        species[5] =
            +1.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.68362875e+03 * invT;
        /*species 6: HO2 */
        species[6] =
            +3.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 7: H2O2 */
        species[7] =
            +3.57316685e+00
            +2.16806820e-03 * tc[1]
            -4.91562940e-07 * tc[2]
            +5.87225893e-11 * tc[3]
            -2.86330712e-15 * tc[4]
            -1.80069609e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +1.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
    }
    return;
}


/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesEnthalpy(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=300 kelvin */
    if (T < 300) {
        /*species 0: H2 */
        species[0] =
            +3.29812431e+00
            +4.12472087e-04 * tc[1]
            -2.71433843e-07 * tc[2]
            -2.36885858e-11 * tc[3]
            +8.26974448e-14 * tc[4]
            -1.01252087e+03 * invT;
        /*species 1: O2 */
        species[1] =
            +3.21293640e+00
            +5.63743175e-04 * tc[1]
            -1.91871682e-07 * tc[2]
            +3.28469308e-10 * tc[3]
            -1.75371078e-13 * tc[4]
            -1.00524902e+03 * invT;
        /*species 2: H2O */
        species[2] =
            +3.38684249e+00
            +1.73749123e-03 * tc[1]
            -2.11823211e-06 * tc[2]
            +1.74214532e-09 * tc[3]
            -5.01317694e-13 * tc[4]
            -3.02081133e+04 * invT;
        /*species 3: H */
        species[3] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716270e+04 * invT;
        /*species 4: O */
        species[4] =
            +2.94642878e+00
            -8.19083245e-04 * tc[1]
            +8.07010567e-07 * tc[2]
            -4.00710797e-10 * tc[3]
            +7.78139272e-14 * tc[4]
            +2.91476445e+04 * invT;
        /*species 5: OH */
        species[5] =
            +4.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.34630913e+03 * invT;
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 7: H2O2 */
        species[7] =
            +3.38875365e+00
            +3.28461290e-03 * tc[1]
            -4.95004193e-08 * tc[2]
            -1.15645138e-09 * tc[3]
            +4.94302950e-13 * tc[4]
            -1.76631465e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
    } else {
        /*species 0: H2 */
        species[0] =
            +2.99142337e+00
            +3.50032206e-04 * tc[1]
            -1.87794290e-08 * tc[2]
            -2.30789455e-12 * tc[3]
            +3.16550358e-16 * tc[4]
            -8.35033997e+02 * invT;
        /*species 1: O2 */
        species[1] =
            +3.69757819e+00
            +3.06759845e-04 * tc[1]
            -4.19613997e-08 * tc[2]
            +4.43820370e-12 * tc[3]
            -2.27287062e-16 * tc[4]
            -1.23393018e+03 * invT;
        /*species 2: H2O */
        species[2] =
            +2.67214561e+00
            +1.52814644e-03 * tc[1]
            -2.91008670e-07 * tc[2]
            +3.00249098e-11 * tc[3]
            -1.27832357e-15 * tc[4]
            -2.98992090e+04 * invT;
        /*species 3: H */
        species[3] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716270e+04 * invT;
        /*species 4: O */
        species[4] =
            +2.54205966e+00
            -1.37753096e-05 * tc[1]
            -1.03426778e-09 * tc[2]
            +1.13776685e-12 * tc[3]
            -8.73610300e-17 * tc[4]
            +2.92308027e+04 * invT;
        /*species 5: OH */
        species[5] =
            +2.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.68362875e+03 * invT;
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 7: H2O2 */
        species[7] =
            +4.57316685e+00
            +2.16806820e-03 * tc[1]
            -4.91562940e-07 * tc[2]
            +5.87225893e-11 * tc[3]
            -2.86330712e-15 * tc[4]
            -1.80069609e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
    }
    return;
}


/*compute the S/R at the given temperature (Eq 21) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesEntropy(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=300 kelvin */
    if (T < 300) {
        /*species 0: H2 */
        species[0] =
            +3.29812431e+00 * tc[0]
            +8.24944174e-04 * tc[1]
            -4.07150765e-07 * tc[2]
            -3.15847811e-11 * tc[3]
            +1.03371806e-13 * tc[4]
            -3.29409409e+00 ;
        /*species 1: O2 */
        species[1] =
            +3.21293640e+00 * tc[0]
            +1.12748635e-03 * tc[1]
            -2.87807523e-07 * tc[2]
            +4.37959077e-10 * tc[3]
            -2.19213848e-13 * tc[4]
            +6.03473759e+00 ;
        /*species 2: H2O */
        species[2] =
            +3.38684249e+00 * tc[0]
            +3.47498246e-03 * tc[1]
            -3.17734817e-06 * tc[2]
            +2.32286042e-09 * tc[3]
            -6.26647117e-13 * tc[4]
            +2.59023285e+00 ;
        /*species 3: H */
        species[3] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.60117608e-01 ;
        /*species 4: O */
        species[4] =
            +2.94642878e+00 * tc[0]
            -1.63816649e-03 * tc[1]
            +1.21051585e-06 * tc[2]
            -5.34281063e-10 * tc[3]
            +9.72674090e-14 * tc[4]
            +2.96399498e+00 ;
        /*species 5: OH */
        species[5] =
            +4.12530561e+00 * tc[0]
            -3.22544939e-03 * tc[1]
            +3.26382346e-06 * tc[2]
            -1.93284548e-09 * tc[3]
            +5.15593447e-13 * tc[4]
            -6.90432960e-01 ;
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00 * tc[0]
            -4.74912051e-03 * tc[1]
            +1.05791445e-05 * tc[2]
            -8.09212980e-09 * tc[3]
            +2.32306281e-12 * tc[4]
            +3.71666245e+00 ;
        /*species 7: H2O2 */
        species[7] =
            +3.38875365e+00 * tc[0]
            +6.56922581e-03 * tc[1]
            -7.42506290e-08 * tc[2]
            -1.54193517e-09 * tc[3]
            +6.17878688e-13 * tc[4]
            +6.78536320e+00 ;
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00 * tc[0]
            +1.40824000e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213750e-13 * tc[4]
            +3.95037200e+00 ;
    } else {
        /*species 0: H2 */
        species[0] =
            +2.99142337e+00 * tc[0]
            +7.00064411e-04 * tc[1]
            -2.81691434e-08 * tc[2]
            -3.07719273e-12 * tc[3]
            +3.95687948e-16 * tc[4]
            -1.35511017e+00 ;
        /*species 1: O2 */
        species[1] =
            +3.69757819e+00 * tc[0]
            +6.13519689e-04 * tc[1]
            -6.29420995e-08 * tc[2]
            +5.91760493e-12 * tc[3]
            -2.84108828e-16 * tc[4]
            +3.18916559e+00 ;
        /*species 2: H2O */
        species[2] =
            +2.67214561e+00 * tc[0]
            +3.05629289e-03 * tc[1]
            -4.36513005e-07 * tc[2]
            +4.00332130e-11 * tc[3]
            -1.59790447e-15 * tc[4]
            +6.86281681e+00 ;
        /*species 3: H */
        species[3] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.60117638e-01 ;
        /*species 4: O */
        species[4] =
            +2.54205966e+00 * tc[0]
            -2.75506191e-05 * tc[1]
            -1.55140167e-09 * tc[2]
            +1.51702247e-12 * tc[3]
            -1.09201287e-16 * tc[4]
            +4.92030811e+00 ;
        /*species 5: OH */
        species[5] =
            +2.86472886e+00 * tc[0]
            +1.05650448e-03 * tc[1]
            -1.29541379e-07 * tc[2]
            +1.01739558e-11 * tc[3]
            -3.32989690e-16 * tc[4]
            +5.70164073e+00 ;
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00 * tc[0]
            +2.23982013e-03 * tc[1]
            -3.16829075e-07 * tc[2]
            +3.80821233e-11 * tc[3]
            -2.69771337e-15 * tc[4]
            +3.78510215e+00 ;
        /*species 7: H2O2 */
        species[7] =
            +4.57316685e+00 * tc[0]
            +4.33613639e-03 * tc[1]
            -7.37344410e-07 * tc[2]
            +7.82967857e-11 * tc[3]
            -3.57913390e-15 * tc[4]
            +5.01136959e-01 ;
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00 * tc[0]
            +1.48797700e-03 * tc[1]
            -2.84238050e-07 * tc[2]
            +3.36568000e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
    }
    return;
}


/*save molecular weights into array */
void molecularWeight(double * restrict wt)
{
    wt[0] = 2.015940; /*H2 */
    wt[1] = 31.998800; /*O2 */
    wt[2] = 18.015340; /*H2O */
    wt[3] = 1.007970; /*H */
    wt[4] = 15.999400; /*O */
    wt[5] = 17.007370; /*OH */
    wt[6] = 33.006770; /*HO2 */
    wt[7] = 34.014740; /*H2O2 */
    wt[8] = 28.013400; /*N2 */

    return;
}


/*save atomic weights into array */
void atomicWeight(double * restrict awt)
{
    awt[0] = 1.007970; /*H */
    awt[1] = 15.999400; /*O */
    awt[2] = 14.006700; /*N */

    return;
}
/* get temperature given internal energy in mass units and mass fracs */
void GET_T_GIVEN_EY(double * restrict e, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict t, int * ierr)
{
#ifdef CONVERGENCE
    const int maxiter = 5000;
    const double tol  = 1.e-12;
#else
    const int maxiter = 200;
    const double tol  = 1.e-6;
#endif
    double ein  = *e;
    double tmin = 250;/*max lower bound for thermo def */
    double tmax = 4000;/*min upper bound for thermo def */
    double e1,emin,emax,cv,t1,dt;
    int i;/* loop counter */
    CKUBMS(&tmin, y, iwrk, rwrk, &emin);
    CKUBMS(&tmax, y, iwrk, rwrk, &emax);
    if (ein < emin) {
        /*Linear Extrapolation below tmin */
        CKCVBS(&tmin, y, iwrk, rwrk, &cv);
        *t = tmin - (emin-ein)/cv;
        *ierr = 1;
        return;
    }
    if (ein > emax) {
        /*Linear Extrapolation above tmax */
        CKCVBS(&tmax, y, iwrk, rwrk, &cv);
        *t = tmax - (emax-ein)/cv;
        *ierr = 1;
        return;
    }
    t1 = *t;
    if (t1 < tmin || t1 > tmax) {
        t1 = tmin + (tmax-tmin)/(emax-emin)*(ein-emin);
    }
    for (i = 0; i < maxiter; ++i) {
        CKUBMS(&t1,y,iwrk,rwrk,&e1);
        CKCVBS(&t1,y,iwrk,rwrk,&cv);
        dt = (ein - e1) / cv;
        if (dt > 100.) { dt = 100.; }
        else if (dt < -100.) { dt = -100.; }
        else if (fabs(dt) < tol) break;
        else if (t1+dt == t1) break;
        t1 += dt;
    }
    *t = t1;
    *ierr = 0;
    return;
}
/* get temperature given enthalpy in mass units and mass fracs */
void GET_T_GIVEN_HY(double * restrict h, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict t, int * ierr)
{
#ifdef CONVERGENCE
    const int maxiter = 5000;
    const double tol  = 1.e-12;
#else
    const int maxiter = 200;
    const double tol  = 1.e-6;
#endif
    double hin  = *h;
    double tmin = 250;/*max lower bound for thermo def */
    double tmax = 4000;/*min upper bound for thermo def */
    double h1,hmin,hmax,cp,t1,dt;
    int i;/* loop counter */
    CKHBMS(&tmin, y, iwrk, rwrk, &hmin);
    CKHBMS(&tmax, y, iwrk, rwrk, &hmax);
    if (hin < hmin) {
        /*Linear Extrapolation below tmin */
        CKCPBS(&tmin, y, iwrk, rwrk, &cp);
        *t = tmin - (hmin-hin)/cp;
        *ierr = 1;
        return;
    }
    if (hin > hmax) {
        /*Linear Extrapolation above tmax */
        CKCPBS(&tmax, y, iwrk, rwrk, &cp);
        *t = tmax - (hmax-hin)/cp;
        *ierr = 1;
        return;
    }
    t1 = *t;
    if (t1 < tmin || t1 > tmax) {
        t1 = tmin + (tmax-tmin)/(hmax-hmin)*(hin-hmin);
    }
    for (i = 0; i < maxiter; ++i) {
        CKHBMS(&t1,y,iwrk,rwrk,&h1);
        CKCPBS(&t1,y,iwrk,rwrk,&cp);
        dt = (hin - h1) / cp;
        if (dt > 100.) { dt = 100.; }
        else if (dt < -100.) { dt = -100.; }
        else if (fabs(dt) < tol) break;
        else if (t1+dt == t1) break;
        t1 += dt;
    }
    *t = t1;
    *ierr = 0;
    return;
}

/* End of file  */
