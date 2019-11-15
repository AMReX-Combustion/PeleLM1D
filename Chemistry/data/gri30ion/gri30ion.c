
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
#define CKCHRG CKCHRG
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
#define CKCHRG ckchrg
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
#define CKCHRG ckchrg_
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
void CKCHRG(int * restrict iwrk, double * restrict rwrk, int * restrict kcharge);
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
static const double imw[20] = {
    1.0 / 15.035060,  /*CH3 */
    1.0 / 16.043030,  /*CH4 */
    1.0 / 1.007970,  /*H */
    1.0 / 17.007370,  /*OH */
    1.0 / 18.015340,  /*H2O */
    1.0 / 15.999400,  /*O */
    1.0 / 29.017971,  /*HCOp */
    1.0 / 30.006100,  /*NO */
    1.0 / 14.006700,  /*N */
    1.0 / 2.015940,  /*H2 */
    1.0 / 31.998800,  /*O2 */
    1.0 / 30.026490,  /*CH2O */
    1.0 / 28.010550,  /*CO */
    1.0 / 44.009950,  /*CO2 */
    1.0 / 13.019120,  /*CH */
    1.0 / 19.022761,  /*H3Op */
    1.0 / 28.013400,  /*N2 */
    1.0 / 46.005500,  /*NO2 */
    1.0 / 0.000549,  /*E */
    1.0 / 14.027090};  /*CH2 */



static double fwd_A[31], fwd_beta[31], fwd_Ea[31];
static double low_A[31], low_beta[31], low_Ea[31];
static double rev_A[31], rev_beta[31], rev_Ea[31];
static double troe_a[31],troe_Ts[31], troe_Tss[31], troe_Tsss[31];
static double sri_a[31], sri_b[31], sri_c[31], sri_d[31], sri_e[31];
static double activation_units[31], prefactor_units[31], phase_units[31];
static int is_PD[31], troe_len[31], sri_len[31], nTB[31], *TBid[31];
static double *TB[31];

static double fwd_A_DEF[31], fwd_beta_DEF[31], fwd_Ea_DEF[31];
static double low_A_DEF[31], low_beta_DEF[31], low_Ea_DEF[31];
static double rev_A_DEF[31], rev_beta_DEF[31], rev_Ea_DEF[31];
static double troe_a_DEF[31],troe_Ts_DEF[31], troe_Tss_DEF[31], troe_Tsss_DEF[31];
static double sri_a_DEF[31], sri_b_DEF[31], sri_c_DEF[31], sri_d_DEF[31], sri_e_DEF[31];
static double activation_units_DEF[31], prefactor_units_DEF[31], phase_units_DEF[31];
static int is_PD_DEF[31], troe_len_DEF[31], sri_len_DEF[31], nTB_DEF[31], *TBid_DEF[31];
static double *TB_DEF[31];
static int rxn_map[31] = {4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,0,1,2,20,21,22,23,24,25,26,27,28,3,29,30};

void GET_REACTION_MAP(int *rmap)
{
    for (int i=0; i<31; ++i) {
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
  if (reaction_id<0 || reaction_id>=31) {
    printf("Bad reaction id = %d",reaction_id);
    abort();
  };
  int mrid = rxn_map[reaction_id];

  if (param_id == THIRD_BODY) {
    if (species_id<0 || species_id>=20) {
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
    for (int i=0; i<31; i++) {
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
    for (int i=0; i<31; i++) {
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
  for (int i=0; i<31; ++i) {
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
    // (0):  CH4 <=> CH3 + H
    fwd_A[4]     = 4230000000000000;
    fwd_beta[4]  = 0;
    fwd_Ea[4]    = 108690;
    prefactor_units[4]  = 1;
    activation_units[4] = 0.50321666580471969;
    phase_units[4]      = 1e-6;
    is_PD[4] = 0;
    nTB[4] = 0;

    // (1):  CH4 + OH <=> CH3 + H2O
    fwd_A[5]     = 200000000000000;
    fwd_beta[5]  = 0;
    fwd_Ea[5]    = 8410;
    prefactor_units[5]  = 1.0000000000000002e-06;
    activation_units[5] = 0.50321666580471969;
    phase_units[5]      = 1e-12;
    is_PD[5] = 0;
    nTB[5] = 0;

    // (2):  CH4 + O <=> CH3 + OH
    fwd_A[6]     = 34800000000000;
    fwd_beta[6]  = 0;
    fwd_Ea[6]    = 8380;
    prefactor_units[6]  = 1.0000000000000002e-06;
    activation_units[6] = 0.50321666580471969;
    phase_units[6]      = 1e-12;
    is_PD[6] = 0;
    nTB[6] = 0;

    // (3):  CH4 + H <=> CH3 + H2
    fwd_A[7]     = 435000000000000;
    fwd_beta[7]  = 0;
    fwd_Ea[7]    = 13740;
    prefactor_units[7]  = 1.0000000000000002e-06;
    activation_units[7] = 0.50321666580471969;
    phase_units[7]      = 1e-12;
    is_PD[7] = 0;
    nTB[7] = 0;

    // (4):  CH3 + O2 <=> CH2O + OH
    fwd_A[8]     = 529000000000;
    fwd_beta[8]  = 0;
    fwd_Ea[8]    = 1700;
    prefactor_units[8]  = 1.0000000000000002e-06;
    activation_units[8] = 0.50321666580471969;
    phase_units[8]      = 1e-12;
    is_PD[8] = 0;
    nTB[8] = 0;

    // (5):  CH2O + OH <=> CO + H2O + H
    fwd_A[9]     = 587000000000000;
    fwd_beta[9]  = 0;
    fwd_Ea[9]    = 4880;
    prefactor_units[9]  = 1.0000000000000002e-06;
    activation_units[9] = 0.50321666580471969;
    phase_units[9]      = 1e-12;
    is_PD[9] = 0;
    nTB[9] = 0;

    // (6):  CO + OH => CO2 + H
    fwd_A[10]     = 1300000000000;
    fwd_beta[10]  = 0;
    fwd_Ea[10]    = 1530;
    prefactor_units[10]  = 1.0000000000000002e-06;
    activation_units[10] = 0.50321666580471969;
    phase_units[10]      = 1e-12;
    is_PD[10] = 0;
    nTB[10] = 0;

    // (7):  CO2 + H => CO + OH
    fwd_A[11]     = 145000000000000;
    fwd_beta[11]  = 0;
    fwd_Ea[11]    = 23760;
    prefactor_units[11]  = 1.0000000000000002e-06;
    activation_units[11] = 0.50321666580471969;
    phase_units[11]      = 1e-12;
    is_PD[11] = 0;
    nTB[11] = 0;

    // (8):  O2 + H => OH + O
    fwd_A[12]     = 224000000000000;
    fwd_beta[12]  = 0;
    fwd_Ea[12]    = 16800;
    prefactor_units[12]  = 1.0000000000000002e-06;
    activation_units[12] = 0.50321666580471969;
    phase_units[12]      = 1e-12;
    is_PD[12] = 0;
    nTB[12] = 0;

    // (9):  OH + O => O2 + H
    fwd_A[13]     = 17100000000000;
    fwd_beta[13]  = 0;
    fwd_Ea[13]    = 870;
    prefactor_units[13]  = 1.0000000000000002e-06;
    activation_units[13] = 0.50321666580471969;
    phase_units[13]      = 1e-12;
    is_PD[13] = 0;
    nTB[13] = 0;

    // (10):  O + H2 => OH + H
    fwd_A[14]     = 17400000000000;
    fwd_beta[14]  = 0;
    fwd_Ea[14]    = 9450;
    prefactor_units[14]  = 1.0000000000000002e-06;
    activation_units[14] = 0.50321666580471969;
    phase_units[14]      = 1e-12;
    is_PD[14] = 0;
    nTB[14] = 0;

    // (11):  OH + H => O + H2
    fwd_A[15]     = 7700000000000;
    fwd_beta[15]  = 0;
    fwd_Ea[15]    = 7580;
    prefactor_units[15]  = 1.0000000000000002e-06;
    activation_units[15] = 0.50321666580471969;
    phase_units[15]      = 1e-12;
    is_PD[15] = 0;
    nTB[15] = 0;

    // (12):  O + H2O => 2 OH
    fwd_A[16]     = 57500000000000;
    fwd_beta[16]  = 0;
    fwd_Ea[16]    = 18100;
    prefactor_units[16]  = 1.0000000000000002e-06;
    activation_units[16] = 0.50321666580471969;
    phase_units[16]      = 1e-12;
    is_PD[16] = 0;
    nTB[16] = 0;

    // (13):  2 OH => O + H2O
    fwd_A[17]     = 5380000000000;
    fwd_beta[17]  = 0;
    fwd_Ea[17]    = 1050;
    prefactor_units[17]  = 1.0000000000000002e-06;
    activation_units[17] = 0.50321666580471969;
    phase_units[17]      = 1e-12;
    is_PD[17] = 0;
    nTB[17] = 0;

    // (14):  OH + H2 => H2O + H
    fwd_A[18]     = 21900000000000;
    fwd_beta[18]  = 0;
    fwd_Ea[18]    = 5150;
    prefactor_units[18]  = 1.0000000000000002e-06;
    activation_units[18] = 0.50321666580471969;
    phase_units[18]      = 1e-12;
    is_PD[18] = 0;
    nTB[18] = 0;

    // (15):  H2O + H => H2 + OH
    fwd_A[19]     = 84100000000000;
    fwd_beta[19]  = 0;
    fwd_Ea[19]    = 20100;
    prefactor_units[19]  = 1.0000000000000002e-06;
    activation_units[19] = 0.50321666580471969;
    phase_units[19]      = 1e-12;
    is_PD[19] = 0;
    nTB[19] = 0;

    // (16):  H + OH + M <=> H2O + M
    fwd_A[0]     = 2e+19;
    fwd_beta[0]  = -1;
    fwd_Ea[0]    = 0;
    prefactor_units[0]  = 1.0000000000000002e-12;
    activation_units[0] = 0.50321666580471969;
    phase_units[0]      = 1e-12;
    is_PD[0] = 0;
    nTB[0] = 3;
    TB[0] = (double *) malloc(3 * sizeof(double));
    TBid[0] = (int *) malloc(3 * sizeof(int));
    TBid[0][0] = 9; TB[0][0] = 0.72999999999999998; // H2
    TBid[0][1] = 4; TB[0][1] = 3.6499999999999999; // H2O
    TBid[0][2] = 1; TB[0][2] = 2; // CH4

    // (17):  O + O + M <=> O2 + M
    fwd_A[1]     = 890000000000000;
    fwd_beta[1]  = -0.5;
    fwd_Ea[1]    = 0;
    prefactor_units[1]  = 1.0000000000000002e-12;
    activation_units[1] = 0.50321666580471969;
    phase_units[1]      = 1e-12;
    is_PD[1] = 0;
    nTB[1] = 5;
    TB[1] = (double *) malloc(5 * sizeof(double));
    TBid[1] = (int *) malloc(5 * sizeof(int));
    TBid[1][0] = 9; TB[1][0] = 2.3999999999999999; // H2
    TBid[1][1] = 4; TB[1][1] = 15.4; // H2O
    TBid[1][2] = 1; TB[1][2] = 2; // CH4
    TBid[1][3] = 12; TB[1][3] = 1.75; // CO
    TBid[1][4] = 13; TB[1][4] = 3.6000000000000001; // CO2

    // (18):  H + H + M <=> H2 + M
    fwd_A[2]     = 1e+18;
    fwd_beta[2]  = -1;
    fwd_Ea[2]    = 0;
    prefactor_units[2]  = 1.0000000000000002e-12;
    activation_units[2] = 0.50321666580471969;
    phase_units[2]      = 1e-12;
    is_PD[2] = 0;
    nTB[2] = 4;
    TB[2] = (double *) malloc(4 * sizeof(double));
    TBid[2] = (int *) malloc(4 * sizeof(int));
    TBid[2][0] = 9; TB[2][0] = 0; // H2
    TBid[2][1] = 4; TB[2][1] = 0; // H2O
    TBid[2][2] = 1; TB[2][2] = 2; // CH4
    TBid[2][3] = 13; TB[2][3] = 0; // CO2

    // (19):  CH3 + O <=> CH + H2O
    fwd_A[20]     = 280000000;
    fwd_beta[20]  = 0;
    fwd_Ea[20]    = 0;
    prefactor_units[20]  = 1.0000000000000002e-06;
    activation_units[20] = 0.50321666580471969;
    phase_units[20]      = 1e-12;
    is_PD[20] = 0;
    nTB[20] = 0;

    // (20):  CH + O <=> HCOp + E
    fwd_A[21]     = 575000000000;
    fwd_beta[21]  = 0;
    fwd_Ea[21]    = 6000;
    prefactor_units[21]  = 1.0000000000000002e-06;
    activation_units[21] = 0.50321666580471969;
    phase_units[21]      = 1e-12;
    is_PD[21] = 0;
    nTB[21] = 0;

    // (21):  HCOp + H2O <=> CO + H3Op
    fwd_A[22]     = 5.02e+17;
    fwd_beta[22]  = 0;
    fwd_Ea[22]    = 24000;
    prefactor_units[22]  = 1.0000000000000002e-06;
    activation_units[22] = 0.50321666580471969;
    phase_units[22]      = 1e-12;
    is_PD[22] = 0;
    nTB[22] = 0;

    // (22):  H3Op + E <=> H2O + H
    fwd_A[23]     = 1.44e+17;
    fwd_beta[23]  = 0;
    fwd_Ea[23]    = 0;
    prefactor_units[23]  = 1.0000000000000002e-06;
    activation_units[23] = 0.50321666580471969;
    phase_units[23]      = 1e-12;
    is_PD[23] = 0;
    nTB[23] = 0;

    // (23):  CH + O2 <=> CO + OH
    fwd_A[24]     = 60000000000;
    fwd_beta[24]  = 0;
    fwd_Ea[24]    = 0;
    prefactor_units[24]  = 1.0000000000000002e-06;
    activation_units[24] = 0.50321666580471969;
    phase_units[24]      = 1e-12;
    is_PD[24] = 0;
    nTB[24] = 0;

    // (24):  O + N2 => NO + N
    fwd_A[25]     = 136000000000000;
    fwd_beta[25]  = 0;
    fwd_Ea[25]    = 75400;
    prefactor_units[25]  = 1.0000000000000002e-06;
    activation_units[25] = 0.50321666580471969;
    phase_units[25]      = 1e-12;
    is_PD[25] = 0;
    nTB[25] = 0;

    // (25):  NO + N => O + N2
    fwd_A[26]     = 31000000000000;
    fwd_beta[26]  = 0;
    fwd_Ea[26]    = 330;
    prefactor_units[26]  = 1.0000000000000002e-06;
    activation_units[26] = 0.50321666580471969;
    phase_units[26]      = 1e-12;
    is_PD[26] = 0;
    nTB[26] = 0;

    // (26):  N + O2 => NO + O
    fwd_A[27]     = 6430000000;
    fwd_beta[27]  = 1;
    fwd_Ea[27]    = 6250;
    prefactor_units[27]  = 1.0000000000000002e-06;
    activation_units[27] = 0.50321666580471969;
    phase_units[27]      = 1e-12;
    is_PD[27] = 0;
    nTB[27] = 0;

    // (27):  NO + O => N + O2
    fwd_A[28]     = 1550000000;
    fwd_beta[28]  = 1;
    fwd_Ea[28]    = 38640;
    prefactor_units[28]  = 1.0000000000000002e-06;
    activation_units[28] = 0.50321666580471969;
    phase_units[28]      = 1e-12;
    is_PD[28] = 0;
    nTB[28] = 0;

    // (28):  NO + O + M <=> NO2 + M
    fwd_A[3]     = 1050000000000000;
    fwd_beta[3]  = 0;
    fwd_Ea[3]    = -1870;
    prefactor_units[3]  = 1.0000000000000002e-12;
    activation_units[3] = 0.50321666580471969;
    phase_units[3]      = 1e-12;
    is_PD[3] = 0;
    nTB[3] = 5;
    TB[3] = (double *) malloc(5 * sizeof(double));
    TBid[3] = (int *) malloc(5 * sizeof(int));
    TBid[3][0] = 9; TB[3][0] = 2; // H2
    TBid[3][1] = 4; TB[3][1] = 6; // H2O
    TBid[3][2] = 1; TB[3][2] = 2; // CH4
    TBid[3][3] = 12; TB[3][3] = 1.5; // CO
    TBid[3][4] = 13; TB[3][4] = 2; // CO2

    // (29):  NO2 + O <=> NO + O2
    fwd_A[29]     = 2100000000000;
    fwd_beta[29]  = 0;
    fwd_Ea[29]    = 0;
    prefactor_units[29]  = 1.0000000000000002e-06;
    activation_units[29] = 0.50321666580471969;
    phase_units[29]      = 1e-12;
    is_PD[29] = 0;
    nTB[29] = 0;

    // (30):  NO2 + H <=> NO + OH
    fwd_A[30]     = 300000000000000;
    fwd_beta[30]  = 0;
    fwd_Ea[30]    = 0;
    prefactor_units[30]  = 1.0000000000000002e-06;
    activation_units[30] = 0.50321666580471969;
    phase_units[30]      = 1e-12;
    is_PD[30] = 0;
    nTB[30] = 0;

    SetAllDefaults();
}



/*A few mechanism parameters */
void CKINDX(int * iwrk, double * restrict rwrk, int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 5;
    *kk = 20;
    *ii = 31;
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
    for (i=0; i<lenkname*5; i++) {
        kname[i] = ' ';
    }

    /* O  */
    kname[ 0*lenkname + 0 ] = 'O';
    kname[ 0*lenkname + 1 ] = ' ';

    /* H  */
    kname[ 1*lenkname + 0 ] = 'H';
    kname[ 1*lenkname + 1 ] = ' ';

    /* C  */
    kname[ 2*lenkname + 0 ] = 'C';
    kname[ 2*lenkname + 1 ] = ' ';

    /* N  */
    kname[ 3*lenkname + 0 ] = 'N';
    kname[ 3*lenkname + 1 ] = ' ';

    /* E  */
    kname[ 4*lenkname + 0 ] = 'E';
    kname[ 4*lenkname + 1 ] = ' ';

}


/* Returns the char strings of species names */
void CKSYMS(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*20; i++) {
        kname[i] = ' ';
    }

    /* CH3  */
    kname[ 0*lenkname + 0 ] = 'C';
    kname[ 0*lenkname + 1 ] = 'H';
    kname[ 0*lenkname + 2 ] = '3';
    kname[ 0*lenkname + 3 ] = ' ';

    /* CH4  */
    kname[ 1*lenkname + 0 ] = 'C';
    kname[ 1*lenkname + 1 ] = 'H';
    kname[ 1*lenkname + 2 ] = '4';
    kname[ 1*lenkname + 3 ] = ' ';

    /* H  */
    kname[ 2*lenkname + 0 ] = 'H';
    kname[ 2*lenkname + 1 ] = ' ';

    /* OH  */
    kname[ 3*lenkname + 0 ] = 'O';
    kname[ 3*lenkname + 1 ] = 'H';
    kname[ 3*lenkname + 2 ] = ' ';

    /* H2O  */
    kname[ 4*lenkname + 0 ] = 'H';
    kname[ 4*lenkname + 1 ] = '2';
    kname[ 4*lenkname + 2 ] = 'O';
    kname[ 4*lenkname + 3 ] = ' ';

    /* O  */
    kname[ 5*lenkname + 0 ] = 'O';
    kname[ 5*lenkname + 1 ] = ' ';

    /* HCOp  */
    kname[ 6*lenkname + 0 ] = 'H';
    kname[ 6*lenkname + 1 ] = 'C';
    kname[ 6*lenkname + 2 ] = 'O';
    kname[ 6*lenkname + 3 ] = 'P';
    kname[ 6*lenkname + 4 ] = ' ';

    /* NO  */
    kname[ 7*lenkname + 0 ] = 'N';
    kname[ 7*lenkname + 1 ] = 'O';
    kname[ 7*lenkname + 2 ] = ' ';

    /* N  */
    kname[ 8*lenkname + 0 ] = 'N';
    kname[ 8*lenkname + 1 ] = ' ';

    /* H2  */
    kname[ 9*lenkname + 0 ] = 'H';
    kname[ 9*lenkname + 1 ] = '2';
    kname[ 9*lenkname + 2 ] = ' ';

    /* O2  */
    kname[ 10*lenkname + 0 ] = 'O';
    kname[ 10*lenkname + 1 ] = '2';
    kname[ 10*lenkname + 2 ] = ' ';

    /* CH2O  */
    kname[ 11*lenkname + 0 ] = 'C';
    kname[ 11*lenkname + 1 ] = 'H';
    kname[ 11*lenkname + 2 ] = '2';
    kname[ 11*lenkname + 3 ] = 'O';
    kname[ 11*lenkname + 4 ] = ' ';

    /* CO  */
    kname[ 12*lenkname + 0 ] = 'C';
    kname[ 12*lenkname + 1 ] = 'O';
    kname[ 12*lenkname + 2 ] = ' ';

    /* CO2  */
    kname[ 13*lenkname + 0 ] = 'C';
    kname[ 13*lenkname + 1 ] = 'O';
    kname[ 13*lenkname + 2 ] = '2';
    kname[ 13*lenkname + 3 ] = ' ';

    /* CH  */
    kname[ 14*lenkname + 0 ] = 'C';
    kname[ 14*lenkname + 1 ] = 'H';
    kname[ 14*lenkname + 2 ] = ' ';

    /* H3Op  */
    kname[ 15*lenkname + 0 ] = 'H';
    kname[ 15*lenkname + 1 ] = '3';
    kname[ 15*lenkname + 2 ] = 'O';
    kname[ 15*lenkname + 3 ] = 'P';
    kname[ 15*lenkname + 4 ] = ' ';

    /* N2  */
    kname[ 16*lenkname + 0 ] = 'N';
    kname[ 16*lenkname + 1 ] = '2';
    kname[ 16*lenkname + 2 ] = ' ';

    /* NO2  */
    kname[ 17*lenkname + 0 ] = 'N';
    kname[ 17*lenkname + 1 ] = 'O';
    kname[ 17*lenkname + 2 ] = '2';
    kname[ 17*lenkname + 3 ] = ' ';

    /* E  */
    kname[ 18*lenkname + 0 ] = 'E';
    kname[ 18*lenkname + 1 ] = ' ';

    /* CH2  */
    kname[ 19*lenkname + 0 ] = 'C';
    kname[ 19*lenkname + 1 ] = 'H';
    kname[ 19*lenkname + 2 ] = '2';
    kname[ 19*lenkname + 3 ] = ' ';

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
    XW += x[0]*15.035060; /*CH3 */
    XW += x[1]*16.043030; /*CH4 */
    XW += x[2]*1.007970; /*H */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*15.999400; /*O */
    XW += x[6]*29.017971; /*HCOp */
    XW += x[7]*30.006100; /*NO */
    XW += x[8]*14.006700; /*N */
    XW += x[9]*2.015940; /*H2 */
    XW += x[10]*31.998800; /*O2 */
    XW += x[11]*30.026490; /*CH2O */
    XW += x[12]*28.010550; /*CO */
    XW += x[13]*44.009950; /*CO2 */
    XW += x[14]*13.019120; /*CH */
    XW += x[15]*19.022761; /*H3Op */
    XW += x[16]*28.013400; /*N2 */
    XW += x[17]*46.005500; /*NO2 */
    XW += x[18]*0.000549; /*E */
    XW += x[19]*14.027090; /*CH2 */
    *P = *rho * 8.31451e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void CKPY(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]*imw[0]; /*CH3 */
    YOW += y[1]*imw[1]; /*CH4 */
    YOW += y[2]*imw[2]; /*H */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*O */
    YOW += y[6]*imw[6]; /*HCOp */
    YOW += y[7]*imw[7]; /*NO */
    YOW += y[8]*imw[8]; /*N */
    YOW += y[9]*imw[9]; /*H2 */
    YOW += y[10]*imw[10]; /*O2 */
    YOW += y[11]*imw[11]; /*CH2O */
    YOW += y[12]*imw[12]; /*CO */
    YOW += y[13]*imw[13]; /*CO2 */
    YOW += y[14]*imw[14]; /*CH */
    YOW += y[15]*imw[15]; /*H3Op */
    YOW += y[16]*imw[16]; /*N2 */
    YOW += y[17]*imw[17]; /*NO2 */
    YOW += y[18]*imw[18]; /*E */
    YOW += y[19]*imw[19]; /*CH2 */
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

    for (int n=0; n<20; n++) {
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
    W += c[0]*15.035060; /*CH3 */
    W += c[1]*16.043030; /*CH4 */
    W += c[2]*1.007970; /*H */
    W += c[3]*17.007370; /*OH */
    W += c[4]*18.015340; /*H2O */
    W += c[5]*15.999400; /*O */
    W += c[6]*29.017971; /*HCOp */
    W += c[7]*30.006100; /*NO */
    W += c[8]*14.006700; /*N */
    W += c[9]*2.015940; /*H2 */
    W += c[10]*31.998800; /*O2 */
    W += c[11]*30.026490; /*CH2O */
    W += c[12]*28.010550; /*CO */
    W += c[13]*44.009950; /*CO2 */
    W += c[14]*13.019120; /*CH */
    W += c[15]*19.022761; /*H3Op */
    W += c[16]*28.013400; /*N2 */
    W += c[17]*46.005500; /*NO2 */
    W += c[18]*0.000549; /*E */
    W += c[19]*14.027090; /*CH2 */

    for (id = 0; id < 20; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.31451e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
void CKRHOX(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict rho)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*15.035060; /*CH3 */
    XW += x[1]*16.043030; /*CH4 */
    XW += x[2]*1.007970; /*H */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*15.999400; /*O */
    XW += x[6]*29.017971; /*HCOp */
    XW += x[7]*30.006100; /*NO */
    XW += x[8]*14.006700; /*N */
    XW += x[9]*2.015940; /*H2 */
    XW += x[10]*31.998800; /*O2 */
    XW += x[11]*30.026490; /*CH2O */
    XW += x[12]*28.010550; /*CO */
    XW += x[13]*44.009950; /*CO2 */
    XW += x[14]*13.019120; /*CH */
    XW += x[15]*19.022761; /*H3Op */
    XW += x[16]*28.013400; /*N2 */
    XW += x[17]*46.005500; /*NO2 */
    XW += x[18]*0.000549; /*E */
    XW += x[19]*14.027090; /*CH2 */
    *rho = *P * XW / (8.31451e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
void CKRHOY(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict rho)
{
    double YOW = 0;
    double tmp[20];

    for (int i = 0; i < 20; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 20; i++)
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
    W += c[0]*15.035060; /*CH3 */
    W += c[1]*16.043030; /*CH4 */
    W += c[2]*1.007970; /*H */
    W += c[3]*17.007370; /*OH */
    W += c[4]*18.015340; /*H2O */
    W += c[5]*15.999400; /*O */
    W += c[6]*29.017971; /*HCOp */
    W += c[7]*30.006100; /*NO */
    W += c[8]*14.006700; /*N */
    W += c[9]*2.015940; /*H2 */
    W += c[10]*31.998800; /*O2 */
    W += c[11]*30.026490; /*CH2O */
    W += c[12]*28.010550; /*CO */
    W += c[13]*44.009950; /*CO2 */
    W += c[14]*13.019120; /*CH */
    W += c[15]*19.022761; /*H3Op */
    W += c[16]*28.013400; /*N2 */
    W += c[17]*46.005500; /*NO2 */
    W += c[18]*0.000549; /*E */
    W += c[19]*14.027090; /*CH2 */

    for (id = 0; id < 20; ++id) {
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
    double tmp[20];

    for (int i = 0; i < 20; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 20; i++)
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
    XW += x[0]*15.035060; /*CH3 */
    XW += x[1]*16.043030; /*CH4 */
    XW += x[2]*1.007970; /*H */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*15.999400; /*O */
    XW += x[6]*29.017971; /*HCOp */
    XW += x[7]*30.006100; /*NO */
    XW += x[8]*14.006700; /*N */
    XW += x[9]*2.015940; /*H2 */
    XW += x[10]*31.998800; /*O2 */
    XW += x[11]*30.026490; /*CH2O */
    XW += x[12]*28.010550; /*CO */
    XW += x[13]*44.009950; /*CO2 */
    XW += x[14]*13.019120; /*CH */
    XW += x[15]*19.022761; /*H3Op */
    XW += x[16]*28.013400; /*N2 */
    XW += x[17]*46.005500; /*NO2 */
    XW += x[18]*0.000549; /*E */
    XW += x[19]*14.027090; /*CH2 */
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
    W += c[0]*15.035060; /*CH3 */
    W += c[1]*16.043030; /*CH4 */
    W += c[2]*1.007970; /*H */
    W += c[3]*17.007370; /*OH */
    W += c[4]*18.015340; /*H2O */
    W += c[5]*15.999400; /*O */
    W += c[6]*29.017971; /*HCOp */
    W += c[7]*30.006100; /*NO */
    W += c[8]*14.006700; /*N */
    W += c[9]*2.015940; /*H2 */
    W += c[10]*31.998800; /*O2 */
    W += c[11]*30.026490; /*CH2O */
    W += c[12]*28.010550; /*CO */
    W += c[13]*44.009950; /*CO2 */
    W += c[14]*13.019120; /*CH */
    W += c[15]*19.022761; /*H3Op */
    W += c[16]*28.013400; /*N2 */
    W += c[17]*46.005500; /*NO2 */
    W += c[18]*0.000549; /*E */
    W += c[19]*14.027090; /*CH2 */

    for (id = 0; id < 20; ++id) {
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
    double tmp[20];

    for (int i = 0; i < 20; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 20; i++)
    {
        YOW += tmp[i];
    }

    double YOWINV = 1.0/YOW;

    for (int i = 0; i < 20; i++)
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

    for (int n=0; n<20; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] = y[n*(*np)+i] * imw[n];
            YOW[i] += x[n*(*np)+i];
        }
    }

    for (int i=0; i<(*np); i++) {
        YOW[i] = 1.0/YOW[i];
    }

    for (int n=0; n<20; n++) {
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
    for (int i = 0; i < 20; i++)
    {
        c[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 20; i++)
    {
        YOW += c[i];
    }

    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*Now compute conversion */

    for (int i = 0; i < 20; i++)
    {
        c[i] = PWORT * y[i] * imw[i];
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict c)
{
    for (int i = 0; i < 20; i++)
    {
        c[i] = (*rho)  * y[i] * imw[i];
    }
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
void CKXTY(double * restrict x, int * iwrk, double * restrict rwrk, double * restrict y)
{
    double XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*15.035060; /*CH3 */
    XW += x[1]*16.043030; /*CH4 */
    XW += x[2]*1.007970; /*H */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*15.999400; /*O */
    XW += x[6]*29.017971; /*HCOp */
    XW += x[7]*30.006100; /*NO */
    XW += x[8]*14.006700; /*N */
    XW += x[9]*2.015940; /*H2 */
    XW += x[10]*31.998800; /*O2 */
    XW += x[11]*30.026490; /*CH2O */
    XW += x[12]*28.010550; /*CO */
    XW += x[13]*44.009950; /*CO2 */
    XW += x[14]*13.019120; /*CH */
    XW += x[15]*19.022761; /*H3Op */
    XW += x[16]*28.013400; /*N2 */
    XW += x[17]*46.005500; /*NO2 */
    XW += x[18]*0.000549; /*E */
    XW += x[19]*14.027090; /*CH2 */
    /*Now compute conversion */
    double XWinv = 1.0/XW;
    y[0] = x[0]*15.035060*XWinv; 
    y[1] = x[1]*16.043030*XWinv; 
    y[2] = x[2]*1.007970*XWinv; 
    y[3] = x[3]*17.007370*XWinv; 
    y[4] = x[4]*18.015340*XWinv; 
    y[5] = x[5]*15.999400*XWinv; 
    y[6] = x[6]*29.017971*XWinv; 
    y[7] = x[7]*30.006100*XWinv; 
    y[8] = x[8]*14.006700*XWinv; 
    y[9] = x[9]*2.015940*XWinv; 
    y[10] = x[10]*31.998800*XWinv; 
    y[11] = x[11]*30.026490*XWinv; 
    y[12] = x[12]*28.010550*XWinv; 
    y[13] = x[13]*44.009950*XWinv; 
    y[14] = x[14]*13.019120*XWinv; 
    y[15] = x[15]*19.022761*XWinv; 
    y[16] = x[16]*28.013400*XWinv; 
    y[17] = x[17]*46.005500*XWinv; 
    y[18] = x[18]*0.000549*XWinv; 
    y[19] = x[19]*14.027090*XWinv; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.31451e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 20; ++id) {
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
    XW += x[0]*15.035060; /*CH3 */
    XW += x[1]*16.043030; /*CH4 */
    XW += x[2]*1.007970; /*H */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*15.999400; /*O */
    XW += x[6]*29.017971; /*HCOp */
    XW += x[7]*30.006100; /*NO */
    XW += x[8]*14.006700; /*N */
    XW += x[9]*2.015940; /*H2 */
    XW += x[10]*31.998800; /*O2 */
    XW += x[11]*30.026490; /*CH2O */
    XW += x[12]*28.010550; /*CO */
    XW += x[13]*44.009950; /*CO2 */
    XW += x[14]*13.019120; /*CH */
    XW += x[15]*19.022761; /*H3Op */
    XW += x[16]*28.013400; /*N2 */
    XW += x[17]*46.005500; /*NO2 */
    XW += x[18]*0.000549; /*E */
    XW += x[19]*14.027090; /*CH2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 20; ++id) {
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
    for (id = 0; id < 20; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    double sumCinv = 1.0/sumC;
    for (id = 0; id < 20; ++id) {
        x[id] = c[id]*sumCinv;
    }

    return;
}


/*convert c[species] (molar conc) to y[species] (mass fracs) */
void CKCTY(double * restrict c, int * iwrk, double * restrict rwrk, double * restrict y)
{
    double CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*15.035060; /*CH3 */
    CW += c[1]*16.043030; /*CH4 */
    CW += c[2]*1.007970; /*H */
    CW += c[3]*17.007370; /*OH */
    CW += c[4]*18.015340; /*H2O */
    CW += c[5]*15.999400; /*O */
    CW += c[6]*29.017971; /*HCOp */
    CW += c[7]*30.006100; /*NO */
    CW += c[8]*14.006700; /*N */
    CW += c[9]*2.015940; /*H2 */
    CW += c[10]*31.998800; /*O2 */
    CW += c[11]*30.026490; /*CH2O */
    CW += c[12]*28.010550; /*CO */
    CW += c[13]*44.009950; /*CO2 */
    CW += c[14]*13.019120; /*CH */
    CW += c[15]*19.022761; /*H3Op */
    CW += c[16]*28.013400; /*N2 */
    CW += c[17]*46.005500; /*NO2 */
    CW += c[18]*0.000549; /*E */
    CW += c[19]*14.027090; /*CH2 */
    /*Now compute conversion */
    double CWinv = 1.0/CW;
    y[0] = c[0]*15.035060*CWinv; 
    y[1] = c[1]*16.043030*CWinv; 
    y[2] = c[2]*1.007970*CWinv; 
    y[3] = c[3]*17.007370*CWinv; 
    y[4] = c[4]*18.015340*CWinv; 
    y[5] = c[5]*15.999400*CWinv; 
    y[6] = c[6]*29.017971*CWinv; 
    y[7] = c[7]*30.006100*CWinv; 
    y[8] = c[8]*14.006700*CWinv; 
    y[9] = c[9]*2.015940*CWinv; 
    y[10] = c[10]*31.998800*CWinv; 
    y[11] = c[11]*30.026490*CWinv; 
    y[12] = c[12]*28.010550*CWinv; 
    y[13] = c[13]*44.009950*CWinv; 
    y[14] = c[14]*13.019120*CWinv; 
    y[15] = c[15]*19.022761*CWinv; 
    y[16] = c[16]*28.013400*CWinv; 
    y[17] = c[17]*46.005500*CWinv; 
    y[18] = c[18]*0.000549*CWinv; 
    y[19] = c[19]*14.027090*CWinv; 

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
    for (id = 0; id < 20; ++id) {
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
    for (id = 0; id < 20; ++id) {
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
    for (id = 0; id < 20; ++id) {
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
    for (id = 0; id < 20; ++id) {
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
    for (id = 0; id < 20; ++id) {
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
    for (id = 0; id < 20; ++id) {
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
    for (id = 0; id < 20; ++id) {
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
    cvms[0] *= 5.530081023953346e+06; /*CH3 */
    cvms[1] *= 5.182630712527496e+06; /*CH4 */
    cvms[2] *= 8.248767324424338e+07; /*H */
    cvms[3] *= 4.888768810227566e+06; /*OH */
    cvms[4] *= 4.615239012974499e+06; /*H2O */
    cvms[5] *= 5.196763628636074e+06; /*O */
    cvms[6] *= 2.865296777533546e+06; /*HCOp */
    cvms[7] *= 2.770939908885194e+06; /*NO */
    cvms[8] *= 5.936094868884177e+06; /*N */
    cvms[9] *= 4.124383662212169e+07; /*H2 */
    cvms[10] *= 2.598381814318037e+06; /*O2 */
    cvms[11] *= 2.769058254894261e+06; /*CH2O */
    cvms[12] *= 2.968349425484326e+06; /*CO */
    cvms[13] *= 1.889234139098090e+06; /*CO2 */
    cvms[14] *= 6.386384025955671e+06; /*CH */
    cvms[15] *= 4.370821783688074e+06; /*H3Op */
    cvms[16] *= 2.968047434442088e+06; /*N2 */
    cvms[17] *= 1.807286085359359e+06; /*NO2 */
    cvms[18] *= 1.515641927222661e+11; /*E */
    cvms[19] *= 5.927466067445207e+06; /*CH2 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
void CKCPMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict cpms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 5.530081023953346e+06; /*CH3 */
    cpms[1] *= 5.182630712527496e+06; /*CH4 */
    cpms[2] *= 8.248767324424338e+07; /*H */
    cpms[3] *= 4.888768810227566e+06; /*OH */
    cpms[4] *= 4.615239012974499e+06; /*H2O */
    cpms[5] *= 5.196763628636074e+06; /*O */
    cpms[6] *= 2.865296777533546e+06; /*HCOp */
    cpms[7] *= 2.770939908885194e+06; /*NO */
    cpms[8] *= 5.936094868884177e+06; /*N */
    cpms[9] *= 4.124383662212169e+07; /*H2 */
    cpms[10] *= 2.598381814318037e+06; /*O2 */
    cpms[11] *= 2.769058254894261e+06; /*CH2O */
    cpms[12] *= 2.968349425484326e+06; /*CO */
    cpms[13] *= 1.889234139098090e+06; /*CO2 */
    cpms[14] *= 6.386384025955671e+06; /*CH */
    cpms[15] *= 4.370821783688074e+06; /*H3Op */
    cpms[16] *= 2.968047434442088e+06; /*N2 */
    cpms[17] *= 1.807286085359359e+06; /*NO2 */
    cpms[18] *= 1.515641927222661e+11; /*E */
    cpms[19] *= 5.927466067445207e+06; /*CH2 */
}


/*Returns internal energy in mass units (Eq 30.) */
void CKUMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    for (int i = 0; i < 20; i++)
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
    for (int i = 0; i < 20; i++)
    {
        hms[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int * restrict np, double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hms)
{
    double tc[5], h[20];

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
        hms[9*(*np)+i] = h[9];
        hms[10*(*np)+i] = h[10];
        hms[11*(*np)+i] = h[11];
        hms[12*(*np)+i] = h[12];
        hms[13*(*np)+i] = h[13];
        hms[14*(*np)+i] = h[14];
        hms[15*(*np)+i] = h[15];
        hms[16*(*np)+i] = h[16];
        hms[17*(*np)+i] = h[17];
        hms[18*(*np)+i] = h[18];
        hms[19*(*np)+i] = h[19];
    }

    for (int n=0; n<20; n++) {
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
    for (int i = 0; i < 20; i++)
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
    for (int i = 0; i < 20; i++)
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
    sms[0] *= 5.530081023953346e+06; /*CH3 */
    sms[1] *= 5.182630712527496e+06; /*CH4 */
    sms[2] *= 8.248767324424338e+07; /*H */
    sms[3] *= 4.888768810227566e+06; /*OH */
    sms[4] *= 4.615239012974499e+06; /*H2O */
    sms[5] *= 5.196763628636074e+06; /*O */
    sms[6] *= 2.865296777533546e+06; /*HCOp */
    sms[7] *= 2.770939908885194e+06; /*NO */
    sms[8] *= 5.936094868884177e+06; /*N */
    sms[9] *= 4.124383662212169e+07; /*H2 */
    sms[10] *= 2.598381814318037e+06; /*O2 */
    sms[11] *= 2.769058254894261e+06; /*CH2O */
    sms[12] *= 2.968349425484326e+06; /*CO */
    sms[13] *= 1.889234139098090e+06; /*CO2 */
    sms[14] *= 6.386384025955671e+06; /*CH */
    sms[15] *= 4.370821783688074e+06; /*H3Op */
    sms[16] *= 2.968047434442088e+06; /*N2 */
    sms[17] *= 1.807286085359359e+06; /*NO2 */
    sms[18] *= 1.515641927222661e+11; /*E */
    sms[19] *= 5.927466067445207e+06; /*CH2 */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[20]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 20; ++id) {
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
    double cpor[20], tresult[20]; /* temporary storage */
    cp_R(cpor, tc);
    for (int i = 0; i < 20; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 20; i++)
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
    double cvor[20]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 20; ++id) {
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
    double cvor[20]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*CH3 */
    result += cvor[1]*y[1]*imw[1]; /*CH4 */
    result += cvor[2]*y[2]*imw[2]; /*H */
    result += cvor[3]*y[3]*imw[3]; /*OH */
    result += cvor[4]*y[4]*imw[4]; /*H2O */
    result += cvor[5]*y[5]*imw[5]; /*O */
    result += cvor[6]*y[6]*imw[6]; /*HCOp */
    result += cvor[7]*y[7]*imw[7]; /*NO */
    result += cvor[8]*y[8]*imw[8]; /*N */
    result += cvor[9]*y[9]*imw[9]; /*H2 */
    result += cvor[10]*y[10]*imw[10]; /*O2 */
    result += cvor[11]*y[11]*imw[11]; /*CH2O */
    result += cvor[12]*y[12]*imw[12]; /*CO */
    result += cvor[13]*y[13]*imw[13]; /*CO2 */
    result += cvor[14]*y[14]*imw[14]; /*CH */
    result += cvor[15]*y[15]*imw[15]; /*H3Op */
    result += cvor[16]*y[16]*imw[16]; /*N2 */
    result += cvor[17]*y[17]*imw[17]; /*NO2 */
    result += cvor[18]*y[18]*imw[18]; /*E */
    result += cvor[19]*y[19]*imw[19]; /*CH2 */

    *cvbs = result * 8.31451e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[20]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 20; ++id) {
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
    double hml[20], tmp[20]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);
    int id;
    for (id = 0; id < 20; ++id) {
        tmp[id] = y[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 20; ++id) {
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
    double uml[20]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 20; ++id) {
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
    double ums[20]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*CH3 */
    result += y[1]*ums[1]*imw[1]; /*CH4 */
    result += y[2]*ums[2]*imw[2]; /*H */
    result += y[3]*ums[3]*imw[3]; /*OH */
    result += y[4]*ums[4]*imw[4]; /*H2O */
    result += y[5]*ums[5]*imw[5]; /*O */
    result += y[6]*ums[6]*imw[6]; /*HCOp */
    result += y[7]*ums[7]*imw[7]; /*NO */
    result += y[8]*ums[8]*imw[8]; /*N */
    result += y[9]*ums[9]*imw[9]; /*H2 */
    result += y[10]*ums[10]*imw[10]; /*O2 */
    result += y[11]*ums[11]*imw[11]; /*CH2O */
    result += y[12]*ums[12]*imw[12]; /*CO */
    result += y[13]*ums[13]*imw[13]; /*CO2 */
    result += y[14]*ums[14]*imw[14]; /*CH */
    result += y[15]*ums[15]*imw[15]; /*H3Op */
    result += y[16]*ums[16]*imw[16]; /*N2 */
    result += y[17]*ums[17]*imw[17]; /*NO2 */
    result += y[18]*ums[18]*imw[18]; /*E */
    result += y[19]*ums[19]*imw[19]; /*CH2 */

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
    double sor[20]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 20; ++id) {
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
    double sor[20]; /* temporary storage */
    double x[20]; /* need a ytx conversion */
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*CH3 */
    YOW += y[1]*imw[1]; /*CH4 */
    YOW += y[2]*imw[2]; /*H */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*O */
    YOW += y[6]*imw[6]; /*HCOp */
    YOW += y[7]*imw[7]; /*NO */
    YOW += y[8]*imw[8]; /*N */
    YOW += y[9]*imw[9]; /*H2 */
    YOW += y[10]*imw[10]; /*O2 */
    YOW += y[11]*imw[11]; /*CH2O */
    YOW += y[12]*imw[12]; /*CO */
    YOW += y[13]*imw[13]; /*CO2 */
    YOW += y[14]*imw[14]; /*CH */
    YOW += y[15]*imw[15]; /*H3Op */
    YOW += y[16]*imw[16]; /*N2 */
    YOW += y[17]*imw[17]; /*NO2 */
    YOW += y[18]*imw[18]; /*E */
    YOW += y[19]*imw[19]; /*CH2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(15.035060*YOW); 
    x[1] = y[1]/(16.043030*YOW); 
    x[2] = y[2]/(1.007970*YOW); 
    x[3] = y[3]/(17.007370*YOW); 
    x[4] = y[4]/(18.015340*YOW); 
    x[5] = y[5]/(15.999400*YOW); 
    x[6] = y[6]/(29.017971*YOW); 
    x[7] = y[7]/(30.006100*YOW); 
    x[8] = y[8]/(14.006700*YOW); 
    x[9] = y[9]/(2.015940*YOW); 
    x[10] = y[10]/(31.998800*YOW); 
    x[11] = y[11]/(30.026490*YOW); 
    x[12] = y[12]/(28.010550*YOW); 
    x[13] = y[13]/(44.009950*YOW); 
    x[14] = y[14]/(13.019120*YOW); 
    x[15] = y[15]/(19.022761*YOW); 
    x[16] = y[16]/(28.013400*YOW); 
    x[17] = y[17]/(46.005500*YOW); 
    x[18] = y[18]/(0.000549*YOW); 
    x[19] = y[19]/(14.027090*YOW); 
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
    result += x[9]*(sor[9]-log((x[9]+1e-100))-logPratio);
    result += x[10]*(sor[10]-log((x[10]+1e-100))-logPratio);
    result += x[11]*(sor[11]-log((x[11]+1e-100))-logPratio);
    result += x[12]*(sor[12]-log((x[12]+1e-100))-logPratio);
    result += x[13]*(sor[13]-log((x[13]+1e-100))-logPratio);
    result += x[14]*(sor[14]-log((x[14]+1e-100))-logPratio);
    result += x[15]*(sor[15]-log((x[15]+1e-100))-logPratio);
    result += x[16]*(sor[16]-log((x[16]+1e-100))-logPratio);
    result += x[17]*(sor[17]-log((x[17]+1e-100))-logPratio);
    result += x[18]*(sor[18]-log((x[18]+1e-100))-logPratio);
    result += x[19]*(sor[19]-log((x[19]+1e-100))-logPratio);
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
    double gort[20]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 20; ++id) {
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
    double gort[20]; /* temporary storage */
    double x[20]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*CH3 */
    YOW += y[1]*imw[1]; /*CH4 */
    YOW += y[2]*imw[2]; /*H */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*O */
    YOW += y[6]*imw[6]; /*HCOp */
    YOW += y[7]*imw[7]; /*NO */
    YOW += y[8]*imw[8]; /*N */
    YOW += y[9]*imw[9]; /*H2 */
    YOW += y[10]*imw[10]; /*O2 */
    YOW += y[11]*imw[11]; /*CH2O */
    YOW += y[12]*imw[12]; /*CO */
    YOW += y[13]*imw[13]; /*CO2 */
    YOW += y[14]*imw[14]; /*CH */
    YOW += y[15]*imw[15]; /*H3Op */
    YOW += y[16]*imw[16]; /*N2 */
    YOW += y[17]*imw[17]; /*NO2 */
    YOW += y[18]*imw[18]; /*E */
    YOW += y[19]*imw[19]; /*CH2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(15.035060*YOW); 
    x[1] = y[1]/(16.043030*YOW); 
    x[2] = y[2]/(1.007970*YOW); 
    x[3] = y[3]/(17.007370*YOW); 
    x[4] = y[4]/(18.015340*YOW); 
    x[5] = y[5]/(15.999400*YOW); 
    x[6] = y[6]/(29.017971*YOW); 
    x[7] = y[7]/(30.006100*YOW); 
    x[8] = y[8]/(14.006700*YOW); 
    x[9] = y[9]/(2.015940*YOW); 
    x[10] = y[10]/(31.998800*YOW); 
    x[11] = y[11]/(30.026490*YOW); 
    x[12] = y[12]/(28.010550*YOW); 
    x[13] = y[13]/(44.009950*YOW); 
    x[14] = y[14]/(13.019120*YOW); 
    x[15] = y[15]/(19.022761*YOW); 
    x[16] = y[16]/(28.013400*YOW); 
    x[17] = y[17]/(46.005500*YOW); 
    x[18] = y[18]/(0.000549*YOW); 
    x[19] = y[19]/(14.027090*YOW); 
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
    result += x[9]*(gort[9]+log((x[9]+1e-100))+logPratio);
    result += x[10]*(gort[10]+log((x[10]+1e-100))+logPratio);
    result += x[11]*(gort[11]+log((x[11]+1e-100))+logPratio);
    result += x[12]*(gort[12]+log((x[12]+1e-100))+logPratio);
    result += x[13]*(gort[13]+log((x[13]+1e-100))+logPratio);
    result += x[14]*(gort[14]+log((x[14]+1e-100))+logPratio);
    result += x[15]*(gort[15]+log((x[15]+1e-100))+logPratio);
    result += x[16]*(gort[16]+log((x[16]+1e-100))+logPratio);
    result += x[17]*(gort[17]+log((x[17]+1e-100))+logPratio);
    result += x[18]*(gort[18]+log((x[18]+1e-100))+logPratio);
    result += x[19]*(gort[19]+log((x[19]+1e-100))+logPratio);
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
    double aort[20]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 20; ++id) {
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
    double aort[20]; /* temporary storage */
    double x[20]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*CH3 */
    YOW += y[1]*imw[1]; /*CH4 */
    YOW += y[2]*imw[2]; /*H */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*O */
    YOW += y[6]*imw[6]; /*HCOp */
    YOW += y[7]*imw[7]; /*NO */
    YOW += y[8]*imw[8]; /*N */
    YOW += y[9]*imw[9]; /*H2 */
    YOW += y[10]*imw[10]; /*O2 */
    YOW += y[11]*imw[11]; /*CH2O */
    YOW += y[12]*imw[12]; /*CO */
    YOW += y[13]*imw[13]; /*CO2 */
    YOW += y[14]*imw[14]; /*CH */
    YOW += y[15]*imw[15]; /*H3Op */
    YOW += y[16]*imw[16]; /*N2 */
    YOW += y[17]*imw[17]; /*NO2 */
    YOW += y[18]*imw[18]; /*E */
    YOW += y[19]*imw[19]; /*CH2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(15.035060*YOW); 
    x[1] = y[1]/(16.043030*YOW); 
    x[2] = y[2]/(1.007970*YOW); 
    x[3] = y[3]/(17.007370*YOW); 
    x[4] = y[4]/(18.015340*YOW); 
    x[5] = y[5]/(15.999400*YOW); 
    x[6] = y[6]/(29.017971*YOW); 
    x[7] = y[7]/(30.006100*YOW); 
    x[8] = y[8]/(14.006700*YOW); 
    x[9] = y[9]/(2.015940*YOW); 
    x[10] = y[10]/(31.998800*YOW); 
    x[11] = y[11]/(30.026490*YOW); 
    x[12] = y[12]/(28.010550*YOW); 
    x[13] = y[13]/(44.009950*YOW); 
    x[14] = y[14]/(13.019120*YOW); 
    x[15] = y[15]/(19.022761*YOW); 
    x[16] = y[16]/(28.013400*YOW); 
    x[17] = y[17]/(46.005500*YOW); 
    x[18] = y[18]/(0.000549*YOW); 
    x[19] = y[19]/(14.027090*YOW); 
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
    result += x[9]*(aort[9]+log((x[9]+1e-100))+logPratio);
    result += x[10]*(aort[10]+log((x[10]+1e-100))+logPratio);
    result += x[11]*(aort[11]+log((x[11]+1e-100))+logPratio);
    result += x[12]*(aort[12]+log((x[12]+1e-100))+logPratio);
    result += x[13]*(aort[13]+log((x[13]+1e-100))+logPratio);
    result += x[14]*(aort[14]+log((x[14]+1e-100))+logPratio);
    result += x[15]*(aort[15]+log((x[15]+1e-100))+logPratio);
    result += x[16]*(aort[16]+log((x[16]+1e-100))+logPratio);
    result += x[17]*(aort[17]+log((x[17]+1e-100))+logPratio);
    result += x[18]*(aort[18]+log((x[18]+1e-100))+logPratio);
    result += x[19]*(aort[19]+log((x[19]+1e-100))+logPratio);
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
void CKWC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 20; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 20; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[20]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*CH3 */
    YOW += y[1]*imw[1]; /*CH4 */
    YOW += y[2]*imw[2]; /*H */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*O */
    YOW += y[6]*imw[6]; /*HCOp */
    YOW += y[7]*imw[7]; /*NO */
    YOW += y[8]*imw[8]; /*N */
    YOW += y[9]*imw[9]; /*H2 */
    YOW += y[10]*imw[10]; /*O2 */
    YOW += y[11]*imw[11]; /*CH2O */
    YOW += y[12]*imw[12]; /*CO */
    YOW += y[13]*imw[13]; /*CO2 */
    YOW += y[14]*imw[14]; /*CH */
    YOW += y[15]*imw[15]; /*H3Op */
    YOW += y[16]*imw[16]; /*N2 */
    YOW += y[17]*imw[17]; /*NO2 */
    YOW += y[18]*imw[18]; /*E */
    YOW += y[19]*imw[19]; /*CH2 */
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
    c[9] = PWORT * y[9]*imw[9]; 
    c[10] = PWORT * y[10]*imw[10]; 
    c[11] = PWORT * y[11]*imw[11]; 
    c[12] = PWORT * y[12]*imw[12]; 
    c[13] = PWORT * y[13]*imw[13]; 
    c[14] = PWORT * y[14]*imw[14]; 
    c[15] = PWORT * y[15]*imw[15]; 
    c[16] = PWORT * y[16]*imw[16]; 
    c[17] = PWORT * y[17]*imw[17]; 
    c[18] = PWORT * y[18]*imw[18]; 
    c[19] = PWORT * y[19]*imw[19]; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 20; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[20]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 20; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 20; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void CKWYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[20]; /*temporary storage */
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
    c[9] = 1e6 * (*rho) * y[9]*imw[9]; 
    c[10] = 1e6 * (*rho) * y[10]*imw[10]; 
    c[11] = 1e6 * (*rho) * y[11]*imw[11]; 
    c[12] = 1e6 * (*rho) * y[12]*imw[12]; 
    c[13] = 1e6 * (*rho) * y[13]*imw[13]; 
    c[14] = 1e6 * (*rho) * y[14]*imw[14]; 
    c[15] = 1e6 * (*rho) * y[15]*imw[15]; 
    c[16] = 1e6 * (*rho) * y[16]*imw[16]; 
    c[17] = 1e6 * (*rho) * y[17]*imw[17]; 
    c[18] = 1e6 * (*rho) * y[18]*imw[18]; 
    c[19] = 1e6 * (*rho) * y[19]*imw[19]; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 20; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void VCKWYR(int * restrict np, double * restrict rho, double * restrict T,
	    double * restrict y, int * restrict iwrk, double * restrict rwrk,
	    double * restrict wdot)
{
    double c[20*(*np)]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    for (int n=0; n<20; n++) {
        for (int i=0; i<(*np); i++) {
            c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];
        }
    }

    /*call productionRate */
    vproductionRate(*np, wdot, c, T);

    /*convert to chemkin units */
    for (int i=0; i<20*(*np); i++) {
        wdot[i] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict wdot)
{
    int id; /*loop counter */
    double c[20]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*15.035060; /*CH3 */
    XW += x[1]*16.043030; /*CH4 */
    XW += x[2]*1.007970; /*H */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*15.999400; /*O */
    XW += x[6]*29.017971; /*HCOp */
    XW += x[7]*30.006100; /*NO */
    XW += x[8]*14.006700; /*N */
    XW += x[9]*2.015940; /*H2 */
    XW += x[10]*31.998800; /*O2 */
    XW += x[11]*30.026490; /*CH2O */
    XW += x[12]*28.010550; /*CO */
    XW += x[13]*44.009950; /*CO2 */
    XW += x[14]*13.019120; /*CH */
    XW += x[15]*19.022761; /*H3Op */
    XW += x[16]*28.013400; /*N2 */
    XW += x[17]*46.005500; /*NO2 */
    XW += x[18]*0.000549; /*E */
    XW += x[19]*14.027090; /*CH2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 20; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 20; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double * restrict T, double * restrict C, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 20; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 20; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 31; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKKFKR(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict q_f, double * restrict q_r)
{
    int id; /*loop counter */
    double c[20]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 20; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 31; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[20]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*CH3 */
    YOW += y[1]*imw[1]; /*CH4 */
    YOW += y[2]*imw[2]; /*H */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*O */
    YOW += y[6]*imw[6]; /*HCOp */
    YOW += y[7]*imw[7]; /*NO */
    YOW += y[8]*imw[8]; /*N */
    YOW += y[9]*imw[9]; /*H2 */
    YOW += y[10]*imw[10]; /*O2 */
    YOW += y[11]*imw[11]; /*CH2O */
    YOW += y[12]*imw[12]; /*CO */
    YOW += y[13]*imw[13]; /*CO2 */
    YOW += y[14]*imw[14]; /*CH */
    YOW += y[15]*imw[15]; /*H3Op */
    YOW += y[16]*imw[16]; /*N2 */
    YOW += y[17]*imw[17]; /*NO2 */
    YOW += y[18]*imw[18]; /*E */
    YOW += y[19]*imw[19]; /*CH2 */
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
    c[9] = PWORT * y[9]*imw[9]; 
    c[10] = PWORT * y[10]*imw[10]; 
    c[11] = PWORT * y[11]*imw[11]; 
    c[12] = PWORT * y[12]*imw[12]; 
    c[13] = PWORT * y[13]*imw[13]; 
    c[14] = PWORT * y[14]*imw[14]; 
    c[15] = PWORT * y[15]*imw[15]; 
    c[16] = PWORT * y[16]*imw[16]; 
    c[17] = PWORT * y[17]*imw[17]; 
    c[18] = PWORT * y[18]*imw[18]; 
    c[19] = PWORT * y[19]*imw[19]; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 31; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[20]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 20; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 31; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[20]; /*temporary storage */
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
    c[9] = 1e6 * (*rho) * y[9]*imw[9]; 
    c[10] = 1e6 * (*rho) * y[10]*imw[10]; 
    c[11] = 1e6 * (*rho) * y[11]*imw[11]; 
    c[12] = 1e6 * (*rho) * y[12]*imw[12]; 
    c[13] = 1e6 * (*rho) * y[13]*imw[13]; 
    c[14] = 1e6 * (*rho) * y[14]*imw[14]; 
    c[15] = 1e6 * (*rho) * y[15]*imw[15]; 
    c[16] = 1e6 * (*rho) * y[16]*imw[16]; 
    c[17] = 1e6 * (*rho) * y[17]*imw[17]; 
    c[18] = 1e6 * (*rho) * y[18]*imw[18]; 
    c[19] = 1e6 * (*rho) * y[19]*imw[19]; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 31; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict qdot)
{
    int id; /*loop counter */
    double c[20]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*15.035060; /*CH3 */
    XW += x[1]*16.043030; /*CH4 */
    XW += x[2]*1.007970; /*H */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*15.999400; /*O */
    XW += x[6]*29.017971; /*HCOp */
    XW += x[7]*30.006100; /*NO */
    XW += x[8]*14.006700; /*N */
    XW += x[9]*2.015940; /*H2 */
    XW += x[10]*31.998800; /*O2 */
    XW += x[11]*30.026490; /*CH2O */
    XW += x[12]*28.010550; /*CO */
    XW += x[13]*44.009950; /*CO2 */
    XW += x[14]*13.019120; /*CH */
    XW += x[15]*19.022761; /*H3Op */
    XW += x[16]*28.013400; /*N2 */
    XW += x[17]*46.005500; /*NO2 */
    XW += x[18]*0.000549; /*E */
    XW += x[19]*14.027090; /*CH2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 20; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 31; ++id) {
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
    for (id = 0; id < 20 * kd; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: H + OH + M <=> H2O + M */
    nuki[ 2 * kd + 0 ] += -1 ;
    nuki[ 3 * kd + 0 ] += -1 ;
    nuki[ 4 * kd + 0 ] += +1 ;

    /*reaction 2: O + O + M <=> O2 + M */
    nuki[ 5 * kd + 1 ] += -1 ;
    nuki[ 5 * kd + 1 ] += -1 ;
    nuki[ 10 * kd + 1 ] += +1 ;

    /*reaction 3: H + H + M <=> H2 + M */
    nuki[ 2 * kd + 2 ] += -1 ;
    nuki[ 2 * kd + 2 ] += -1 ;
    nuki[ 9 * kd + 2 ] += +1 ;

    /*reaction 4: NO + O + M <=> NO2 + M */
    nuki[ 7 * kd + 3 ] += -1 ;
    nuki[ 5 * kd + 3 ] += -1 ;
    nuki[ 17 * kd + 3 ] += +1 ;

    /*reaction 5: CH4 <=> CH3 + H */
    nuki[ 1 * kd + 4 ] += -1 ;
    nuki[ 0 * kd + 4 ] += +1 ;
    nuki[ 2 * kd + 4 ] += +1 ;

    /*reaction 6: CH4 + OH <=> CH3 + H2O */
    nuki[ 1 * kd + 5 ] += -1 ;
    nuki[ 3 * kd + 5 ] += -1 ;
    nuki[ 0 * kd + 5 ] += +1 ;
    nuki[ 4 * kd + 5 ] += +1 ;

    /*reaction 7: CH4 + O <=> CH3 + OH */
    nuki[ 1 * kd + 6 ] += -1 ;
    nuki[ 5 * kd + 6 ] += -1 ;
    nuki[ 0 * kd + 6 ] += +1 ;
    nuki[ 3 * kd + 6 ] += +1 ;

    /*reaction 8: CH4 + H <=> CH3 + H2 */
    nuki[ 1 * kd + 7 ] += -1 ;
    nuki[ 2 * kd + 7 ] += -1 ;
    nuki[ 0 * kd + 7 ] += +1 ;
    nuki[ 9 * kd + 7 ] += +1 ;

    /*reaction 9: CH3 + O2 <=> CH2O + OH */
    nuki[ 0 * kd + 8 ] += -1 ;
    nuki[ 10 * kd + 8 ] += -1 ;
    nuki[ 11 * kd + 8 ] += +1 ;
    nuki[ 3 * kd + 8 ] += +1 ;

    /*reaction 10: CH2O + OH <=> CO + H2O + H */
    nuki[ 11 * kd + 9 ] += -1 ;
    nuki[ 3 * kd + 9 ] += -1 ;
    nuki[ 12 * kd + 9 ] += +1 ;
    nuki[ 4 * kd + 9 ] += +1 ;
    nuki[ 2 * kd + 9 ] += +1 ;

    /*reaction 11: CO + OH => CO2 + H */
    nuki[ 12 * kd + 10 ] += -1 ;
    nuki[ 3 * kd + 10 ] += -1 ;
    nuki[ 13 * kd + 10 ] += +1 ;
    nuki[ 2 * kd + 10 ] += +1 ;

    /*reaction 12: CO2 + H => CO + OH */
    nuki[ 13 * kd + 11 ] += -1 ;
    nuki[ 2 * kd + 11 ] += -1 ;
    nuki[ 12 * kd + 11 ] += +1 ;
    nuki[ 3 * kd + 11 ] += +1 ;

    /*reaction 13: O2 + H => OH + O */
    nuki[ 10 * kd + 12 ] += -1 ;
    nuki[ 2 * kd + 12 ] += -1 ;
    nuki[ 3 * kd + 12 ] += +1 ;
    nuki[ 5 * kd + 12 ] += +1 ;

    /*reaction 14: OH + O => O2 + H */
    nuki[ 3 * kd + 13 ] += -1 ;
    nuki[ 5 * kd + 13 ] += -1 ;
    nuki[ 10 * kd + 13 ] += +1 ;
    nuki[ 2 * kd + 13 ] += +1 ;

    /*reaction 15: O + H2 => OH + H */
    nuki[ 5 * kd + 14 ] += -1 ;
    nuki[ 9 * kd + 14 ] += -1 ;
    nuki[ 3 * kd + 14 ] += +1 ;
    nuki[ 2 * kd + 14 ] += +1 ;

    /*reaction 16: OH + H => O + H2 */
    nuki[ 3 * kd + 15 ] += -1 ;
    nuki[ 2 * kd + 15 ] += -1 ;
    nuki[ 5 * kd + 15 ] += +1 ;
    nuki[ 9 * kd + 15 ] += +1 ;

    /*reaction 17: O + H2O => 2 OH */
    nuki[ 5 * kd + 16 ] += -1 ;
    nuki[ 4 * kd + 16 ] += -1 ;
    nuki[ 3 * kd + 16 ] += +2 ;

    /*reaction 18: 2 OH => O + H2O */
    nuki[ 3 * kd + 17 ] += -2 ;
    nuki[ 5 * kd + 17 ] += +1 ;
    nuki[ 4 * kd + 17 ] += +1 ;

    /*reaction 19: OH + H2 => H2O + H */
    nuki[ 3 * kd + 18 ] += -1 ;
    nuki[ 9 * kd + 18 ] += -1 ;
    nuki[ 4 * kd + 18 ] += +1 ;
    nuki[ 2 * kd + 18 ] += +1 ;

    /*reaction 20: H2O + H => H2 + OH */
    nuki[ 4 * kd + 19 ] += -1 ;
    nuki[ 2 * kd + 19 ] += -1 ;
    nuki[ 9 * kd + 19 ] += +1 ;
    nuki[ 3 * kd + 19 ] += +1 ;

    /*reaction 21: CH3 + O <=> CH + H2O */
    nuki[ 0 * kd + 20 ] += -1 ;
    nuki[ 5 * kd + 20 ] += -1 ;
    nuki[ 14 * kd + 20 ] += +1 ;
    nuki[ 4 * kd + 20 ] += +1 ;

    /*reaction 22: CH + O <=> HCOp + E */
    nuki[ 14 * kd + 21 ] += -1 ;
    nuki[ 5 * kd + 21 ] += -1 ;
    nuki[ 6 * kd + 21 ] += +1 ;
    nuki[ 18 * kd + 21 ] += +1 ;

    /*reaction 23: HCOp + H2O <=> CO + H3Op */
    nuki[ 6 * kd + 22 ] += -1 ;
    nuki[ 4 * kd + 22 ] += -1 ;
    nuki[ 12 * kd + 22 ] += +1 ;
    nuki[ 15 * kd + 22 ] += +1 ;

    /*reaction 24: H3Op + E <=> H2O + H */
    nuki[ 15 * kd + 23 ] += -1 ;
    nuki[ 18 * kd + 23 ] += -1 ;
    nuki[ 4 * kd + 23 ] += +1 ;
    nuki[ 2 * kd + 23 ] += +1 ;

    /*reaction 25: CH + O2 <=> CO + OH */
    nuki[ 14 * kd + 24 ] += -1 ;
    nuki[ 10 * kd + 24 ] += -1 ;
    nuki[ 12 * kd + 24 ] += +1 ;
    nuki[ 3 * kd + 24 ] += +1 ;

    /*reaction 26: O + N2 => NO + N */
    nuki[ 5 * kd + 25 ] += -1 ;
    nuki[ 16 * kd + 25 ] += -1 ;
    nuki[ 7 * kd + 25 ] += +1 ;
    nuki[ 8 * kd + 25 ] += +1 ;

    /*reaction 27: NO + N => O + N2 */
    nuki[ 7 * kd + 26 ] += -1 ;
    nuki[ 8 * kd + 26 ] += -1 ;
    nuki[ 5 * kd + 26 ] += +1 ;
    nuki[ 16 * kd + 26 ] += +1 ;

    /*reaction 28: N + O2 => NO + O */
    nuki[ 8 * kd + 27 ] += -1 ;
    nuki[ 10 * kd + 27 ] += -1 ;
    nuki[ 7 * kd + 27 ] += +1 ;
    nuki[ 5 * kd + 27 ] += +1 ;

    /*reaction 29: NO + O => N + O2 */
    nuki[ 7 * kd + 28 ] += -1 ;
    nuki[ 5 * kd + 28 ] += -1 ;
    nuki[ 8 * kd + 28 ] += +1 ;
    nuki[ 10 * kd + 28 ] += +1 ;

    /*reaction 30: NO2 + O <=> NO + O2 */
    nuki[ 17 * kd + 29 ] += -1 ;
    nuki[ 5 * kd + 29 ] += -1 ;
    nuki[ 7 * kd + 29 ] += +1 ;
    nuki[ 10 * kd + 29 ] += +1 ;

    /*reaction 31: NO2 + H <=> NO + OH */
    nuki[ 17 * kd + 30 ] += -1 ;
    nuki[ 2 * kd + 30 ] += -1 ;
    nuki[ 7 * kd + 30 ] += +1 ;
    nuki[ 3 * kd + 30 ] += +1 ;
}


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim, int * iwrk, double * restrict rwrk, int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < kd * 20; ++ id) {
         ncf[id] = 0; 
    }

    /*CH3 */
    ncf[ 0 * kd + 2 ] = 1; /*C */
    ncf[ 0 * kd + 1 ] = 3; /*H */

    /*CH4 */
    ncf[ 1 * kd + 2 ] = 1; /*C */
    ncf[ 1 * kd + 1 ] = 4; /*H */

    /*H */
    ncf[ 2 * kd + 1 ] = 1; /*H */

    /*OH */
    ncf[ 3 * kd + 0 ] = 1; /*O */
    ncf[ 3 * kd + 1 ] = 1; /*H */

    /*H2O */
    ncf[ 4 * kd + 1 ] = 2; /*H */
    ncf[ 4 * kd + 0 ] = 1; /*O */

    /*O */
    ncf[ 5 * kd + 0 ] = 1; /*O */

    /*HCOp */
    ncf[ 6 * kd + 2 ] = 1; /*C */
    ncf[ 6 * kd + 1 ] = 1; /*H */
    ncf[ 6 * kd + 0 ] = 1; /*O */
    ncf[ 6 * kd + 4 ] = -1; /*E */

    /*NO */
    ncf[ 7 * kd + 3 ] = 1; /*N */
    ncf[ 7 * kd + 0 ] = 1; /*O */

    /*N */
    ncf[ 8 * kd + 3 ] = 1; /*N */

    /*H2 */
    ncf[ 9 * kd + 1 ] = 2; /*H */

    /*O2 */
    ncf[ 10 * kd + 0 ] = 2; /*O */

    /*CH2O */
    ncf[ 11 * kd + 1 ] = 2; /*H */
    ncf[ 11 * kd + 2 ] = 1; /*C */
    ncf[ 11 * kd + 0 ] = 1; /*O */

    /*CO */
    ncf[ 12 * kd + 2 ] = 1; /*C */
    ncf[ 12 * kd + 0 ] = 1; /*O */

    /*CO2 */
    ncf[ 13 * kd + 2 ] = 1; /*C */
    ncf[ 13 * kd + 0 ] = 2; /*O */

    /*CH */
    ncf[ 14 * kd + 2 ] = 1; /*C */
    ncf[ 14 * kd + 1 ] = 1; /*H */

    /*H3Op */
    ncf[ 15 * kd + 1 ] = 3; /*H */
    ncf[ 15 * kd + 0 ] = 1; /*O */
    ncf[ 15 * kd + 4 ] = -1; /*E */

    /*N2 */
    ncf[ 16 * kd + 3 ] = 2; /*N */

    /*NO2 */
    ncf[ 17 * kd + 3 ] = 1; /*N */
    ncf[ 17 * kd + 0 ] = 2; /*O */

    /*E */
    ncf[ 18 * kd + 4 ] = 1; /*E */

    /*CH2 */
    ncf[ 19 * kd + 2 ] = 1; /*C */
    ncf[ 19 * kd + 1 ] = 2; /*H */

}


/*Returns the arrehenius coefficients  */
/*for all reactions */
void CKABE(int * iwrk, double * restrict rwrk, double * restrict a, double * restrict b, double * restrict e)
{
    for (int i=0; i<31; ++i) {
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
    double gort[20]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + OH + M <=> H2O + M */
    eqcon[0] *= 1e+06; 

    /*reaction 2: O + O + M <=> O2 + M */
    eqcon[1] *= 1e+06; 

    /*reaction 3: H + H + M <=> H2 + M */
    eqcon[2] *= 1e+06; 

    /*reaction 4: NO + O + M <=> NO2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: CH4 <=> CH3 + H */
    eqcon[4] *= 1e-06; 

    /*reaction 6: CH4 + OH <=> CH3 + H2O */
    /*eqcon[5] *= 1;  */

    /*reaction 7: CH4 + O <=> CH3 + OH */
    /*eqcon[6] *= 1;  */

    /*reaction 8: CH4 + H <=> CH3 + H2 */
    /*eqcon[7] *= 1;  */

    /*reaction 9: CH3 + O2 <=> CH2O + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: CH2O + OH <=> CO + H2O + H */
    eqcon[9] *= 1e-06; 

    /*reaction 11: CO + OH => CO2 + H */
    /*eqcon[10] *= 1;  */

    /*reaction 12: CO2 + H => CO + OH */
    /*eqcon[11] *= 1;  */

    /*reaction 13: O2 + H => OH + O */
    /*eqcon[12] *= 1;  */

    /*reaction 14: OH + O => O2 + H */
    /*eqcon[13] *= 1;  */

    /*reaction 15: O + H2 => OH + H */
    /*eqcon[14] *= 1;  */

    /*reaction 16: OH + H => O + H2 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: O + H2O => 2 OH */
    /*eqcon[16] *= 1;  */

    /*reaction 18: 2 OH => O + H2O */
    /*eqcon[17] *= 1;  */

    /*reaction 19: OH + H2 => H2O + H */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H2O + H => H2 + OH */
    /*eqcon[19] *= 1;  */

    /*reaction 21: CH3 + O <=> CH + H2O */
    /*eqcon[20] *= 1;  */

    /*reaction 22: CH + O <=> HCOp + E */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HCOp + H2O <=> CO + H3Op */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H3Op + E <=> H2O + H */
    /*eqcon[23] *= 1;  */

    /*reaction 25: CH + O2 <=> CO + OH */
    /*eqcon[24] *= 1;  */

    /*reaction 26: O + N2 => NO + N */
    /*eqcon[25] *= 1;  */

    /*reaction 27: NO + N => O + N2 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: N + O2 => NO + O */
    /*eqcon[27] *= 1;  */

    /*reaction 29: NO + O => N + O2 */
    /*eqcon[28] *= 1;  */

    /*reaction 30: NO2 + O <=> NO + O2 */
    /*eqcon[29] *= 1;  */

    /*reaction 31: NO2 + H <=> NO + OH */
    /*eqcon[30] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[20]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + OH + M <=> H2O + M */
    eqcon[0] *= 1e+06; 

    /*reaction 2: O + O + M <=> O2 + M */
    eqcon[1] *= 1e+06; 

    /*reaction 3: H + H + M <=> H2 + M */
    eqcon[2] *= 1e+06; 

    /*reaction 4: NO + O + M <=> NO2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: CH4 <=> CH3 + H */
    eqcon[4] *= 1e-06; 

    /*reaction 6: CH4 + OH <=> CH3 + H2O */
    /*eqcon[5] *= 1;  */

    /*reaction 7: CH4 + O <=> CH3 + OH */
    /*eqcon[6] *= 1;  */

    /*reaction 8: CH4 + H <=> CH3 + H2 */
    /*eqcon[7] *= 1;  */

    /*reaction 9: CH3 + O2 <=> CH2O + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: CH2O + OH <=> CO + H2O + H */
    eqcon[9] *= 1e-06; 

    /*reaction 11: CO + OH => CO2 + H */
    /*eqcon[10] *= 1;  */

    /*reaction 12: CO2 + H => CO + OH */
    /*eqcon[11] *= 1;  */

    /*reaction 13: O2 + H => OH + O */
    /*eqcon[12] *= 1;  */

    /*reaction 14: OH + O => O2 + H */
    /*eqcon[13] *= 1;  */

    /*reaction 15: O + H2 => OH + H */
    /*eqcon[14] *= 1;  */

    /*reaction 16: OH + H => O + H2 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: O + H2O => 2 OH */
    /*eqcon[16] *= 1;  */

    /*reaction 18: 2 OH => O + H2O */
    /*eqcon[17] *= 1;  */

    /*reaction 19: OH + H2 => H2O + H */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H2O + H => H2 + OH */
    /*eqcon[19] *= 1;  */

    /*reaction 21: CH3 + O <=> CH + H2O */
    /*eqcon[20] *= 1;  */

    /*reaction 22: CH + O <=> HCOp + E */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HCOp + H2O <=> CO + H3Op */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H3Op + E <=> H2O + H */
    /*eqcon[23] *= 1;  */

    /*reaction 25: CH + O2 <=> CO + OH */
    /*eqcon[24] *= 1;  */

    /*reaction 26: O + N2 => NO + N */
    /*eqcon[25] *= 1;  */

    /*reaction 27: NO + N => O + N2 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: N + O2 => NO + O */
    /*eqcon[27] *= 1;  */

    /*reaction 29: NO + O => N + O2 */
    /*eqcon[28] *= 1;  */

    /*reaction 30: NO2 + O <=> NO + O2 */
    /*eqcon[29] *= 1;  */

    /*reaction 31: NO2 + H <=> NO + OH */
    /*eqcon[30] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double * restrict P, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[20]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + OH + M <=> H2O + M */
    eqcon[0] *= 1e+06; 

    /*reaction 2: O + O + M <=> O2 + M */
    eqcon[1] *= 1e+06; 

    /*reaction 3: H + H + M <=> H2 + M */
    eqcon[2] *= 1e+06; 

    /*reaction 4: NO + O + M <=> NO2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: CH4 <=> CH3 + H */
    eqcon[4] *= 1e-06; 

    /*reaction 6: CH4 + OH <=> CH3 + H2O */
    /*eqcon[5] *= 1;  */

    /*reaction 7: CH4 + O <=> CH3 + OH */
    /*eqcon[6] *= 1;  */

    /*reaction 8: CH4 + H <=> CH3 + H2 */
    /*eqcon[7] *= 1;  */

    /*reaction 9: CH3 + O2 <=> CH2O + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: CH2O + OH <=> CO + H2O + H */
    eqcon[9] *= 1e-06; 

    /*reaction 11: CO + OH => CO2 + H */
    /*eqcon[10] *= 1;  */

    /*reaction 12: CO2 + H => CO + OH */
    /*eqcon[11] *= 1;  */

    /*reaction 13: O2 + H => OH + O */
    /*eqcon[12] *= 1;  */

    /*reaction 14: OH + O => O2 + H */
    /*eqcon[13] *= 1;  */

    /*reaction 15: O + H2 => OH + H */
    /*eqcon[14] *= 1;  */

    /*reaction 16: OH + H => O + H2 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: O + H2O => 2 OH */
    /*eqcon[16] *= 1;  */

    /*reaction 18: 2 OH => O + H2O */
    /*eqcon[17] *= 1;  */

    /*reaction 19: OH + H2 => H2O + H */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H2O + H => H2 + OH */
    /*eqcon[19] *= 1;  */

    /*reaction 21: CH3 + O <=> CH + H2O */
    /*eqcon[20] *= 1;  */

    /*reaction 22: CH + O <=> HCOp + E */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HCOp + H2O <=> CO + H3Op */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H3Op + E <=> H2O + H */
    /*eqcon[23] *= 1;  */

    /*reaction 25: CH + O2 <=> CO + OH */
    /*eqcon[24] *= 1;  */

    /*reaction 26: O + N2 => NO + N */
    /*eqcon[25] *= 1;  */

    /*reaction 27: NO + N => O + N2 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: N + O2 => NO + O */
    /*eqcon[27] *= 1;  */

    /*reaction 29: NO + O => N + O2 */
    /*eqcon[28] *= 1;  */

    /*reaction 30: NO2 + O <=> NO + O2 */
    /*eqcon[29] *= 1;  */

    /*reaction 31: NO2 + H <=> NO + OH */
    /*eqcon[30] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[20]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + OH + M <=> H2O + M */
    eqcon[0] *= 1e+06; 

    /*reaction 2: O + O + M <=> O2 + M */
    eqcon[1] *= 1e+06; 

    /*reaction 3: H + H + M <=> H2 + M */
    eqcon[2] *= 1e+06; 

    /*reaction 4: NO + O + M <=> NO2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: CH4 <=> CH3 + H */
    eqcon[4] *= 1e-06; 

    /*reaction 6: CH4 + OH <=> CH3 + H2O */
    /*eqcon[5] *= 1;  */

    /*reaction 7: CH4 + O <=> CH3 + OH */
    /*eqcon[6] *= 1;  */

    /*reaction 8: CH4 + H <=> CH3 + H2 */
    /*eqcon[7] *= 1;  */

    /*reaction 9: CH3 + O2 <=> CH2O + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: CH2O + OH <=> CO + H2O + H */
    eqcon[9] *= 1e-06; 

    /*reaction 11: CO + OH => CO2 + H */
    /*eqcon[10] *= 1;  */

    /*reaction 12: CO2 + H => CO + OH */
    /*eqcon[11] *= 1;  */

    /*reaction 13: O2 + H => OH + O */
    /*eqcon[12] *= 1;  */

    /*reaction 14: OH + O => O2 + H */
    /*eqcon[13] *= 1;  */

    /*reaction 15: O + H2 => OH + H */
    /*eqcon[14] *= 1;  */

    /*reaction 16: OH + H => O + H2 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: O + H2O => 2 OH */
    /*eqcon[16] *= 1;  */

    /*reaction 18: 2 OH => O + H2O */
    /*eqcon[17] *= 1;  */

    /*reaction 19: OH + H2 => H2O + H */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H2O + H => H2 + OH */
    /*eqcon[19] *= 1;  */

    /*reaction 21: CH3 + O <=> CH + H2O */
    /*eqcon[20] *= 1;  */

    /*reaction 22: CH + O <=> HCOp + E */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HCOp + H2O <=> CO + H3Op */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H3Op + E <=> H2O + H */
    /*eqcon[23] *= 1;  */

    /*reaction 25: CH + O2 <=> CO + OH */
    /*eqcon[24] *= 1;  */

    /*reaction 26: O + N2 => NO + N */
    /*eqcon[25] *= 1;  */

    /*reaction 27: NO + N => O + N2 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: N + O2 => NO + O */
    /*eqcon[27] *= 1;  */

    /*reaction 29: NO + O => N + O2 */
    /*eqcon[28] *= 1;  */

    /*reaction 30: NO2 + O <=> NO + O2 */
    /*eqcon[29] *= 1;  */

    /*reaction 31: NO2 + H <=> NO + OH */
    /*eqcon[30] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double * restrict rho, double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[20]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + OH + M <=> H2O + M */
    eqcon[0] *= 1e+06; 

    /*reaction 2: O + O + M <=> O2 + M */
    eqcon[1] *= 1e+06; 

    /*reaction 3: H + H + M <=> H2 + M */
    eqcon[2] *= 1e+06; 

    /*reaction 4: NO + O + M <=> NO2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: CH4 <=> CH3 + H */
    eqcon[4] *= 1e-06; 

    /*reaction 6: CH4 + OH <=> CH3 + H2O */
    /*eqcon[5] *= 1;  */

    /*reaction 7: CH4 + O <=> CH3 + OH */
    /*eqcon[6] *= 1;  */

    /*reaction 8: CH4 + H <=> CH3 + H2 */
    /*eqcon[7] *= 1;  */

    /*reaction 9: CH3 + O2 <=> CH2O + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: CH2O + OH <=> CO + H2O + H */
    eqcon[9] *= 1e-06; 

    /*reaction 11: CO + OH => CO2 + H */
    /*eqcon[10] *= 1;  */

    /*reaction 12: CO2 + H => CO + OH */
    /*eqcon[11] *= 1;  */

    /*reaction 13: O2 + H => OH + O */
    /*eqcon[12] *= 1;  */

    /*reaction 14: OH + O => O2 + H */
    /*eqcon[13] *= 1;  */

    /*reaction 15: O + H2 => OH + H */
    /*eqcon[14] *= 1;  */

    /*reaction 16: OH + H => O + H2 */
    /*eqcon[15] *= 1;  */

    /*reaction 17: O + H2O => 2 OH */
    /*eqcon[16] *= 1;  */

    /*reaction 18: 2 OH => O + H2O */
    /*eqcon[17] *= 1;  */

    /*reaction 19: OH + H2 => H2O + H */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H2O + H => H2 + OH */
    /*eqcon[19] *= 1;  */

    /*reaction 21: CH3 + O <=> CH + H2O */
    /*eqcon[20] *= 1;  */

    /*reaction 22: CH + O <=> HCOp + E */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HCOp + H2O <=> CO + H3Op */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H3Op + E <=> H2O + H */
    /*eqcon[23] *= 1;  */

    /*reaction 25: CH + O2 <=> CO + OH */
    /*eqcon[24] *= 1;  */

    /*reaction 26: O + N2 => NO + N */
    /*eqcon[25] *= 1;  */

    /*reaction 27: NO + N => O + N2 */
    /*eqcon[26] *= 1;  */

    /*reaction 28: N + O2 => NO + O */
    /*eqcon[27] *= 1;  */

    /*reaction 29: NO + O => N + O2 */
    /*eqcon[28] *= 1;  */

    /*reaction 30: NO2 + O <=> NO + O2 */
    /*eqcon[29] *= 1;  */

    /*reaction 31: NO2 + H <=> NO + OH */
    /*eqcon[30] *= 1;  */
}


/*Returns the electronic charges of the species */
void CKCHRG(int * restrict iwrk, double * restrict rwrk, int * restrict kcharge)
{
    kcharge[0] = 0; /* CH3 */
    kcharge[1] = 0; /* CH4 */
    kcharge[2] = 0; /* H */
    kcharge[3] = 0; /* OH */
    kcharge[4] = 0; /* H2O */
    kcharge[5] = 0; /* O */
    kcharge[6] = 1; /* HCOp */
    kcharge[7] = 0; /* NO */
    kcharge[8] = 0; /* N */
    kcharge[9] = 0; /* H2 */
    kcharge[10] = 0; /* O2 */
    kcharge[11] = 0; /* CH2O */
    kcharge[12] = 0; /* CO */
    kcharge[13] = 0; /* CO2 */
    kcharge[14] = 0; /* CH */
    kcharge[15] = 1; /* H3Op */
    kcharge[16] = 0; /* N2 */
    kcharge[17] = 0; /* NO2 */
    kcharge[18] = -1; /* E */
    kcharge[19] = 0; /* CH2 */
}

static double T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static double k_f_save[31];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static double Kc_save[31];
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

    double qdot, q_f[31], q_r[31];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 20; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[2] -= qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[1]-q_r[1];
    wdot[5] -= qdot;
    wdot[5] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[2]-q_r[2];
    wdot[2] -= qdot;
    wdot[2] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[3]-q_r[3];
    wdot[5] -= qdot;
    wdot[7] -= qdot;
    wdot[17] += qdot;

    qdot = q_f[4]-q_r[4];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[2] += qdot;

    qdot = q_f[5]-q_r[5];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[6]-q_r[6];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[5] -= qdot;

    qdot = q_f[7]-q_r[7];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[2] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[8]-q_r[8];
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[10] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[9]-q_r[9];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[11] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[10]-q_r[10];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[12] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[11]-q_r[11];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[12] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[12]-q_r[12];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[5] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[13]-q_r[13];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[5] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[14]-q_r[14];
    wdot[2] += qdot;
    wdot[3] += qdot;
    wdot[5] -= qdot;
    wdot[9] -= qdot;

    qdot = q_f[15]-q_r[15];
    wdot[2] -= qdot;
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[9] += qdot;

    qdot = q_f[16]-q_r[16];
    wdot[3] += 2 * qdot;
    wdot[4] -= qdot;
    wdot[5] -= qdot;

    qdot = q_f[17]-q_r[17];
    wdot[3] -= 2 * qdot;
    wdot[4] += qdot;
    wdot[5] += qdot;

    qdot = q_f[18]-q_r[18];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[9] -= qdot;

    qdot = q_f[19]-q_r[19];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[20]-q_r[20];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[5] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[21]-q_r[21];
    wdot[5] -= qdot;
    wdot[6] += qdot;
    wdot[14] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[22]-q_r[22];
    wdot[4] -= qdot;
    wdot[6] -= qdot;
    wdot[12] += qdot;
    wdot[15] += qdot;

    qdot = q_f[23]-q_r[23];
    wdot[2] += qdot;
    wdot[4] += qdot;
    wdot[15] -= qdot;
    wdot[18] -= qdot;

    qdot = q_f[24]-q_r[24];
    wdot[3] += qdot;
    wdot[10] -= qdot;
    wdot[12] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[25]-q_r[25];
    wdot[5] -= qdot;
    wdot[7] += qdot;
    wdot[8] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[26]-q_r[26];
    wdot[5] += qdot;
    wdot[7] -= qdot;
    wdot[8] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[27]-q_r[27];
    wdot[5] += qdot;
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[10] -= qdot;

    qdot = q_f[28]-q_r[28];
    wdot[5] -= qdot;
    wdot[7] -= qdot;
    wdot[8] += qdot;
    wdot[10] += qdot;

    qdot = q_f[29]-q_r[29];
    wdot[5] -= qdot;
    wdot[7] += qdot;
    wdot[10] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[30]-q_r[30];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[7] += qdot;
    wdot[17] -= qdot;

    return;
}

void comp_k_f(double * restrict tc, double invT, double * restrict k_f)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<31; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(double * restrict tc, double invT, double * restrict Kc)
{
    /*compute the Gibbs free energy */
    double g_RT[20];
    gibbs(g_RT, tc);

    Kc[0] = g_RT[2] + g_RT[3] - g_RT[4];
    Kc[1] = g_RT[5] + g_RT[5] - g_RT[10];
    Kc[2] = g_RT[2] + g_RT[2] - g_RT[9];
    Kc[3] = g_RT[5] + g_RT[7] - g_RT[17];
    Kc[4] = -g_RT[0] + g_RT[1] - g_RT[2];
    Kc[5] = -g_RT[0] + g_RT[1] + g_RT[3] - g_RT[4];
    Kc[6] = -g_RT[0] + g_RT[1] - g_RT[3] + g_RT[5];
    Kc[7] = -g_RT[0] + g_RT[1] + g_RT[2] - g_RT[9];
    Kc[8] = g_RT[0] - g_RT[3] + g_RT[10] - g_RT[11];
    Kc[9] = -g_RT[2] + g_RT[3] - g_RT[4] + g_RT[11] - g_RT[12];
    Kc[10] = -g_RT[2] + g_RT[3] + g_RT[12] - g_RT[13];
    Kc[11] = g_RT[2] - g_RT[3] - g_RT[12] + g_RT[13];
    Kc[12] = g_RT[2] - g_RT[3] - g_RT[5] + g_RT[10];
    Kc[13] = -g_RT[2] + g_RT[3] + g_RT[5] - g_RT[10];
    Kc[14] = -g_RT[2] - g_RT[3] + g_RT[5] + g_RT[9];
    Kc[15] = g_RT[2] + g_RT[3] - g_RT[5] - g_RT[9];
    Kc[16] = -2*g_RT[3] + g_RT[4] + g_RT[5];
    Kc[17] = 2*g_RT[3] - g_RT[4] - g_RT[5];
    Kc[18] = -g_RT[2] + g_RT[3] - g_RT[4] + g_RT[9];
    Kc[19] = g_RT[2] - g_RT[3] + g_RT[4] - g_RT[9];
    Kc[20] = g_RT[0] - g_RT[4] + g_RT[5] - g_RT[14];
    Kc[21] = g_RT[5] - g_RT[6] + g_RT[14] - g_RT[18];
    Kc[22] = g_RT[4] + g_RT[6] - g_RT[12] - g_RT[15];
    Kc[23] = -g_RT[2] - g_RT[4] + g_RT[15] + g_RT[18];
    Kc[24] = -g_RT[3] + g_RT[10] - g_RT[12] + g_RT[14];
    Kc[25] = g_RT[5] - g_RT[7] - g_RT[8] + g_RT[16];
    Kc[26] = -g_RT[5] + g_RT[7] + g_RT[8] - g_RT[16];
    Kc[27] = -g_RT[5] - g_RT[7] + g_RT[8] + g_RT[10];
    Kc[28] = g_RT[5] + g_RT[7] - g_RT[8] - g_RT[10];
    Kc[29] = g_RT[5] - g_RT[7] - g_RT[10] + g_RT[17];
    Kc[30] = g_RT[2] - g_RT[3] - g_RT[7] + g_RT[17];

#ifdef __INTEL_COMPILER
     #pragma simd
#endif
    for (int i=0; i<31; ++i) {
        Kc[i] = exp(Kc[i]);
    };

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 * invT;
    double refCinv = 1 / refC;

    Kc[0] *= refCinv;
    Kc[1] *= refCinv;
    Kc[2] *= refCinv;
    Kc[3] *= refCinv;
    Kc[4] *= refC;
    Kc[9] *= refC;

    return;
}

void comp_qfqr(double * restrict qf, double * restrict qr, double * restrict sc, double * restrict tc, double invT)
{

    /*reaction 1: H + OH + M <=> H2O + M */
    qf[0] = sc[2]*sc[3];
    qr[0] = sc[4];

    /*reaction 2: O + O + M <=> O2 + M */
    qf[1] = sc[5]*sc[5];
    qr[1] = sc[10];

    /*reaction 3: H + H + M <=> H2 + M */
    qf[2] = sc[2]*sc[2];
    qr[2] = sc[9];

    /*reaction 4: NO + O + M <=> NO2 + M */
    qf[3] = sc[5]*sc[7];
    qr[3] = sc[17];

    /*reaction 5: CH4 <=> CH3 + H */
    qf[4] = sc[1];
    qr[4] = sc[0]*sc[2];

    /*reaction 6: CH4 + OH <=> CH3 + H2O */
    qf[5] = sc[1]*sc[3];
    qr[5] = sc[0]*sc[4];

    /*reaction 7: CH4 + O <=> CH3 + OH */
    qf[6] = sc[1]*sc[5];
    qr[6] = sc[0]*sc[3];

    /*reaction 8: CH4 + H <=> CH3 + H2 */
    qf[7] = sc[1]*sc[2];
    qr[7] = sc[0]*sc[9];

    /*reaction 9: CH3 + O2 <=> CH2O + OH */
    qf[8] = sc[0]*sc[10];
    qr[8] = sc[3]*sc[11];

    /*reaction 10: CH2O + OH <=> CO + H2O + H */
    qf[9] = sc[3]*sc[11];
    qr[9] = sc[2]*sc[4]*sc[12];

    /*reaction 11: CO + OH => CO2 + H */
    qf[10] = sc[3]*sc[12];
    qr[10] = 0.0;

    /*reaction 12: CO2 + H => CO + OH */
    qf[11] = sc[2]*sc[13];
    qr[11] = 0.0;

    /*reaction 13: O2 + H => OH + O */
    qf[12] = sc[2]*sc[10];
    qr[12] = 0.0;

    /*reaction 14: OH + O => O2 + H */
    qf[13] = sc[3]*sc[5];
    qr[13] = 0.0;

    /*reaction 15: O + H2 => OH + H */
    qf[14] = sc[5]*sc[9];
    qr[14] = 0.0;

    /*reaction 16: OH + H => O + H2 */
    qf[15] = sc[2]*sc[3];
    qr[15] = 0.0;

    /*reaction 17: O + H2O => 2 OH */
    qf[16] = sc[4]*sc[5];
    qr[16] = 0.0;

    /*reaction 18: 2 OH => O + H2O */
    qf[17] = sc[3]*sc[3];
    qr[17] = 0.0;

    /*reaction 19: OH + H2 => H2O + H */
    qf[18] = sc[3]*sc[9];
    qr[18] = 0.0;

    /*reaction 20: H2O + H => H2 + OH */
    qf[19] = sc[2]*sc[4];
    qr[19] = 0.0;

    /*reaction 21: CH3 + O <=> CH + H2O */
    qf[20] = sc[0]*sc[5];
    qr[20] = sc[4]*sc[14];

    /*reaction 22: CH + O <=> HCOp + E */
    qf[21] = sc[5]*sc[14];
    qr[21] = sc[6]*sc[18];

    /*reaction 23: HCOp + H2O <=> CO + H3Op */
    qf[22] = sc[4]*sc[6];
    qr[22] = sc[12]*sc[15];

    /*reaction 24: H3Op + E <=> H2O + H */
    qf[23] = sc[15]*sc[18];
    qr[23] = sc[2]*sc[4];

    /*reaction 25: CH + O2 <=> CO + OH */
    qf[24] = sc[10]*sc[14];
    qr[24] = sc[3]*sc[12];

    /*reaction 26: O + N2 => NO + N */
    qf[25] = sc[5]*sc[16];
    qr[25] = 0.0;

    /*reaction 27: NO + N => O + N2 */
    qf[26] = sc[7]*sc[8];
    qr[26] = 0.0;

    /*reaction 28: N + O2 => NO + O */
    qf[27] = sc[8]*sc[10];
    qr[27] = 0.0;

    /*reaction 29: NO + O => N + O2 */
    qf[28] = sc[5]*sc[7];
    qr[28] = 0.0;

    /*reaction 30: NO2 + O <=> NO + O2 */
    qf[29] = sc[5]*sc[17];
    qr[29] = sc[7]*sc[10];

    /*reaction 31: NO2 + H <=> NO + OH */
    qf[30] = sc[2]*sc[17];
    qr[30] = sc[3]*sc[7];

    double T = tc[1];

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int i = 0; i < 20; ++i) {
        mixture += sc[i];
    }

    double Corr[31];
    for (int i = 0; i < 31; ++i) {
        Corr[i] = 1.0;
    }

    /* simple three-body correction */
    {
        double alpha;
        alpha = mixture + (TB[0][0] - 1)*sc[9] + (TB[0][1] - 1)*sc[4] + (TB[0][2] - 1)*sc[1];
        Corr[0] = alpha;
        alpha = mixture + (TB[1][0] - 1)*sc[9] + (TB[1][1] - 1)*sc[4] + (TB[1][2] - 1)*sc[1] + (TB[1][3] - 1)*sc[12] + (TB[1][4] - 1)*sc[13];
        Corr[1] = alpha;
        alpha = mixture + (TB[2][0] - 1)*sc[9] + (TB[2][1] - 1)*sc[4] + (TB[2][2] - 1)*sc[1] + (TB[2][3] - 1)*sc[13];
        Corr[2] = alpha;
        alpha = mixture + (TB[3][0] - 1)*sc[9] + (TB[3][1] - 1)*sc[4] + (TB[3][2] - 1)*sc[1] + (TB[3][3] - 1)*sc[12] + (TB[3][4] - 1)*sc[13];
        Corr[3] = alpha;
    }

    for (int i=0; i<31; i++)
    {
        qf[i] *= Corr[i] * k_f_save[i];
        qr[i] *= Corr[i] * k_f_save[i] / Kc_save[i];
    }

    return;
}


/*compute the production rate for each species */
void vproductionRate(int npt, double * restrict wdot, double * restrict sc, double * restrict T)
{
    double k_f_s[31*npt], Kc_s[31*npt], mixture[npt], g_RT[20*npt];
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

    for (int n=0; n<20; n++) {
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
        k_f_s[21*npt+i] = prefactor_units[21] * fwd_A[21] * exp(fwd_beta[21] * tc[i] - activation_units[21] * fwd_Ea[21] * invT[i]);
        k_f_s[22*npt+i] = prefactor_units[22] * fwd_A[22] * exp(fwd_beta[22] * tc[i] - activation_units[22] * fwd_Ea[22] * invT[i]);
        k_f_s[23*npt+i] = prefactor_units[23] * fwd_A[23] * exp(fwd_beta[23] * tc[i] - activation_units[23] * fwd_Ea[23] * invT[i]);
        k_f_s[24*npt+i] = prefactor_units[24] * fwd_A[24] * exp(fwd_beta[24] * tc[i] - activation_units[24] * fwd_Ea[24] * invT[i]);
        k_f_s[25*npt+i] = prefactor_units[25] * fwd_A[25] * exp(fwd_beta[25] * tc[i] - activation_units[25] * fwd_Ea[25] * invT[i]);
        k_f_s[26*npt+i] = prefactor_units[26] * fwd_A[26] * exp(fwd_beta[26] * tc[i] - activation_units[26] * fwd_Ea[26] * invT[i]);
        k_f_s[27*npt+i] = prefactor_units[27] * fwd_A[27] * exp(fwd_beta[27] * tc[i] - activation_units[27] * fwd_Ea[27] * invT[i]);
        k_f_s[28*npt+i] = prefactor_units[28] * fwd_A[28] * exp(fwd_beta[28] * tc[i] - activation_units[28] * fwd_Ea[28] * invT[i]);
        k_f_s[29*npt+i] = prefactor_units[29] * fwd_A[29] * exp(fwd_beta[29] * tc[i] - activation_units[29] * fwd_Ea[29] * invT[i]);
        k_f_s[30*npt+i] = prefactor_units[30] * fwd_A[30] * exp(fwd_beta[30] * tc[i] - activation_units[30] * fwd_Ea[30] * invT[i]);
    }
}

void vcomp_gibbs(int npt, double * restrict g_RT, double * restrict tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        double tg[5], g[20];
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
        g_RT[9*npt+i] = g[9];
        g_RT[10*npt+i] = g[10];
        g_RT[11*npt+i] = g[11];
        g_RT[12*npt+i] = g[12];
        g_RT[13*npt+i] = g[13];
        g_RT[14*npt+i] = g[14];
        g_RT[15*npt+i] = g[15];
        g_RT[16*npt+i] = g[16];
        g_RT[17*npt+i] = g[17];
        g_RT[18*npt+i] = g[18];
        g_RT[19*npt+i] = g[19];
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

        Kc_s[0*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[3*npt+i]) - (g_RT[4*npt+i]));
        Kc_s[1*npt+i] = refCinv * exp((g_RT[5*npt+i] + g_RT[5*npt+i]) - (g_RT[10*npt+i]));
        Kc_s[2*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[2*npt+i]) - (g_RT[9*npt+i]));
        Kc_s[3*npt+i] = refCinv * exp((g_RT[5*npt+i] + g_RT[7*npt+i]) - (g_RT[17*npt+i]));
        Kc_s[4*npt+i] = refC * exp((g_RT[1*npt+i]) - (g_RT[0*npt+i] + g_RT[2*npt+i]));
        Kc_s[5*npt+i] = exp((g_RT[1*npt+i] + g_RT[3*npt+i]) - (g_RT[0*npt+i] + g_RT[4*npt+i]));
        Kc_s[6*npt+i] = exp((g_RT[1*npt+i] + g_RT[5*npt+i]) - (g_RT[0*npt+i] + g_RT[3*npt+i]));
        Kc_s[7*npt+i] = exp((g_RT[1*npt+i] + g_RT[2*npt+i]) - (g_RT[0*npt+i] + g_RT[9*npt+i]));
        Kc_s[8*npt+i] = exp((g_RT[0*npt+i] + g_RT[10*npt+i]) - (g_RT[3*npt+i] + g_RT[11*npt+i]));
        Kc_s[9*npt+i] = refC * exp((g_RT[3*npt+i] + g_RT[11*npt+i]) - (g_RT[2*npt+i] + g_RT[4*npt+i] + g_RT[12*npt+i]));
        Kc_s[10*npt+i] = exp((g_RT[3*npt+i] + g_RT[12*npt+i]) - (g_RT[2*npt+i] + g_RT[13*npt+i]));
        Kc_s[11*npt+i] = exp((g_RT[2*npt+i] + g_RT[13*npt+i]) - (g_RT[3*npt+i] + g_RT[12*npt+i]));
        Kc_s[12*npt+i] = exp((g_RT[2*npt+i] + g_RT[10*npt+i]) - (g_RT[3*npt+i] + g_RT[5*npt+i]));
        Kc_s[13*npt+i] = exp((g_RT[3*npt+i] + g_RT[5*npt+i]) - (g_RT[2*npt+i] + g_RT[10*npt+i]));
        Kc_s[14*npt+i] = exp((g_RT[5*npt+i] + g_RT[9*npt+i]) - (g_RT[2*npt+i] + g_RT[3*npt+i]));
        Kc_s[15*npt+i] = exp((g_RT[2*npt+i] + g_RT[3*npt+i]) - (g_RT[5*npt+i] + g_RT[9*npt+i]));
        Kc_s[16*npt+i] = exp((g_RT[4*npt+i] + g_RT[5*npt+i]) - (2 * g_RT[3*npt+i]));
        Kc_s[17*npt+i] = exp((2 * g_RT[3*npt+i]) - (g_RT[4*npt+i] + g_RT[5*npt+i]));
        Kc_s[18*npt+i] = exp((g_RT[3*npt+i] + g_RT[9*npt+i]) - (g_RT[2*npt+i] + g_RT[4*npt+i]));
        Kc_s[19*npt+i] = exp((g_RT[2*npt+i] + g_RT[4*npt+i]) - (g_RT[3*npt+i] + g_RT[9*npt+i]));
        Kc_s[20*npt+i] = exp((g_RT[0*npt+i] + g_RT[5*npt+i]) - (g_RT[4*npt+i] + g_RT[14*npt+i]));
        Kc_s[21*npt+i] = exp((g_RT[5*npt+i] + g_RT[14*npt+i]) - (g_RT[6*npt+i] + g_RT[18*npt+i]));
        Kc_s[22*npt+i] = exp((g_RT[4*npt+i] + g_RT[6*npt+i]) - (g_RT[12*npt+i] + g_RT[15*npt+i]));
        Kc_s[23*npt+i] = exp((g_RT[15*npt+i] + g_RT[18*npt+i]) - (g_RT[2*npt+i] + g_RT[4*npt+i]));
        Kc_s[24*npt+i] = exp((g_RT[10*npt+i] + g_RT[14*npt+i]) - (g_RT[3*npt+i] + g_RT[12*npt+i]));
        Kc_s[25*npt+i] = exp((g_RT[5*npt+i] + g_RT[16*npt+i]) - (g_RT[7*npt+i] + g_RT[8*npt+i]));
        Kc_s[26*npt+i] = exp((g_RT[7*npt+i] + g_RT[8*npt+i]) - (g_RT[5*npt+i] + g_RT[16*npt+i]));
        Kc_s[27*npt+i] = exp((g_RT[8*npt+i] + g_RT[10*npt+i]) - (g_RT[5*npt+i] + g_RT[7*npt+i]));
        Kc_s[28*npt+i] = exp((g_RT[5*npt+i] + g_RT[7*npt+i]) - (g_RT[8*npt+i] + g_RT[10*npt+i]));
        Kc_s[29*npt+i] = exp((g_RT[5*npt+i] + g_RT[17*npt+i]) - (g_RT[7*npt+i] + g_RT[10*npt+i]));
        Kc_s[30*npt+i] = exp((g_RT[2*npt+i] + g_RT[17*npt+i]) - (g_RT[3*npt+i] + g_RT[7*npt+i]));
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

        /*reaction 1: H + OH + M <=> H2O + M */
        phi_f = sc[2*npt+i]*sc[3*npt+i];
        alpha = mixture[i] + (TB[0][0] - 1)*sc[9*npt+i] + (TB[0][1] - 1)*sc[4*npt+i] + (TB[0][2] - 1)*sc[1*npt+i];
        k_f = alpha * k_f_s[0*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i];
        Kc = Kc_s[0*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 2: O + O + M <=> O2 + M */
        phi_f = sc[5*npt+i]*sc[5*npt+i];
        alpha = mixture[i] + (TB[1][0] - 1)*sc[9*npt+i] + (TB[1][1] - 1)*sc[4*npt+i] + (TB[1][2] - 1)*sc[1*npt+i] + (TB[1][3] - 1)*sc[12*npt+i] + (TB[1][4] - 1)*sc[13*npt+i];
        k_f = alpha * k_f_s[1*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[10*npt+i];
        Kc = Kc_s[1*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 3: H + H + M <=> H2 + M */
        phi_f = sc[2*npt+i]*sc[2*npt+i];
        alpha = mixture[i] + (TB[2][0] - 1)*sc[9*npt+i] + (TB[2][1] - 1)*sc[4*npt+i] + (TB[2][2] - 1)*sc[1*npt+i] + (TB[2][3] - 1)*sc[13*npt+i];
        k_f = alpha * k_f_s[2*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i];
        Kc = Kc_s[2*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 4: NO + O + M <=> NO2 + M */
        phi_f = sc[5*npt+i]*sc[7*npt+i];
        alpha = mixture[i] + (TB[3][0] - 1)*sc[9*npt+i] + (TB[3][1] - 1)*sc[4*npt+i] + (TB[3][2] - 1)*sc[1*npt+i] + (TB[3][3] - 1)*sc[12*npt+i] + (TB[3][4] - 1)*sc[13*npt+i];
        k_f = alpha * k_f_s[3*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[17*npt+i];
        Kc = Kc_s[3*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[7*npt+i] -= qdot;
        wdot[17*npt+i] += qdot;

        /*reaction 5: CH4 <=> CH3 + H */
        phi_f = sc[1*npt+i];
        k_f = k_f_s[4*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[2*npt+i];
        Kc = Kc_s[4*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;

        /*reaction 6: CH4 + OH <=> CH3 + H2O */
        phi_f = sc[1*npt+i]*sc[3*npt+i];
        k_f = k_f_s[5*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[4*npt+i];
        Kc = Kc_s[5*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 7: CH4 + O <=> CH3 + OH */
        phi_f = sc[1*npt+i]*sc[5*npt+i];
        k_f = k_f_s[6*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[3*npt+i];
        Kc = Kc_s[6*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;

        /*reaction 8: CH4 + H <=> CH3 + H2 */
        phi_f = sc[1*npt+i]*sc[2*npt+i];
        k_f = k_f_s[7*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[9*npt+i];
        Kc = Kc_s[7*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 9: CH3 + O2 <=> CH2O + OH */
        phi_f = sc[0*npt+i]*sc[10*npt+i];
        k_f = k_f_s[8*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[11*npt+i];
        Kc = Kc_s[8*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 10: CH2O + OH <=> CO + H2O + H */
        phi_f = sc[3*npt+i]*sc[11*npt+i];
        k_f = k_f_s[9*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[4*npt+i]*sc[12*npt+i];
        Kc = Kc_s[9*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 11: CO + OH => CO2 + H */
        phi_f = sc[3*npt+i]*sc[12*npt+i];
        k_f = k_f_s[10*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 12: CO2 + H => CO + OH */
        phi_f = sc[2*npt+i]*sc[13*npt+i];
        k_f = k_f_s[11*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 13: O2 + H => OH + O */
        phi_f = sc[2*npt+i]*sc[10*npt+i];
        k_f = k_f_s[12*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;

        /*reaction 14: OH + O => O2 + H */
        phi_f = sc[3*npt+i]*sc[5*npt+i];
        k_f = k_f_s[13*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 15: O + H2 => OH + H */
        phi_f = sc[5*npt+i]*sc[9*npt+i];
        k_f = k_f_s[14*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;

        /*reaction 16: OH + H => O + H2 */
        phi_f = sc[2*npt+i]*sc[3*npt+i];
        k_f = k_f_s[15*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 17: O + H2O => 2 OH */
        phi_f = sc[4*npt+i]*sc[5*npt+i];
        k_f = k_f_s[16*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] += 2 * qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;

        /*reaction 18: 2 OH => O + H2O */
        phi_f = sc[3*npt+i]*sc[3*npt+i];
        k_f = k_f_s[17*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= 2 * qdot;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 19: OH + H2 => H2O + H */
        phi_f = sc[3*npt+i]*sc[9*npt+i];
        k_f = k_f_s[18*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;

        /*reaction 20: H2O + H => H2 + OH */
        phi_f = sc[2*npt+i]*sc[4*npt+i];
        k_f = k_f_s[19*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 21: CH3 + O <=> CH + H2O */
        phi_f = sc[0*npt+i]*sc[5*npt+i];
        k_f = k_f_s[20*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[14*npt+i];
        Kc = Kc_s[20*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 22: CH + O <=> HCOp + E */
        phi_f = sc[5*npt+i]*sc[14*npt+i];
        k_f = k_f_s[21*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[18*npt+i];
        Kc = Kc_s[21*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;
        wdot[18*npt+i] += qdot;

        /*reaction 23: HCOp + H2O <=> CO + H3Op */
        phi_f = sc[4*npt+i]*sc[6*npt+i];
        k_f = k_f_s[22*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[12*npt+i]*sc[15*npt+i];
        Kc = Kc_s[22*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 24: H3Op + E <=> H2O + H */
        phi_f = sc[15*npt+i]*sc[18*npt+i];
        k_f = k_f_s[23*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[4*npt+i];
        Kc = Kc_s[23*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[4*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 25: CH + O2 <=> CO + OH */
        phi_f = sc[10*npt+i]*sc[14*npt+i];
        k_f = k_f_s[24*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[12*npt+i];
        Kc = Kc_s[24*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 26: O + N2 => NO + N */
        phi_f = sc[5*npt+i]*sc[16*npt+i];
        k_f = k_f_s[25*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 27: NO + N => O + N2 */
        phi_f = sc[7*npt+i]*sc[8*npt+i];
        k_f = k_f_s[26*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 28: N + O2 => NO + O */
        phi_f = sc[8*npt+i]*sc[10*npt+i];
        k_f = k_f_s[27*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;

        /*reaction 29: NO + O => N + O2 */
        phi_f = sc[5*npt+i]*sc[7*npt+i];
        k_f = k_f_s[28*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 30: NO2 + O <=> NO + O2 */
        phi_f = sc[5*npt+i]*sc[17*npt+i];
        k_f = k_f_s[29*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[10*npt+i];
        Kc = Kc_s[29*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 31: NO2 + H <=> NO + OH */
        phi_f = sc[2*npt+i]*sc[17*npt+i];
        k_f = k_f_s[30*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[7*npt+i];
        Kc = Kc_s[30*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[7*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;
    }
}

/*compute the reaction Jacobian */
void DWDOT(double * restrict J, double * restrict sc, double * restrict Tp, int * consP)
{
    double c[20];

    for (int k=0; k<20; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian(J, c, *Tp, *consP);

    /* dwdot[k]/dT */
    for (int k=0; k<20; k++) {
        J[420+k] *= 1.e-6;
    }

    /* dTdot/d[X] */
    for (int k=0; k<20; k++) {
        J[k*21+20] *= 1.e6;
    }

    return;
}

/*compute the reaction Jacobian */
void aJacobian(double * restrict J, double * restrict sc, double T, int consP)
{
    for (int i=0; i<441; i++) {
        J[i] = 0.0;
    }

    double wdot[20];
    for (int k=0; k<20; k++) {
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
    for (int k = 0; k < 20; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[20];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[20];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[20];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    /*reaction 1: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[0][0] - 1)*sc[9] + (TB[0][1] - 1)*sc[4] + (TB[0][2] - 1)*sc[1];
    /* forward */
    phi_f = sc[2]*sc[3];
    k_f = prefactor_units[0] * fwd_A[0]
                * exp(fwd_beta[0] * tc[0] - activation_units[0] * fwd_Ea[0] * invT);
    dlnkfdT = fwd_beta[0] * invT + activation_units[0] * fwd_Ea[0] * invT2;
    /* reverse */
    phi_r = sc[4];
    Kc = refCinv * exp(g_RT[2] + g_RT[3] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[3]) + (h_RT[4]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[CH4] */
        dqdci = (TB[0][2] - 1)*q_nocor;
        J[23] -= dqdci;               /* dwdot[H]/d[CH4] */
        J[24] -= dqdci;               /* dwdot[OH]/d[CH4] */
        J[25] += dqdci;               /* dwdot[H2O]/d[CH4] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[3];
        J[44] -= dqdci;               /* dwdot[H]/d[H] */
        J[45] -= dqdci;               /* dwdot[OH]/d[H] */
        J[46] += dqdci;               /* dwdot[H2O]/d[H] */
        /* d()/d[OH] */
        dqdci =  + k_f*sc[2];
        J[65] -= dqdci;               /* dwdot[H]/d[OH] */
        J[66] -= dqdci;               /* dwdot[OH]/d[OH] */
        J[67] += dqdci;               /* dwdot[H2O]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[0][1] - 1)*q_nocor - k_r;
        J[86] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[87] -= dqdci;               /* dwdot[OH]/d[H2O] */
        J[88] += dqdci;               /* dwdot[H2O]/d[H2O] */
        /* d()/d[H2] */
        dqdci = (TB[0][0] - 1)*q_nocor;
        J[191] -= dqdci;              /* dwdot[H]/d[H2] */
        J[192] -= dqdci;              /* dwdot[OH]/d[H2] */
        J[193] += dqdci;              /* dwdot[H2O]/d[H2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = TB[0][2];
        dqdc[2] = dcdc_fac + k_f*sc[3];
        dqdc[3] = dcdc_fac + k_f*sc[2];
        dqdc[4] = TB[0][1] - k_r;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[0][0];
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        for (int k=0; k<20; k++) {
            J[21*k+2] -= dqdc[k];
            J[21*k+3] -= dqdc[k];
            J[21*k+4] += dqdc[k];
        }
    }
    J[422] -= dqdT; /* dwdot[H]/dT */
    J[423] -= dqdT; /* dwdot[OH]/dT */
    J[424] += dqdT; /* dwdot[H2O]/dT */

    /*reaction 2: O + O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[1][0] - 1)*sc[9] + (TB[1][1] - 1)*sc[4] + (TB[1][2] - 1)*sc[1] + (TB[1][3] - 1)*sc[12] + (TB[1][4] - 1)*sc[13];
    /* forward */
    phi_f = sc[5]*sc[5];
    k_f = prefactor_units[1] * fwd_A[1]
                * exp(fwd_beta[1] * tc[0] - activation_units[1] * fwd_Ea[1] * invT);
    dlnkfdT = fwd_beta[1] * invT + activation_units[1] * fwd_Ea[1] * invT2;
    /* reverse */
    phi_r = sc[10];
    Kc = refCinv * exp(g_RT[5] + g_RT[5] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[5]) + (h_RT[10]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= 2 * q; /* O */
    wdot[10] += q; /* O2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[CH4] */
        dqdci = (TB[1][2] - 1)*q_nocor;
        J[26] += -2 * dqdci;          /* dwdot[O]/d[CH4] */
        J[31] += dqdci;               /* dwdot[O2]/d[CH4] */
        /* d()/d[H2O] */
        dqdci = (TB[1][1] - 1)*q_nocor;
        J[89] += -2 * dqdci;          /* dwdot[O]/d[H2O] */
        J[94] += dqdci;               /* dwdot[O2]/d[H2O] */
        /* d()/d[O] */
        dqdci =  + k_f*2*sc[5];
        J[110] += -2 * dqdci;         /* dwdot[O]/d[O] */
        J[115] += dqdci;              /* dwdot[O2]/d[O] */
        /* d()/d[H2] */
        dqdci = (TB[1][0] - 1)*q_nocor;
        J[194] += -2 * dqdci;         /* dwdot[O]/d[H2] */
        J[199] += dqdci;              /* dwdot[O2]/d[H2] */
        /* d()/d[O2] */
        dqdci =  - k_r;
        J[215] += -2 * dqdci;         /* dwdot[O]/d[O2] */
        J[220] += dqdci;              /* dwdot[O2]/d[O2] */
        /* d()/d[CO] */
        dqdci = (TB[1][3] - 1)*q_nocor;
        J[257] += -2 * dqdci;         /* dwdot[O]/d[CO] */
        J[262] += dqdci;              /* dwdot[O2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[1][4] - 1)*q_nocor;
        J[278] += -2 * dqdci;         /* dwdot[O]/d[CO2] */
        J[283] += dqdci;              /* dwdot[O2]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = TB[1][2];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = TB[1][1];
        dqdc[5] = dcdc_fac + k_f*2*sc[5];
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[1][0];
        dqdc[10] = dcdc_fac - k_r;
        dqdc[11] = dcdc_fac;
        dqdc[12] = TB[1][3];
        dqdc[13] = TB[1][4];
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        for (int k=0; k<20; k++) {
            J[21*k+5] += -2 * dqdc[k];
            J[21*k+10] += dqdc[k];
        }
    }
    J[425] += -2 * dqdT; /* dwdot[O]/dT */
    J[430] += dqdT; /* dwdot[O2]/dT */

    /*reaction 3: H + H + M <=> H2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[2][0] - 1)*sc[9] + (TB[2][1] - 1)*sc[4] + (TB[2][2] - 1)*sc[1] + (TB[2][3] - 1)*sc[13];
    /* forward */
    phi_f = sc[2]*sc[2];
    k_f = prefactor_units[2] * fwd_A[2]
                * exp(fwd_beta[2] * tc[0] - activation_units[2] * fwd_Ea[2] * invT);
    dlnkfdT = fwd_beta[2] * invT + activation_units[2] * fwd_Ea[2] * invT2;
    /* reverse */
    phi_r = sc[9];
    Kc = refCinv * exp(g_RT[2] + g_RT[2] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[2]) + (h_RT[9]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= 2 * q; /* H */
    wdot[9] += q; /* H2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[CH4] */
        dqdci = (TB[2][2] - 1)*q_nocor;
        J[23] += -2 * dqdci;          /* dwdot[H]/d[CH4] */
        J[30] += dqdci;               /* dwdot[H2]/d[CH4] */
        /* d()/d[H] */
        dqdci =  + k_f*2*sc[2];
        J[44] += -2 * dqdci;          /* dwdot[H]/d[H] */
        J[51] += dqdci;               /* dwdot[H2]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[2][1] - 1)*q_nocor;
        J[86] += -2 * dqdci;          /* dwdot[H]/d[H2O] */
        J[93] += dqdci;               /* dwdot[H2]/d[H2O] */
        /* d()/d[H2] */
        dqdci = (TB[2][0] - 1)*q_nocor - k_r;
        J[191] += -2 * dqdci;         /* dwdot[H]/d[H2] */
        J[198] += dqdci;              /* dwdot[H2]/d[H2] */
        /* d()/d[CO2] */
        dqdci = (TB[2][3] - 1)*q_nocor;
        J[275] += -2 * dqdci;         /* dwdot[H]/d[CO2] */
        J[282] += dqdci;              /* dwdot[H2]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = TB[2][2];
        dqdc[2] = dcdc_fac + k_f*2*sc[2];
        dqdc[3] = dcdc_fac;
        dqdc[4] = TB[2][1];
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[2][0] - k_r;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = TB[2][3];
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        for (int k=0; k<20; k++) {
            J[21*k+2] += -2 * dqdc[k];
            J[21*k+9] += dqdc[k];
        }
    }
    J[422] += -2 * dqdT; /* dwdot[H]/dT */
    J[429] += dqdT; /* dwdot[H2]/dT */

    /*reaction 4: NO + O + M <=> NO2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[3][0] - 1)*sc[9] + (TB[3][1] - 1)*sc[4] + (TB[3][2] - 1)*sc[1] + (TB[3][3] - 1)*sc[12] + (TB[3][4] - 1)*sc[13];
    /* forward */
    phi_f = sc[5]*sc[7];
    k_f = prefactor_units[3] * fwd_A[3]
                * exp(fwd_beta[3] * tc[0] - activation_units[3] * fwd_Ea[3] * invT);
    dlnkfdT = fwd_beta[3] * invT + activation_units[3] * fwd_Ea[3] * invT2;
    /* reverse */
    phi_r = sc[17];
    Kc = refCinv * exp(g_RT[5] + g_RT[7] - g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[17]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[7] -= q; /* NO */
    wdot[17] += q; /* NO2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[CH4] */
        dqdci = (TB[3][2] - 1)*q_nocor;
        J[26] -= dqdci;               /* dwdot[O]/d[CH4] */
        J[28] -= dqdci;               /* dwdot[NO]/d[CH4] */
        J[38] += dqdci;               /* dwdot[NO2]/d[CH4] */
        /* d()/d[H2O] */
        dqdci = (TB[3][1] - 1)*q_nocor;
        J[89] -= dqdci;               /* dwdot[O]/d[H2O] */
        J[91] -= dqdci;               /* dwdot[NO]/d[H2O] */
        J[101] += dqdci;              /* dwdot[NO2]/d[H2O] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[7];
        J[110] -= dqdci;              /* dwdot[O]/d[O] */
        J[112] -= dqdci;              /* dwdot[NO]/d[O] */
        J[122] += dqdci;              /* dwdot[NO2]/d[O] */
        /* d()/d[NO] */
        dqdci =  + k_f*sc[5];
        J[152] -= dqdci;              /* dwdot[O]/d[NO] */
        J[154] -= dqdci;              /* dwdot[NO]/d[NO] */
        J[164] += dqdci;              /* dwdot[NO2]/d[NO] */
        /* d()/d[H2] */
        dqdci = (TB[3][0] - 1)*q_nocor;
        J[194] -= dqdci;              /* dwdot[O]/d[H2] */
        J[196] -= dqdci;              /* dwdot[NO]/d[H2] */
        J[206] += dqdci;              /* dwdot[NO2]/d[H2] */
        /* d()/d[CO] */
        dqdci = (TB[3][3] - 1)*q_nocor;
        J[257] -= dqdci;              /* dwdot[O]/d[CO] */
        J[259] -= dqdci;              /* dwdot[NO]/d[CO] */
        J[269] += dqdci;              /* dwdot[NO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[3][4] - 1)*q_nocor;
        J[278] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[280] -= dqdci;              /* dwdot[NO]/d[CO2] */
        J[290] += dqdci;              /* dwdot[NO2]/d[CO2] */
        /* d()/d[NO2] */
        dqdci =  - k_r;
        J[362] -= dqdci;              /* dwdot[O]/d[NO2] */
        J[364] -= dqdci;              /* dwdot[NO]/d[NO2] */
        J[374] += dqdci;              /* dwdot[NO2]/d[NO2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = TB[3][2];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = TB[3][1];
        dqdc[5] = dcdc_fac + k_f*sc[7];
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac + k_f*sc[5];
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[3][0];
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = TB[3][3];
        dqdc[13] = TB[3][4];
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac - k_r;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        for (int k=0; k<20; k++) {
            J[21*k+5] -= dqdc[k];
            J[21*k+7] -= dqdc[k];
            J[21*k+17] += dqdc[k];
        }
    }
    J[425] -= dqdT; /* dwdot[O]/dT */
    J[427] -= dqdT; /* dwdot[NO]/dT */
    J[437] += dqdT; /* dwdot[NO2]/dT */

    /*reaction 5: CH4 <=> CH3 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1];
    k_f = prefactor_units[4] * fwd_A[4]
                * exp(fwd_beta[4] * tc[0] - activation_units[4] * fwd_Ea[4] * invT);
    dlnkfdT = fwd_beta[4] * invT + activation_units[4] * fwd_Ea[4] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[2];
    Kc = refC * exp(-g_RT[0] + g_RT[1] - g_RT[2]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1]) + (h_RT[0] + h_RT[2]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* CH3 */
    wdot[1] -= q; /* CH4 */
    wdot[2] += q; /* H */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[2];
    J[0] += dqdci;                /* dwdot[CH3]/d[CH3] */
    J[1] -= dqdci;                /* dwdot[CH4]/d[CH3] */
    J[2] += dqdci;                /* dwdot[H]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f;
    J[21] += dqdci;               /* dwdot[CH3]/d[CH4] */
    J[22] -= dqdci;               /* dwdot[CH4]/d[CH4] */
    J[23] += dqdci;               /* dwdot[H]/d[CH4] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[0];
    J[42] += dqdci;               /* dwdot[CH3]/d[H] */
    J[43] -= dqdci;               /* dwdot[CH4]/d[H] */
    J[44] += dqdci;               /* dwdot[H]/d[H] */
    /* d()/dT */
    J[420] += dqdT;               /* dwdot[CH3]/dT */
    J[421] -= dqdT;               /* dwdot[CH4]/dT */
    J[422] += dqdT;               /* dwdot[H]/dT */

    /*reaction 6: CH4 + OH <=> CH3 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = prefactor_units[5] * fwd_A[5]
                * exp(fwd_beta[5] * tc[0] - activation_units[5] * fwd_Ea[5] * invT);
    dlnkfdT = fwd_beta[5] * invT + activation_units[5] * fwd_Ea[5] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[4];
    Kc = exp(-g_RT[0] + g_RT[1] + g_RT[3] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[0] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* CH3 */
    wdot[1] -= q; /* CH4 */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[4];
    J[0] += dqdci;                /* dwdot[CH3]/d[CH3] */
    J[1] -= dqdci;                /* dwdot[CH4]/d[CH3] */
    J[3] -= dqdci;                /* dwdot[OH]/d[CH3] */
    J[4] += dqdci;                /* dwdot[H2O]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[3];
    J[21] += dqdci;               /* dwdot[CH3]/d[CH4] */
    J[22] -= dqdci;               /* dwdot[CH4]/d[CH4] */
    J[24] -= dqdci;               /* dwdot[OH]/d[CH4] */
    J[25] += dqdci;               /* dwdot[H2O]/d[CH4] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[1];
    J[63] += dqdci;               /* dwdot[CH3]/d[OH] */
    J[64] -= dqdci;               /* dwdot[CH4]/d[OH] */
    J[66] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[67] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[0];
    J[84] += dqdci;               /* dwdot[CH3]/d[H2O] */
    J[85] -= dqdci;               /* dwdot[CH4]/d[H2O] */
    J[87] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[88] += dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[420] += dqdT;               /* dwdot[CH3]/dT */
    J[421] -= dqdT;               /* dwdot[CH4]/dT */
    J[423] -= dqdT;               /* dwdot[OH]/dT */
    J[424] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 7: CH4 + O <=> CH3 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[5];
    k_f = prefactor_units[6] * fwd_A[6]
                * exp(fwd_beta[6] * tc[0] - activation_units[6] * fwd_Ea[6] * invT);
    dlnkfdT = fwd_beta[6] * invT + activation_units[6] * fwd_Ea[6] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[3] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[5]) + (h_RT[0] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* CH3 */
    wdot[1] -= q; /* CH4 */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[3];
    J[0] += dqdci;                /* dwdot[CH3]/d[CH3] */
    J[1] -= dqdci;                /* dwdot[CH4]/d[CH3] */
    J[3] += dqdci;                /* dwdot[OH]/d[CH3] */
    J[5] -= dqdci;                /* dwdot[O]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[5];
    J[21] += dqdci;               /* dwdot[CH3]/d[CH4] */
    J[22] -= dqdci;               /* dwdot[CH4]/d[CH4] */
    J[24] += dqdci;               /* dwdot[OH]/d[CH4] */
    J[26] -= dqdci;               /* dwdot[O]/d[CH4] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[0];
    J[63] += dqdci;               /* dwdot[CH3]/d[OH] */
    J[64] -= dqdci;               /* dwdot[CH4]/d[OH] */
    J[66] += dqdci;               /* dwdot[OH]/d[OH] */
    J[68] -= dqdci;               /* dwdot[O]/d[OH] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[1];
    J[105] += dqdci;              /* dwdot[CH3]/d[O] */
    J[106] -= dqdci;              /* dwdot[CH4]/d[O] */
    J[108] += dqdci;              /* dwdot[OH]/d[O] */
    J[110] -= dqdci;              /* dwdot[O]/d[O] */
    /* d()/dT */
    J[420] += dqdT;               /* dwdot[CH3]/dT */
    J[421] -= dqdT;               /* dwdot[CH4]/dT */
    J[423] += dqdT;               /* dwdot[OH]/dT */
    J[425] -= dqdT;               /* dwdot[O]/dT */

    /*reaction 8: CH4 + H <=> CH3 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = prefactor_units[7] * fwd_A[7]
                * exp(fwd_beta[7] * tc[0] - activation_units[7] * fwd_Ea[7] * invT);
    dlnkfdT = fwd_beta[7] * invT + activation_units[7] * fwd_Ea[7] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[9];
    Kc = exp(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[0] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* CH3 */
    wdot[1] -= q; /* CH4 */
    wdot[2] -= q; /* H */
    wdot[9] += q; /* H2 */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[9];
    J[0] += dqdci;                /* dwdot[CH3]/d[CH3] */
    J[1] -= dqdci;                /* dwdot[CH4]/d[CH3] */
    J[2] -= dqdci;                /* dwdot[H]/d[CH3] */
    J[9] += dqdci;                /* dwdot[H2]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[2];
    J[21] += dqdci;               /* dwdot[CH3]/d[CH4] */
    J[22] -= dqdci;               /* dwdot[CH4]/d[CH4] */
    J[23] -= dqdci;               /* dwdot[H]/d[CH4] */
    J[30] += dqdci;               /* dwdot[H2]/d[CH4] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[1];
    J[42] += dqdci;               /* dwdot[CH3]/d[H] */
    J[43] -= dqdci;               /* dwdot[CH4]/d[H] */
    J[44] -= dqdci;               /* dwdot[H]/d[H] */
    J[51] += dqdci;               /* dwdot[H2]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[0];
    J[189] += dqdci;              /* dwdot[CH3]/d[H2] */
    J[190] -= dqdci;              /* dwdot[CH4]/d[H2] */
    J[191] -= dqdci;              /* dwdot[H]/d[H2] */
    J[198] += dqdci;              /* dwdot[H2]/d[H2] */
    /* d()/dT */
    J[420] += dqdT;               /* dwdot[CH3]/dT */
    J[421] -= dqdT;               /* dwdot[CH4]/dT */
    J[422] -= dqdT;               /* dwdot[H]/dT */
    J[429] += dqdT;               /* dwdot[H2]/dT */

    /*reaction 9: CH3 + O2 <=> CH2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[10];
    k_f = prefactor_units[8] * fwd_A[8]
                * exp(fwd_beta[8] * tc[0] - activation_units[8] * fwd_Ea[8] * invT);
    dlnkfdT = fwd_beta[8] * invT + activation_units[8] * fwd_Ea[8] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[11];
    Kc = exp(g_RT[0] - g_RT[3] + g_RT[10] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[10]) + (h_RT[3] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* CH3 */
    wdot[3] += q; /* OH */
    wdot[10] -= q; /* O2 */
    wdot[11] += q; /* CH2O */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[10];
    J[0] -= dqdci;                /* dwdot[CH3]/d[CH3] */
    J[3] += dqdci;                /* dwdot[OH]/d[CH3] */
    J[10] -= dqdci;               /* dwdot[O2]/d[CH3] */
    J[11] += dqdci;               /* dwdot[CH2O]/d[CH3] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[11];
    J[63] -= dqdci;               /* dwdot[CH3]/d[OH] */
    J[66] += dqdci;               /* dwdot[OH]/d[OH] */
    J[73] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[74] += dqdci;               /* dwdot[CH2O]/d[OH] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[0];
    J[210] -= dqdci;              /* dwdot[CH3]/d[O2] */
    J[213] += dqdci;              /* dwdot[OH]/d[O2] */
    J[220] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[221] += dqdci;              /* dwdot[CH2O]/d[O2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[3];
    J[231] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[234] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[241] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[242] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[420] -= dqdT;               /* dwdot[CH3]/dT */
    J[423] += dqdT;               /* dwdot[OH]/dT */
    J[430] -= dqdT;               /* dwdot[O2]/dT */
    J[431] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 10: CH2O + OH <=> CO + H2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[11];
    k_f = prefactor_units[9] * fwd_A[9]
                * exp(fwd_beta[9] * tc[0] - activation_units[9] * fwd_Ea[9] * invT);
    dlnkfdT = fwd_beta[9] * invT + activation_units[9] * fwd_Ea[9] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[4]*sc[12];
    Kc = refC * exp(-g_RT[2] + g_RT[3] - g_RT[4] + g_RT[11] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[11]) + (h_RT[2] + h_RT[4] + h_RT[12]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[11] -= q; /* CH2O */
    wdot[12] += q; /* CO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4]*sc[12];
    J[44] += dqdci;               /* dwdot[H]/d[H] */
    J[45] -= dqdci;               /* dwdot[OH]/d[H] */
    J[46] += dqdci;               /* dwdot[H2O]/d[H] */
    J[53] -= dqdci;               /* dwdot[CH2O]/d[H] */
    J[54] += dqdci;               /* dwdot[CO]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[11];
    J[65] += dqdci;               /* dwdot[H]/d[OH] */
    J[66] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[67] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[74] -= dqdci;               /* dwdot[CH2O]/d[OH] */
    J[75] += dqdci;               /* dwdot[CO]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2]*sc[12];
    J[86] += dqdci;               /* dwdot[H]/d[H2O] */
    J[87] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[88] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[95] -= dqdci;               /* dwdot[CH2O]/d[H2O] */
    J[96] += dqdci;               /* dwdot[CO]/d[H2O] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[3];
    J[233] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[234] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[235] += dqdci;              /* dwdot[H2O]/d[CH2O] */
    J[242] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[243] += dqdci;              /* dwdot[CO]/d[CH2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[2]*sc[4];
    J[254] += dqdci;              /* dwdot[H]/d[CO] */
    J[255] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[256] += dqdci;              /* dwdot[H2O]/d[CO] */
    J[263] -= dqdci;              /* dwdot[CH2O]/d[CO] */
    J[264] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[422] += dqdT;               /* dwdot[H]/dT */
    J[423] -= dqdT;               /* dwdot[OH]/dT */
    J[424] += dqdT;               /* dwdot[H2O]/dT */
    J[431] -= dqdT;               /* dwdot[CH2O]/dT */
    J[432] += dqdT;               /* dwdot[CO]/dT */

    /*reaction 11: CO + OH => CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[12];
    k_f = prefactor_units[10] * fwd_A[10]
                * exp(fwd_beta[10] * tc[0] - activation_units[10] * fwd_Ea[10] * invT);
    dlnkfdT = fwd_beta[10] * invT + activation_units[10] * fwd_Ea[10] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[3] -= q; /* OH */
    wdot[12] -= q; /* CO */
    wdot[13] += q; /* CO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[12];
    J[65] += dqdci;               /* dwdot[H]/d[OH] */
    J[66] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[75] -= dqdci;               /* dwdot[CO]/d[OH] */
    J[76] += dqdci;               /* dwdot[CO2]/d[OH] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[3];
    J[254] += dqdci;              /* dwdot[H]/d[CO] */
    J[255] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[264] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[265] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/dT */
    J[422] += dqdT;               /* dwdot[H]/dT */
    J[423] -= dqdT;               /* dwdot[OH]/dT */
    J[432] -= dqdT;               /* dwdot[CO]/dT */
    J[433] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 12: CO2 + H => CO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[13];
    k_f = prefactor_units[11] * fwd_A[11]
                * exp(fwd_beta[11] * tc[0] - activation_units[11] * fwd_Ea[11] * invT);
    dlnkfdT = fwd_beta[11] * invT + activation_units[11] * fwd_Ea[11] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[3] += q; /* OH */
    wdot[12] += q; /* CO */
    wdot[13] -= q; /* CO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[13];
    J[44] -= dqdci;               /* dwdot[H]/d[H] */
    J[45] += dqdci;               /* dwdot[OH]/d[H] */
    J[54] += dqdci;               /* dwdot[CO]/d[H] */
    J[55] -= dqdci;               /* dwdot[CO2]/d[H] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[2];
    J[275] -= dqdci;              /* dwdot[H]/d[CO2] */
    J[276] += dqdci;              /* dwdot[OH]/d[CO2] */
    J[285] += dqdci;              /* dwdot[CO]/d[CO2] */
    J[286] -= dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[422] -= dqdT;               /* dwdot[H]/dT */
    J[423] += dqdT;               /* dwdot[OH]/dT */
    J[432] += dqdT;               /* dwdot[CO]/dT */
    J[433] -= dqdT;               /* dwdot[CO2]/dT */

    /*reaction 13: O2 + H => OH + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[10];
    k_f = prefactor_units[12] * fwd_A[12]
                * exp(fwd_beta[12] * tc[0] - activation_units[12] * fwd_Ea[12] * invT);
    dlnkfdT = fwd_beta[12] * invT + activation_units[12] * fwd_Ea[12] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[3] += q; /* OH */
    wdot[5] += q; /* O */
    wdot[10] -= q; /* O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[10];
    J[44] -= dqdci;               /* dwdot[H]/d[H] */
    J[45] += dqdci;               /* dwdot[OH]/d[H] */
    J[47] += dqdci;               /* dwdot[O]/d[H] */
    J[52] -= dqdci;               /* dwdot[O2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[2];
    J[212] -= dqdci;              /* dwdot[H]/d[O2] */
    J[213] += dqdci;              /* dwdot[OH]/d[O2] */
    J[215] += dqdci;              /* dwdot[O]/d[O2] */
    J[220] -= dqdci;              /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[422] -= dqdT;               /* dwdot[H]/dT */
    J[423] += dqdT;               /* dwdot[OH]/dT */
    J[425] += dqdT;               /* dwdot[O]/dT */
    J[430] -= dqdT;               /* dwdot[O2]/dT */

    /*reaction 14: OH + O => O2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[5];
    k_f = prefactor_units[13] * fwd_A[13]
                * exp(fwd_beta[13] * tc[0] - activation_units[13] * fwd_Ea[13] * invT);
    dlnkfdT = fwd_beta[13] * invT + activation_units[13] * fwd_Ea[13] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[3] -= q; /* OH */
    wdot[5] -= q; /* O */
    wdot[10] += q; /* O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[5];
    J[65] += dqdci;               /* dwdot[H]/d[OH] */
    J[66] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[68] -= dqdci;               /* dwdot[O]/d[OH] */
    J[73] += dqdci;               /* dwdot[O2]/d[OH] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[3];
    J[107] += dqdci;              /* dwdot[H]/d[O] */
    J[108] -= dqdci;              /* dwdot[OH]/d[O] */
    J[110] -= dqdci;              /* dwdot[O]/d[O] */
    J[115] += dqdci;              /* dwdot[O2]/d[O] */
    /* d()/dT */
    J[422] += dqdT;               /* dwdot[H]/dT */
    J[423] -= dqdT;               /* dwdot[OH]/dT */
    J[425] -= dqdT;               /* dwdot[O]/dT */
    J[430] += dqdT;               /* dwdot[O2]/dT */

    /*reaction 15: O + H2 => OH + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[9];
    k_f = prefactor_units[14] * fwd_A[14]
                * exp(fwd_beta[14] * tc[0] - activation_units[14] * fwd_Ea[14] * invT);
    dlnkfdT = fwd_beta[14] * invT + activation_units[14] * fwd_Ea[14] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O */
    wdot[9] -= q; /* H2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[9];
    J[107] += dqdci;              /* dwdot[H]/d[O] */
    J[108] += dqdci;              /* dwdot[OH]/d[O] */
    J[110] -= dqdci;              /* dwdot[O]/d[O] */
    J[114] -= dqdci;              /* dwdot[H2]/d[O] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[5];
    J[191] += dqdci;              /* dwdot[H]/d[H2] */
    J[192] += dqdci;              /* dwdot[OH]/d[H2] */
    J[194] -= dqdci;              /* dwdot[O]/d[H2] */
    J[198] -= dqdci;              /* dwdot[H2]/d[H2] */
    /* d()/dT */
    J[422] += dqdT;               /* dwdot[H]/dT */
    J[423] += dqdT;               /* dwdot[OH]/dT */
    J[425] -= dqdT;               /* dwdot[O]/dT */
    J[429] -= dqdT;               /* dwdot[H2]/dT */

    /*reaction 16: OH + H => O + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[3];
    k_f = prefactor_units[15] * fwd_A[15]
                * exp(fwd_beta[15] * tc[0] - activation_units[15] * fwd_Ea[15] * invT);
    dlnkfdT = fwd_beta[15] * invT + activation_units[15] * fwd_Ea[15] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[3] -= q; /* OH */
    wdot[5] += q; /* O */
    wdot[9] += q; /* H2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3];
    J[44] -= dqdci;               /* dwdot[H]/d[H] */
    J[45] -= dqdci;               /* dwdot[OH]/d[H] */
    J[47] += dqdci;               /* dwdot[O]/d[H] */
    J[51] += dqdci;               /* dwdot[H2]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[2];
    J[65] -= dqdci;               /* dwdot[H]/d[OH] */
    J[66] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[68] += dqdci;               /* dwdot[O]/d[OH] */
    J[72] += dqdci;               /* dwdot[H2]/d[OH] */
    /* d()/dT */
    J[422] -= dqdT;               /* dwdot[H]/dT */
    J[423] -= dqdT;               /* dwdot[OH]/dT */
    J[425] += dqdT;               /* dwdot[O]/dT */
    J[429] += dqdT;               /* dwdot[H2]/dT */

    /*reaction 17: O + H2O => 2 OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[5];
    k_f = prefactor_units[16] * fwd_A[16]
                * exp(fwd_beta[16] * tc[0] - activation_units[16] * fwd_Ea[16] * invT);
    dlnkfdT = fwd_beta[16] * invT + activation_units[16] * fwd_Ea[16] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += 2 * q; /* OH */
    wdot[4] -= q; /* H2O */
    wdot[5] -= q; /* O */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[5];
    J[87] += 2 * dqdci;           /* dwdot[OH]/d[H2O] */
    J[88] -= dqdci;               /* dwdot[H2O]/d[H2O] */
    J[89] -= dqdci;               /* dwdot[O]/d[H2O] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[4];
    J[108] += 2 * dqdci;          /* dwdot[OH]/d[O] */
    J[109] -= dqdci;              /* dwdot[H2O]/d[O] */
    J[110] -= dqdci;              /* dwdot[O]/d[O] */
    /* d()/dT */
    J[423] += 2 * dqdT;           /* dwdot[OH]/dT */
    J[424] -= dqdT;               /* dwdot[H2O]/dT */
    J[425] -= dqdT;               /* dwdot[O]/dT */

    /*reaction 18: 2 OH => O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[3];
    k_f = prefactor_units[17] * fwd_A[17]
                * exp(fwd_beta[17] * tc[0] - activation_units[17] * fwd_Ea[17] * invT);
    dlnkfdT = fwd_beta[17] * invT + activation_units[17] * fwd_Ea[17] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= 2 * q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[5] += q; /* O */
    /* d()/d[OH] */
    dqdci =  + k_f*2*sc[3];
    J[66] += -2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[67] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[68] += dqdci;               /* dwdot[O]/d[OH] */
    /* d()/dT */
    J[423] += -2 * dqdT;          /* dwdot[OH]/dT */
    J[424] += dqdT;               /* dwdot[H2O]/dT */
    J[425] += dqdT;               /* dwdot[O]/dT */

    /*reaction 19: OH + H2 => H2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[9];
    k_f = prefactor_units[18] * fwd_A[18]
                * exp(fwd_beta[18] * tc[0] - activation_units[18] * fwd_Ea[18] * invT);
    dlnkfdT = fwd_beta[18] * invT + activation_units[18] * fwd_Ea[18] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[9] -= q; /* H2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[9];
    J[65] += dqdci;               /* dwdot[H]/d[OH] */
    J[66] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[67] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[72] -= dqdci;               /* dwdot[H2]/d[OH] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[3];
    J[191] += dqdci;              /* dwdot[H]/d[H2] */
    J[192] -= dqdci;              /* dwdot[OH]/d[H2] */
    J[193] += dqdci;              /* dwdot[H2O]/d[H2] */
    J[198] -= dqdci;              /* dwdot[H2]/d[H2] */
    /* d()/dT */
    J[422] += dqdT;               /* dwdot[H]/dT */
    J[423] -= dqdT;               /* dwdot[OH]/dT */
    J[424] += dqdT;               /* dwdot[H2O]/dT */
    J[429] -= dqdT;               /* dwdot[H2]/dT */

    /*reaction 20: H2O + H => H2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[4];
    k_f = prefactor_units[19] * fwd_A[19]
                * exp(fwd_beta[19] * tc[0] - activation_units[19] * fwd_Ea[19] * invT);
    dlnkfdT = fwd_beta[19] * invT + activation_units[19] * fwd_Ea[19] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[3] += q; /* OH */
    wdot[4] -= q; /* H2O */
    wdot[9] += q; /* H2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[4];
    J[44] -= dqdci;               /* dwdot[H]/d[H] */
    J[45] += dqdci;               /* dwdot[OH]/d[H] */
    J[46] -= dqdci;               /* dwdot[H2O]/d[H] */
    J[51] += dqdci;               /* dwdot[H2]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[2];
    J[86] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[87] += dqdci;               /* dwdot[OH]/d[H2O] */
    J[88] -= dqdci;               /* dwdot[H2O]/d[H2O] */
    J[93] += dqdci;               /* dwdot[H2]/d[H2O] */
    /* d()/dT */
    J[422] -= dqdT;               /* dwdot[H]/dT */
    J[423] += dqdT;               /* dwdot[OH]/dT */
    J[424] -= dqdT;               /* dwdot[H2O]/dT */
    J[429] += dqdT;               /* dwdot[H2]/dT */

    /*reaction 21: CH3 + O <=> CH + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[5];
    k_f = prefactor_units[20] * fwd_A[20]
                * exp(fwd_beta[20] * tc[0] - activation_units[20] * fwd_Ea[20] * invT);
    dlnkfdT = fwd_beta[20] * invT + activation_units[20] * fwd_Ea[20] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[14];
    Kc = exp(g_RT[0] - g_RT[4] + g_RT[5] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[4] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* CH3 */
    wdot[4] += q; /* H2O */
    wdot[5] -= q; /* O */
    wdot[14] += q; /* CH */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[5];
    J[0] -= dqdci;                /* dwdot[CH3]/d[CH3] */
    J[4] += dqdci;                /* dwdot[H2O]/d[CH3] */
    J[5] -= dqdci;                /* dwdot[O]/d[CH3] */
    J[14] += dqdci;               /* dwdot[CH]/d[CH3] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[14];
    J[84] -= dqdci;               /* dwdot[CH3]/d[H2O] */
    J[88] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[89] -= dqdci;               /* dwdot[O]/d[H2O] */
    J[98] += dqdci;               /* dwdot[CH]/d[H2O] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[0];
    J[105] -= dqdci;              /* dwdot[CH3]/d[O] */
    J[109] += dqdci;              /* dwdot[H2O]/d[O] */
    J[110] -= dqdci;              /* dwdot[O]/d[O] */
    J[119] += dqdci;              /* dwdot[CH]/d[O] */
    /* d()/d[CH] */
    dqdci =  - k_r*sc[4];
    J[294] -= dqdci;              /* dwdot[CH3]/d[CH] */
    J[298] += dqdci;              /* dwdot[H2O]/d[CH] */
    J[299] -= dqdci;              /* dwdot[O]/d[CH] */
    J[308] += dqdci;              /* dwdot[CH]/d[CH] */
    /* d()/dT */
    J[420] -= dqdT;               /* dwdot[CH3]/dT */
    J[424] += dqdT;               /* dwdot[H2O]/dT */
    J[425] -= dqdT;               /* dwdot[O]/dT */
    J[434] += dqdT;               /* dwdot[CH]/dT */

    /*reaction 22: CH + O <=> HCOp + E */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[14];
    k_f = prefactor_units[21] * fwd_A[21]
                * exp(fwd_beta[21] * tc[0] - activation_units[21] * fwd_Ea[21] * invT);
    dlnkfdT = fwd_beta[21] * invT + activation_units[21] * fwd_Ea[21] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[18];
    Kc = exp(g_RT[5] - g_RT[6] + g_RT[14] - g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[14]) + (h_RT[6] + h_RT[18]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[6] += q; /* HCOp */
    wdot[14] -= q; /* CH */
    wdot[18] += q; /* E */
    /* d()/d[O] */
    dqdci =  + k_f*sc[14];
    J[110] -= dqdci;              /* dwdot[O]/d[O] */
    J[111] += dqdci;              /* dwdot[HCOp]/d[O] */
    J[119] -= dqdci;              /* dwdot[CH]/d[O] */
    J[123] += dqdci;              /* dwdot[E]/d[O] */
    /* d()/d[HCOp] */
    dqdci =  - k_r*sc[18];
    J[131] -= dqdci;              /* dwdot[O]/d[HCOp] */
    J[132] += dqdci;              /* dwdot[HCOp]/d[HCOp] */
    J[140] -= dqdci;              /* dwdot[CH]/d[HCOp] */
    J[144] += dqdci;              /* dwdot[E]/d[HCOp] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[5];
    J[299] -= dqdci;              /* dwdot[O]/d[CH] */
    J[300] += dqdci;              /* dwdot[HCOp]/d[CH] */
    J[308] -= dqdci;              /* dwdot[CH]/d[CH] */
    J[312] += dqdci;              /* dwdot[E]/d[CH] */
    /* d()/d[E] */
    dqdci =  - k_r*sc[6];
    J[383] -= dqdci;              /* dwdot[O]/d[E] */
    J[384] += dqdci;              /* dwdot[HCOp]/d[E] */
    J[392] -= dqdci;              /* dwdot[CH]/d[E] */
    J[396] += dqdci;              /* dwdot[E]/d[E] */
    /* d()/dT */
    J[425] -= dqdT;               /* dwdot[O]/dT */
    J[426] += dqdT;               /* dwdot[HCOp]/dT */
    J[434] -= dqdT;               /* dwdot[CH]/dT */
    J[438] += dqdT;               /* dwdot[E]/dT */

    /*reaction 23: HCOp + H2O <=> CO + H3Op */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[6];
    k_f = prefactor_units[22] * fwd_A[22]
                * exp(fwd_beta[22] * tc[0] - activation_units[22] * fwd_Ea[22] * invT);
    dlnkfdT = fwd_beta[22] * invT + activation_units[22] * fwd_Ea[22] * invT2;
    /* reverse */
    phi_r = sc[12]*sc[15];
    Kc = exp(g_RT[4] + g_RT[6] - g_RT[12] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[6]) + (h_RT[12] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* H2O */
    wdot[6] -= q; /* HCOp */
    wdot[12] += q; /* CO */
    wdot[15] += q; /* H3Op */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[6];
    J[88] -= dqdci;               /* dwdot[H2O]/d[H2O] */
    J[90] -= dqdci;               /* dwdot[HCOp]/d[H2O] */
    J[96] += dqdci;               /* dwdot[CO]/d[H2O] */
    J[99] += dqdci;               /* dwdot[H3Op]/d[H2O] */
    /* d()/d[HCOp] */
    dqdci =  + k_f*sc[4];
    J[130] -= dqdci;              /* dwdot[H2O]/d[HCOp] */
    J[132] -= dqdci;              /* dwdot[HCOp]/d[HCOp] */
    J[138] += dqdci;              /* dwdot[CO]/d[HCOp] */
    J[141] += dqdci;              /* dwdot[H3Op]/d[HCOp] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[15];
    J[256] -= dqdci;              /* dwdot[H2O]/d[CO] */
    J[258] -= dqdci;              /* dwdot[HCOp]/d[CO] */
    J[264] += dqdci;              /* dwdot[CO]/d[CO] */
    J[267] += dqdci;              /* dwdot[H3Op]/d[CO] */
    /* d()/d[H3Op] */
    dqdci =  - k_r*sc[12];
    J[319] -= dqdci;              /* dwdot[H2O]/d[H3Op] */
    J[321] -= dqdci;              /* dwdot[HCOp]/d[H3Op] */
    J[327] += dqdci;              /* dwdot[CO]/d[H3Op] */
    J[330] += dqdci;              /* dwdot[H3Op]/d[H3Op] */
    /* d()/dT */
    J[424] -= dqdT;               /* dwdot[H2O]/dT */
    J[426] -= dqdT;               /* dwdot[HCOp]/dT */
    J[432] += dqdT;               /* dwdot[CO]/dT */
    J[435] += dqdT;               /* dwdot[H3Op]/dT */

    /*reaction 24: H3Op + E <=> H2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[15]*sc[18];
    k_f = prefactor_units[23] * fwd_A[23]
                * exp(fwd_beta[23] * tc[0] - activation_units[23] * fwd_Ea[23] * invT);
    dlnkfdT = fwd_beta[23] * invT + activation_units[23] * fwd_Ea[23] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[4];
    Kc = exp(-g_RT[2] - g_RT[4] + g_RT[15] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[15] + h_RT[18]) + (h_RT[2] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[4] += q; /* H2O */
    wdot[15] -= q; /* H3Op */
    wdot[18] -= q; /* E */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4];
    J[44] += dqdci;               /* dwdot[H]/d[H] */
    J[46] += dqdci;               /* dwdot[H2O]/d[H] */
    J[57] -= dqdci;               /* dwdot[H3Op]/d[H] */
    J[60] -= dqdci;               /* dwdot[E]/d[H] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[86] += dqdci;               /* dwdot[H]/d[H2O] */
    J[88] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[99] -= dqdci;               /* dwdot[H3Op]/d[H2O] */
    J[102] -= dqdci;              /* dwdot[E]/d[H2O] */
    /* d()/d[H3Op] */
    dqdci =  + k_f*sc[18];
    J[317] += dqdci;              /* dwdot[H]/d[H3Op] */
    J[319] += dqdci;              /* dwdot[H2O]/d[H3Op] */
    J[330] -= dqdci;              /* dwdot[H3Op]/d[H3Op] */
    J[333] -= dqdci;              /* dwdot[E]/d[H3Op] */
    /* d()/d[E] */
    dqdci =  + k_f*sc[15];
    J[380] += dqdci;              /* dwdot[H]/d[E] */
    J[382] += dqdci;              /* dwdot[H2O]/d[E] */
    J[393] -= dqdci;              /* dwdot[H3Op]/d[E] */
    J[396] -= dqdci;              /* dwdot[E]/d[E] */
    /* d()/dT */
    J[422] += dqdT;               /* dwdot[H]/dT */
    J[424] += dqdT;               /* dwdot[H2O]/dT */
    J[435] -= dqdT;               /* dwdot[H3Op]/dT */
    J[438] -= dqdT;               /* dwdot[E]/dT */

    /*reaction 25: CH + O2 <=> CO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[14];
    k_f = prefactor_units[24] * fwd_A[24]
                * exp(fwd_beta[24] * tc[0] - activation_units[24] * fwd_Ea[24] * invT);
    dlnkfdT = fwd_beta[24] * invT + activation_units[24] * fwd_Ea[24] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[12];
    Kc = exp(-g_RT[3] + g_RT[10] - g_RT[12] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[14]) + (h_RT[3] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[10] -= q; /* O2 */
    wdot[12] += q; /* CO */
    wdot[14] -= q; /* CH */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[12];
    J[66] += dqdci;               /* dwdot[OH]/d[OH] */
    J[73] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[75] += dqdci;               /* dwdot[CO]/d[OH] */
    J[77] -= dqdci;               /* dwdot[CH]/d[OH] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[14];
    J[213] += dqdci;              /* dwdot[OH]/d[O2] */
    J[220] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[222] += dqdci;              /* dwdot[CO]/d[O2] */
    J[224] -= dqdci;              /* dwdot[CH]/d[O2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[3];
    J[255] += dqdci;              /* dwdot[OH]/d[CO] */
    J[262] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[264] += dqdci;              /* dwdot[CO]/d[CO] */
    J[266] -= dqdci;              /* dwdot[CH]/d[CO] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[10];
    J[297] += dqdci;              /* dwdot[OH]/d[CH] */
    J[304] -= dqdci;              /* dwdot[O2]/d[CH] */
    J[306] += dqdci;              /* dwdot[CO]/d[CH] */
    J[308] -= dqdci;              /* dwdot[CH]/d[CH] */
    /* d()/dT */
    J[423] += dqdT;               /* dwdot[OH]/dT */
    J[430] -= dqdT;               /* dwdot[O2]/dT */
    J[432] += dqdT;               /* dwdot[CO]/dT */
    J[434] -= dqdT;               /* dwdot[CH]/dT */

    /*reaction 26: O + N2 => NO + N */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[16];
    k_f = prefactor_units[25] * fwd_A[25]
                * exp(fwd_beta[25] * tc[0] - activation_units[25] * fwd_Ea[25] * invT);
    dlnkfdT = fwd_beta[25] * invT + activation_units[25] * fwd_Ea[25] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[7] += q; /* NO */
    wdot[8] += q; /* N */
    wdot[16] -= q; /* N2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[16];
    J[110] -= dqdci;              /* dwdot[O]/d[O] */
    J[112] += dqdci;              /* dwdot[NO]/d[O] */
    J[113] += dqdci;              /* dwdot[N]/d[O] */
    J[121] -= dqdci;              /* dwdot[N2]/d[O] */
    /* d()/d[N2] */
    dqdci =  + k_f*sc[5];
    J[341] -= dqdci;              /* dwdot[O]/d[N2] */
    J[343] += dqdci;              /* dwdot[NO]/d[N2] */
    J[344] += dqdci;              /* dwdot[N]/d[N2] */
    J[352] -= dqdci;              /* dwdot[N2]/d[N2] */
    /* d()/dT */
    J[425] -= dqdT;               /* dwdot[O]/dT */
    J[427] += dqdT;               /* dwdot[NO]/dT */
    J[428] += dqdT;               /* dwdot[N]/dT */
    J[436] -= dqdT;               /* dwdot[N2]/dT */

    /*reaction 27: NO + N => O + N2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[8];
    k_f = prefactor_units[26] * fwd_A[26]
                * exp(fwd_beta[26] * tc[0] - activation_units[26] * fwd_Ea[26] * invT);
    dlnkfdT = fwd_beta[26] * invT + activation_units[26] * fwd_Ea[26] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O */
    wdot[7] -= q; /* NO */
    wdot[8] -= q; /* N */
    wdot[16] += q; /* N2 */
    /* d()/d[NO] */
    dqdci =  + k_f*sc[8];
    J[152] += dqdci;              /* dwdot[O]/d[NO] */
    J[154] -= dqdci;              /* dwdot[NO]/d[NO] */
    J[155] -= dqdci;              /* dwdot[N]/d[NO] */
    J[163] += dqdci;              /* dwdot[N2]/d[NO] */
    /* d()/d[N] */
    dqdci =  + k_f*sc[7];
    J[173] += dqdci;              /* dwdot[O]/d[N] */
    J[175] -= dqdci;              /* dwdot[NO]/d[N] */
    J[176] -= dqdci;              /* dwdot[N]/d[N] */
    J[184] += dqdci;              /* dwdot[N2]/d[N] */
    /* d()/dT */
    J[425] += dqdT;               /* dwdot[O]/dT */
    J[427] -= dqdT;               /* dwdot[NO]/dT */
    J[428] -= dqdT;               /* dwdot[N]/dT */
    J[436] += dqdT;               /* dwdot[N2]/dT */

    /*reaction 28: N + O2 => NO + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[10];
    k_f = prefactor_units[27] * fwd_A[27]
                * exp(fwd_beta[27] * tc[0] - activation_units[27] * fwd_Ea[27] * invT);
    dlnkfdT = fwd_beta[27] * invT + activation_units[27] * fwd_Ea[27] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O */
    wdot[7] += q; /* NO */
    wdot[8] -= q; /* N */
    wdot[10] -= q; /* O2 */
    /* d()/d[N] */
    dqdci =  + k_f*sc[10];
    J[173] += dqdci;              /* dwdot[O]/d[N] */
    J[175] += dqdci;              /* dwdot[NO]/d[N] */
    J[176] -= dqdci;              /* dwdot[N]/d[N] */
    J[178] -= dqdci;              /* dwdot[O2]/d[N] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[8];
    J[215] += dqdci;              /* dwdot[O]/d[O2] */
    J[217] += dqdci;              /* dwdot[NO]/d[O2] */
    J[218] -= dqdci;              /* dwdot[N]/d[O2] */
    J[220] -= dqdci;              /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[425] += dqdT;               /* dwdot[O]/dT */
    J[427] += dqdT;               /* dwdot[NO]/dT */
    J[428] -= dqdT;               /* dwdot[N]/dT */
    J[430] -= dqdT;               /* dwdot[O2]/dT */

    /*reaction 29: NO + O => N + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[7];
    k_f = prefactor_units[28] * fwd_A[28]
                * exp(fwd_beta[28] * tc[0] - activation_units[28] * fwd_Ea[28] * invT);
    dlnkfdT = fwd_beta[28] * invT + activation_units[28] * fwd_Ea[28] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[7] -= q; /* NO */
    wdot[8] += q; /* N */
    wdot[10] += q; /* O2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[7];
    J[110] -= dqdci;              /* dwdot[O]/d[O] */
    J[112] -= dqdci;              /* dwdot[NO]/d[O] */
    J[113] += dqdci;              /* dwdot[N]/d[O] */
    J[115] += dqdci;              /* dwdot[O2]/d[O] */
    /* d()/d[NO] */
    dqdci =  + k_f*sc[5];
    J[152] -= dqdci;              /* dwdot[O]/d[NO] */
    J[154] -= dqdci;              /* dwdot[NO]/d[NO] */
    J[155] += dqdci;              /* dwdot[N]/d[NO] */
    J[157] += dqdci;              /* dwdot[O2]/d[NO] */
    /* d()/dT */
    J[425] -= dqdT;               /* dwdot[O]/dT */
    J[427] -= dqdT;               /* dwdot[NO]/dT */
    J[428] += dqdT;               /* dwdot[N]/dT */
    J[430] += dqdT;               /* dwdot[O2]/dT */

    /*reaction 30: NO2 + O <=> NO + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[17];
    k_f = prefactor_units[29] * fwd_A[29]
                * exp(fwd_beta[29] * tc[0] - activation_units[29] * fwd_Ea[29] * invT);
    dlnkfdT = fwd_beta[29] * invT + activation_units[29] * fwd_Ea[29] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[10];
    Kc = exp(g_RT[5] - g_RT[7] - g_RT[10] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[17]) + (h_RT[7] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[7] += q; /* NO */
    wdot[10] += q; /* O2 */
    wdot[17] -= q; /* NO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[17];
    J[110] -= dqdci;              /* dwdot[O]/d[O] */
    J[112] += dqdci;              /* dwdot[NO]/d[O] */
    J[115] += dqdci;              /* dwdot[O2]/d[O] */
    J[122] -= dqdci;              /* dwdot[NO2]/d[O] */
    /* d()/d[NO] */
    dqdci =  - k_r*sc[10];
    J[152] -= dqdci;              /* dwdot[O]/d[NO] */
    J[154] += dqdci;              /* dwdot[NO]/d[NO] */
    J[157] += dqdci;              /* dwdot[O2]/d[NO] */
    J[164] -= dqdci;              /* dwdot[NO2]/d[NO] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[215] -= dqdci;              /* dwdot[O]/d[O2] */
    J[217] += dqdci;              /* dwdot[NO]/d[O2] */
    J[220] += dqdci;              /* dwdot[O2]/d[O2] */
    J[227] -= dqdci;              /* dwdot[NO2]/d[O2] */
    /* d()/d[NO2] */
    dqdci =  + k_f*sc[5];
    J[362] -= dqdci;              /* dwdot[O]/d[NO2] */
    J[364] += dqdci;              /* dwdot[NO]/d[NO2] */
    J[367] += dqdci;              /* dwdot[O2]/d[NO2] */
    J[374] -= dqdci;              /* dwdot[NO2]/d[NO2] */
    /* d()/dT */
    J[425] -= dqdT;               /* dwdot[O]/dT */
    J[427] += dqdT;               /* dwdot[NO]/dT */
    J[430] += dqdT;               /* dwdot[O2]/dT */
    J[437] -= dqdT;               /* dwdot[NO2]/dT */

    /*reaction 31: NO2 + H <=> NO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[17];
    k_f = prefactor_units[30] * fwd_A[30]
                * exp(fwd_beta[30] * tc[0] - activation_units[30] * fwd_Ea[30] * invT);
    dlnkfdT = fwd_beta[30] * invT + activation_units[30] * fwd_Ea[30] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[7];
    Kc = exp(g_RT[2] - g_RT[3] - g_RT[7] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[17]) + (h_RT[3] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[3] += q; /* OH */
    wdot[7] += q; /* NO */
    wdot[17] -= q; /* NO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[17];
    J[44] -= dqdci;               /* dwdot[H]/d[H] */
    J[45] += dqdci;               /* dwdot[OH]/d[H] */
    J[49] += dqdci;               /* dwdot[NO]/d[H] */
    J[59] -= dqdci;               /* dwdot[NO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[7];
    J[65] -= dqdci;               /* dwdot[H]/d[OH] */
    J[66] += dqdci;               /* dwdot[OH]/d[OH] */
    J[70] += dqdci;               /* dwdot[NO]/d[OH] */
    J[80] -= dqdci;               /* dwdot[NO2]/d[OH] */
    /* d()/d[NO] */
    dqdci =  - k_r*sc[3];
    J[149] -= dqdci;              /* dwdot[H]/d[NO] */
    J[150] += dqdci;              /* dwdot[OH]/d[NO] */
    J[154] += dqdci;              /* dwdot[NO]/d[NO] */
    J[164] -= dqdci;              /* dwdot[NO2]/d[NO] */
    /* d()/d[NO2] */
    dqdci =  + k_f*sc[2];
    J[359] -= dqdci;              /* dwdot[H]/d[NO2] */
    J[360] += dqdci;              /* dwdot[OH]/d[NO2] */
    J[364] += dqdci;              /* dwdot[NO]/d[NO2] */
    J[374] -= dqdci;              /* dwdot[NO2]/d[NO2] */
    /* d()/dT */
    J[422] -= dqdT;               /* dwdot[H]/dT */
    J[423] += dqdT;               /* dwdot[OH]/dT */
    J[427] += dqdT;               /* dwdot[NO]/dT */
    J[437] -= dqdT;               /* dwdot[NO2]/dT */

    double c_R[20], dcRdT[20], e_RT[20];
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
    for (int k = 0; k < 20; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[420+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 20; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 20; ++m) {
            dehmixdc += eh_RT[m]*J[k*21+m];
        }
        J[k*21+20] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[440] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}


/*compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void dcvpRdT(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: CH3 */
        species[0] =
            +2.01095175e-03
            +1.14604371e-05 * tc[1]
            -2.06135228e-08 * tc[2]
            +1.01754294e-11 * tc[3];
        /*species 1: CH4 */
        species[1] =
            -1.36709788e-02
            +9.83601198e-05 * tc[1]
            -1.45422908e-07 * tc[2]
            +6.66775824e-11 * tc[3];
        /*species 2: H */
        species[2] =
            +7.05332819e-13
            -3.99183928e-15 * tc[1]
            +6.90244896e-18 * tc[2]
            -3.71092933e-21 * tc[3];
        /*species 3: OH */
        species[3] =
            -2.40131752e-03
            +9.23587682e-06 * tc[1]
            -1.16434000e-08 * tc[2]
            +5.45645880e-12 * tc[3];
        /*species 4: H2O */
        species[4] =
            -2.03643410e-03
            +1.30408042e-05 * tc[1]
            -1.64639119e-08 * tc[2]
            +7.08791268e-12 * tc[3];
        /*species 5: O */
        species[5] =
            -3.27931884e-03
            +1.32861279e-05 * tc[1]
            -1.83841987e-08 * tc[2]
            +8.45063884e-12 * tc[3];
        /*species 7: NO */
        species[7] =
            -4.63897600e-03
            +2.20820440e-05 * tc[1]
            -2.80084062e-08 * tc[2]
            +1.12143080e-11 * tc[3];
        /*species 8: N */
        species[8] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 9: H2 */
        species[9] =
            +7.98052075e-03
            -3.89563020e-05 * tc[1]
            +6.04716282e-08 * tc[2]
            -2.95044704e-11 * tc[3];
        /*species 10: O2 */
        species[10] =
            -2.99673416e-03
            +1.96946040e-05 * tc[1]
            -2.90438853e-08 * tc[2]
            +1.29749135e-11 * tc[3];
        /*species 11: CH2O */
        species[11] =
            -9.90833369e-03
            +7.46440016e-05 * tc[1]
            -1.13785578e-07 * tc[2]
            +5.27090608e-11 * tc[3];
        /*species 12: CO */
        species[12] =
            -6.10353680e-04
            +2.03362866e-06 * tc[1]
            +2.72101765e-09 * tc[2]
            -3.61769800e-12 * tc[3];
        /*species 13: CO2 */
        species[13] =
            +8.98459677e-03
            -1.42471254e-05 * tc[1]
            +7.37757066e-09 * tc[2]
            -5.74798192e-13 * tc[3];
        /*species 14: CH */
        species[14] =
            +3.23835541e-04
            -3.37798130e-06 * tc[1]
            +9.48651981e-09 * tc[2]
            -5.62436268e-12 * tc[3];
        /*species 15: H3Op */
        species[15] =
            -9.10852723e-04
            +2.32727042e-05 * tc[1]
            -3.64094595e-08 * tc[2]
            +1.70463850e-11 * tc[3];
        /*species 16: N2 */
        species[16] =
            +1.40824040e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77941600e-12 * tc[3];
        /*species 17: NO2 */
        species[17] =
            -1.58542900e-03
            +3.33156240e-05 * tc[1]
            -6.14262780e-08 * tc[2]
            +3.13402256e-11 * tc[3];
        /*species 18: E */
        species[18] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 19: CH2 */
        species[19] =
            +9.68872143e-04
            +5.58979682e-06 * tc[1]
            -1.15527346e-08 * tc[2]
            +6.74966876e-12 * tc[3];
    } else {
        /*species 0: CH3 */
        species[0] =
            +7.23990037e-03
            -5.97428696e-06 * tc[1]
            +1.78705393e-09 * tc[2]
            -1.86861758e-13 * tc[3];
        /*species 1: CH4 */
        species[1] =
            +1.33909467e-02
            -1.14657162e-05 * tc[1]
            +3.66877605e-09 * tc[2]
            -4.07260920e-13 * tc[3];
        /*species 2: H */
        species[2] =
            -2.30842973e-11
            +3.23123896e-14 * tc[1]
            -1.42054571e-17 * tc[2]
            +1.99278943e-21 * tc[3];
        /*species 3: OH */
        species[3] =
            +5.48429716e-04
            +2.53010456e-07 * tc[1]
            -2.63838467e-10 * tc[2]
            +4.69649504e-14 * tc[3];
        /*species 4: H2O */
        species[4] =
            +2.17691804e-03
            -3.28145036e-07 * tc[1]
            -2.91125961e-10 * tc[2]
            +6.72803968e-14 * tc[3];
        /*species 5: O */
        species[5] =
            -8.59741137e-05
            +8.38969178e-08 * tc[1]
            -3.00533397e-11 * tc[2]
            +4.91334764e-15 * tc[3];
        /*species 7: NO */
        species[7] =
            +1.19110430e-03
            -8.58340960e-07 * tc[1]
            +2.08373007e-10 * tc[2]
            -1.61344396e-14 * tc[3];
        /*species 8: N */
        species[8] =
            +1.74890650e-04
            -2.38047380e-07 * tc[1]
            +9.06787350e-11 * tc[2]
            -8.14439280e-15 * tc[3];
        /*species 9: H2 */
        species[9] =
            -4.94024731e-05
            +9.98913556e-07 * tc[1]
            -5.38699182e-10 * tc[2]
            +8.01021504e-14 * tc[3];
        /*species 10: O2 */
        species[10] =
            +1.48308754e-03
            -1.51593334e-06 * tc[1]
            +6.28411665e-10 * tc[2]
            -8.66871176e-14 * tc[3];
        /*species 11: CH2O */
        species[11] =
            +9.20000082e-03
            -8.84517626e-06 * tc[1]
            +3.01923636e-09 * tc[2]
            -3.53542256e-13 * tc[3];
        /*species 12: CO */
        species[12] =
            +2.06252743e-03
            -1.99765154e-06 * tc[1]
            +6.90159024e-10 * tc[2]
            -8.14590864e-14 * tc[3];
        /*species 13: CO2 */
        species[13] =
            +4.41437026e-03
            -4.42962808e-06 * tc[1]
            +1.57047056e-09 * tc[2]
            -1.88833666e-13 * tc[3];
        /*species 14: CH */
        species[14] =
            +9.70913681e-04
            +2.88891310e-07 * tc[1]
            -3.92063547e-10 * tc[2]
            +7.04317532e-14 * tc[3];
        /*species 15: H3Op */
        species[15] =
            +5.72844840e-03
            -3.67906478e-06 * tc[1]
            +8.20732044e-10 * tc[2]
            -6.16375668e-14 * tc[3];
        /*species 16: N2 */
        species[16] =
            +1.48797680e-03
            -1.13695200e-06 * tc[1]
            +3.02911140e-10 * tc[2]
            -2.70134040e-14 * tc[3];
        /*species 17: NO2 */
        species[17] =
            +2.17239560e-03
            -1.65613812e-06 * tc[1]
            +4.72425300e-10 * tc[2]
            -4.20435800e-14 * tc[3];
        /*species 18: E */
        species[18] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 19: CH2 */
        species[19] =
            +3.65639292e-03
            -2.81789194e-06 * tc[1]
            +7.80538647e-10 * tc[2]
            -7.50910268e-14 * tc[3];
    }

    /*species with midpoint at T=3654.18 kelvin */
    if (T < 3654.18) {
        /*species 6: HCOp */
        species[6] =
            +4.94860200e-03
            -4.61706120e-06 * tc[1]
            +1.53365256e-09 * tc[2]
            -1.74082524e-13 * tc[3];
    } else {
        /*species 6: HCOp */
        species[6] =
            +1.66311070e-04
            -3.69526860e-08 * tc[1]
            +2.72569287e-12 * tc[2]
            -6.58177640e-17 * tc[3];
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

    double q_f[31], q_r[31];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 31; ++i) {
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

    /*reaction 1: H + OH + M <=> H2O + M */
    kc[0] = 1.0 / (refC) * exp((g_RT[2] + g_RT[3]) - (g_RT[4]));

    /*reaction 2: O + O + M <=> O2 + M */
    kc[1] = 1.0 / (refC) * exp((g_RT[5] + g_RT[5]) - (g_RT[10]));

    /*reaction 3: H + H + M <=> H2 + M */
    kc[2] = 1.0 / (refC) * exp((g_RT[2] + g_RT[2]) - (g_RT[9]));

    /*reaction 4: NO + O + M <=> NO2 + M */
    kc[3] = 1.0 / (refC) * exp((g_RT[7] + g_RT[5]) - (g_RT[17]));

    /*reaction 5: CH4 <=> CH3 + H */
    kc[4] = refC * exp((g_RT[1]) - (g_RT[0] + g_RT[2]));

    /*reaction 6: CH4 + OH <=> CH3 + H2O */
    kc[5] = exp((g_RT[1] + g_RT[3]) - (g_RT[0] + g_RT[4]));

    /*reaction 7: CH4 + O <=> CH3 + OH */
    kc[6] = exp((g_RT[1] + g_RT[5]) - (g_RT[0] + g_RT[3]));

    /*reaction 8: CH4 + H <=> CH3 + H2 */
    kc[7] = exp((g_RT[1] + g_RT[2]) - (g_RT[0] + g_RT[9]));

    /*reaction 9: CH3 + O2 <=> CH2O + OH */
    kc[8] = exp((g_RT[0] + g_RT[10]) - (g_RT[11] + g_RT[3]));

    /*reaction 10: CH2O + OH <=> CO + H2O + H */
    kc[9] = refC * exp((g_RT[11] + g_RT[3]) - (g_RT[12] + g_RT[4] + g_RT[2]));

    /*reaction 11: CO + OH => CO2 + H */
    kc[10] = exp((g_RT[12] + g_RT[3]) - (g_RT[13] + g_RT[2]));

    /*reaction 12: CO2 + H => CO + OH */
    kc[11] = exp((g_RT[13] + g_RT[2]) - (g_RT[12] + g_RT[3]));

    /*reaction 13: O2 + H => OH + O */
    kc[12] = exp((g_RT[10] + g_RT[2]) - (g_RT[3] + g_RT[5]));

    /*reaction 14: OH + O => O2 + H */
    kc[13] = exp((g_RT[3] + g_RT[5]) - (g_RT[10] + g_RT[2]));

    /*reaction 15: O + H2 => OH + H */
    kc[14] = exp((g_RT[5] + g_RT[9]) - (g_RT[3] + g_RT[2]));

    /*reaction 16: OH + H => O + H2 */
    kc[15] = exp((g_RT[3] + g_RT[2]) - (g_RT[5] + g_RT[9]));

    /*reaction 17: O + H2O => 2 OH */
    kc[16] = exp((g_RT[5] + g_RT[4]) - (2 * g_RT[3]));

    /*reaction 18: 2 OH => O + H2O */
    kc[17] = exp((2 * g_RT[3]) - (g_RT[5] + g_RT[4]));

    /*reaction 19: OH + H2 => H2O + H */
    kc[18] = exp((g_RT[3] + g_RT[9]) - (g_RT[4] + g_RT[2]));

    /*reaction 20: H2O + H => H2 + OH */
    kc[19] = exp((g_RT[4] + g_RT[2]) - (g_RT[9] + g_RT[3]));

    /*reaction 21: CH3 + O <=> CH + H2O */
    kc[20] = exp((g_RT[0] + g_RT[5]) - (g_RT[14] + g_RT[4]));

    /*reaction 22: CH + O <=> HCOp + E */
    kc[21] = exp((g_RT[14] + g_RT[5]) - (g_RT[6] + g_RT[18]));

    /*reaction 23: HCOp + H2O <=> CO + H3Op */
    kc[22] = exp((g_RT[6] + g_RT[4]) - (g_RT[12] + g_RT[15]));

    /*reaction 24: H3Op + E <=> H2O + H */
    kc[23] = exp((g_RT[15] + g_RT[18]) - (g_RT[4] + g_RT[2]));

    /*reaction 25: CH + O2 <=> CO + OH */
    kc[24] = exp((g_RT[14] + g_RT[10]) - (g_RT[12] + g_RT[3]));

    /*reaction 26: O + N2 => NO + N */
    kc[25] = exp((g_RT[5] + g_RT[16]) - (g_RT[7] + g_RT[8]));

    /*reaction 27: NO + N => O + N2 */
    kc[26] = exp((g_RT[7] + g_RT[8]) - (g_RT[5] + g_RT[16]));

    /*reaction 28: N + O2 => NO + O */
    kc[27] = exp((g_RT[8] + g_RT[10]) - (g_RT[7] + g_RT[5]));

    /*reaction 29: NO + O => N + O2 */
    kc[28] = exp((g_RT[7] + g_RT[5]) - (g_RT[8] + g_RT[10]));

    /*reaction 30: NO2 + O <=> NO + O2 */
    kc[29] = exp((g_RT[17] + g_RT[5]) - (g_RT[7] + g_RT[10]));

    /*reaction 31: NO2 + H <=> NO + OH */
    kc[30] = exp((g_RT[17] + g_RT[2]) - (g_RT[7] + g_RT[3]));

    return;
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void gibbs(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: CH3 */
        species[0] =
            +1.644499880000000e+04 * invT
            +2.069026070000000e+00
            -3.673590400000000e+00 * tc[0]
            -1.005475875000000e-03 * tc[1]
            -9.550364266666668e-07 * tc[2]
            +5.725978541666666e-10 * tc[3]
            -1.271928670000000e-13 * tc[4];
        /*species 1: CH4 */
        species[1] =
            -1.024664760000000e+04 * invT
            +9.791179889999999e+00
            -5.149876130000000e+00 * tc[0]
            +6.835489400000000e-03 * tc[1]
            -8.196676650000000e-06 * tc[2]
            +4.039525216666667e-09 * tc[3]
            -8.334697800000000e-13 * tc[4];
        /*species 2: H */
        species[2] =
            +2.547365990000000e+04 * invT
            +2.946682853000000e+00
            -2.500000000000000e+00 * tc[0]
            -3.526664095000000e-13 * tc[1]
            +3.326532733333333e-16 * tc[2]
            -1.917346933333333e-19 * tc[3]
            +4.638661660000000e-23 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.615080560000000e+03 * invT
            +4.095940888000000e+00
            -3.992015430000000e+00 * tc[0]
            +1.200658760000000e-03 * tc[1]
            -7.696564016666666e-07 * tc[2]
            +3.234277775000000e-10 * tc[3]
            -6.820573500000000e-14 * tc[4];
        /*species 4: H2O */
        species[4] =
            -3.029372670000000e+04 * invT
            +5.047672768000000e+00
            -4.198640560000000e+00 * tc[0]
            +1.018217050000000e-03 * tc[1]
            -1.086733685000000e-06 * tc[2]
            +4.573308850000000e-10 * tc[3]
            -8.859890850000000e-14 * tc[4];
        /*species 5: O */
        species[5] =
            +2.912225920000000e+04 * invT
            +1.116333640000000e+00
            -3.168267100000000e+00 * tc[0]
            +1.639659420000000e-03 * tc[1]
            -1.107177326666667e-06 * tc[2]
            +5.106721866666666e-10 * tc[3]
            -1.056329855000000e-13 * tc[4];
        /*species 7: NO */
        species[7] =
            +9.844623000000000e+03 * invT
            +1.937629900000000e+00
            -4.218476300000000e+00 * tc[0]
            +2.319488000000000e-03 * tc[1]
            -1.840170333333333e-06 * tc[2]
            +7.780112833333333e-10 * tc[3]
            -1.401788500000000e-13 * tc[4];
        /*species 8: N */
        species[8] =
            +5.610463700000000e+04 * invT
            -1.693908700000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 9: H2 */
        species[9] =
            -9.179351730000000e+02 * invT
            +1.661320882000000e+00
            -2.344331120000000e+00 * tc[0]
            -3.990260375000000e-03 * tc[1]
            +3.246358500000000e-06 * tc[2]
            -1.679767450000000e-09 * tc[3]
            +3.688058805000000e-13 * tc[4];
        /*species 10: O2 */
        species[10] =
            -1.063943560000000e+03 * invT
            +1.247806300000001e-01
            -3.782456360000000e+00 * tc[0]
            +1.498367080000000e-03 * tc[1]
            -1.641217001666667e-06 * tc[2]
            +8.067745908333334e-10 * tc[3]
            -1.621864185000000e-13 * tc[4];
        /*species 11: CH2O */
        species[11] =
            -1.430895670000000e+04 * invT
            +4.190910250000000e+00
            -4.793723150000000e+00 * tc[0]
            +4.954166845000000e-03 * tc[1]
            -6.220333466666666e-06 * tc[2]
            +3.160710508333333e-09 * tc[3]
            -6.588632600000000e-13 * tc[4];
        /*species 12: CO */
        species[12] =
            -1.434408600000000e+04 * invT
            +7.112418999999992e-02
            -3.579533470000000e+00 * tc[0]
            +3.051768400000000e-04 * tc[1]
            -1.694690550000000e-07 * tc[2]
            -7.558382366666667e-11 * tc[3]
            +4.522122495000000e-14 * tc[4];
        /*species 13: CO2 */
        species[13] =
            -4.837196970000000e+04 * invT
            -7.544278700000000e+00
            -2.356773520000000e+00 * tc[0]
            -4.492298385000000e-03 * tc[1]
            +1.187260448333333e-06 * tc[2]
            -2.049325183333333e-10 * tc[3]
            +7.184977399999999e-15 * tc[4];
        /*species 14: CH */
        species[14] =
            +7.079729340000000e+04 * invT
            +1.405805570000000e+00
            -3.489816650000000e+00 * tc[0]
            -1.619177705000000e-04 * tc[1]
            +2.814984416666667e-07 * tc[2]
            -2.635144391666666e-10 * tc[3]
            +7.030453350000001e-14 * tc[4];
        /*species 15: H3Op */
        species[15] =
            +7.140275180000000e+04 * invT
            +2.321383240000000e+00
            -3.792952510000000e+00 * tc[0]
            +4.554263615000000e-04 * tc[1]
            -1.939392016666666e-06 * tc[2]
            +1.011373875000000e-09 * tc[3]
            -2.130798120000000e-13 * tc[4];
        /*species 16: N2 */
        species[16] =
            -1.020899900000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041202000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427000000000e-13 * tc[4];
        /*species 17: NO2 */
        species[17] =
            +2.896617900000000e+03 * invT
            -2.367960500000000e+00
            -3.944031200000000e+00 * tc[0]
            +7.927145000000000e-04 * tc[1]
            -2.776302000000000e-06 * tc[2]
            +1.706285500000000e-09 * tc[3]
            -3.917528200000000e-13 * tc[4];
        /*species 18: E */
        species[18] =
            -7.453750000000000e+02 * invT
            +1.422081220000000e+01
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 19: CH2 */
        species[19] =
            +4.600404010000000e+04 * invT
            +2.200146820000000e+00
            -3.762678670000000e+00 * tc[0]
            -4.844360715000000e-04 * tc[1]
            -4.658164016666667e-07 * tc[2]
            +3.209092941666667e-10 * tc[3]
            -8.437085950000000e-14 * tc[4];
    } else {
        /*species 0: CH3 */
        species[0] =
            +1.677558430000000e+04 * invT
            -6.194354070000000e+00
            -2.285717720000000e+00 * tc[0]
            -3.619950185000000e-03 * tc[1]
            +4.978572466666667e-07 * tc[2]
            -4.964038700000000e-11 * tc[3]
            +2.335771970000000e-15 * tc[4];
        /*species 1: CH4 */
        species[1] =
            -9.468344590000001e+03 * invT
            -1.836246650500000e+01
            -7.485149500000000e-02 * tc[0]
            -6.695473350000000e-03 * tc[1]
            +9.554763483333333e-07 * tc[2]
            -1.019104458333333e-10 * tc[3]
            +5.090761500000000e-15 * tc[4];
        /*species 2: H */
        species[2] =
            +2.547365990000000e+04 * invT
            +2.946682924000000e+00
            -2.500000010000000e+00 * tc[0]
            +1.154214865000000e-11 * tc[1]
            -2.692699133333334e-15 * tc[2]
            +3.945960291666667e-19 * tc[3]
            -2.490986785000000e-23 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.858657000000000e+03 * invT
            -1.383808430000000e+00
            -3.092887670000000e+00 * tc[0]
            -2.742148580000000e-04 * tc[1]
            -2.108420466666667e-08 * tc[2]
            +7.328846300000000e-12 * tc[3]
            -5.870618800000000e-16 * tc[4];
        /*species 4: H2O */
        species[4] =
            -3.000429710000000e+04 * invT
            -1.932777610000000e+00
            -3.033992490000000e+00 * tc[0]
            -1.088459020000000e-03 * tc[1]
            +2.734541966666666e-08 * tc[2]
            +8.086832250000000e-12 * tc[3]
            -8.410049600000000e-16 * tc[4];
        /*species 5: O */
        species[5] =
            +2.921757910000000e+04 * invT
            -2.214917859999999e+00
            -2.569420780000000e+00 * tc[0]
            +4.298705685000000e-05 * tc[1]
            -6.991409816666667e-09 * tc[2]
            +8.348149916666666e-13 * tc[3]
            -6.141684549999999e-17 * tc[4];
        /*species 7: NO */
        species[7] =
            +9.920974600000000e+03 * invT
            -3.108697100000001e+00
            -3.260605600000000e+00 * tc[0]
            -5.955521500000000e-04 * tc[1]
            +7.152841333333333e-08 * tc[2]
            -5.788139083333334e-12 * tc[3]
            +2.016804950000000e-16 * tc[4];
        /*species 8: N */
        species[8] =
            +5.613377300000000e+04 * invT
            -2.233666700000000e+00
            -2.415942900000000e+00 * tc[0]
            -8.744532500000000e-05 * tc[1]
            +1.983728166666667e-08 * tc[2]
            -2.518853750000000e-12 * tc[3]
            +1.018049100000000e-16 * tc[4];
        /*species 9: H2 */
        species[9] =
            -9.501589220000000e+02 * invT
            +6.542302510000000e+00
            -3.337279200000000e+00 * tc[0]
            +2.470123655000000e-05 * tc[1]
            -8.324279633333333e-08 * tc[2]
            +1.496386616666667e-11 * tc[3]
            -1.001276880000000e-15 * tc[4];
        /*species 10: O2 */
        species[10] =
            -1.088457720000000e+03 * invT
            -2.170693450000000e+00
            -3.282537840000000e+00 * tc[0]
            -7.415437700000000e-04 * tc[1]
            +1.263277781666667e-07 * tc[2]
            -1.745587958333333e-11 * tc[3]
            +1.083588970000000e-15 * tc[4];
        /*species 11: CH2O */
        species[11] =
            -1.399583230000000e+04 * invT
            -1.189563292000000e+01
            -1.760690080000000e+00 * tc[0]
            -4.600000410000000e-03 * tc[1]
            +7.370980216666666e-07 * tc[2]
            -8.386767666666666e-11 * tc[3]
            +4.419278200000001e-15 * tc[4];
        /*species 12: CO */
        species[12] =
            -1.415187240000000e+04 * invT
            -5.103502110000000e+00
            -2.715185610000000e+00 * tc[0]
            -1.031263715000000e-03 * tc[1]
            +1.664709618333334e-07 * tc[2]
            -1.917108400000000e-11 * tc[3]
            +1.018238580000000e-15 * tc[4];
        /*species 13: CO2 */
        species[13] =
            -4.875916600000000e+04 * invT
            +1.585822230000000e+00
            -3.857460290000000e+00 * tc[0]
            -2.207185130000000e-03 * tc[1]
            +3.691356733333334e-07 * tc[2]
            -4.362418233333334e-11 * tc[3]
            +2.360420820000000e-15 * tc[4];
        /*species 14: CH */
        species[14] =
            +7.101243640000001e+04 * invT
            -2.606515260000000e+00
            -2.878464730000000e+00 * tc[0]
            -4.854568405000000e-04 * tc[1]
            -2.407427583333333e-08 * tc[2]
            +1.089065408333333e-11 * tc[3]
            -8.803969149999999e-16 * tc[4];
        /*species 15: H3Op */
        species[15] =
            +7.162442270000000e+04 * invT
            -4.962027280000000e+00
            -2.496477650000000e+00 * tc[0]
            -2.864224200000000e-03 * tc[1]
            +3.065887316666667e-07 * tc[2]
            -2.279811233333333e-11 * tc[3]
            +7.704695850000000e-16 * tc[4];
        /*species 16: N2 */
        species[16] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439884000000000e-04 * tc[1]
            +9.474600000000001e-08 * tc[2]
            -8.414198333333333e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
        /*species 17: NO2 */
        species[17] =
            +2.316498300000000e+03 * invT
            +5.002171150000000e+00
            -4.884754200000000e+00 * tc[0]
            -1.086197800000000e-03 * tc[1]
            +1.380115100000000e-07 * tc[2]
            -1.312292500000000e-11 * tc[3]
            +5.255447500000000e-16 * tc[4];
        /*species 18: E */
        species[18] =
            -7.453750000000000e+02 * invT
            +1.422081220000000e+01
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 19: CH2 */
        species[19] =
            +4.626360400000000e+04 * invT
            -3.297092110000000e+00
            -2.874101130000000e+00 * tc[0]
            -1.828196460000000e-03 * tc[1]
            +2.348243283333333e-07 * tc[2]
            -2.168162908333333e-11 * tc[3]
            +9.386378350000000e-16 * tc[4];
    }

    /*species with midpoint at T=3654.18 kelvin */
    if (T < 3654.18) {
        /*species 6: HCOp */
        species[6] =
            +9.913002200000000e+04 * invT
            -3.625749800000000e+00
            -2.876276800000000e+00 * tc[0]
            -2.474301000000000e-03 * tc[1]
            +3.847551000000000e-07 * tc[2]
            -4.260146000000000e-11 * tc[3]
            +2.176031550000000e-15 * tc[4];
    } else {
        /*species 6: HCOp */
        species[6] =
            +9.612933300000000e+04 * invT
            +2.500523650000000e+01
            -6.915726500000000e+00 * tc[0]
            -8.315553500000000e-05 * tc[1]
            +3.079390500000000e-09 * tc[2]
            -7.571369083333333e-14 * tc[3]
            +8.227220500000000e-19 * tc[4];
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

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: CH3 */
        species[0] =
            +1.64449988e+04 * invT
            +1.06902607e+00
            -3.67359040e+00 * tc[0]
            -1.00547588e-03 * tc[1]
            -9.55036427e-07 * tc[2]
            +5.72597854e-10 * tc[3]
            -1.27192867e-13 * tc[4];
        /*species 1: CH4 */
        species[1] =
            -1.02466476e+04 * invT
            +8.79117989e+00
            -5.14987613e+00 * tc[0]
            +6.83548940e-03 * tc[1]
            -8.19667665e-06 * tc[2]
            +4.03952522e-09 * tc[3]
            -8.33469780e-13 * tc[4];
        /*species 2: H */
        species[2] =
            +2.54736599e+04 * invT
            +1.94668285e+00
            -2.50000000e+00 * tc[0]
            -3.52666409e-13 * tc[1]
            +3.32653273e-16 * tc[2]
            -1.91734693e-19 * tc[3]
            +4.63866166e-23 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.61508056e+03 * invT
            +3.09594089e+00
            -3.99201543e+00 * tc[0]
            +1.20065876e-03 * tc[1]
            -7.69656402e-07 * tc[2]
            +3.23427778e-10 * tc[3]
            -6.82057350e-14 * tc[4];
        /*species 4: H2O */
        species[4] =
            -3.02937267e+04 * invT
            +4.04767277e+00
            -4.19864056e+00 * tc[0]
            +1.01821705e-03 * tc[1]
            -1.08673369e-06 * tc[2]
            +4.57330885e-10 * tc[3]
            -8.85989085e-14 * tc[4];
        /*species 5: O */
        species[5] =
            +2.91222592e+04 * invT
            +1.16333640e-01
            -3.16826710e+00 * tc[0]
            +1.63965942e-03 * tc[1]
            -1.10717733e-06 * tc[2]
            +5.10672187e-10 * tc[3]
            -1.05632985e-13 * tc[4];
        /*species 7: NO */
        species[7] =
            +9.84462300e+03 * invT
            +9.37629900e-01
            -4.21847630e+00 * tc[0]
            +2.31948800e-03 * tc[1]
            -1.84017033e-06 * tc[2]
            +7.78011283e-10 * tc[3]
            -1.40178850e-13 * tc[4];
        /*species 8: N */
        species[8] =
            +5.61046370e+04 * invT
            -2.69390870e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 9: H2 */
        species[9] =
            -9.17935173e+02 * invT
            +6.61320882e-01
            -2.34433112e+00 * tc[0]
            -3.99026037e-03 * tc[1]
            +3.24635850e-06 * tc[2]
            -1.67976745e-09 * tc[3]
            +3.68805881e-13 * tc[4];
        /*species 10: O2 */
        species[10] =
            -1.06394356e+03 * invT
            -8.75219370e-01
            -3.78245636e+00 * tc[0]
            +1.49836708e-03 * tc[1]
            -1.64121700e-06 * tc[2]
            +8.06774591e-10 * tc[3]
            -1.62186418e-13 * tc[4];
        /*species 11: CH2O */
        species[11] =
            -1.43089567e+04 * invT
            +3.19091025e+00
            -4.79372315e+00 * tc[0]
            +4.95416684e-03 * tc[1]
            -6.22033347e-06 * tc[2]
            +3.16071051e-09 * tc[3]
            -6.58863260e-13 * tc[4];
        /*species 12: CO */
        species[12] =
            -1.43440860e+04 * invT
            -9.28875810e-01
            -3.57953347e+00 * tc[0]
            +3.05176840e-04 * tc[1]
            -1.69469055e-07 * tc[2]
            -7.55838237e-11 * tc[3]
            +4.52212249e-14 * tc[4];
        /*species 13: CO2 */
        species[13] =
            -4.83719697e+04 * invT
            -8.54427870e+00
            -2.35677352e+00 * tc[0]
            -4.49229839e-03 * tc[1]
            +1.18726045e-06 * tc[2]
            -2.04932518e-10 * tc[3]
            +7.18497740e-15 * tc[4];
        /*species 14: CH */
        species[14] =
            +7.07972934e+04 * invT
            +4.05805570e-01
            -3.48981665e+00 * tc[0]
            -1.61917771e-04 * tc[1]
            +2.81498442e-07 * tc[2]
            -2.63514439e-10 * tc[3]
            +7.03045335e-14 * tc[4];
        /*species 15: H3Op */
        species[15] =
            +7.14027518e+04 * invT
            +1.32138324e+00
            -3.79295251e+00 * tc[0]
            +4.55426362e-04 * tc[1]
            -1.93939202e-06 * tc[2]
            +1.01137387e-09 * tc[3]
            -2.13079812e-13 * tc[4];
        /*species 16: N2 */
        species[16] =
            -1.02089990e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
        /*species 17: NO2 */
        species[17] =
            +2.89661790e+03 * invT
            -3.36796050e+00
            -3.94403120e+00 * tc[0]
            +7.92714500e-04 * tc[1]
            -2.77630200e-06 * tc[2]
            +1.70628550e-09 * tc[3]
            -3.91752820e-13 * tc[4];
        /*species 18: E */
        species[18] =
            -7.45375000e+02 * invT
            +1.32208122e+01
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 19: CH2 */
        species[19] =
            +4.60040401e+04 * invT
            +1.20014682e+00
            -3.76267867e+00 * tc[0]
            -4.84436072e-04 * tc[1]
            -4.65816402e-07 * tc[2]
            +3.20909294e-10 * tc[3]
            -8.43708595e-14 * tc[4];
    } else {
        /*species 0: CH3 */
        species[0] =
            +1.67755843e+04 * invT
            -7.19435407e+00
            -2.28571772e+00 * tc[0]
            -3.61995018e-03 * tc[1]
            +4.97857247e-07 * tc[2]
            -4.96403870e-11 * tc[3]
            +2.33577197e-15 * tc[4];
        /*species 1: CH4 */
        species[1] =
            -9.46834459e+03 * invT
            -1.93624665e+01
            -7.48514950e-02 * tc[0]
            -6.69547335e-03 * tc[1]
            +9.55476348e-07 * tc[2]
            -1.01910446e-10 * tc[3]
            +5.09076150e-15 * tc[4];
        /*species 2: H */
        species[2] =
            +2.54736599e+04 * invT
            +1.94668292e+00
            -2.50000001e+00 * tc[0]
            +1.15421486e-11 * tc[1]
            -2.69269913e-15 * tc[2]
            +3.94596029e-19 * tc[3]
            -2.49098679e-23 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.85865700e+03 * invT
            -2.38380843e+00
            -3.09288767e+00 * tc[0]
            -2.74214858e-04 * tc[1]
            -2.10842047e-08 * tc[2]
            +7.32884630e-12 * tc[3]
            -5.87061880e-16 * tc[4];
        /*species 4: H2O */
        species[4] =
            -3.00042971e+04 * invT
            -2.93277761e+00
            -3.03399249e+00 * tc[0]
            -1.08845902e-03 * tc[1]
            +2.73454197e-08 * tc[2]
            +8.08683225e-12 * tc[3]
            -8.41004960e-16 * tc[4];
        /*species 5: O */
        species[5] =
            +2.92175791e+04 * invT
            -3.21491786e+00
            -2.56942078e+00 * tc[0]
            +4.29870569e-05 * tc[1]
            -6.99140982e-09 * tc[2]
            +8.34814992e-13 * tc[3]
            -6.14168455e-17 * tc[4];
        /*species 7: NO */
        species[7] =
            +9.92097460e+03 * invT
            -4.10869710e+00
            -3.26060560e+00 * tc[0]
            -5.95552150e-04 * tc[1]
            +7.15284133e-08 * tc[2]
            -5.78813908e-12 * tc[3]
            +2.01680495e-16 * tc[4];
        /*species 8: N */
        species[8] =
            +5.61337730e+04 * invT
            -3.23366670e+00
            -2.41594290e+00 * tc[0]
            -8.74453250e-05 * tc[1]
            +1.98372817e-08 * tc[2]
            -2.51885375e-12 * tc[3]
            +1.01804910e-16 * tc[4];
        /*species 9: H2 */
        species[9] =
            -9.50158922e+02 * invT
            +5.54230251e+00
            -3.33727920e+00 * tc[0]
            +2.47012365e-05 * tc[1]
            -8.32427963e-08 * tc[2]
            +1.49638662e-11 * tc[3]
            -1.00127688e-15 * tc[4];
        /*species 10: O2 */
        species[10] =
            -1.08845772e+03 * invT
            -3.17069345e+00
            -3.28253784e+00 * tc[0]
            -7.41543770e-04 * tc[1]
            +1.26327778e-07 * tc[2]
            -1.74558796e-11 * tc[3]
            +1.08358897e-15 * tc[4];
        /*species 11: CH2O */
        species[11] =
            -1.39958323e+04 * invT
            -1.28956329e+01
            -1.76069008e+00 * tc[0]
            -4.60000041e-03 * tc[1]
            +7.37098022e-07 * tc[2]
            -8.38676767e-11 * tc[3]
            +4.41927820e-15 * tc[4];
        /*species 12: CO */
        species[12] =
            -1.41518724e+04 * invT
            -6.10350211e+00
            -2.71518561e+00 * tc[0]
            -1.03126372e-03 * tc[1]
            +1.66470962e-07 * tc[2]
            -1.91710840e-11 * tc[3]
            +1.01823858e-15 * tc[4];
        /*species 13: CO2 */
        species[13] =
            -4.87591660e+04 * invT
            +5.85822230e-01
            -3.85746029e+00 * tc[0]
            -2.20718513e-03 * tc[1]
            +3.69135673e-07 * tc[2]
            -4.36241823e-11 * tc[3]
            +2.36042082e-15 * tc[4];
        /*species 14: CH */
        species[14] =
            +7.10124364e+04 * invT
            -3.60651526e+00
            -2.87846473e+00 * tc[0]
            -4.85456840e-04 * tc[1]
            -2.40742758e-08 * tc[2]
            +1.08906541e-11 * tc[3]
            -8.80396915e-16 * tc[4];
        /*species 15: H3Op */
        species[15] =
            +7.16244227e+04 * invT
            -5.96202728e+00
            -2.49647765e+00 * tc[0]
            -2.86422420e-03 * tc[1]
            +3.06588732e-07 * tc[2]
            -2.27981123e-11 * tc[3]
            +7.70469585e-16 * tc[4];
        /*species 16: N2 */
        species[16] =
            -9.22797700e+02 * invT
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988400e-04 * tc[1]
            +9.47460000e-08 * tc[2]
            -8.41419833e-12 * tc[3]
            +3.37667550e-16 * tc[4];
        /*species 17: NO2 */
        species[17] =
            +2.31649830e+03 * invT
            +4.00217115e+00
            -4.88475420e+00 * tc[0]
            -1.08619780e-03 * tc[1]
            +1.38011510e-07 * tc[2]
            -1.31229250e-11 * tc[3]
            +5.25544750e-16 * tc[4];
        /*species 18: E */
        species[18] =
            -7.45375000e+02 * invT
            +1.32208122e+01
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 19: CH2 */
        species[19] =
            +4.62636040e+04 * invT
            -4.29709211e+00
            -2.87410113e+00 * tc[0]
            -1.82819646e-03 * tc[1]
            +2.34824328e-07 * tc[2]
            -2.16816291e-11 * tc[3]
            +9.38637835e-16 * tc[4];
    }

    /*species with midpoint at T=3654.18 kelvin */
    if (T < 3654.18) {
        /*species 6: HCOp */
        species[6] =
            +9.91300220e+04 * invT
            -4.62574980e+00
            -2.87627680e+00 * tc[0]
            -2.47430100e-03 * tc[1]
            +3.84755100e-07 * tc[2]
            -4.26014600e-11 * tc[3]
            +2.17603155e-15 * tc[4];
    } else {
        /*species 6: HCOp */
        species[6] =
            +9.61293330e+04 * invT
            +2.40052365e+01
            -6.91572650e+00 * tc[0]
            -8.31555350e-05 * tc[1]
            +3.07939050e-09 * tc[2]
            -7.57136908e-14 * tc[3]
            +8.22722050e-19 * tc[4];
    }
    return;
}


/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void cv_R(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: CH3 */
        species[0] =
            +2.67359040e+00
            +2.01095175e-03 * tc[1]
            +5.73021856e-06 * tc[2]
            -6.87117425e-09 * tc[3]
            +2.54385734e-12 * tc[4];
        /*species 1: CH4 */
        species[1] =
            +4.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 2: H */
        species[2] =
            +1.50000000e+00
            +7.05332819e-13 * tc[1]
            -1.99591964e-15 * tc[2]
            +2.30081632e-18 * tc[3]
            -9.27732332e-22 * tc[4];
        /*species 3: OH */
        species[3] =
            +2.99201543e+00
            -2.40131752e-03 * tc[1]
            +4.61793841e-06 * tc[2]
            -3.88113333e-09 * tc[3]
            +1.36411470e-12 * tc[4];
        /*species 4: H2O */
        species[4] =
            +3.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 5: O */
        species[5] =
            +2.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 7: NO */
        species[7] =
            +3.21847630e+00
            -4.63897600e-03 * tc[1]
            +1.10410220e-05 * tc[2]
            -9.33613540e-09 * tc[3]
            +2.80357700e-12 * tc[4];
        /*species 8: N */
        species[8] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 9: H2 */
        species[9] =
            +1.34433112e+00
            +7.98052075e-03 * tc[1]
            -1.94781510e-05 * tc[2]
            +2.01572094e-08 * tc[3]
            -7.37611761e-12 * tc[4];
        /*species 10: O2 */
        species[10] =
            +2.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 11: CH2O */
        species[11] =
            +3.79372315e+00
            -9.90833369e-03 * tc[1]
            +3.73220008e-05 * tc[2]
            -3.79285261e-08 * tc[3]
            +1.31772652e-11 * tc[4];
        /*species 12: CO */
        species[12] =
            +2.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 13: CO2 */
        species[13] =
            +1.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 14: CH */
        species[14] =
            +2.48981665e+00
            +3.23835541e-04 * tc[1]
            -1.68899065e-06 * tc[2]
            +3.16217327e-09 * tc[3]
            -1.40609067e-12 * tc[4];
        /*species 15: H3Op */
        species[15] =
            +2.79295251e+00
            -9.10852723e-04 * tc[1]
            +1.16363521e-05 * tc[2]
            -1.21364865e-08 * tc[3]
            +4.26159624e-12 * tc[4];
        /*species 16: N2 */
        species[16] =
            +2.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
        /*species 17: NO2 */
        species[17] =
            +2.94403120e+00
            -1.58542900e-03 * tc[1]
            +1.66578120e-05 * tc[2]
            -2.04754260e-08 * tc[3]
            +7.83505640e-12 * tc[4];
        /*species 18: E */
        species[18] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 19: CH2 */
        species[19] =
            +2.76267867e+00
            +9.68872143e-04 * tc[1]
            +2.79489841e-06 * tc[2]
            -3.85091153e-09 * tc[3]
            +1.68741719e-12 * tc[4];
    } else {
        /*species 0: CH3 */
        species[0] =
            +1.28571772e+00
            +7.23990037e-03 * tc[1]
            -2.98714348e-06 * tc[2]
            +5.95684644e-10 * tc[3]
            -4.67154394e-14 * tc[4];
        /*species 1: CH4 */
        species[1] =
            -9.25148505e-01
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 2: H */
        species[2] =
            +1.50000001e+00
            -2.30842973e-11 * tc[1]
            +1.61561948e-14 * tc[2]
            -4.73515235e-18 * tc[3]
            +4.98197357e-22 * tc[4];
        /*species 3: OH */
        species[3] =
            +2.09288767e+00
            +5.48429716e-04 * tc[1]
            +1.26505228e-07 * tc[2]
            -8.79461556e-11 * tc[3]
            +1.17412376e-14 * tc[4];
        /*species 4: H2O */
        species[4] =
            +2.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 5: O */
        species[5] =
            +1.56942078e+00
            -8.59741137e-05 * tc[1]
            +4.19484589e-08 * tc[2]
            -1.00177799e-11 * tc[3]
            +1.22833691e-15 * tc[4];
        /*species 7: NO */
        species[7] =
            +2.26060560e+00
            +1.19110430e-03 * tc[1]
            -4.29170480e-07 * tc[2]
            +6.94576690e-11 * tc[3]
            -4.03360990e-15 * tc[4];
        /*species 8: N */
        species[8] =
            +1.41594290e+00
            +1.74890650e-04 * tc[1]
            -1.19023690e-07 * tc[2]
            +3.02262450e-11 * tc[3]
            -2.03609820e-15 * tc[4];
        /*species 9: H2 */
        species[9] =
            +2.33727920e+00
            -4.94024731e-05 * tc[1]
            +4.99456778e-07 * tc[2]
            -1.79566394e-10 * tc[3]
            +2.00255376e-14 * tc[4];
        /*species 10: O2 */
        species[10] =
            +2.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 11: CH2O */
        species[11] =
            +7.60690080e-01
            +9.20000082e-03 * tc[1]
            -4.42258813e-06 * tc[2]
            +1.00641212e-09 * tc[3]
            -8.83855640e-14 * tc[4];
        /*species 12: CO */
        species[12] =
            +1.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 13: CO2 */
        species[13] =
            +2.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 14: CH */
        species[14] =
            +1.87846473e+00
            +9.70913681e-04 * tc[1]
            +1.44445655e-07 * tc[2]
            -1.30687849e-10 * tc[3]
            +1.76079383e-14 * tc[4];
        /*species 15: H3Op */
        species[15] =
            +1.49647765e+00
            +5.72844840e-03 * tc[1]
            -1.83953239e-06 * tc[2]
            +2.73577348e-10 * tc[3]
            -1.54093917e-14 * tc[4];
        /*species 16: N2 */
        species[16] =
            +1.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 17: NO2 */
        species[17] =
            +3.88475420e+00
            +2.17239560e-03 * tc[1]
            -8.28069060e-07 * tc[2]
            +1.57475100e-10 * tc[3]
            -1.05108950e-14 * tc[4];
        /*species 18: E */
        species[18] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 19: CH2 */
        species[19] =
            +1.87410113e+00
            +3.65639292e-03 * tc[1]
            -1.40894597e-06 * tc[2]
            +2.60179549e-10 * tc[3]
            -1.87727567e-14 * tc[4];
    }

    /*species with midpoint at T=3654.18 kelvin */
    if (T < 3654.18) {
        /*species 6: HCOp */
        species[6] =
            +1.87627680e+00
            +4.94860200e-03 * tc[1]
            -2.30853060e-06 * tc[2]
            +5.11217520e-10 * tc[3]
            -4.35206310e-14 * tc[4];
    } else {
        /*species 6: HCOp */
        species[6] =
            +5.91572650e+00
            +1.66311070e-04 * tc[1]
            -1.84763430e-08 * tc[2]
            +9.08564290e-13 * tc[3]
            -1.64544410e-17 * tc[4];
    }
    return;
}


/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void cp_R(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: CH3 */
        species[0] =
            +3.67359040e+00
            +2.01095175e-03 * tc[1]
            +5.73021856e-06 * tc[2]
            -6.87117425e-09 * tc[3]
            +2.54385734e-12 * tc[4];
        /*species 1: CH4 */
        species[1] =
            +5.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 2: H */
        species[2] =
            +2.50000000e+00
            +7.05332819e-13 * tc[1]
            -1.99591964e-15 * tc[2]
            +2.30081632e-18 * tc[3]
            -9.27732332e-22 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.99201543e+00
            -2.40131752e-03 * tc[1]
            +4.61793841e-06 * tc[2]
            -3.88113333e-09 * tc[3]
            +1.36411470e-12 * tc[4];
        /*species 4: H2O */
        species[4] =
            +4.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 5: O */
        species[5] =
            +3.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 7: NO */
        species[7] =
            +4.21847630e+00
            -4.63897600e-03 * tc[1]
            +1.10410220e-05 * tc[2]
            -9.33613540e-09 * tc[3]
            +2.80357700e-12 * tc[4];
        /*species 8: N */
        species[8] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 9: H2 */
        species[9] =
            +2.34433112e+00
            +7.98052075e-03 * tc[1]
            -1.94781510e-05 * tc[2]
            +2.01572094e-08 * tc[3]
            -7.37611761e-12 * tc[4];
        /*species 10: O2 */
        species[10] =
            +3.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 11: CH2O */
        species[11] =
            +4.79372315e+00
            -9.90833369e-03 * tc[1]
            +3.73220008e-05 * tc[2]
            -3.79285261e-08 * tc[3]
            +1.31772652e-11 * tc[4];
        /*species 12: CO */
        species[12] =
            +3.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 13: CO2 */
        species[13] =
            +2.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 14: CH */
        species[14] =
            +3.48981665e+00
            +3.23835541e-04 * tc[1]
            -1.68899065e-06 * tc[2]
            +3.16217327e-09 * tc[3]
            -1.40609067e-12 * tc[4];
        /*species 15: H3Op */
        species[15] =
            +3.79295251e+00
            -9.10852723e-04 * tc[1]
            +1.16363521e-05 * tc[2]
            -1.21364865e-08 * tc[3]
            +4.26159624e-12 * tc[4];
        /*species 16: N2 */
        species[16] =
            +3.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
        /*species 17: NO2 */
        species[17] =
            +3.94403120e+00
            -1.58542900e-03 * tc[1]
            +1.66578120e-05 * tc[2]
            -2.04754260e-08 * tc[3]
            +7.83505640e-12 * tc[4];
        /*species 18: E */
        species[18] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 19: CH2 */
        species[19] =
            +3.76267867e+00
            +9.68872143e-04 * tc[1]
            +2.79489841e-06 * tc[2]
            -3.85091153e-09 * tc[3]
            +1.68741719e-12 * tc[4];
    } else {
        /*species 0: CH3 */
        species[0] =
            +2.28571772e+00
            +7.23990037e-03 * tc[1]
            -2.98714348e-06 * tc[2]
            +5.95684644e-10 * tc[3]
            -4.67154394e-14 * tc[4];
        /*species 1: CH4 */
        species[1] =
            +7.48514950e-02
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 2: H */
        species[2] =
            +2.50000001e+00
            -2.30842973e-11 * tc[1]
            +1.61561948e-14 * tc[2]
            -4.73515235e-18 * tc[3]
            +4.98197357e-22 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.09288767e+00
            +5.48429716e-04 * tc[1]
            +1.26505228e-07 * tc[2]
            -8.79461556e-11 * tc[3]
            +1.17412376e-14 * tc[4];
        /*species 4: H2O */
        species[4] =
            +3.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 5: O */
        species[5] =
            +2.56942078e+00
            -8.59741137e-05 * tc[1]
            +4.19484589e-08 * tc[2]
            -1.00177799e-11 * tc[3]
            +1.22833691e-15 * tc[4];
        /*species 7: NO */
        species[7] =
            +3.26060560e+00
            +1.19110430e-03 * tc[1]
            -4.29170480e-07 * tc[2]
            +6.94576690e-11 * tc[3]
            -4.03360990e-15 * tc[4];
        /*species 8: N */
        species[8] =
            +2.41594290e+00
            +1.74890650e-04 * tc[1]
            -1.19023690e-07 * tc[2]
            +3.02262450e-11 * tc[3]
            -2.03609820e-15 * tc[4];
        /*species 9: H2 */
        species[9] =
            +3.33727920e+00
            -4.94024731e-05 * tc[1]
            +4.99456778e-07 * tc[2]
            -1.79566394e-10 * tc[3]
            +2.00255376e-14 * tc[4];
        /*species 10: O2 */
        species[10] =
            +3.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 11: CH2O */
        species[11] =
            +1.76069008e+00
            +9.20000082e-03 * tc[1]
            -4.42258813e-06 * tc[2]
            +1.00641212e-09 * tc[3]
            -8.83855640e-14 * tc[4];
        /*species 12: CO */
        species[12] =
            +2.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 13: CO2 */
        species[13] =
            +3.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 14: CH */
        species[14] =
            +2.87846473e+00
            +9.70913681e-04 * tc[1]
            +1.44445655e-07 * tc[2]
            -1.30687849e-10 * tc[3]
            +1.76079383e-14 * tc[4];
        /*species 15: H3Op */
        species[15] =
            +2.49647765e+00
            +5.72844840e-03 * tc[1]
            -1.83953239e-06 * tc[2]
            +2.73577348e-10 * tc[3]
            -1.54093917e-14 * tc[4];
        /*species 16: N2 */
        species[16] =
            +2.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 17: NO2 */
        species[17] =
            +4.88475420e+00
            +2.17239560e-03 * tc[1]
            -8.28069060e-07 * tc[2]
            +1.57475100e-10 * tc[3]
            -1.05108950e-14 * tc[4];
        /*species 18: E */
        species[18] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 19: CH2 */
        species[19] =
            +2.87410113e+00
            +3.65639292e-03 * tc[1]
            -1.40894597e-06 * tc[2]
            +2.60179549e-10 * tc[3]
            -1.87727567e-14 * tc[4];
    }

    /*species with midpoint at T=3654.18 kelvin */
    if (T < 3654.18) {
        /*species 6: HCOp */
        species[6] =
            +2.87627680e+00
            +4.94860200e-03 * tc[1]
            -2.30853060e-06 * tc[2]
            +5.11217520e-10 * tc[3]
            -4.35206310e-14 * tc[4];
    } else {
        /*species 6: HCOp */
        species[6] =
            +6.91572650e+00
            +1.66311070e-04 * tc[1]
            -1.84763430e-08 * tc[2]
            +9.08564290e-13 * tc[3]
            -1.64544410e-17 * tc[4];
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

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: CH3 */
        species[0] =
            +2.67359040e+00
            +1.00547588e-03 * tc[1]
            +1.91007285e-06 * tc[2]
            -1.71779356e-09 * tc[3]
            +5.08771468e-13 * tc[4]
            +1.64449988e+04 * invT;
        /*species 1: CH4 */
        species[1] =
            +4.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 * invT;
        /*species 2: H */
        species[2] =
            +1.50000000e+00
            +3.52666409e-13 * tc[1]
            -6.65306547e-16 * tc[2]
            +5.75204080e-19 * tc[3]
            -1.85546466e-22 * tc[4]
            +2.54736599e+04 * invT;
        /*species 3: OH */
        species[3] =
            +2.99201543e+00
            -1.20065876e-03 * tc[1]
            +1.53931280e-06 * tc[2]
            -9.70283332e-10 * tc[3]
            +2.72822940e-13 * tc[4]
            +3.61508056e+03 * invT;
        /*species 4: H2O */
        species[4] =
            +3.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
        /*species 5: O */
        species[5] =
            +2.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 * invT;
        /*species 7: NO */
        species[7] =
            +3.21847630e+00
            -2.31948800e-03 * tc[1]
            +3.68034067e-06 * tc[2]
            -2.33403385e-09 * tc[3]
            +5.60715400e-13 * tc[4]
            +9.84462300e+03 * invT;
        /*species 8: N */
        species[8] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +5.61046370e+04 * invT;
        /*species 9: H2 */
        species[9] =
            +1.34433112e+00
            +3.99026037e-03 * tc[1]
            -6.49271700e-06 * tc[2]
            +5.03930235e-09 * tc[3]
            -1.47522352e-12 * tc[4]
            -9.17935173e+02 * invT;
        /*species 10: O2 */
        species[10] =
            +2.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 11: CH2O */
        species[11] =
            +3.79372315e+00
            -4.95416684e-03 * tc[1]
            +1.24406669e-05 * tc[2]
            -9.48213152e-09 * tc[3]
            +2.63545304e-12 * tc[4]
            -1.43089567e+04 * invT;
        /*species 12: CO */
        species[12] =
            +2.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 13: CO2 */
        species[13] =
            +1.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
        /*species 14: CH */
        species[14] =
            +2.48981665e+00
            +1.61917771e-04 * tc[1]
            -5.62996883e-07 * tc[2]
            +7.90543317e-10 * tc[3]
            -2.81218134e-13 * tc[4]
            +7.07972934e+04 * invT;
        /*species 15: H3Op */
        species[15] =
            +2.79295251e+00
            -4.55426362e-04 * tc[1]
            +3.87878403e-06 * tc[2]
            -3.03412162e-09 * tc[3]
            +8.52319248e-13 * tc[4]
            +7.14027518e+04 * invT;
        /*species 16: N2 */
        species[16] =
            +2.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
        /*species 17: NO2 */
        species[17] =
            +2.94403120e+00
            -7.92714500e-04 * tc[1]
            +5.55260400e-06 * tc[2]
            -5.11885650e-09 * tc[3]
            +1.56701128e-12 * tc[4]
            +2.89661790e+03 * invT;
        /*species 18: E */
        species[18] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 19: CH2 */
        species[19] =
            +2.76267867e+00
            +4.84436072e-04 * tc[1]
            +9.31632803e-07 * tc[2]
            -9.62727883e-10 * tc[3]
            +3.37483438e-13 * tc[4]
            +4.60040401e+04 * invT;
    } else {
        /*species 0: CH3 */
        species[0] =
            +1.28571772e+00
            +3.61995018e-03 * tc[1]
            -9.95714493e-07 * tc[2]
            +1.48921161e-10 * tc[3]
            -9.34308788e-15 * tc[4]
            +1.67755843e+04 * invT;
        /*species 1: CH4 */
        species[1] =
            -9.25148505e-01
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 * invT;
        /*species 2: H */
        species[2] =
            +1.50000001e+00
            -1.15421486e-11 * tc[1]
            +5.38539827e-15 * tc[2]
            -1.18378809e-18 * tc[3]
            +9.96394714e-23 * tc[4]
            +2.54736599e+04 * invT;
        /*species 3: OH */
        species[3] =
            +2.09288767e+00
            +2.74214858e-04 * tc[1]
            +4.21684093e-08 * tc[2]
            -2.19865389e-11 * tc[3]
            +2.34824752e-15 * tc[4]
            +3.85865700e+03 * invT;
        /*species 4: H2O */
        species[4] =
            +2.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
        /*species 5: O */
        species[5] =
            +1.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 * invT;
        /*species 7: NO */
        species[7] =
            +2.26060560e+00
            +5.95552150e-04 * tc[1]
            -1.43056827e-07 * tc[2]
            +1.73644173e-11 * tc[3]
            -8.06721980e-16 * tc[4]
            +9.92097460e+03 * invT;
        /*species 8: N */
        species[8] =
            +1.41594290e+00
            +8.74453250e-05 * tc[1]
            -3.96745633e-08 * tc[2]
            +7.55656125e-12 * tc[3]
            -4.07219640e-16 * tc[4]
            +5.61337730e+04 * invT;
        /*species 9: H2 */
        species[9] =
            +2.33727920e+00
            -2.47012365e-05 * tc[1]
            +1.66485593e-07 * tc[2]
            -4.48915985e-11 * tc[3]
            +4.00510752e-15 * tc[4]
            -9.50158922e+02 * invT;
        /*species 10: O2 */
        species[10] =
            +2.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 11: CH2O */
        species[11] =
            +7.60690080e-01
            +4.60000041e-03 * tc[1]
            -1.47419604e-06 * tc[2]
            +2.51603030e-10 * tc[3]
            -1.76771128e-14 * tc[4]
            -1.39958323e+04 * invT;
        /*species 12: CO */
        species[12] =
            +1.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 13: CO2 */
        species[13] =
            +2.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
        /*species 14: CH */
        species[14] =
            +1.87846473e+00
            +4.85456840e-04 * tc[1]
            +4.81485517e-08 * tc[2]
            -3.26719623e-11 * tc[3]
            +3.52158766e-15 * tc[4]
            +7.10124364e+04 * invT;
        /*species 15: H3Op */
        species[15] =
            +1.49647765e+00
            +2.86422420e-03 * tc[1]
            -6.13177463e-07 * tc[2]
            +6.83943370e-11 * tc[3]
            -3.08187834e-15 * tc[4]
            +7.16244227e+04 * invT;
        /*species 16: N2 */
        species[16] =
            +1.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 17: NO2 */
        species[17] =
            +3.88475420e+00
            +1.08619780e-03 * tc[1]
            -2.76023020e-07 * tc[2]
            +3.93687750e-11 * tc[3]
            -2.10217900e-15 * tc[4]
            +2.31649830e+03 * invT;
        /*species 18: E */
        species[18] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 19: CH2 */
        species[19] =
            +1.87410113e+00
            +1.82819646e-03 * tc[1]
            -4.69648657e-07 * tc[2]
            +6.50448872e-11 * tc[3]
            -3.75455134e-15 * tc[4]
            +4.62636040e+04 * invT;
    }

    /*species with midpoint at T=3654.18 kelvin */
    if (T < 3654.18) {
        /*species 6: HCOp */
        species[6] =
            +1.87627680e+00
            +2.47430100e-03 * tc[1]
            -7.69510200e-07 * tc[2]
            +1.27804380e-10 * tc[3]
            -8.70412620e-15 * tc[4]
            +9.91300220e+04 * invT;
    } else {
        /*species 6: HCOp */
        species[6] =
            +5.91572650e+00
            +8.31555350e-05 * tc[1]
            -6.15878100e-09 * tc[2]
            +2.27141072e-13 * tc[3]
            -3.29088820e-18 * tc[4]
            +9.61293330e+04 * invT;
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

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: CH3 */
        species[0] =
            +3.67359040e+00
            +1.00547588e-03 * tc[1]
            +1.91007285e-06 * tc[2]
            -1.71779356e-09 * tc[3]
            +5.08771468e-13 * tc[4]
            +1.64449988e+04 * invT;
        /*species 1: CH4 */
        species[1] =
            +5.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 * invT;
        /*species 2: H */
        species[2] =
            +2.50000000e+00
            +3.52666409e-13 * tc[1]
            -6.65306547e-16 * tc[2]
            +5.75204080e-19 * tc[3]
            -1.85546466e-22 * tc[4]
            +2.54736599e+04 * invT;
        /*species 3: OH */
        species[3] =
            +3.99201543e+00
            -1.20065876e-03 * tc[1]
            +1.53931280e-06 * tc[2]
            -9.70283332e-10 * tc[3]
            +2.72822940e-13 * tc[4]
            +3.61508056e+03 * invT;
        /*species 4: H2O */
        species[4] =
            +4.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
        /*species 5: O */
        species[5] =
            +3.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 * invT;
        /*species 7: NO */
        species[7] =
            +4.21847630e+00
            -2.31948800e-03 * tc[1]
            +3.68034067e-06 * tc[2]
            -2.33403385e-09 * tc[3]
            +5.60715400e-13 * tc[4]
            +9.84462300e+03 * invT;
        /*species 8: N */
        species[8] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +5.61046370e+04 * invT;
        /*species 9: H2 */
        species[9] =
            +2.34433112e+00
            +3.99026037e-03 * tc[1]
            -6.49271700e-06 * tc[2]
            +5.03930235e-09 * tc[3]
            -1.47522352e-12 * tc[4]
            -9.17935173e+02 * invT;
        /*species 10: O2 */
        species[10] =
            +3.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 11: CH2O */
        species[11] =
            +4.79372315e+00
            -4.95416684e-03 * tc[1]
            +1.24406669e-05 * tc[2]
            -9.48213152e-09 * tc[3]
            +2.63545304e-12 * tc[4]
            -1.43089567e+04 * invT;
        /*species 12: CO */
        species[12] =
            +3.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 13: CO2 */
        species[13] =
            +2.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
        /*species 14: CH */
        species[14] =
            +3.48981665e+00
            +1.61917771e-04 * tc[1]
            -5.62996883e-07 * tc[2]
            +7.90543317e-10 * tc[3]
            -2.81218134e-13 * tc[4]
            +7.07972934e+04 * invT;
        /*species 15: H3Op */
        species[15] =
            +3.79295251e+00
            -4.55426362e-04 * tc[1]
            +3.87878403e-06 * tc[2]
            -3.03412162e-09 * tc[3]
            +8.52319248e-13 * tc[4]
            +7.14027518e+04 * invT;
        /*species 16: N2 */
        species[16] =
            +3.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
        /*species 17: NO2 */
        species[17] =
            +3.94403120e+00
            -7.92714500e-04 * tc[1]
            +5.55260400e-06 * tc[2]
            -5.11885650e-09 * tc[3]
            +1.56701128e-12 * tc[4]
            +2.89661790e+03 * invT;
        /*species 18: E */
        species[18] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 19: CH2 */
        species[19] =
            +3.76267867e+00
            +4.84436072e-04 * tc[1]
            +9.31632803e-07 * tc[2]
            -9.62727883e-10 * tc[3]
            +3.37483438e-13 * tc[4]
            +4.60040401e+04 * invT;
    } else {
        /*species 0: CH3 */
        species[0] =
            +2.28571772e+00
            +3.61995018e-03 * tc[1]
            -9.95714493e-07 * tc[2]
            +1.48921161e-10 * tc[3]
            -9.34308788e-15 * tc[4]
            +1.67755843e+04 * invT;
        /*species 1: CH4 */
        species[1] =
            +7.48514950e-02
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 * invT;
        /*species 2: H */
        species[2] =
            +2.50000001e+00
            -1.15421486e-11 * tc[1]
            +5.38539827e-15 * tc[2]
            -1.18378809e-18 * tc[3]
            +9.96394714e-23 * tc[4]
            +2.54736599e+04 * invT;
        /*species 3: OH */
        species[3] =
            +3.09288767e+00
            +2.74214858e-04 * tc[1]
            +4.21684093e-08 * tc[2]
            -2.19865389e-11 * tc[3]
            +2.34824752e-15 * tc[4]
            +3.85865700e+03 * invT;
        /*species 4: H2O */
        species[4] =
            +3.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
        /*species 5: O */
        species[5] =
            +2.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 * invT;
        /*species 7: NO */
        species[7] =
            +3.26060560e+00
            +5.95552150e-04 * tc[1]
            -1.43056827e-07 * tc[2]
            +1.73644173e-11 * tc[3]
            -8.06721980e-16 * tc[4]
            +9.92097460e+03 * invT;
        /*species 8: N */
        species[8] =
            +2.41594290e+00
            +8.74453250e-05 * tc[1]
            -3.96745633e-08 * tc[2]
            +7.55656125e-12 * tc[3]
            -4.07219640e-16 * tc[4]
            +5.61337730e+04 * invT;
        /*species 9: H2 */
        species[9] =
            +3.33727920e+00
            -2.47012365e-05 * tc[1]
            +1.66485593e-07 * tc[2]
            -4.48915985e-11 * tc[3]
            +4.00510752e-15 * tc[4]
            -9.50158922e+02 * invT;
        /*species 10: O2 */
        species[10] =
            +3.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 11: CH2O */
        species[11] =
            +1.76069008e+00
            +4.60000041e-03 * tc[1]
            -1.47419604e-06 * tc[2]
            +2.51603030e-10 * tc[3]
            -1.76771128e-14 * tc[4]
            -1.39958323e+04 * invT;
        /*species 12: CO */
        species[12] =
            +2.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 13: CO2 */
        species[13] =
            +3.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
        /*species 14: CH */
        species[14] =
            +2.87846473e+00
            +4.85456840e-04 * tc[1]
            +4.81485517e-08 * tc[2]
            -3.26719623e-11 * tc[3]
            +3.52158766e-15 * tc[4]
            +7.10124364e+04 * invT;
        /*species 15: H3Op */
        species[15] =
            +2.49647765e+00
            +2.86422420e-03 * tc[1]
            -6.13177463e-07 * tc[2]
            +6.83943370e-11 * tc[3]
            -3.08187834e-15 * tc[4]
            +7.16244227e+04 * invT;
        /*species 16: N2 */
        species[16] =
            +2.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 17: NO2 */
        species[17] =
            +4.88475420e+00
            +1.08619780e-03 * tc[1]
            -2.76023020e-07 * tc[2]
            +3.93687750e-11 * tc[3]
            -2.10217900e-15 * tc[4]
            +2.31649830e+03 * invT;
        /*species 18: E */
        species[18] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 19: CH2 */
        species[19] =
            +2.87410113e+00
            +1.82819646e-03 * tc[1]
            -4.69648657e-07 * tc[2]
            +6.50448872e-11 * tc[3]
            -3.75455134e-15 * tc[4]
            +4.62636040e+04 * invT;
    }

    /*species with midpoint at T=3654.18 kelvin */
    if (T < 3654.18) {
        /*species 6: HCOp */
        species[6] =
            +2.87627680e+00
            +2.47430100e-03 * tc[1]
            -7.69510200e-07 * tc[2]
            +1.27804380e-10 * tc[3]
            -8.70412620e-15 * tc[4]
            +9.91300220e+04 * invT;
    } else {
        /*species 6: HCOp */
        species[6] =
            +6.91572650e+00
            +8.31555350e-05 * tc[1]
            -6.15878100e-09 * tc[2]
            +2.27141072e-13 * tc[3]
            -3.29088820e-18 * tc[4]
            +9.61293330e+04 * invT;
    }
    return;
}


/*compute the S/R at the given temperature (Eq 21) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesEntropy(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: CH3 */
        species[0] =
            +3.67359040e+00 * tc[0]
            +2.01095175e-03 * tc[1]
            +2.86510928e-06 * tc[2]
            -2.29039142e-09 * tc[3]
            +6.35964335e-13 * tc[4]
            +1.60456433e+00 ;
        /*species 1: CH4 */
        species[1] =
            +5.14987613e+00 * tc[0]
            -1.36709788e-02 * tc[1]
            +2.45900299e-05 * tc[2]
            -1.61581009e-08 * tc[3]
            +4.16734890e-12 * tc[4]
            -4.64130376e+00 ;
        /*species 2: H */
        species[2] =
            +2.50000000e+00 * tc[0]
            +7.05332819e-13 * tc[1]
            -9.97959820e-16 * tc[2]
            +7.66938773e-19 * tc[3]
            -2.31933083e-22 * tc[4]
            -4.46682853e-01 ;
        /*species 3: OH */
        species[3] =
            +3.99201543e+00 * tc[0]
            -2.40131752e-03 * tc[1]
            +2.30896920e-06 * tc[2]
            -1.29371111e-09 * tc[3]
            +3.41028675e-13 * tc[4]
            -1.03925458e-01 ;
        /*species 4: H2O */
        species[4] =
            +4.19864056e+00 * tc[0]
            -2.03643410e-03 * tc[1]
            +3.26020105e-06 * tc[2]
            -1.82932354e-09 * tc[3]
            +4.42994543e-13 * tc[4]
            -8.49032208e-01 ;
        /*species 5: O */
        species[5] =
            +3.16826710e+00 * tc[0]
            -3.27931884e-03 * tc[1]
            +3.32153198e-06 * tc[2]
            -2.04268875e-09 * tc[3]
            +5.28164927e-13 * tc[4]
            +2.05193346e+00 ;
        /*species 7: NO */
        species[7] =
            +4.21847630e+00 * tc[0]
            -4.63897600e-03 * tc[1]
            +5.52051100e-06 * tc[2]
            -3.11204513e-09 * tc[3]
            +7.00894250e-13 * tc[4]
            +2.28084640e+00 ;
        /*species 8: N */
        species[8] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +4.19390870e+00 ;
        /*species 9: H2 */
        species[9] =
            +2.34433112e+00 * tc[0]
            +7.98052075e-03 * tc[1]
            -9.73907550e-06 * tc[2]
            +6.71906980e-09 * tc[3]
            -1.84402940e-12 * tc[4]
            +6.83010238e-01 ;
        /*species 10: O2 */
        species[10] =
            +3.78245636e+00 * tc[0]
            -2.99673416e-03 * tc[1]
            +4.92365101e-06 * tc[2]
            -3.22709836e-09 * tc[3]
            +8.10932092e-13 * tc[4]
            +3.65767573e+00 ;
        /*species 11: CH2O */
        species[11] =
            +4.79372315e+00 * tc[0]
            -9.90833369e-03 * tc[1]
            +1.86610004e-05 * tc[2]
            -1.26428420e-08 * tc[3]
            +3.29431630e-12 * tc[4]
            +6.02812900e-01 ;
        /*species 12: CO */
        species[12] =
            +3.57953347e+00 * tc[0]
            -6.10353680e-04 * tc[1]
            +5.08407165e-07 * tc[2]
            +3.02335295e-10 * tc[3]
            -2.26106125e-13 * tc[4]
            +3.50840928e+00 ;
        /*species 13: CO2 */
        species[13] =
            +2.35677352e+00 * tc[0]
            +8.98459677e-03 * tc[1]
            -3.56178134e-06 * tc[2]
            +8.19730073e-10 * tc[3]
            -3.59248870e-14 * tc[4]
            +9.90105222e+00 ;
        /*species 14: CH */
        species[14] =
            +3.48981665e+00 * tc[0]
            +3.23835541e-04 * tc[1]
            -8.44495325e-07 * tc[2]
            +1.05405776e-09 * tc[3]
            -3.51522668e-13 * tc[4]
            +2.08401108e+00 ;
        /*species 15: H3Op */
        species[15] =
            +3.79295251e+00 * tc[0]
            -9.10852723e-04 * tc[1]
            +5.81817605e-06 * tc[2]
            -4.04549550e-09 * tc[3]
            +1.06539906e-12 * tc[4]
            +1.47156927e+00 ;
        /*species 16: N2 */
        species[16] =
            +3.29867700e+00 * tc[0]
            +1.40824040e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213500e-13 * tc[4]
            +3.95037200e+00 ;
        /*species 17: NO2 */
        species[17] =
            +3.94403120e+00 * tc[0]
            -1.58542900e-03 * tc[1]
            +8.32890600e-06 * tc[2]
            -6.82514200e-09 * tc[3]
            +1.95876410e-12 * tc[4]
            +6.31199170e+00 ;
        /*species 18: E */
        species[18] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -1.17208122e+01 ;
        /*species 19: CH2 */
        species[19] =
            +3.76267867e+00 * tc[0]
            +9.68872143e-04 * tc[1]
            +1.39744921e-06 * tc[2]
            -1.28363718e-09 * tc[3]
            +4.21854298e-13 * tc[4]
            +1.56253185e+00 ;
    } else {
        /*species 0: CH3 */
        species[0] =
            +2.28571772e+00 * tc[0]
            +7.23990037e-03 * tc[1]
            -1.49357174e-06 * tc[2]
            +1.98561548e-10 * tc[3]
            -1.16788599e-14 * tc[4]
            +8.48007179e+00 ;
        /*species 1: CH4 */
        species[1] =
            +7.48514950e-02 * tc[0]
            +1.33909467e-02 * tc[1]
            -2.86642905e-06 * tc[2]
            +4.07641783e-10 * tc[3]
            -2.54538075e-14 * tc[4]
            +1.84373180e+01 ;
        /*species 2: H */
        species[2] =
            +2.50000001e+00 * tc[0]
            -2.30842973e-11 * tc[1]
            +8.07809740e-15 * tc[2]
            -1.57838412e-18 * tc[3]
            +1.24549339e-22 * tc[4]
            -4.46682914e-01 ;
        /*species 3: OH */
        species[3] =
            +3.09288767e+00 * tc[0]
            +5.48429716e-04 * tc[1]
            +6.32526140e-08 * tc[2]
            -2.93153852e-11 * tc[3]
            +2.93530940e-15 * tc[4]
            +4.47669610e+00 ;
        /*species 4: H2O */
        species[4] =
            +3.03399249e+00 * tc[0]
            +2.17691804e-03 * tc[1]
            -8.20362590e-08 * tc[2]
            -3.23473290e-11 * tc[3]
            +4.20502480e-15 * tc[4]
            +4.96677010e+00 ;
        /*species 5: O */
        species[5] =
            +2.56942078e+00 * tc[0]
            -8.59741137e-05 * tc[1]
            +2.09742295e-08 * tc[2]
            -3.33925997e-12 * tc[3]
            +3.07084227e-16 * tc[4]
            +4.78433864e+00 ;
        /*species 7: NO */
        species[7] =
            +3.26060560e+00 * tc[0]
            +1.19110430e-03 * tc[1]
            -2.14585240e-07 * tc[2]
            +2.31525563e-11 * tc[3]
            -1.00840247e-15 * tc[4]
            +6.36930270e+00 ;
        /*species 8: N */
        species[8] =
            +2.41594290e+00 * tc[0]
            +1.74890650e-04 * tc[1]
            -5.95118450e-08 * tc[2]
            +1.00754150e-11 * tc[3]
            -5.09024550e-16 * tc[4]
            +4.64960960e+00 ;
        /*species 9: H2 */
        species[9] =
            +3.33727920e+00 * tc[0]
            -4.94024731e-05 * tc[1]
            +2.49728389e-07 * tc[2]
            -5.98554647e-11 * tc[3]
            +5.00638440e-15 * tc[4]
            -3.20502331e+00 ;
        /*species 10: O2 */
        species[10] =
            +3.28253784e+00 * tc[0]
            +1.48308754e-03 * tc[1]
            -3.78983334e-07 * tc[2]
            +6.98235183e-11 * tc[3]
            -5.41794485e-15 * tc[4]
            +5.45323129e+00 ;
        /*species 11: CH2O */
        species[11] =
            +1.76069008e+00 * tc[0]
            +9.20000082e-03 * tc[1]
            -2.21129406e-06 * tc[2]
            +3.35470707e-10 * tc[3]
            -2.20963910e-14 * tc[4]
            +1.36563230e+01 ;
        /*species 12: CO */
        species[12] =
            +2.71518561e+00 * tc[0]
            +2.06252743e-03 * tc[1]
            -4.99412886e-07 * tc[2]
            +7.66843360e-11 * tc[3]
            -5.09119290e-15 * tc[4]
            +7.81868772e+00 ;
        /*species 13: CO2 */
        species[13] =
            +3.85746029e+00 * tc[0]
            +4.41437026e-03 * tc[1]
            -1.10740702e-06 * tc[2]
            +1.74496729e-10 * tc[3]
            -1.18021041e-14 * tc[4]
            +2.27163806e+00 ;
        /*species 14: CH */
        species[14] =
            +2.87846473e+00 * tc[0]
            +9.70913681e-04 * tc[1]
            +7.22228275e-08 * tc[2]
            -4.35626163e-11 * tc[3]
            +4.40198457e-15 * tc[4]
            +5.48497999e+00 ;
        /*species 15: H3Op */
        species[15] =
            +2.49647765e+00 * tc[0]
            +5.72844840e-03 * tc[1]
            -9.19766195e-07 * tc[2]
            +9.11924493e-11 * tc[3]
            -3.85234792e-15 * tc[4]
            +7.45850493e+00 ;
        /*species 16: N2 */
        species[16] =
            +2.92664000e+00 * tc[0]
            +1.48797680e-03 * tc[1]
            -2.84238000e-07 * tc[2]
            +3.36567933e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
        /*species 17: NO2 */
        species[17] =
            +4.88475420e+00 * tc[0]
            +2.17239560e-03 * tc[1]
            -4.14034530e-07 * tc[2]
            +5.24917000e-11 * tc[3]
            -2.62772375e-15 * tc[4]
            -1.17416950e-01 ;
        /*species 18: E */
        species[18] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -1.17208122e+01 ;
        /*species 19: CH2 */
        species[19] =
            +2.87410113e+00 * tc[0]
            +3.65639292e-03 * tc[1]
            -7.04472985e-07 * tc[2]
            +8.67265163e-11 * tc[3]
            -4.69318918e-15 * tc[4]
            +6.17119324e+00 ;
    }

    /*species with midpoint at T=3654.18 kelvin */
    if (T < 3654.18) {
        /*species 6: HCOp */
        species[6] =
            +2.87627680e+00 * tc[0]
            +4.94860200e-03 * tc[1]
            -1.15426530e-06 * tc[2]
            +1.70405840e-10 * tc[3]
            -1.08801578e-14 * tc[4]
            +6.50202660e+00 ;
    } else {
        /*species 6: HCOp */
        species[6] =
            +6.91572650e+00 * tc[0]
            +1.66311070e-04 * tc[1]
            -9.23817150e-09 * tc[2]
            +3.02854763e-13 * tc[3]
            -4.11361025e-18 * tc[4]
            -1.80895100e+01 ;
    }
    return;
}


/*save molecular weights into array */
void molecularWeight(double * restrict wt)
{
    wt[0] = 15.035060; /*CH3 */
    wt[1] = 16.043030; /*CH4 */
    wt[2] = 1.007970; /*H */
    wt[3] = 17.007370; /*OH */
    wt[4] = 18.015340; /*H2O */
    wt[5] = 15.999400; /*O */
    wt[6] = 29.017971; /*HCOp */
    wt[7] = 30.006100; /*NO */
    wt[8] = 14.006700; /*N */
    wt[9] = 2.015940; /*H2 */
    wt[10] = 31.998800; /*O2 */
    wt[11] = 30.026490; /*CH2O */
    wt[12] = 28.010550; /*CO */
    wt[13] = 44.009950; /*CO2 */
    wt[14] = 13.019120; /*CH */
    wt[15] = 19.022761; /*H3Op */
    wt[16] = 28.013400; /*N2 */
    wt[17] = 46.005500; /*NO2 */
    wt[18] = 0.000549; /*E */
    wt[19] = 14.027090; /*CH2 */

    return;
}


/*save atomic weights into array */
void atomicWeight(double * restrict awt)
{
    awt[0] = 15.999400; /*O */
    awt[1] = 1.007970; /*H */
    awt[2] = 12.011150; /*C */
    awt[3] = 14.006700; /*N */
    awt[4] = 0.000549; /*E */

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
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetLENIMC EGTRANSETLENIMC
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetLENIMC egtransetlenimc
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetLENIMC egtransetlenimc_
#endif
void egtransetLENIMC(int* LENIMC) {
  *LENIMC =           83;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetLENRMC EGTRANSETLENRMC
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetLENRMC egtransetlenrmc
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetLENRMC egtransetlenrmc_
#endif
void egtransetLENRMC(int* LENRMC) {
  *LENRMC =         8380;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNO EGTRANSETNO
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNO egtransetno
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNO egtransetno_
#endif
void egtransetNO(int* NO) {
  *NO =            4;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetKK EGTRANSETKK
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetKK egtransetkk
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetKK egtransetkk_
#endif
void egtransetKK(int* KK) {
  *KK =           20;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNLITE EGTRANSETNLITE
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNLITE egtransetnlite
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNLITE egtransetnlite_
#endif
void egtransetNLITE(int* NLITE) {
  *NLITE =            3;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetPATM EGTRANSETPATM
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetPATM egtransetpatm
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetPATM egtransetpatm_
#endif
void egtransetPATM(double* PATM) {
  *PATM =   0.1013250000000000E+07;}
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetWT EGTRANSETWT
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetWT egtransetwt
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetWT egtransetwt_
#endif
void egtransetWT(double* WT) {
  WT[           0] =   0.1503506028652191E+02;
  WT[           1] =   0.1604303026199341E+02;
  WT[           2] =   0.1007969975471497E+01;
  WT[           3] =   0.1700737011432648E+02;
  WT[           4] =   0.1801534008979797E+02;
  WT[           5] =   0.1599940013885498E+02;
  WT[           6] =   0.2901797189645004E+02;
  WT[           7] =   0.3000609970092773E+02;
  WT[           8] =   0.1400669956207275E+02;
  WT[           9] =   0.2015939950942993E+01;
  WT[          10] =   0.3199880027770996E+02;
  WT[          11] =   0.3002649044990540E+02;
  WT[          12] =   0.2801055049896240E+02;
  WT[          13] =   0.4400995063781738E+02;
  WT[          14] =   0.1301912033557892E+02;
  WT[          15] =   0.1902276148728561E+02;
  WT[          16] =   0.2801339912414551E+02;
  WT[          17] =   0.4600549983978271E+02;
  WT[          18] =   0.5485779838636518E-03;
  WT[          19] =   0.1402709031105042E+02;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetEPS EGTRANSETEPS
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetEPS egtranseteps
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetEPS egtranseteps_
#endif
void egtransetEPS(double* EPS) {
  EPS[           0] =   0.1440000000000000E+03;
  EPS[           1] =   0.1414000000000000E+03;
  EPS[           2] =   0.1450000000000000E+03;
  EPS[           3] =   0.8000000000000000E+02;
  EPS[           4] =   0.5724000000000000E+03;
  EPS[           5] =   0.8000000000000000E+02;
  EPS[           6] =   0.4980000000000000E+03;
  EPS[           7] =   0.9753000000000000E+02;
  EPS[           8] =   0.7140000000000001E+02;
  EPS[           9] =   0.3800000000000000E+02;
  EPS[          10] =   0.1074000000000000E+03;
  EPS[          11] =   0.4980000000000000E+03;
  EPS[          12] =   0.9809999999999999E+02;
  EPS[          13] =   0.2440000000000000E+03;
  EPS[          14] =   0.8000000000000000E+02;
  EPS[          15] =   0.5724000000000000E+03;
  EPS[          16] =   0.9753000000000000E+02;
  EPS[          17] =   0.2000000000000000E+03;
  EPS[          18] =   0.8500000000000000E+03;
  EPS[          19] =   0.1440000000000000E+03;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetSIG EGTRANSETSIG
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetSIG egtransetsig
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetSIG egtransetsig_
#endif
void egtransetSIG(double* SIG) {
  SIG[           0] =   0.3800000000000000E+01;
  SIG[           1] =   0.3746000000000000E+01;
  SIG[           2] =   0.2050000000000000E+01;
  SIG[           3] =   0.2750000000000000E+01;
  SIG[           4] =   0.2605000000000000E+01;
  SIG[           5] =   0.2750000000000000E+01;
  SIG[           6] =   0.3590000000000000E+01;
  SIG[           7] =   0.3621000000000000E+01;
  SIG[           8] =   0.3298000000000000E+01;
  SIG[           9] =   0.2920000000000000E+01;
  SIG[          10] =   0.3458000000000000E+01;
  SIG[          11] =   0.3590000000000000E+01;
  SIG[          12] =   0.3650000000000000E+01;
  SIG[          13] =   0.3763000000000000E+01;
  SIG[          14] =   0.2750000000000000E+01;
  SIG[          15] =   0.2605000000000000E+01;
  SIG[          16] =   0.3621000000000000E+01;
  SIG[          17] =   0.3500000000000000E+01;
  SIG[          18] =   0.4250000000000000E+03;
  SIG[          19] =   0.3800000000000000E+01;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetDIP EGTRANSETDIP
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetDIP egtransetdip
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetDIP egtransetdip_
#endif
void egtransetDIP(double* DIP) {
  DIP[           0] =   0.0000000000000000E+00;
  DIP[           1] =   0.0000000000000000E+00;
  DIP[           2] =   0.0000000000000000E+00;
  DIP[           3] =   0.0000000000000000E+00;
  DIP[           4] =   0.1844000000000000E+01;
  DIP[           5] =   0.0000000000000000E+00;
  DIP[           6] =   0.0000000000000000E+00;
  DIP[           7] =   0.0000000000000000E+00;
  DIP[           8] =   0.0000000000000000E+00;
  DIP[           9] =   0.0000000000000000E+00;
  DIP[          10] =   0.0000000000000000E+00;
  DIP[          11] =   0.0000000000000000E+00;
  DIP[          12] =   0.0000000000000000E+00;
  DIP[          13] =   0.0000000000000000E+00;
  DIP[          14] =   0.0000000000000000E+00;
  DIP[          15] =   0.1844000000000000E+01;
  DIP[          16] =   0.0000000000000000E+00;
  DIP[          17] =   0.0000000000000000E+00;
  DIP[          18] =   0.0000000000000000E+00;
  DIP[          19] =   0.0000000000000000E+00;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetPOL EGTRANSETPOL
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetPOL egtransetpol
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetPOL egtransetpol_
#endif
void egtransetPOL(double* POL) {
  POL[           0] =   0.0000000000000000E+00;
  POL[           1] =   0.2600000000000000E+01;
  POL[           2] =   0.0000000000000000E+00;
  POL[           3] =   0.0000000000000000E+00;
  POL[           4] =   0.0000000000000000E+00;
  POL[           5] =   0.0000000000000000E+00;
  POL[           6] =   0.0000000000000000E+00;
  POL[           7] =   0.1760000000000000E+01;
  POL[           8] =   0.0000000000000000E+00;
  POL[           9] =   0.7900000000000000E+00;
  POL[          10] =   0.1600000000000000E+01;
  POL[          11] =   0.0000000000000000E+00;
  POL[          12] =   0.1950000000000000E+01;
  POL[          13] =   0.2650000000000000E+01;
  POL[          14] =   0.0000000000000000E+00;
  POL[          15] =   0.0000000000000000E+00;
  POL[          16] =   0.1760000000000000E+01;
  POL[          17] =   0.0000000000000000E+00;
  POL[          18] =   0.0000000000000000E+00;
  POL[          19] =   0.0000000000000000E+00;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetZROT EGTRANSETZROT
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetZROT egtransetzrot
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetZROT egtransetzrot_
#endif
void egtransetZROT(double* ZROT) {
  ZROT[           0] =   0.0000000000000000E+00;
  ZROT[           1] =   0.1300000000000000E+02;
  ZROT[           2] =   0.0000000000000000E+00;
  ZROT[           3] =   0.0000000000000000E+00;
  ZROT[           4] =   0.4000000000000000E+01;
  ZROT[           5] =   0.0000000000000000E+00;
  ZROT[           6] =   0.0000000000000000E+00;
  ZROT[           7] =   0.4000000000000000E+01;
  ZROT[           8] =   0.0000000000000000E+00;
  ZROT[           9] =   0.2800000000000000E+03;
  ZROT[          10] =   0.3800000000000000E+01;
  ZROT[          11] =   0.2000000000000000E+01;
  ZROT[          12] =   0.1800000000000000E+01;
  ZROT[          13] =   0.2100000000000000E+01;
  ZROT[          14] =   0.0000000000000000E+00;
  ZROT[          15] =   0.4000000000000000E+01;
  ZROT[          16] =   0.4000000000000000E+01;
  ZROT[          17] =   0.1000000000000000E+01;
  ZROT[          18] =   0.1000000000000000E+01;
  ZROT[          19] =   0.0000000000000000E+00;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNLIN EGTRANSETNLIN
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNLIN egtransetnlin
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNLIN egtransetnlin_
#endif
void egtransetNLIN(int* NLIN) {
  NLIN[           0] =            1;
  NLIN[           1] =            2;
  NLIN[           2] =            0;
  NLIN[           3] =            1;
  NLIN[           4] =            2;
  NLIN[           5] =            0;
  NLIN[           6] =            2;
  NLIN[           7] =            1;
  NLIN[           8] =            0;
  NLIN[           9] =            1;
  NLIN[          10] =            1;
  NLIN[          11] =            2;
  NLIN[          12] =            1;
  NLIN[          13] =            1;
  NLIN[          14] =            1;
  NLIN[          15] =            2;
  NLIN[          16] =            1;
  NLIN[          17] =            2;
  NLIN[          18] =            0;
  NLIN[          19] =            1;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFLAM EGTRANSETCOFLAM
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFLAM egtransetcoflam
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFLAM egtransetcoflam_
#endif
void egtransetCOFLAM(double* COFLAM) {
  COFLAM[           0] =   0.1208656985380625E+02;
  COFLAM[           1] =  -0.3791666220447708E+01;
  COFLAM[           2] =   0.7819100491045913E+00;
  COFLAM[           3] =  -0.4156129826797447E-01;
  COFLAM[           4] =   0.9559526846070922E+01;
  COFLAM[           5] =  -0.3251648121323282E+01;
  COFLAM[           6] =   0.7750954738197482E+00;
  COFLAM[           7] =  -0.4347541297630472E-01;
  COFLAM[           8] =  -0.3766745631253565E+00;
  COFLAM[           9] =   0.3437535150998352E+01;
  COFLAM[          10] =  -0.3660041411242859E+00;
  COFLAM[          11] =   0.1599148758178118E-01;
  COFLAM[          12] =   0.1475874620420910E+02;
  COFLAM[          13] =  -0.3531252599661289E+01;
  COFLAM[          14] =   0.5788603821266735E+00;
  COFLAM[          15] =  -0.2562509217392182E-01;
  COFLAM[          16] =   0.2251607076002226E+02;
  COFLAM[          17] =  -0.8597061155164372E+01;
  COFLAM[          18] =   0.1475703352610380E+01;
  COFLAM[          19] =  -0.7336194664755748E-01;
  COFLAM[          20] =   0.1901192499833660E+01;
  COFLAM[          21] =   0.1830277475782733E+01;
  COFLAM[          22] =  -0.1589498673817636E+00;
  COFLAM[          23] =   0.7096188225056049E-02;
  COFLAM[          24] =  -0.1717011951962953E+00;
  COFLAM[          25] =   0.6146277139624015E+00;
  COFLAM[          26] =   0.2241740897542817E+00;
  COFLAM[          27] =  -0.1810845630249059E-01;
  COFLAM[          28] =   0.8432249653115694E+01;
  COFLAM[          29] =  -0.1618954562122087E+01;
  COFLAM[          30] =   0.3778158534240268E+00;
  COFLAM[          31] =  -0.1955692806758166E-01;
  COFLAM[          32] =   0.2221385418645368E+01;
  COFLAM[          33] =   0.1590510240655957E+01;
  COFLAM[          34] =  -0.1270537390077182E+00;
  COFLAM[          35] =   0.5682425577944382E-02;
  COFLAM[          36] =   0.1083144193584373E+02;
  COFLAM[          37] =  -0.1205052794384021E+01;
  COFLAM[          38] =   0.2282713175734062E+00;
  COFLAM[          39] =  -0.8274911289650255E-02;
  COFLAM[          40] =  -0.2588874986012296E+01;
  COFLAM[          41] =   0.3196698444251531E+01;
  COFLAM[          42] =  -0.3183369839404002E+00;
  COFLAM[          43] =   0.1394344715230612E-01;
  COFLAM[          44] =   0.1489363739789972E+01;
  COFLAM[          45] =  -0.6039120155555791E+00;
  COFLAM[          46] =   0.4684882817936775E+00;
  COFLAM[          47] =  -0.3225988482135683E-01;
  COFLAM[          48] =   0.9900423525072517E+01;
  COFLAM[          49] =  -0.2273352721794634E+01;
  COFLAM[          50] =   0.4717809389728932E+00;
  COFLAM[          51] =  -0.2394578139317899E-01;
  COFLAM[          52] =  -0.1211717236713152E+02;
  COFLAM[          53] =   0.6227082396229542E+01;
  COFLAM[          54] =  -0.6210515389533838E+00;
  COFLAM[          55] =   0.2298522004617450E-01;
  COFLAM[          56] =   0.2161890379462439E+02;
  COFLAM[          57] =  -0.6630557890027998E+01;
  COFLAM[          58] =   0.1043654141215424E+01;
  COFLAM[          59] =  -0.4821124906708143E-01;
  COFLAM[          60] =   0.1864807242121871E+02;
  COFLAM[          61] =  -0.7525097955112172E+01;
  COFLAM[          62] =   0.1409276672416255E+01;
  COFLAM[          63] =  -0.7373222670706672E-01;
  COFLAM[          64] =   0.1131914099090367E+02;
  COFLAM[          65] =  -0.2817308169021882E+01;
  COFLAM[          66] =   0.5415014032317143E+00;
  COFLAM[          67] =  -0.2689153797571326E-01;
  COFLAM[          68] =  -0.1450044718997646E+02;
  COFLAM[          69] =   0.7556916313990506E+01;
  COFLAM[          70] =  -0.8438490520213084E+00;
  COFLAM[          71] =   0.3495614481844105E-01;
  COFLAM[          72] =   0.5585161074259044E+01;
  COFLAM[          73] =  -0.3194509231847466E+01;
  COFLAM[          74] =   0.6682611688624166E+00;
  COFLAM[          75] =  -0.3538987109886362E-01;
  COFLAM[          76] =   0.1205209988598598E+02;
  COFLAM[          77] =  -0.3362497748037753E+01;
  COFLAM[          78] =   0.6619316650714779E+00;
  COFLAM[          79] =  -0.3381355791959224E-01;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFETA EGTRANSETCOFETA
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFETA egtransetcofeta
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFETA egtransetcofeta_
#endif
void egtransetCOFETA(double* COFETA) {
  COFETA[           0] =  -0.1975834070168401E+02;
  COFETA[           1] =   0.3419070933965452E+01;
  COFETA[           2] =  -0.3637502184371358E+00;
  COFETA[           3] =   0.1589993485342581E-01;
  COFETA[           4] =  -0.1957393167710191E+02;
  COFETA[           5] =   0.3375443809041139E+01;
  COFETA[           6] =  -0.3585148910195592E+00;
  COFETA[           7] =   0.1569150857325605E-01;
  COFETA[           8] =  -0.1992658985236223E+02;
  COFETA[           9] =   0.3437535150998333E+01;
  COFETA[          10] =  -0.3660041411242830E+00;
  COFETA[          11] =   0.1599148758178105E-01;
  COFETA[          12] =  -0.1485356216414066E+02;
  COFETA[          13] =   0.1830277475782665E+01;
  COFETA[          14] =  -0.1589498673817540E+00;
  COFETA[          15] =   0.7096188225055598E-02;
  COFETA[          16] =  -0.1167680175985407E+02;
  COFETA[          17] =  -0.8724679593380245E+00;
  COFETA[          18] =   0.3458379658811396E+00;
  COFETA[          19] =  -0.2040291235495196E-01;
  COFETA[          20] =  -0.1488410994215470E+02;
  COFETA[          21] =   0.1830277475782681E+01;
  COFETA[          22] =  -0.1589498673817564E+00;
  COFETA[          23] =   0.7096188225055718E-02;
  COFETA[          24] =  -0.2072209556320896E+02;
  COFETA[          25] =   0.3098279967435479E+01;
  COFETA[          26] =  -0.2275331939119109E+00;
  COFETA[          27] =   0.6257207134157165E-02;
  COFETA[          28] =  -0.1632308560239169E+02;
  COFETA[          29] =   0.2292176009378712E+01;
  COFETA[          30] =  -0.2194972386392208E+00;
  COFETA[          31] =   0.9740858903808247E-02;
  COFETA[          32] =  -0.1469693249825934E+02;
  COFETA[          33] =   0.1590510240656090E+01;
  COFETA[          34] =  -0.1270537390077376E+00;
  COFETA[          35] =   0.5682425577945314E-02;
  COFETA[          36] =  -0.1378707509843740E+02;
  COFETA[          37] =   0.9792780569374036E+00;
  COFETA[          38] =  -0.4652432169221873E-01;
  COFETA[          39] =   0.2150486800510954E-02;
  COFETA[          40] =  -0.1690583585391785E+02;
  COFETA[          41] =   0.2562428627672089E+01;
  COFETA[          42] =  -0.2546284432807669E+00;
  COFETA[          43] =   0.1126298284349495E-01;
  COFETA[          44] =  -0.2070501323899667E+02;
  COFETA[          45] =   0.3098279967435433E+01;
  COFETA[          46] =  -0.2275331939119042E+00;
  COFETA[          47] =   0.6257207134156837E-02;
  COFETA[          48] =  -0.1641639830677721E+02;
  COFETA[          49] =   0.2308637683893279E+01;
  COFETA[          50] =  -0.2216421664642076E+00;
  COFETA[          51] =   0.9833984530793158E-02;
  COFETA[          52] =  -0.2363833542953974E+02;
  COFETA[          53] =   0.4982494306355357E+01;
  COFETA[          54] =  -0.5505909825468872E+00;
  COFETA[          55] =   0.2334007366151606E-01;
  COFETA[          56] =  -0.1498717602128921E+02;
  COFETA[          57] =   0.1830277475782689E+01;
  COFETA[          58] =  -0.1589498673817574E+00;
  COFETA[          59] =   0.7096188225055758E-02;
  COFETA[          60] =  -0.1164959545334551E+02;
  COFETA[          61] =  -0.8724679593381346E+00;
  COFETA[          62] =   0.3458379658811553E+00;
  COFETA[          63] =  -0.2040291235495271E-01;
  COFETA[          64] =  -0.1635744447664137E+02;
  COFETA[          65] =   0.2292176009378658E+01;
  COFETA[          66] =  -0.2194972386392133E+00;
  COFETA[          67] =   0.9740858903807891E-02;
  COFETA[          68] =  -0.2198022474125375E+02;
  COFETA[          69] =   0.4488526152631485E+01;
  COFETA[          70] =  -0.4955170006329650E+00;
  COFETA[          71] =   0.2131963014236630E-01;
  COFETA[          72] =  -0.2148087370970207E+02;
  COFETA[          73] =  -0.3194509231847701E+01;
  COFETA[          74] =   0.6682611688624503E+00;
  COFETA[          75] =  -0.3538987109886522E-01;
  COFETA[          76] =  -0.1979303787372817E+02;
  COFETA[          77] =   0.3419070933965388E+01;
  COFETA[          78] =  -0.3637502184371259E+00;
  COFETA[          79] =   0.1589993485342531E-01;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFD EGTRANSETCOFD
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFD egtransetcofd
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFD egtransetcofd_
#endif
void egtransetCOFD(double* COFD) {
  COFD[           0] =  -0.1731477887424627E+02;
  COFD[           1] =   0.4182999326187860E+01;
  COFD[           2] =  -0.3264888667764371E+00;
  COFD[           3] =   0.1406581894032641E-01;
  COFD[           4] =  -0.1726709159699041E+02;
  COFD[           5] =   0.4166249447916941E+01;
  COFD[           6] =  -0.3245523257185624E+00;
  COFD[           7] =   0.1399183491478787E-01;
  COFD[           8] =  -0.1752194960164599E+02;
  COFD[           9] =   0.4875287931841000E+01;
  COFD[          10] =  -0.4147968034669867E+00;
  COFD[          11] =   0.1781089573575908E-01;
  COFD[          12] =  -0.1543681755468930E+02;
  COFD[          13] =   0.3615102026808611E+01;
  COFD[          14] =  -0.2578306288187833E+00;
  COFD[          15] =   0.1130180672438872E-01;
  COFD[          16] =  -0.2027289153258947E+02;
  COFD[          17] =   0.5157651529525397E+01;
  COFD[          18] =  -0.4233239455087539E+00;
  COFD[          19] =   0.1708864601311192E-01;
  COFD[          20] =  -0.1543125672017102E+02;
  COFD[          21] =   0.3618855115638650E+01;
  COFD[          22] =  -0.2583419121546929E+00;
  COFD[          23] =   0.1132494930338730E-01;
  COFD[          24] =  -0.2048209697647840E+02;
  COFD[          25] =   0.5127870094995610E+01;
  COFD[          26] =  -0.4231024044426618E+00;
  COFD[          27] =   0.1721759814463075E-01;
  COFD[          28] =  -0.1645223742454043E+02;
  COFD[          29] =   0.3854697730165210E+01;
  COFD[          30] =  -0.2878037882685343E+00;
  COFD[          31] =   0.1255204150178310E-01;
  COFD[          32] =  -0.1514651789245071E+02;
  COFD[          33] =   0.3462312981288904E+01;
  COFD[          34] =  -0.2380726416670257E+00;
  COFD[          35] =   0.1044895488982180E-01;
  COFD[          36] =  -0.1322163692827858E+02;
  COFD[          37] =   0.3075522930065609E+01;
  COFD[          38] =  -0.1904123874279204E+00;
  COFD[          39] =   0.8492834273761591E-02;
  COFD[          40] =  -0.1672933652118024E+02;
  COFD[          41] =   0.3965785816566581E+01;
  COFD[          42] =  -0.3013281840320347E+00;
  COFD[          43] =   0.1310001478947402E-01;
  COFD[          44] =  -0.2047207935890673E+02;
  COFD[          45] =   0.5120549000988871E+01;
  COFD[          46] =  -0.4219907391201317E+00;
  COFD[          47] =   0.1716268084967035E-01;
  COFD[          48] =  -0.1645180642806876E+02;
  COFD[          49] =   0.3855410394996957E+01;
  COFD[          50] =  -0.2878581523898601E+00;
  COFD[          51] =   0.1255291913453395E-01;
  COFD[          52] =  -0.1923215106925674E+02;
  COFD[          53] =   0.4769770645950153E+01;
  COFD[          54] =  -0.3939017067778038E+00;
  COFD[          55] =   0.1662898220909356E-01;
  COFD[          56] =  -0.1543306786434106E+02;
  COFD[          57] =   0.3641484216742679E+01;
  COFD[          58] =  -0.2614192836082300E+00;
  COFD[          59] =   0.1146400341320006E-01;
  COFD[          60] =  -0.2028829562643207E+02;
  COFD[          61] =   0.5159058519743692E+01;
  COFD[          62] =  -0.4235262462637739E+00;
  COFD[          63] =   0.1709816078216507E-01;
  COFD[          64] =  -0.1641389233309519E+02;
  COFD[          65] =   0.3843808661762178E+01;
  COFD[          66] =  -0.2863465660410651E+00;
  COFD[          67] =   0.1248720098494784E-01;
  COFD[          68] =  -0.1868821899264003E+02;
  COFD[          69] =   0.4623091902441412E+01;
  COFD[          70] =  -0.3789790809161654E+00;
  COFD[          71] =   0.1615045683658846E-01;
  COFD[          72] =  -0.2380880735019050E+02;
  COFD[          73] =   0.4985872967066852E+01;
  COFD[          74] =  -0.3856996034210413E+00;
  COFD[          75] =   0.1488146671087461E-01;
  COFD[          76] =  -0.1729871483444523E+02;
  COFD[          77] =   0.4183622454402964E+01;
  COFD[          78] =  -0.3265686107794125E+00;
  COFD[          79] =   0.1406921012900881E-01;
  COFD[          80] =  -0.1726709159699041E+02;
  COFD[          81] =   0.4166249447916941E+01;
  COFD[          82] =  -0.3245523257185624E+00;
  COFD[          83] =   0.1399183491478787E-01;
  COFD[          84] =  -0.1722618348214878E+02;
  COFD[          85] =   0.4152147706270270E+01;
  COFD[          86] =  -0.3229851063924484E+00;
  COFD[          87] =   0.1393497521609246E-01;
  COFD[          88] =  -0.1742491429482865E+02;
  COFD[          89] =   0.4846655191576849E+01;
  COFD[          90] =  -0.4113732013955531E+00;
  COFD[          91] =   0.1767520282723403E-01;
  COFD[          92] =  -0.1538242046452401E+02;
  COFD[          93] =   0.3594755693746865E+01;
  COFD[          94] =  -0.2552359406239311E+00;
  COFD[          95] =   0.1119129111015014E-01;
  COFD[          96] =  -0.2046809044123080E+02;
  COFD[          97] =   0.5168356480779646E+01;
  COFD[          98] =  -0.4170015852313753E+00;
  COFD[          99] =   0.1651418922298926E-01;
  COFD[         100] =  -0.1537895325066121E+02;
  COFD[         101] =   0.3599582752195144E+01;
  COFD[         102] =  -0.2558927079381389E+00;
  COFD[         103] =   0.1122098249033523E-01;
  COFD[         104] =  -0.2046876559898188E+02;
  COFD[         105] =   0.5125566020964355E+01;
  COFD[         106] =  -0.4233467066361414E+00;
  COFD[         107] =   0.1724767006424224E-01;
  COFD[         108] =  -0.1635833918443450E+02;
  COFD[         109] =   0.3815691616616636E+01;
  COFD[         110] =  -0.2827061894445888E+00;
  COFD[         111] =   0.1232981052514379E-01;
  COFD[         112] =  -0.1508728524143146E+02;
  COFD[         113] =   0.3439843254137499E+01;
  COFD[         114] =  -0.2351547205941974E+00;
  COFD[         115] =   0.1032244834588984E-01;
  COFD[         116] =  -0.1315730506576576E+02;
  COFD[         117] =   0.3055295555288239E+01;
  COFD[         118] =  -0.1877084814632856E+00;
  COFD[         119] =   0.8372316362831445E-02;
  COFD[         120] =  -0.1677759619477616E+02;
  COFD[         121] =   0.3988901785374923E+01;
  COFD[         122] =  -0.3051808273299473E+00;
  COFD[         123] =   0.1330327276290500E-01;
  COFD[         124] =  -0.2045967452432749E+02;
  COFD[         125] =   0.5118525449350731E+01;
  COFD[         126] =  -0.4222729551558472E+00;
  COFD[         127] =   0.1719441142995254E-01;
  COFD[         128] =  -0.1635709831001740E+02;
  COFD[         129] =   0.3816293164170347E+01;
  COFD[         130] =  -0.2827489538440444E+00;
  COFD[         131] =   0.1233030225007116E-01;
  COFD[         132] =  -0.1917783133928890E+02;
  COFD[         133] =   0.4748871842879942E+01;
  COFD[         134] =  -0.3916351371001520E+00;
  COFD[         135] =   0.1654935017103951E-01;
  COFD[         136] =  -0.1538700627222775E+02;
  COFD[         137] =   0.3625494228341230E+01;
  COFD[         138] =  -0.2594144135440631E+00;
  COFD[         139] =   0.1138002455729723E-01;
  COFD[         140] =  -0.2048707524059223E+02;
  COFD[         141] =   0.5171197873032807E+01;
  COFD[         142] =  -0.4174185529276651E+00;
  COFD[         143] =   0.1653418373058960E-01;
  COFD[         144] =  -0.1632103748685788E+02;
  COFD[         145] =   0.3805435800233515E+01;
  COFD[         146] =  -0.2813339083459423E+00;
  COFD[         147] =   0.1226875862421682E-01;
  COFD[         148] =  -0.1862959184503705E+02;
  COFD[         149] =   0.4599775876853535E+01;
  COFD[         150] =  -0.3762918861611520E+00;
  COFD[         151] =   0.1604808326150353E-01;
  COFD[         152] =  -0.2384916272500252E+02;
  COFD[         153] =   0.5010115663668988E+01;
  COFD[         154] =  -0.3897718053261842E+00;
  COFD[         155] =   0.1509403961940095E-01;
  COFD[         156] =  -0.1725286194864467E+02;
  COFD[         157] =   0.4167814533234269E+01;
  COFD[         158] =  -0.3247525848550428E+00;
  COFD[         159] =   0.1400034963858472E-01;
  COFD[         160] =  -0.1752194960164599E+02;
  COFD[         161] =   0.4875287931841000E+01;
  COFD[         162] =  -0.4147968034669867E+00;
  COFD[         163] =   0.1781089573575908E-01;
  COFD[         164] =  -0.1742491429482865E+02;
  COFD[         165] =   0.4846655191576849E+01;
  COFD[         166] =  -0.4113732013955531E+00;
  COFD[         167] =   0.1767520282723403E-01;
  COFD[         168] =  -0.1476789275356774E+02;
  COFD[         169] =   0.4196229061773622E+01;
  COFD[         170] =  -0.3280338158230447E+00;
  COFD[         171] =   0.1412562078478874E-01;
  COFD[         172] =  -0.1509416127412756E+02;
  COFD[         173] =   0.4170272724773333E+01;
  COFD[         174] =  -0.3329631459117541E+00;
  COFD[         175] =   0.1468212402294482E-01;
  COFD[         176] =  -0.1684749465833484E+02;
  COFD[         177] =   0.4372293660809547E+01;
  COFD[         178] =  -0.3053505024899502E+00;
  COFD[         179] =   0.1131162956306083E-01;
  COFD[         180] =  -0.1507688104939430E+02;
  COFD[         181] =   0.4163842058809374E+01;
  COFD[         182] =  -0.3320920473747643E+00;
  COFD[         183] =   0.1464289366683665E-01;
  COFD[         184] =  -0.1636810982086219E+02;
  COFD[         185] =   0.4006746340787458E+01;
  COFD[         186] =  -0.2528221027540429E+00;
  COFD[         187] =   0.8807021201168675E-02;
  COFD[         188] =  -0.1689232568564453E+02;
  COFD[         189] =   0.4730273556795838E+01;
  COFD[         190] =  -0.4048906742596171E+00;
  COFD[         191] =   0.1776022334223405E-01;
  COFD[         192] =  -0.1508990560668830E+02;
  COFD[         193] =   0.4095531915411983E+01;
  COFD[         194] =  -0.3236090887729158E+00;
  COFD[         195] =   0.1429063082534494E-01;
  COFD[         196] =  -0.1174012693931505E+02;
  COFD[         197] =   0.2906468221415321E+01;
  COFD[         198] =  -0.1669926932934059E+00;
  COFD[         199] =   0.7416448355870972E-02;
  COFD[         200] =  -0.1669804215307749E+02;
  COFD[         201] =   0.4646006533188189E+01;
  COFD[         202] =  -0.3904610949204326E+00;
  COFD[         203] =   0.1698646309375116E-01;
  COFD[         204] =  -0.1635692863517560E+02;
  COFD[         205] =   0.4001315639862050E+01;
  COFD[         206] =  -0.2520000132543742E+00;
  COFD[         207] =   0.8766525251705412E-02;
  COFD[         208] =  -0.1694687300815454E+02;
  COFD[         209] =   0.4748148455016244E+01;
  COFD[         210] =  -0.4072480572350639E+00;
  COFD[         211] =   0.1786385780735676E-01;
  COFD[         212] =  -0.1782081542338524E+02;
  COFD[         213] =   0.4822281209606961E+01;
  COFD[         214] =  -0.3887473651235728E+00;
  COFD[         215] =   0.1588256266086469E-01;
  COFD[         216] =  -0.1501308370949002E+02;
  COFD[         217] =   0.4140201601820930E+01;
  COFD[         218] =  -0.3288897802797416E+00;
  COFD[         219] =   0.1449868103483304E-01;
  COFD[         220] =  -0.1683309905666032E+02;
  COFD[         221] =   0.4365079762720271E+01;
  COFD[         222] =  -0.3042705287369754E+00;
  COFD[         223] =   0.1125889127314920E-01;
  COFD[         224] =  -0.1686978847126216E+02;
  COFD[         225] =   0.4721521983629889E+01;
  COFD[         226] =  -0.4037192718464361E+00;
  COFD[         227] =   0.1770807186574834E-01;
  COFD[         228] =  -0.1772949944660123E+02;
  COFD[         229] =   0.4898286384206769E+01;
  COFD[         230] =  -0.4076359849112164E+00;
  COFD[         231] =   0.1707858545899634E-01;
  COFD[         232] =  -0.2378143430900525E+02;
  COFD[         233] =   0.4975082596789566E+01;
  COFD[         234] =  -0.3839197228971605E+00;
  COFD[         235] =   0.1478949218777245E-01;
  COFD[         236] =  -0.1749395439646448E+02;
  COFD[         237] =   0.4865192761077812E+01;
  COFD[         238] =  -0.4135081828421872E+00;
  COFD[         239] =   0.1775619961599839E-01;
  COFD[         240] =  -0.1543681755468930E+02;
  COFD[         241] =   0.3615102026808611E+01;
  COFD[         242] =  -0.2578306288187833E+00;
  COFD[         243] =   0.1130180672438872E-01;
  COFD[         244] =  -0.1538242046452401E+02;
  COFD[         245] =   0.3594755693746865E+01;
  COFD[         246] =  -0.2552359406239311E+00;
  COFD[         247] =   0.1119129111015014E-01;
  COFD[         248] =  -0.1509416127412756E+02;
  COFD[         249] =   0.4170272724773333E+01;
  COFD[         250] =  -0.3329631459117541E+00;
  COFD[         251] =   0.1468212402294482E-01;
  COFD[         252] =  -0.1338328959053655E+02;
  COFD[         253] =   0.2963583718457605E+01;
  COFD[         254] =  -0.1740734963972286E+00;
  COFD[         255] =   0.7709816773593613E-02;
  COFD[         256] =  -0.1902589893327734E+02;
  COFD[         257] =   0.4963044595189791E+01;
  COFD[         258] =  -0.4152448917094257E+00;
  COFD[         259] =   0.1741189546558487E-01;
  COFD[         260] =  -0.1336878175481828E+02;
  COFD[         261] =   0.2964009975320682E+01;
  COFD[         262] =  -0.1741323405483551E+00;
  COFD[         263] =   0.7712512499580482E-02;
  COFD[         264] =  -0.1910208778121467E+02;
  COFD[         265] =   0.4839856760115432E+01;
  COFD[         266] =  -0.4005219275419654E+00;
  COFD[         267] =   0.1681707424541106E-01;
  COFD[         268] =  -0.1454363744692276E+02;
  COFD[         269] =   0.3252775360002836E+01;
  COFD[         270] =  -0.2120947525124356E+00;
  COFD[         271] =   0.9375917476367752E-02;
  COFD[         272] =  -0.1323521591814711E+02;
  COFD[         273] =   0.2853857957507343E+01;
  COFD[         274] =  -0.1596303547846372E+00;
  COFD[         275] =   0.7075449458871568E-02;
  COFD[         276] =  -0.1064556451628366E+02;
  COFD[         277] =   0.2174543761808112E+01;
  COFD[         278] =  -0.6767589520388503E-01;
  COFD[         279] =   0.2921968469178159E-02;
  COFD[         280] =  -0.1480131803810992E+02;
  COFD[         281] =   0.3366342560298504E+01;
  COFD[         282] =  -0.2268953355973816E+00;
  COFD[         283] =   0.1001940799548309E-01;
  COFD[         284] =  -0.1910754101049230E+02;
  COFD[         285] =   0.4838702285957792E+01;
  COFD[         286] =  -0.4002430660796238E+00;
  COFD[         287] =   0.1679936376165156E-01;
  COFD[         288] =  -0.1453137901013359E+02;
  COFD[         289] =   0.3248474872050811E+01;
  COFD[         290] =  -0.2114664210884821E+00;
  COFD[         291] =   0.9345802638823836E-02;
  COFD[         292] =  -0.1755935993388272E+02;
  COFD[         293] =   0.4314361964780352E+01;
  COFD[         294] =  -0.3441940729036914E+00;
  COFD[         295] =   0.1485616821003848E-01;
  COFD[         296] =  -0.1333077163202364E+02;
  COFD[         297] =   0.2971576640967610E+01;
  COFD[         298] =  -0.1751768817382390E+00;
  COFD[         299] =   0.7760363602679025E-02;
  COFD[         300] =  -0.1903987320392389E+02;
  COFD[         301] =   0.4963240312982511E+01;
  COFD[         302] =  -0.4152420737994104E+00;
  COFD[         303] =   0.1741048089845645E-01;
  COFD[         304] =  -0.1450174966203873E+02;
  COFD[         305] =   0.3240507460697169E+01;
  COFD[         306] =  -0.2104181641466212E+00;
  COFD[         307] =   0.9299811294433670E-02;
  COFD[         308] =  -0.1697689861681293E+02;
  COFD[         309] =   0.4143440294873640E+01;
  COFD[         310] =  -0.3245402843230118E+00;
  COFD[         311] =   0.1411020067125902E-01;
  COFD[         312] =  -0.2349867378552522E+02;
  COFD[         313] =   0.5062620989049528E+01;
  COFD[         314] =  -0.4159679664653171E+00;
  COFD[         315] =   0.1696493513352471E-01;
  COFD[         316] =  -0.1541179728688726E+02;
  COFD[         317] =   0.3612468014590696E+01;
  COFD[         318] =  -0.2574709286409149E+00;
  COFD[         319] =   0.1128548676893656E-01;
  COFD[         320] =  -0.2027289153258947E+02;
  COFD[         321] =   0.5157651529525397E+01;
  COFD[         322] =  -0.4233239455087539E+00;
  COFD[         323] =   0.1708864601311192E-01;
  COFD[         324] =  -0.2046809044123080E+02;
  COFD[         325] =   0.5168356480779646E+01;
  COFD[         326] =  -0.4170015852313753E+00;
  COFD[         327] =   0.1651418922298926E-01;
  COFD[         328] =  -0.1684749465833484E+02;
  COFD[         329] =   0.4372293660809547E+01;
  COFD[         330] =  -0.3053505024899502E+00;
  COFD[         331] =   0.1131162956306083E-01;
  COFD[         332] =  -0.1902589893327734E+02;
  COFD[         333] =   0.4963044595189791E+01;
  COFD[         334] =  -0.4152448917094257E+00;
  COFD[         335] =   0.1741189546558487E-01;
  COFD[         336] =  -0.1276497097890644E+02;
  COFD[         337] =   0.1325671264699943E+01;
  COFD[         338] =   0.1805397215957542E+00;
  COFD[         339] =  -0.1280666613762271E-01;
  COFD[         340] =  -0.1901095003690279E+02;
  COFD[         341] =   0.4963229311521320E+01;
  COFD[         342] =  -0.4152358422934130E+00;
  COFD[         343] =   0.1740999487266200E-01;
  COFD[         344] =  -0.1803570966595964E+02;
  COFD[         345] =   0.3728478005860654E+01;
  COFD[         346] =  -0.1834167152948895E+00;
  COFD[         347] =   0.4682246201731059E-02;
  COFD[         348] =  -0.2023623629656903E+02;
  COFD[         349] =   0.5159004982731338E+01;
  COFD[         350] =  -0.4282251257434644E+00;
  COFD[         351] =   0.1747859336162540E-01;
  COFD[         352] =  -0.1881852755665777E+02;
  COFD[         353] =   0.4848454576458683E+01;
  COFD[         354] =  -0.4024998441095268E+00;
  COFD[         355] =   0.1694398072710132E-01;
  COFD[         356] =  -0.1788429952817511E+02;
  COFD[         357] =   0.4910613916875719E+01;
  COFD[         358] =  -0.4181982422933196E+00;
  COFD[         359] =   0.1792763300105206E-01;
  COFD[         360] =  -0.2026577457804009E+02;
  COFD[         361] =   0.5155490277841742E+01;
  COFD[         362] =  -0.4244639218946512E+00;
  COFD[         363] =   0.1718635599262250E-01;
  COFD[         364] =  -0.1803249510513255E+02;
  COFD[         365] =   0.3724553191438368E+01;
  COFD[         366] =  -0.1829095413108150E+00;
  COFD[         367] =   0.4661837926329214E-02;
  COFD[         368] =  -0.2030409270584218E+02;
  COFD[         369] =   0.5183432603071942E+01;
  COFD[         370] =  -0.4310452456553683E+00;
  COFD[         371] =   0.1758768191971549E-01;
  COFD[         372] =  -0.1943698466014210E+02;
  COFD[         373] =   0.4447274740330095E+01;
  COFD[         374] =  -0.2978200583800502E+00;
  COFD[         375] =   0.1038260495688982E-01;
  COFD[         376] =  -0.1894845529969471E+02;
  COFD[         377] =   0.4959287092588942E+01;
  COFD[         378] =  -0.4144391273613986E+00;
  COFD[         379] =   0.1736295122665337E-01;
  COFD[         380] =  -0.1277807183938872E+02;
  COFD[         381] =   0.1325555296302682E+01;
  COFD[         382] =   0.1805522647353651E+00;
  COFD[         383] =  -0.1280700633693305E-01;
  COFD[         384] =  -0.2025241161910038E+02;
  COFD[         385] =   0.5172982615634614E+01;
  COFD[         386] =  -0.4303950755374367E+00;
  COFD[         387] =   0.1758759466858050E-01;
  COFD[         388] =  -0.1994374085526542E+02;
  COFD[         389] =   0.4830081980319046E+01;
  COFD[         390] =  -0.3644950600153112E+00;
  COFD[         391] =   0.1391192544233011E-01;
  COFD[         392] =  -0.1726613964269955E+02;
  COFD[         393] =   0.1754487534188974E+01;
  COFD[         394] =   0.1044739320789194E+00;
  COFD[         395] =  -0.8933007198223725E-02;
  COFD[         396] =  -0.2025720433242960E+02;
  COFD[         397] =   0.5159185593402682E+01;
  COFD[         398] =  -0.4235430393599882E+00;
  COFD[         399] =   0.1709888299272256E-01;
  COFD[         400] =  -0.1543125672017102E+02;
  COFD[         401] =   0.3618855115638650E+01;
  COFD[         402] =  -0.2583419121546929E+00;
  COFD[         403] =   0.1132494930338730E-01;
  COFD[         404] =  -0.1537895325066121E+02;
  COFD[         405] =   0.3599582752195144E+01;
  COFD[         406] =  -0.2558927079381389E+00;
  COFD[         407] =   0.1122098249033523E-01;
  COFD[         408] =  -0.1507688104939430E+02;
  COFD[         409] =   0.4163842058809374E+01;
  COFD[         410] =  -0.3320920473747643E+00;
  COFD[         411] =   0.1464289366683665E-01;
  COFD[         412] =  -0.1336878175481828E+02;
  COFD[         413] =   0.2964009975320682E+01;
  COFD[         414] =  -0.1741323405483551E+00;
  COFD[         415] =   0.7712512499580482E-02;
  COFD[         416] =  -0.1901095003690279E+02;
  COFD[         417] =   0.4963229311521320E+01;
  COFD[         418] =  -0.4152358422934130E+00;
  COFD[         419] =   0.1740999487266200E-01;
  COFD[         420] =  -0.1335274181252259E+02;
  COFD[         421] =   0.2963583718457631E+01;
  COFD[         422] =  -0.1740734963972324E+00;
  COFD[         423] =   0.7709816773593801E-02;
  COFD[         424] =  -0.1908113950961190E+02;
  COFD[         425] =   0.4837754694510448E+01;
  COFD[         426] =  -0.4000161110029841E+00;
  COFD[         427] =   0.1678499662472105E-01;
  COFD[         428] =  -0.1455112852542608E+02;
  COFD[         429] =   0.3264209638419735E+01;
  COFD[         430] =  -0.2136573392296478E+00;
  COFD[         431] =   0.9446846235261253E-02;
  COFD[         432] =  -0.1322032806769774E+02;
  COFD[         433] =   0.2853486528451670E+01;
  COFD[         434] =  -0.1595787301314707E+00;
  COFD[         435] =   0.7073074545706415E-02;
  COFD[         436] =  -0.1064253861628107E+02;
  COFD[         437] =   0.2174550227335954E+01;
  COFD[         438] =  -0.6767643372502988E-01;
  COFD[         439] =   0.2921979246312271E-02;
  COFD[         440] =  -0.1481108040195451E+02;
  COFD[         441] =   0.3378918244574303E+01;
  COFD[         442] =  -0.2286174610007034E+00;
  COFD[         443] =   0.1009773226034280E-01;
  COFD[         444] =  -0.1908621188318923E+02;
  COFD[         445] =   0.4836508433338681E+01;
  COFD[         446] =  -0.3997198537984301E+00;
  COFD[         447] =   0.1676629685125576E-01;
  COFD[         448] =  -0.1453849109247393E+02;
  COFD[         449] =   0.3259546448493950E+01;
  COFD[         450] =  -0.2129797622055096E+00;
  COFD[         451] =   0.9414509352545126E-02;
  COFD[         452] =  -0.1757016199925516E+02;
  COFD[         453] =   0.4327375250167051E+01;
  COFD[         454] =  -0.3458635010619581E+00;
  COFD[         455] =   0.1492736901823257E-01;
  COFD[         456] =  -0.1330980310661613E+02;
  COFD[         457] =   0.2968380352605879E+01;
  COFD[         458] =  -0.1747356558119914E+00;
  COFD[         459] =   0.7740150913675598E-02;
  COFD[         460] =  -0.1902343852801645E+02;
  COFD[         461] =   0.4962858381054779E+01;
  COFD[         462] =  -0.4151378687378806E+00;
  COFD[         463] =   0.1740353197744324E-01;
  COFD[         464] =  -0.1450839204389819E+02;
  COFD[         465] =   0.3251381640250417E+01;
  COFD[         466] =  -0.2119042844536192E+00;
  COFD[         467] =   0.9367271609371054E-02;
  COFD[         468] =  -0.1698858354738882E+02;
  COFD[         469] =   0.4157216748608500E+01;
  COFD[         470] =  -0.3263536052448225E+00;
  COFD[         471] =   0.1418959320103276E-01;
  COFD[         472] =  -0.2349858517247871E+02;
  COFD[         473] =   0.5062581035220192E+01;
  COFD[         474] =  -0.4159619955775484E+00;
  COFD[         475] =   0.1696464437402649E-01;
  COFD[         476] =  -0.1540309090870326E+02;
  COFD[         477] =   0.3614696445371578E+01;
  COFD[         478] =  -0.2577753207878676E+00;
  COFD[         479] =   0.1129930080443808E-01;
  COFD[         480] =  -0.2048209697647840E+02;
  COFD[         481] =   0.5127870094995610E+01;
  COFD[         482] =  -0.4231024044426618E+00;
  COFD[         483] =   0.1721759814463075E-01;
  COFD[         484] =  -0.2046876559898188E+02;
  COFD[         485] =   0.5125566020964355E+01;
  COFD[         486] =  -0.4233467066361414E+00;
  COFD[         487] =   0.1724767006424224E-01;
  COFD[         488] =  -0.1636810982086219E+02;
  COFD[         489] =   0.4006746340787458E+01;
  COFD[         490] =  -0.2528221027540429E+00;
  COFD[         491] =   0.8807021201168675E-02;
  COFD[         492] =  -0.1910208778121467E+02;
  COFD[         493] =   0.4839856760115432E+01;
  COFD[         494] =  -0.4005219275419654E+00;
  COFD[         495] =   0.1681707424541106E-01;
  COFD[         496] =  -0.1803570966595964E+02;
  COFD[         497] =   0.3728478005860654E+01;
  COFD[         498] =  -0.1834167152948895E+00;
  COFD[         499] =   0.4682246201731059E-02;
  COFD[         500] =  -0.1908113950961190E+02;
  COFD[         501] =   0.4837754694510448E+01;
  COFD[         502] =  -0.4000161110029841E+00;
  COFD[         503] =   0.1678499662472105E-01;
  COFD[         504] =  -0.1926515259058734E+02;
  COFD[         505] =   0.4113883642651055E+01;
  COFD[         506] =  -0.2408304609559599E+00;
  COFD[         507] =   0.7421376940711312E-02;
  COFD[         508] =  -0.2001841325156043E+02;
  COFD[         509] =   0.5003088215312937E+01;
  COFD[         510] =  -0.4191428557016870E+00;
  COFD[         511] =   0.1753003548356206E-01;
  COFD[         512] =  -0.1899155557536321E+02;
  COFD[         513] =   0.4783494119468799E+01;
  COFD[         514] =  -0.3960093718772315E+00;
  COFD[         515] =   0.1673518007514816E-01;
  COFD[         516] =  -0.1711220661498793E+02;
  COFD[         517] =   0.4512570318598010E+01;
  COFD[         518] =  -0.3700405486787200E+00;
  COFD[         519] =   0.1597569624722800E-01;
  COFD[         520] =  -0.2016028662112187E+02;
  COFD[         521] =   0.5045245430044817E+01;
  COFD[         522] =  -0.4220163495492107E+00;
  COFD[         523] =   0.1755145130634340E-01;
  COFD[         524] =  -0.1927323685125299E+02;
  COFD[         525] =   0.4113721827894565E+01;
  COFD[         526] =  -0.2408083280145364E+00;
  COFD[         527] =   0.7420403940387676E-02;
  COFD[         528] =  -0.2002404489922154E+02;
  COFD[         529] =   0.5007446109627781E+01;
  COFD[         530] =  -0.4195576074513644E+00;
  COFD[         531] =   0.1754187817233131E-01;
  COFD[         532] =  -0.2106889603355429E+02;
  COFD[         533] =   0.5078234125926331E+01;
  COFD[         534] =  -0.3990596716857998E+00;
  COFD[         535] =   0.1550927632108276E-01;
  COFD[         536] =  -0.1900572705128738E+02;
  COFD[         537] =   0.4829629200885549E+01;
  COFD[         538] =  -0.3981162237001974E+00;
  COFD[         539] =   0.1666585588428635E-01;
  COFD[         540] =  -0.1806710517381950E+02;
  COFD[         541] =   0.3734466949622932E+01;
  COFD[         542] =  -0.1841905263564501E+00;
  COFD[         543] =   0.4713377953567380E-02;
  COFD[         544] =  -0.2000273357738514E+02;
  COFD[         545] =   0.5003613670288710E+01;
  COFD[         546] =  -0.4191992889668390E+00;
  COFD[         547] =   0.1753189913176521E-01;
  COFD[         548] =  -0.2090688825017299E+02;
  COFD[         549] =   0.5106481690478173E+01;
  COFD[         550] =  -0.4094432279607686E+00;
  COFD[         551] =   0.1620751709314507E-01;
  COFD[         552] =  -0.1871985392716656E+02;
  COFD[         553] =   0.2406373385171115E+01;
  COFD[         554] =   0.1073057268211112E-01;
  COFD[         555] =  -0.4539455264233341E-02;
  COFD[         556] =  -0.2042648432871773E+02;
  COFD[         557] =   0.5112800688990101E+01;
  COFD[         558] =  -0.4208144038227885E+00;
  COFD[         559] =   0.1710457743576905E-01;
  COFD[         560] =  -0.1645223742454043E+02;
  COFD[         561] =   0.3854697730165210E+01;
  COFD[         562] =  -0.2878037882685343E+00;
  COFD[         563] =   0.1255204150178310E-01;
  COFD[         564] =  -0.1635833918443450E+02;
  COFD[         565] =   0.3815691616616636E+01;
  COFD[         566] =  -0.2827061894445888E+00;
  COFD[         567] =   0.1232981052514379E-01;
  COFD[         568] =  -0.1689232568564453E+02;
  COFD[         569] =   0.4730273556795838E+01;
  COFD[         570] =  -0.4048906742596171E+00;
  COFD[         571] =   0.1776022334223405E-01;
  COFD[         572] =  -0.1454363744692276E+02;
  COFD[         573] =   0.3252775360002836E+01;
  COFD[         574] =  -0.2120947525124356E+00;
  COFD[         575] =   0.9375917476367752E-02;
  COFD[         576] =  -0.2023623629656903E+02;
  COFD[         577] =   0.5159004982731338E+01;
  COFD[         578] =  -0.4282251257434644E+00;
  COFD[         579] =   0.1747859336162540E-01;
  COFD[         580] =  -0.1455112852542608E+02;
  COFD[         581] =   0.3264209638419735E+01;
  COFD[         582] =  -0.2136573392296478E+00;
  COFD[         583] =   0.9446846235261253E-02;
  COFD[         584] =  -0.2001841325156043E+02;
  COFD[         585] =   0.5003088215312937E+01;
  COFD[         586] =  -0.4191428557016870E+00;
  COFD[         587] =   0.1753003548356206E-01;
  COFD[         588] =  -0.1525210940750062E+02;
  COFD[         589] =   0.3347991675957177E+01;
  COFD[         590] =  -0.2231616752417136E+00;
  COFD[         591] =   0.9800264415023702E-02;
  COFD[         592] =  -0.1430858947735535E+02;
  COFD[         593] =   0.3127427060682429E+01;
  COFD[         594] =  -0.1960980103738432E+00;
  COFD[         595] =   0.8695656429692035E-02;
  COFD[         596] =  -0.1216112311189428E+02;
  COFD[         597] =   0.2673271002528260E+01;
  COFD[         598] =  -0.1366481193799611E+00;
  COFD[         599] =   0.6094303407106706E-02;
  COFD[         600] =  -0.1557170923538268E+02;
  COFD[         601] =   0.3481269416370885E+01;
  COFD[         602] =  -0.2405002957692937E+00;
  COFD[         603] =   0.1055277270605560E-01;
  COFD[         604] =  -0.2002802788514385E+02;
  COFD[         605] =   0.5003477827938354E+01;
  COFD[         606] =  -0.4191917872787976E+00;
  COFD[         607] =   0.1753204616014205E-01;
  COFD[         608] =  -0.1526513030514044E+02;
  COFD[         609] =   0.3356640370381600E+01;
  COFD[         610] =  -0.2242905925759323E+00;
  COFD[         611] =   0.9849412614345979E-02;
  COFD[         612] =  -0.1818484070796676E+02;
  COFD[         613] =   0.4358837296217939E+01;
  COFD[         614] =  -0.3476944141512220E+00;
  COFD[         615] =   0.1491847840720065E-01;
  COFD[         616] =  -0.1457931870566920E+02;
  COFD[         617] =   0.3305301710138485E+01;
  COFD[         618] =  -0.2192722946715860E+00;
  COFD[         619] =   0.9701700683944362E-02;
  COFD[         620] =  -0.2027645055820890E+02;
  COFD[         621] =   0.5170147189288856E+01;
  COFD[         622] =  -0.4299550299922875E+00;
  COFD[         623] =   0.1756549562746075E-01;
  COFD[         624] =  -0.1523620653214519E+02;
  COFD[         625] =   0.3348643253645854E+01;
  COFD[         626] =  -0.2232496547504182E+00;
  COFD[         627] =   0.9804212145714645E-02;
  COFD[         628] =  -0.1753994518424459E+02;
  COFD[         629] =   0.4158828151149182E+01;
  COFD[         630] =  -0.3241998792087618E+00;
  COFD[         631] =   0.1400163924604856E-01;
  COFD[         632] =  -0.2363569469280461E+02;
  COFD[         633] =   0.5050661168714743E+01;
  COFD[         634] =  -0.4076790606763266E+00;
  COFD[         635] =   0.1634373025964126E-01;
  COFD[         636] =  -0.1645726609264868E+02;
  COFD[         637] =   0.3866322991209679E+01;
  COFD[         638] =  -0.2893595968515034E+00;
  COFD[         639] =   0.1262127247206273E-01;
  COFD[         640] =  -0.1514651789245071E+02;
  COFD[         641] =   0.3462312981288904E+01;
  COFD[         642] =  -0.2380726416670257E+00;
  COFD[         643] =   0.1044895488982180E-01;
  COFD[         644] =  -0.1508728524143146E+02;
  COFD[         645] =   0.3439843254137499E+01;
  COFD[         646] =  -0.2351547205941974E+00;
  COFD[         647] =   0.1032244834588984E-01;
  COFD[         648] =  -0.1508990560668830E+02;
  COFD[         649] =   0.4095531915411983E+01;
  COFD[         650] =  -0.3236090887729158E+00;
  COFD[         651] =   0.1429063082534494E-01;
  COFD[         652] =  -0.1323521591814711E+02;
  COFD[         653] =   0.2853857957507343E+01;
  COFD[         654] =  -0.1596303547846372E+00;
  COFD[         655] =   0.7075449458871568E-02;
  COFD[         656] =  -0.1881852755665777E+02;
  COFD[         657] =   0.4848454576458683E+01;
  COFD[         658] =  -0.4024998441095268E+00;
  COFD[         659] =   0.1694398072710132E-01;
  COFD[         660] =  -0.1322032806769774E+02;
  COFD[         661] =   0.2853486528451670E+01;
  COFD[         662] =  -0.1595787301314707E+00;
  COFD[         663] =   0.7073074545706415E-02;
  COFD[         664] =  -0.1899155557536321E+02;
  COFD[         665] =   0.4783494119468799E+01;
  COFD[         666] =  -0.3960093718772315E+00;
  COFD[         667] =   0.1673518007514816E-01;
  COFD[         668] =  -0.1430858947735535E+02;
  COFD[         669] =   0.3127427060682429E+01;
  COFD[         670] =  -0.1960980103738432E+00;
  COFD[         671] =   0.8695656429692035E-02;
  COFD[         672] =  -0.1307749522108035E+02;
  COFD[         673] =   0.2744620023573023E+01;
  COFD[         674] =  -0.1452254920719674E+00;
  COFD[         675] =   0.6441841776411875E-02;
  COFD[         676] =  -0.1066740661833492E+02;
  COFD[         677] =   0.2118871949060134E+01;
  COFD[         678] =  -0.6030629519370559E-01;
  COFD[         679] =   0.2595581600631073E-02;
  COFD[         680] =  -0.1455521688916326E+02;
  COFD[         681] =   0.3233064865160284E+01;
  COFD[         682] =  -0.2095166917359136E+00;
  COFD[         683] =   0.9263348706734287E-02;
  COFD[         684] =  -0.1900014767186941E+02;
  COFD[         685] =   0.4784112582239475E+01;
  COFD[         686] =  -0.3959959506254504E+00;
  COFD[         687] =   0.1673047086070456E-01;
  COFD[         688] =  -0.1430384868728555E+02;
  COFD[         689] =   0.3125877575358964E+01;
  COFD[         690] =  -0.1958410518042518E+00;
  COFD[         691] =   0.8682193205361071E-02;
  COFD[         692] =  -0.1746769854797995E+02;
  COFD[         693] =   0.4264008519853375E+01;
  COFD[         694] =  -0.3398132864002393E+00;
  COFD[         695] =   0.1475670473360388E-01;
  COFD[         696] =  -0.1318009297471136E+02;
  COFD[         697] =   0.2857967371156021E+01;
  COFD[         698] =  -0.1601949749791629E+00;
  COFD[         699] =   0.7101221676809524E-02;
  COFD[         700] =  -0.1883093215750172E+02;
  COFD[         701] =   0.4848462442058104E+01;
  COFD[         702] =  -0.4024597684582722E+00;
  COFD[         703] =   0.1694041152472868E-01;
  COFD[         704] =  -0.1427294640051003E+02;
  COFD[         705] =   0.3117012895539554E+01;
  COFD[         706] =  -0.1946596484063036E+00;
  COFD[         707] =   0.8629721402195967E-02;
  COFD[         708] =  -0.1672667778972758E+02;
  COFD[         709] =   0.4023102094303094E+01;
  COFD[         710] =  -0.3103163375902530E+00;
  COFD[         711] =   0.1355332875778464E-01;
  COFD[         712] =  -0.2328225161733864E+02;
  COFD[         713] =   0.5008847120014931E+01;
  COFD[         714] =  -0.4122620203443967E+00;
  COFD[         715] =   0.1693220599831681E-01;
  COFD[         716] =  -0.1512111139020410E+02;
  COFD[         717] =   0.3458819295912228E+01;
  COFD[         718] =  -0.2375970826594765E+00;
  COFD[         719] =   0.1042745102315496E-01;
  COFD[         720] =  -0.1322163692827858E+02;
  COFD[         721] =   0.3075522930065609E+01;
  COFD[         722] =  -0.1904123874279204E+00;
  COFD[         723] =   0.8492834273761591E-02;
  COFD[         724] =  -0.1315730506576576E+02;
  COFD[         725] =   0.3055295555288239E+01;
  COFD[         726] =  -0.1877084814632856E+00;
  COFD[         727] =   0.8372316362831445E-02;
  COFD[         728] =  -0.1174012693931505E+02;
  COFD[         729] =   0.2906468221415321E+01;
  COFD[         730] =  -0.1669926932934059E+00;
  COFD[         731] =   0.7416448355870972E-02;
  COFD[         732] =  -0.1064556451628366E+02;
  COFD[         733] =   0.2174543761808112E+01;
  COFD[         734] =  -0.6767589520388503E-01;
  COFD[         735] =   0.2921968469178159E-02;
  COFD[         736] =  -0.1788429952817511E+02;
  COFD[         737] =   0.4910613916875719E+01;
  COFD[         738] =  -0.4181982422933196E+00;
  COFD[         739] =   0.1792763300105206E-01;
  COFD[         740] =  -0.1064253861628107E+02;
  COFD[         741] =   0.2174550227335954E+01;
  COFD[         742] =  -0.6767643372502988E-01;
  COFD[         743] =   0.2921979246312271E-02;
  COFD[         744] =  -0.1711220661498793E+02;
  COFD[         745] =   0.4512570318598010E+01;
  COFD[         746] =  -0.3700405486787200E+00;
  COFD[         747] =   0.1597569624722800E-01;
  COFD[         748] =  -0.1216112311189428E+02;
  COFD[         749] =   0.2673271002528260E+01;
  COFD[         750] =  -0.1366481193799611E+00;
  COFD[         751] =   0.6094303407106706E-02;
  COFD[         752] =  -0.1066740661833492E+02;
  COFD[         753] =   0.2118871949060134E+01;
  COFD[         754] =  -0.6030629519370559E-01;
  COFD[         755] =   0.2595581600631073E-02;
  COFD[         756] =  -0.1024490542151301E+02;
  COFD[         757] =   0.2159570672972225E+01;
  COFD[         758] =  -0.7052547816977366E-01;
  COFD[         759] =   0.3272744024019213E-02;
  COFD[         760] =  -0.1227253578601752E+02;
  COFD[         761] =   0.2734362066169806E+01;
  COFD[         762] =  -0.1449054835127276E+00;
  COFD[         763] =   0.6468227545096617E-02;
  COFD[         764] =  -0.1711971921675324E+02;
  COFD[         765] =   0.4515093947815878E+01;
  COFD[         766] =  -0.3703631440941066E+00;
  COFD[         767] =   0.1598940936575525E-01;
  COFD[         768] =  -0.1217676276160025E+02;
  COFD[         769] =   0.2676702831726427E+01;
  COFD[         770] =  -0.1371110647811714E+00;
  COFD[         771] =   0.6115228536440196E-02;
  COFD[         772] =  -0.1537900621577918E+02;
  COFD[         773] =   0.3905666070712026E+01;
  COFD[         774] =  -0.3001957706103549E+00;
  COFD[         775] =   0.1333185506910619E-01;
  COFD[         776] =  -0.1063104236071595E+02;
  COFD[         777] =   0.2174582020370832E+01;
  COFD[         778] =  -0.6767931641091007E-01;
  COFD[         779] =   0.2922054655651823E-02;
  COFD[         780] =  -0.1789738457973026E+02;
  COFD[         781] =   0.4914689305590001E+01;
  COFD[         782] =  -0.4187174010126591E+00;
  COFD[         783] =   0.1794961882909394E-01;
  COFD[         784] =  -0.1215152187647601E+02;
  COFD[         785] =   0.2670065262560549E+01;
  COFD[         786] =  -0.1362026991746787E+00;
  COFD[         787] =   0.6073761265771448E-02;
  COFD[         788] =  -0.1437677954608083E+02;
  COFD[         789] =   0.3542340190226259E+01;
  COFD[         790] =  -0.2518343923430894E+00;
  COFD[         791] =   0.1118661472548827E-01;
  COFD[         792] =  -0.2177638593119188E+02;
  COFD[         793] =   0.4563067342827384E+01;
  COFD[         794] =  -0.3686945303206421E+00;
  COFD[         795] =   0.1561188821007837E-01;
  COFD[         796] =  -0.1320038983011787E+02;
  COFD[         797] =   0.3068261259040047E+01;
  COFD[         798] =  -0.1894129351750649E+00;
  COFD[         799] =   0.8447150657362242E-02;
  COFD[         800] =  -0.1672933652118024E+02;
  COFD[         801] =   0.3965785816566581E+01;
  COFD[         802] =  -0.3013281840320347E+00;
  COFD[         803] =   0.1310001478947402E-01;
  COFD[         804] =  -0.1677759619477616E+02;
  COFD[         805] =   0.3988901785374923E+01;
  COFD[         806] =  -0.3051808273299473E+00;
  COFD[         807] =   0.1330327276290500E-01;
  COFD[         808] =  -0.1669804215307749E+02;
  COFD[         809] =   0.4646006533188189E+01;
  COFD[         810] =  -0.3904610949204326E+00;
  COFD[         811] =   0.1698646309375116E-01;
  COFD[         812] =  -0.1480131803810992E+02;
  COFD[         813] =   0.3366342560298504E+01;
  COFD[         814] =  -0.2268953355973816E+00;
  COFD[         815] =   0.1001940799548309E-01;
  COFD[         816] =  -0.2026577457804009E+02;
  COFD[         817] =   0.5155490277841742E+01;
  COFD[         818] =  -0.4244639218946512E+00;
  COFD[         819] =   0.1718635599262250E-01;
  COFD[         820] =  -0.1481108040195451E+02;
  COFD[         821] =   0.3378918244574303E+01;
  COFD[         822] =  -0.2286174610007034E+00;
  COFD[         823] =   0.1009773226034280E-01;
  COFD[         824] =  -0.2016028662112187E+02;
  COFD[         825] =   0.5045245430044817E+01;
  COFD[         826] =  -0.4220163495492107E+00;
  COFD[         827] =   0.1755145130634340E-01;
  COFD[         828] =  -0.1557170923538268E+02;
  COFD[         829] =   0.3481269416370885E+01;
  COFD[         830] =  -0.2405002957692937E+00;
  COFD[         831] =   0.1055277270605560E-01;
  COFD[         832] =  -0.1455521688916326E+02;
  COFD[         833] =   0.3233064865160284E+01;
  COFD[         834] =  -0.2095166917359136E+00;
  COFD[         835] =   0.9263348706734287E-02;
  COFD[         836] =  -0.1227253578601752E+02;
  COFD[         837] =   0.2734362066169806E+01;
  COFD[         838] =  -0.1449054835127276E+00;
  COFD[         839] =   0.6468227545096617E-02;
  COFD[         840] =  -0.1586646081706729E+02;
  COFD[         841] =   0.3603776072047724E+01;
  COFD[         842] =  -0.2562827075890430E+00;
  COFD[         843] =   0.1123156340288876E-01;
  COFD[         844] =  -0.2017053688099255E+02;
  COFD[         845] =   0.5045824644771611E+01;
  COFD[         846] =  -0.4220969778346808E+00;
  COFD[         847] =   0.1755513792183601E-01;
  COFD[         848] =  -0.1558609099838845E+02;
  COFD[         849] =   0.3490617714246639E+01;
  COFD[         850] =  -0.2417221951734723E+00;
  COFD[         851] =   0.1060604952078433E-01;
  COFD[         852] =  -0.1854821020259210E+02;
  COFD[         853] =   0.4501466293539689E+01;
  COFD[         854] =  -0.3655565309363614E+00;
  COFD[         855] =   0.1566780821214448E-01;
  COFD[         856] =  -0.1484612543141704E+02;
  COFD[         857] =   0.3423492643386516E+01;
  COFD[         858] =  -0.2347207973903428E+00;
  COFD[         859] =   0.1037529625376389E-01;
  COFD[         860] =  -0.2031055243789758E+02;
  COFD[         861] =   0.5168389671669797E+01;
  COFD[         862] =  -0.4264383492219733E+00;
  COFD[         863] =   0.1728448151830319E-01;
  COFD[         864] =  -0.1555717721575701E+02;
  COFD[         865] =   0.3482723600957817E+01;
  COFD[         866] =  -0.2406978172870443E+00;
  COFD[         867] =   0.1056168468865037E-01;
  COFD[         868] =  -0.1777649791730375E+02;
  COFD[         869] =   0.4248843195751910E+01;
  COFD[         870] =  -0.3346930401290839E+00;
  COFD[         871] =   0.1440673440796217E-01;
  COFD[         872] =  -0.2370839244124438E+02;
  COFD[         873] =   0.5048027935648303E+01;
  COFD[         874] =  -0.4040835196333057E+00;
  COFD[         875] =   0.1606332363762197E-01;
  COFD[         876] =  -0.1672939192904382E+02;
  COFD[         877] =   0.3975404184226489E+01;
  COFD[         878] =  -0.3025898020717506E+00;
  COFD[         879] =   0.1315505412218577E-01;
  COFD[         880] =  -0.2047207935890673E+02;
  COFD[         881] =   0.5120549000988871E+01;
  COFD[         882] =  -0.4219907391201317E+00;
  COFD[         883] =   0.1716268084967035E-01;
  COFD[         884] =  -0.2045967452432749E+02;
  COFD[         885] =   0.5118525449350731E+01;
  COFD[         886] =  -0.4222729551558472E+00;
  COFD[         887] =   0.1719441142995254E-01;
  COFD[         888] =  -0.1635692863517560E+02;
  COFD[         889] =   0.4001315639862050E+01;
  COFD[         890] =  -0.2520000132543742E+00;
  COFD[         891] =   0.8766525251705412E-02;
  COFD[         892] =  -0.1910754101049230E+02;
  COFD[         893] =   0.4838702285957792E+01;
  COFD[         894] =  -0.4002430660796238E+00;
  COFD[         895] =   0.1679936376165156E-01;
  COFD[         896] =  -0.1803249510513255E+02;
  COFD[         897] =   0.3724553191438368E+01;
  COFD[         898] =  -0.1829095413108150E+00;
  COFD[         899] =   0.4661837926329214E-02;
  COFD[         900] =  -0.1908621188318923E+02;
  COFD[         901] =   0.4836508433338681E+01;
  COFD[         902] =  -0.3997198537984301E+00;
  COFD[         903] =   0.1676629685125576E-01;
  COFD[         904] =  -0.1927323685125299E+02;
  COFD[         905] =   0.4113721827894565E+01;
  COFD[         906] =  -0.2408083280145364E+00;
  COFD[         907] =   0.7420403940387676E-02;
  COFD[         908] =  -0.2002802788514385E+02;
  COFD[         909] =   0.5003477827938354E+01;
  COFD[         910] =  -0.4191917872787976E+00;
  COFD[         911] =   0.1753204616014205E-01;
  COFD[         912] =  -0.1900014767186941E+02;
  COFD[         913] =   0.4784112582239475E+01;
  COFD[         914] =  -0.3959959506254504E+00;
  COFD[         915] =   0.1673047086070456E-01;
  COFD[         916] =  -0.1711971921675324E+02;
  COFD[         917] =   0.4515093947815878E+01;
  COFD[         918] =  -0.3703631440941066E+00;
  COFD[         919] =   0.1598940936575525E-01;
  COFD[         920] =  -0.2017053688099255E+02;
  COFD[         921] =   0.5045824644771611E+01;
  COFD[         922] =  -0.4220969778346808E+00;
  COFD[         923] =   0.1755513792183601E-01;
  COFD[         924] =  -0.1928223491479955E+02;
  COFD[         925] =   0.4113883642651062E+01;
  COFD[         926] =  -0.2408304609559608E+00;
  COFD[         927] =   0.7421376940711353E-02;
  COFD[         928] =  -0.2003241114299805E+02;
  COFD[         929] =   0.5007355703272904E+01;
  COFD[         930] =  -0.4195287613916270E+00;
  COFD[         931] =   0.1753985908026623E-01;
  COFD[         932] =  -0.2109531964949764E+02;
  COFD[         933] =   0.5085430028540380E+01;
  COFD[         934] =  -0.4001088392443208E+00;
  COFD[         935] =   0.1555927647634490E-01;
  COFD[         936] =  -0.1900967481367958E+02;
  COFD[         937] =   0.4828170520528887E+01;
  COFD[         938] =  -0.3977810908822604E+00;
  COFD[         939] =   0.1664498875757776E-01;
  COFD[         940] =  -0.1806461227977115E+02;
  COFD[         941] =   0.3730745706118281E+01;
  COFD[         942] =  -0.1837097321884064E+00;
  COFD[         943] =   0.4694035704567372E-02;
  COFD[         944] =  -0.2001104497452950E+02;
  COFD[         945] =   0.5003493097552491E+01;
  COFD[         946] =  -0.4191652678118554E+00;
  COFD[         947] =   0.1752960197527472E-01;
  COFD[         948] =  -0.2093240510980354E+02;
  COFD[         949] =   0.5113264124268066E+01;
  COFD[         950] =  -0.4104425915775513E+00;
  COFD[         951] =   0.1625565099449545E-01;
  COFD[         952] =  -0.1871985825163239E+02;
  COFD[         953] =   0.2406372669599885E+01;
  COFD[         954] =   0.1073111423051258E-01;
  COFD[         955] =  -0.4539506957064400E-02;
  COFD[         956] =  -0.2041541490829319E+02;
  COFD[         957] =   0.5105115933032721E+01;
  COFD[         958] =  -0.4196478966186132E+00;
  COFD[         959] =   0.1704696721515930E-01;
  COFD[         960] =  -0.1645180642806876E+02;
  COFD[         961] =   0.3855410394996957E+01;
  COFD[         962] =  -0.2878581523898601E+00;
  COFD[         963] =   0.1255291913453395E-01;
  COFD[         964] =  -0.1635709831001740E+02;
  COFD[         965] =   0.3816293164170347E+01;
  COFD[         966] =  -0.2827489538440444E+00;
  COFD[         967] =   0.1233030225007116E-01;
  COFD[         968] =  -0.1694687300815454E+02;
  COFD[         969] =   0.4748148455016244E+01;
  COFD[         970] =  -0.4072480572350639E+00;
  COFD[         971] =   0.1786385780735676E-01;
  COFD[         972] =  -0.1453137901013359E+02;
  COFD[         973] =   0.3248474872050811E+01;
  COFD[         974] =  -0.2114664210884821E+00;
  COFD[         975] =   0.9345802638823836E-02;
  COFD[         976] =  -0.2030409270584218E+02;
  COFD[         977] =   0.5183432603071942E+01;
  COFD[         978] =  -0.4310452456553683E+00;
  COFD[         979] =   0.1758768191971549E-01;
  COFD[         980] =  -0.1453849109247393E+02;
  COFD[         981] =   0.3259546448493950E+01;
  COFD[         982] =  -0.2129797622055096E+00;
  COFD[         983] =   0.9414509352545126E-02;
  COFD[         984] =  -0.2002404489922154E+02;
  COFD[         985] =   0.5007446109627781E+01;
  COFD[         986] =  -0.4195576074513644E+00;
  COFD[         987] =   0.1754187817233131E-01;
  COFD[         988] =  -0.1526513030514044E+02;
  COFD[         989] =   0.3356640370381600E+01;
  COFD[         990] =  -0.2242905925759323E+00;
  COFD[         991] =   0.9849412614345979E-02;
  COFD[         992] =  -0.1430384868728555E+02;
  COFD[         993] =   0.3125877575358964E+01;
  COFD[         994] =  -0.1958410518042518E+00;
  COFD[         995] =   0.8682193205361071E-02;
  COFD[         996] =  -0.1217676276160025E+02;
  COFD[         997] =   0.2676702831726427E+01;
  COFD[         998] =  -0.1371110647811714E+00;
  COFD[         999] =   0.6115228536440196E-02;
  COFD[        1000] =  -0.1558609099838845E+02;
  COFD[        1001] =   0.3490617714246639E+01;
  COFD[        1002] =  -0.2417221951734723E+00;
  COFD[        1003] =   0.1060604952078433E-01;
  COFD[        1004] =  -0.2003241114299805E+02;
  COFD[        1005] =   0.5007355703272904E+01;
  COFD[        1006] =  -0.4195287613916270E+00;
  COFD[        1007] =   0.1753985908026623E-01;
  COFD[        1008] =  -0.1527554240109983E+02;
  COFD[        1009] =   0.3363967114483633E+01;
  COFD[        1010] =  -0.2252394474488534E+00;
  COFD[        1011] =   0.9890416446510186E-02;
  COFD[        1012] =  -0.1820884838587697E+02;
  COFD[        1013] =   0.4372531756977406E+01;
  COFD[        1014] =  -0.3493925192307190E+00;
  COFD[        1015] =   0.1498866670036591E-01;
  COFD[        1016] =  -0.1456704681065080E+02;
  COFD[        1017] =   0.3300134269607024E+01;
  COFD[        1018] =  -0.2185270169493516E+00;
  COFD[        1019] =   0.9666339881768727E-02;
  COFD[        1020] =  -0.2034278270814375E+02;
  COFD[        1021] =   0.5194044053733033E+01;
  COFD[        1022] =  -0.4326906002865860E+00;
  COFD[        1023] =   0.1767025742696755E-01;
  COFD[        1024] =  -0.1524693443631203E+02;
  COFD[        1025] =   0.3356096228699615E+01;
  COFD[        1026] =  -0.2242170999002023E+00;
  COFD[        1027] =   0.9846114264084087E-02;
  COFD[        1028] =  -0.1755616817780421E+02;
  COFD[        1029] =   0.4169500620680274E+01;
  COFD[        1030] =  -0.3255016307159463E+00;
  COFD[        1031] =   0.1405440364613088E-01;
  COFD[        1032] =  -0.2363377195653007E+02;
  COFD[        1033] =   0.5047702307335338E+01;
  COFD[        1034] =  -0.4070594187739618E+00;
  COFD[        1035] =   0.1630755566219766E-01;
  COFD[        1036] =  -0.1645625614130163E+02;
  COFD[        1037] =   0.3866580965342516E+01;
  COFD[        1038] =  -0.2893531456511056E+00;
  COFD[        1039] =   0.1261944474539078E-01;
  COFD[        1040] =  -0.1923215106925674E+02;
  COFD[        1041] =   0.4769770645950153E+01;
  COFD[        1042] =  -0.3939017067778038E+00;
  COFD[        1043] =   0.1662898220909356E-01;
  COFD[        1044] =  -0.1917783133928890E+02;
  COFD[        1045] =   0.4748871842879942E+01;
  COFD[        1046] =  -0.3916351371001520E+00;
  COFD[        1047] =   0.1654935017103951E-01;
  COFD[        1048] =  -0.1782081542338524E+02;
  COFD[        1049] =   0.4822281209606961E+01;
  COFD[        1050] =  -0.3887473651235728E+00;
  COFD[        1051] =   0.1588256266086469E-01;
  COFD[        1052] =  -0.1755935993388272E+02;
  COFD[        1053] =   0.4314361964780352E+01;
  COFD[        1054] =  -0.3441940729036914E+00;
  COFD[        1055] =   0.1485616821003848E-01;
  COFD[        1056] =  -0.1943698466014210E+02;
  COFD[        1057] =   0.4447274740330095E+01;
  COFD[        1058] =  -0.2978200583800502E+00;
  COFD[        1059] =   0.1038260495688982E-01;
  COFD[        1060] =  -0.1757016199925516E+02;
  COFD[        1061] =   0.4327375250167051E+01;
  COFD[        1062] =  -0.3458635010619581E+00;
  COFD[        1063] =   0.1492736901823257E-01;
  COFD[        1064] =  -0.2106889603355429E+02;
  COFD[        1065] =   0.5078234125926331E+01;
  COFD[        1066] =  -0.3990596716857998E+00;
  COFD[        1067] =   0.1550927632108276E-01;
  COFD[        1068] =  -0.1818484070796676E+02;
  COFD[        1069] =   0.4358837296217939E+01;
  COFD[        1070] =  -0.3476944141512220E+00;
  COFD[        1071] =   0.1491847840720065E-01;
  COFD[        1072] =  -0.1746769854797995E+02;
  COFD[        1073] =   0.4264008519853375E+01;
  COFD[        1074] =  -0.3398132864002393E+00;
  COFD[        1075] =   0.1475670473360388E-01;
  COFD[        1076] =  -0.1537900621577918E+02;
  COFD[        1077] =   0.3905666070712026E+01;
  COFD[        1078] =  -0.3001957706103549E+00;
  COFD[        1079] =   0.1333185506910619E-01;
  COFD[        1080] =  -0.1854821020259210E+02;
  COFD[        1081] =   0.4501466293539689E+01;
  COFD[        1082] =  -0.3655565309363614E+00;
  COFD[        1083] =   0.1566780821214448E-01;
  COFD[        1084] =  -0.2109531964949764E+02;
  COFD[        1085] =   0.5085430028540380E+01;
  COFD[        1086] =  -0.4001088392443208E+00;
  COFD[        1087] =   0.1555927647634490E-01;
  COFD[        1088] =  -0.1820884838587697E+02;
  COFD[        1089] =   0.4372531756977406E+01;
  COFD[        1090] =  -0.3493925192307190E+00;
  COFD[        1091] =   0.1498866670036591E-01;
  COFD[        1092] =  -0.2067984186725416E+02;
  COFD[        1093] =   0.5097112061739998E+01;
  COFD[        1094] =  -0.4258921464419217E+00;
  COFD[        1095] =   0.1760689931137487E-01;
  COFD[        1096] =  -0.1760160974021058E+02;
  COFD[        1097] =   0.4370332628371653E+01;
  COFD[        1098] =  -0.3513719691423631E+00;
  COFD[        1099] =   0.1516221066908801E-01;
  COFD[        1100] =  -0.1951118321112224E+02;
  COFD[        1101] =   0.4471569363537695E+01;
  COFD[        1102] =  -0.3013374208218036E+00;
  COFD[        1103] =   0.1054895382382678E-01;
  COFD[        1104] =  -0.1817998188473039E+02;
  COFD[        1105] =   0.4364878900220307E+01;
  COFD[        1106] =  -0.3484381669569521E+00;
  COFD[        1107] =   0.1494878692100033E-01;
  COFD[        1108] =  -0.2026340088492126E+02;
  COFD[        1109] =   0.5009715621864962E+01;
  COFD[        1110] =  -0.4198847647706185E+00;
  COFD[        1111] =   0.1755705740818908E-01;
  COFD[        1112] =  -0.2289906109626653E+02;
  COFD[        1113] =   0.4409354958846299E+01;
  COFD[        1114] =  -0.2878825260320446E+00;
  COFD[        1115] =   0.9793621866625414E-02;
  COFD[        1116] =  -0.1921093934794400E+02;
  COFD[        1117] =   0.4770560824446395E+01;
  COFD[        1118] =  -0.3938285475558091E+00;
  COFD[        1119] =   0.1661810137086212E-01;
  COFD[        1120] =  -0.1543306786434106E+02;
  COFD[        1121] =   0.3641484216742679E+01;
  COFD[        1122] =  -0.2614192836082300E+00;
  COFD[        1123] =   0.1146400341320006E-01;
  COFD[        1124] =  -0.1538700627222775E+02;
  COFD[        1125] =   0.3625494228341230E+01;
  COFD[        1126] =  -0.2594144135440631E+00;
  COFD[        1127] =   0.1138002455729723E-01;
  COFD[        1128] =  -0.1501308370949002E+02;
  COFD[        1129] =   0.4140201601820930E+01;
  COFD[        1130] =  -0.3288897802797416E+00;
  COFD[        1131] =   0.1449868103483304E-01;
  COFD[        1132] =  -0.1333077163202364E+02;
  COFD[        1133] =   0.2971576640967610E+01;
  COFD[        1134] =  -0.1751768817382390E+00;
  COFD[        1135] =   0.7760363602679025E-02;
  COFD[        1136] =  -0.1894845529969471E+02;
  COFD[        1137] =   0.4959287092588942E+01;
  COFD[        1138] =  -0.4144391273613986E+00;
  COFD[        1139] =   0.1736295122665337E-01;
  COFD[        1140] =  -0.1330980310661613E+02;
  COFD[        1141] =   0.2968380352605879E+01;
  COFD[        1142] =  -0.1747356558119914E+00;
  COFD[        1143] =   0.7740150913675598E-02;
  COFD[        1144] =  -0.1900572705128738E+02;
  COFD[        1145] =   0.4829629200885549E+01;
  COFD[        1146] =  -0.3981162237001974E+00;
  COFD[        1147] =   0.1666585588428635E-01;
  COFD[        1148] =  -0.1457931870566920E+02;
  COFD[        1149] =   0.3305301710138485E+01;
  COFD[        1150] =  -0.2192722946715860E+00;
  COFD[        1151] =   0.9701700683944362E-02;
  COFD[        1152] =  -0.1318009297471136E+02;
  COFD[        1153] =   0.2857967371156021E+01;
  COFD[        1154] =  -0.1601949749791629E+00;
  COFD[        1155] =   0.7101221676809524E-02;
  COFD[        1156] =  -0.1063104236071595E+02;
  COFD[        1157] =   0.2174582020370832E+01;
  COFD[        1158] =  -0.6767931641091007E-01;
  COFD[        1159] =   0.2922054655651823E-02;
  COFD[        1160] =  -0.1484612543141704E+02;
  COFD[        1161] =   0.3423492643386516E+01;
  COFD[        1162] =  -0.2347207973903428E+00;
  COFD[        1163] =   0.1037529625376389E-01;
  COFD[        1164] =  -0.1900967481367958E+02;
  COFD[        1165] =   0.4828170520528887E+01;
  COFD[        1166] =  -0.3977810908822604E+00;
  COFD[        1167] =   0.1664498875757776E-01;
  COFD[        1168] =  -0.1456704681065080E+02;
  COFD[        1169] =   0.3300134269607024E+01;
  COFD[        1170] =  -0.2185270169493516E+00;
  COFD[        1171] =   0.9666339881768727E-02;
  COFD[        1172] =  -0.1760160974021058E+02;
  COFD[        1173] =   0.4370332628371653E+01;
  COFD[        1174] =  -0.3513719691423631E+00;
  COFD[        1175] =   0.1516221066908801E-01;
  COFD[        1176] =  -0.1324967573338813E+02;
  COFD[        1177] =   0.2963583718457639E+01;
  COFD[        1178] =  -0.1740734963972337E+00;
  COFD[        1179] =   0.7709816773593857E-02;
  COFD[        1180] =  -0.1895623695499981E+02;
  COFD[        1181] =   0.4957160985956500E+01;
  COFD[        1182] =  -0.4140469150912638E+00;
  COFD[        1183] =   0.1734089632858991E-01;
  COFD[        1184] =  -0.1453532015237173E+02;
  COFD[        1185] =   0.3291285383780231E+01;
  COFD[        1186] =  -0.2173571659474091E+00;
  COFD[        1187] =   0.9614779127177037E-02;
  COFD[        1188] =  -0.1702257488378831E+02;
  COFD[        1189] =   0.4202534175316763E+01;
  COFD[        1190] =  -0.3323169248630543E+00;
  COFD[        1191] =   0.1445062729566104E-01;
  COFD[        1192] =  -0.2349830619682465E+02;
  COFD[        1193] =   0.5062455579552836E+01;
  COFD[        1194] =  -0.4159432465895745E+00;
  COFD[        1195] =   0.1696373135599558E-01;
  COFD[        1196] =  -0.1539426704211137E+02;
  COFD[        1197] =   0.3632170199052923E+01;
  COFD[        1198] =  -0.2601531779068780E+00;
  COFD[        1199] =   0.1140681634760632E-01;
  COFD[        1200] =  -0.2028829562643207E+02;
  COFD[        1201] =   0.5159058519743692E+01;
  COFD[        1202] =  -0.4235262462637739E+00;
  COFD[        1203] =   0.1709816078216507E-01;
  COFD[        1204] =  -0.2048707524059223E+02;
  COFD[        1205] =   0.5171197873032807E+01;
  COFD[        1206] =  -0.4174185529276651E+00;
  COFD[        1207] =   0.1653418373058960E-01;
  COFD[        1208] =  -0.1683309905666032E+02;
  COFD[        1209] =   0.4365079762720271E+01;
  COFD[        1210] =  -0.3042705287369754E+00;
  COFD[        1211] =   0.1125889127314920E-01;
  COFD[        1212] =  -0.1903987320392389E+02;
  COFD[        1213] =   0.4963240312982511E+01;
  COFD[        1214] =  -0.4152420737994104E+00;
  COFD[        1215] =   0.1741048089845645E-01;
  COFD[        1216] =  -0.1277807183938872E+02;
  COFD[        1217] =   0.1325555296302682E+01;
  COFD[        1218] =   0.1805522647353651E+00;
  COFD[        1219] =  -0.1280700633693305E-01;
  COFD[        1220] =  -0.1902343852801645E+02;
  COFD[        1221] =   0.4962858381054779E+01;
  COFD[        1222] =  -0.4151378687378806E+00;
  COFD[        1223] =   0.1740353197744324E-01;
  COFD[        1224] =  -0.1806710517381950E+02;
  COFD[        1225] =   0.3734466949622932E+01;
  COFD[        1226] =  -0.1841905263564501E+00;
  COFD[        1227] =   0.4713377953567380E-02;
  COFD[        1228] =  -0.2027645055820890E+02;
  COFD[        1229] =   0.5170147189288856E+01;
  COFD[        1230] =  -0.4299550299922875E+00;
  COFD[        1231] =   0.1756549562746075E-01;
  COFD[        1232] =  -0.1883093215750172E+02;
  COFD[        1233] =   0.4848462442058104E+01;
  COFD[        1234] =  -0.4024597684582722E+00;
  COFD[        1235] =   0.1694041152472868E-01;
  COFD[        1236] =  -0.1789738457973026E+02;
  COFD[        1237] =   0.4914689305590001E+01;
  COFD[        1238] =  -0.4187174010126591E+00;
  COFD[        1239] =   0.1794961882909394E-01;
  COFD[        1240] =  -0.2031055243789758E+02;
  COFD[        1241] =   0.5168389671669797E+01;
  COFD[        1242] =  -0.4264383492219733E+00;
  COFD[        1243] =   0.1728448151830319E-01;
  COFD[        1244] =  -0.1806461227977115E+02;
  COFD[        1245] =   0.3730745706118281E+01;
  COFD[        1246] =  -0.1837097321884064E+00;
  COFD[        1247] =   0.4694035704567372E-02;
  COFD[        1248] =  -0.2034278270814375E+02;
  COFD[        1249] =   0.5194044053733033E+01;
  COFD[        1250] =  -0.4326906002865860E+00;
  COFD[        1251] =   0.1767025742696755E-01;
  COFD[        1252] =  -0.1951118321112224E+02;
  COFD[        1253] =   0.4471569363537695E+01;
  COFD[        1254] =  -0.3013374208218036E+00;
  COFD[        1255] =   0.1054895382382678E-01;
  COFD[        1256] =  -0.1895623695499981E+02;
  COFD[        1257] =   0.4957160985956500E+01;
  COFD[        1258] =  -0.4140469150912638E+00;
  COFD[        1259] =   0.1734089632858991E-01;
  COFD[        1260] =  -0.1279217728541477E+02;
  COFD[        1261] =   0.1325671264699948E+01;
  COFD[        1262] =   0.1805397215957534E+00;
  COFD[        1263] =  -0.1280666613762267E-01;
  COFD[        1264] =  -0.2029063758714785E+02;
  COFD[        1265] =   0.5183379095141047E+01;
  COFD[        1266] =  -0.4320078479547174E+00;
  COFD[        1267] =   0.1766855996768827E-01;
  COFD[        1268] =  -0.2001595054016060E+02;
  COFD[        1269] =   0.4853627069870797E+01;
  COFD[        1270] =  -0.3679409707902928E+00;
  COFD[        1271] =   0.1407678359913585E-01;
  COFD[        1272] =  -0.1726614977220465E+02;
  COFD[        1273] =   0.1754490020116512E+01;
  COFD[        1274] =   0.1044738525084192E+00;
  COFD[        1275] =  -0.8933019725221741E-02;
  COFD[        1276] =  -0.2026805290514954E+02;
  COFD[        1277] =   0.5158728558021261E+01;
  COFD[        1278] =  -0.4234667524907835E+00;
  COFD[        1279] =   0.1709481128969107E-01;
  COFD[        1280] =  -0.1641389233309519E+02;
  COFD[        1281] =   0.3843808661762178E+01;
  COFD[        1282] =  -0.2863465660410651E+00;
  COFD[        1283] =   0.1248720098494784E-01;
  COFD[        1284] =  -0.1632103748685788E+02;
  COFD[        1285] =   0.3805435800233515E+01;
  COFD[        1286] =  -0.2813339083459423E+00;
  COFD[        1287] =   0.1226875862421682E-01;
  COFD[        1288] =  -0.1686978847126216E+02;
  COFD[        1289] =   0.4721521983629889E+01;
  COFD[        1290] =  -0.4037192718464361E+00;
  COFD[        1291] =   0.1770807186574834E-01;
  COFD[        1292] =  -0.1450174966203873E+02;
  COFD[        1293] =   0.3240507460697169E+01;
  COFD[        1294] =  -0.2104181641466212E+00;
  COFD[        1295] =   0.9299811294433670E-02;
  COFD[        1296] =  -0.2025241161910038E+02;
  COFD[        1297] =   0.5172982615634614E+01;
  COFD[        1298] =  -0.4303950755374367E+00;
  COFD[        1299] =   0.1758759466858050E-01;
  COFD[        1300] =  -0.1450839204389819E+02;
  COFD[        1301] =   0.3251381640250417E+01;
  COFD[        1302] =  -0.2119042844536192E+00;
  COFD[        1303] =   0.9367271609371054E-02;
  COFD[        1304] =  -0.2000273357738514E+02;
  COFD[        1305] =   0.5003613670288710E+01;
  COFD[        1306] =  -0.4191992889668390E+00;
  COFD[        1307] =   0.1753189913176521E-01;
  COFD[        1308] =  -0.1523620653214519E+02;
  COFD[        1309] =   0.3348643253645854E+01;
  COFD[        1310] =  -0.2232496547504182E+00;
  COFD[        1311] =   0.9804212145714645E-02;
  COFD[        1312] =  -0.1427294640051003E+02;
  COFD[        1313] =   0.3117012895539554E+01;
  COFD[        1314] =  -0.1946596484063036E+00;
  COFD[        1315] =   0.8629721402195967E-02;
  COFD[        1316] =  -0.1215152187647601E+02;
  COFD[        1317] =   0.2670065262560549E+01;
  COFD[        1318] =  -0.1362026991746787E+00;
  COFD[        1319] =   0.6073761265771448E-02;
  COFD[        1320] =  -0.1555717721575701E+02;
  COFD[        1321] =   0.3482723600957817E+01;
  COFD[        1322] =  -0.2406978172870443E+00;
  COFD[        1323] =   0.1056168468865037E-01;
  COFD[        1324] =  -0.2001104497452950E+02;
  COFD[        1325] =   0.5003493097552491E+01;
  COFD[        1326] =  -0.4191652678118554E+00;
  COFD[        1327] =   0.1752960197527472E-01;
  COFD[        1328] =  -0.1524693443631203E+02;
  COFD[        1329] =   0.3356096228699615E+01;
  COFD[        1330] =  -0.2242170999002023E+00;
  COFD[        1331] =   0.9846114264084087E-02;
  COFD[        1332] =  -0.1817998188473039E+02;
  COFD[        1333] =   0.4364878900220307E+01;
  COFD[        1334] =  -0.3484381669569521E+00;
  COFD[        1335] =   0.1494878692100033E-01;
  COFD[        1336] =  -0.1453532015237173E+02;
  COFD[        1337] =   0.3291285383780231E+01;
  COFD[        1338] =  -0.2173571659474091E+00;
  COFD[        1339] =   0.9614779127177037E-02;
  COFD[        1340] =  -0.2029063758714785E+02;
  COFD[        1341] =   0.5183379095141047E+01;
  COFD[        1342] =  -0.4320078479547174E+00;
  COFD[        1343] =   0.1766855996768827E-01;
  COFD[        1344] =  -0.1521775053325088E+02;
  COFD[        1345] =   0.3347991675957199E+01;
  COFD[        1346] =  -0.2231616752417168E+00;
  COFD[        1347] =   0.9800264415023850E-02;
  COFD[        1348] =  -0.1753711041071091E+02;
  COFD[        1349] =   0.4166036237684116E+01;
  COFD[        1350] =  -0.3251274594092163E+00;
  COFD[        1351] =   0.1404131420270148E-01;
  COFD[        1352] =  -0.2363542175234041E+02;
  COFD[        1353] =   0.5050538606714960E+01;
  COFD[        1354] =  -0.4076609632033716E+00;
  COFD[        1355] =   0.1634285774662198E-01;
  COFD[        1356] =  -0.1641792325775208E+02;
  COFD[        1357] =   0.3854808256627624E+01;
  COFD[        1358] =  -0.2878185797680770E+00;
  COFD[        1359] =   0.1255269968409222E-01;
  COFD[        1360] =  -0.1868821899264003E+02;
  COFD[        1361] =   0.4623091902441412E+01;
  COFD[        1362] =  -0.3789790809161654E+00;
  COFD[        1363] =   0.1615045683658846E-01;
  COFD[        1364] =  -0.1862959184503705E+02;
  COFD[        1365] =   0.4599775876853535E+01;
  COFD[        1366] =  -0.3762918861611520E+00;
  COFD[        1367] =   0.1604808326150353E-01;
  COFD[        1368] =  -0.1772949944660123E+02;
  COFD[        1369] =   0.4898286384206769E+01;
  COFD[        1370] =  -0.4076359849112164E+00;
  COFD[        1371] =   0.1707858545899634E-01;
  COFD[        1372] =  -0.1697689861681293E+02;
  COFD[        1373] =   0.4143440294873640E+01;
  COFD[        1374] =  -0.3245402843230118E+00;
  COFD[        1375] =   0.1411020067125902E-01;
  COFD[        1376] =  -0.1994374085526542E+02;
  COFD[        1377] =   0.4830081980319046E+01;
  COFD[        1378] =  -0.3644950600153112E+00;
  COFD[        1379] =   0.1391192544233011E-01;
  COFD[        1380] =  -0.1698858354738882E+02;
  COFD[        1381] =   0.4157216748608500E+01;
  COFD[        1382] =  -0.3263536052448225E+00;
  COFD[        1383] =   0.1418959320103276E-01;
  COFD[        1384] =  -0.2090688825017299E+02;
  COFD[        1385] =   0.5106481690478173E+01;
  COFD[        1386] =  -0.4094432279607686E+00;
  COFD[        1387] =   0.1620751709314507E-01;
  COFD[        1388] =  -0.1753994518424459E+02;
  COFD[        1389] =   0.4158828151149182E+01;
  COFD[        1390] =  -0.3241998792087618E+00;
  COFD[        1391] =   0.1400163924604856E-01;
  COFD[        1392] =  -0.1672667778972758E+02;
  COFD[        1393] =   0.4023102094303094E+01;
  COFD[        1394] =  -0.3103163375902530E+00;
  COFD[        1395] =   0.1355332875778464E-01;
  COFD[        1396] =  -0.1437677954608083E+02;
  COFD[        1397] =   0.3542340190226259E+01;
  COFD[        1398] =  -0.2518343923430894E+00;
  COFD[        1399] =   0.1118661472548827E-01;
  COFD[        1400] =  -0.1777649791730375E+02;
  COFD[        1401] =   0.4248843195751910E+01;
  COFD[        1402] =  -0.3346930401290839E+00;
  COFD[        1403] =   0.1440673440796217E-01;
  COFD[        1404] =  -0.2093240510980354E+02;
  COFD[        1405] =   0.5113264124268066E+01;
  COFD[        1406] =  -0.4104425915775513E+00;
  COFD[        1407] =   0.1625565099449545E-01;
  COFD[        1408] =  -0.1755616817780421E+02;
  COFD[        1409] =   0.4169500620680274E+01;
  COFD[        1410] =  -0.3255016307159463E+00;
  COFD[        1411] =   0.1405440364613088E-01;
  COFD[        1412] =  -0.2026340088492126E+02;
  COFD[        1413] =   0.5009715621864962E+01;
  COFD[        1414] =  -0.4198847647706185E+00;
  COFD[        1415] =   0.1755705740818908E-01;
  COFD[        1416] =  -0.1702257488378831E+02;
  COFD[        1417] =   0.4202534175316763E+01;
  COFD[        1418] =  -0.3323169248630543E+00;
  COFD[        1419] =   0.1445062729566104E-01;
  COFD[        1420] =  -0.2001595054016060E+02;
  COFD[        1421] =   0.4853627069870797E+01;
  COFD[        1422] =  -0.3679409707902928E+00;
  COFD[        1423] =   0.1407678359913585E-01;
  COFD[        1424] =  -0.1753711041071091E+02;
  COFD[        1425] =   0.4166036237684116E+01;
  COFD[        1426] =  -0.3251274594092163E+00;
  COFD[        1427] =   0.1404131420270148E-01;
  COFD[        1428] =  -0.1969673206479491E+02;
  COFD[        1429] =   0.4855662192313339E+01;
  COFD[        1430] =  -0.4041439398948602E+00;
  COFD[        1431] =   0.1704289102344915E-01;
  COFD[        1432] =  -0.2337253695762064E+02;
  COFD[        1433] =   0.4682127905352519E+01;
  COFD[        1434] =  -0.3323278796888292E+00;
  COFD[        1435] =   0.1205004259254705E-01;
  COFD[        1436] =  -0.1867449917403642E+02;
  COFD[        1437] =   0.4627416013756130E+01;
  COFD[        1438] =  -0.3794324844557107E+00;
  COFD[        1439] =   0.1616528198979139E-01;
  COFD[        1440] =  -0.2380880735019050E+02;
  COFD[        1441] =   0.4985872967066852E+01;
  COFD[        1442] =  -0.3856996034210413E+00;
  COFD[        1443] =   0.1488146671087461E-01;
  COFD[        1444] =  -0.2384916272500252E+02;
  COFD[        1445] =   0.5010115663668988E+01;
  COFD[        1446] =  -0.3897718053261842E+00;
  COFD[        1447] =   0.1509403961940095E-01;
  COFD[        1448] =  -0.2378143430900525E+02;
  COFD[        1449] =   0.4975082596789566E+01;
  COFD[        1450] =  -0.3839197228971605E+00;
  COFD[        1451] =   0.1478949218777245E-01;
  COFD[        1452] =  -0.2349867378552522E+02;
  COFD[        1453] =   0.5062620989049528E+01;
  COFD[        1454] =  -0.4159679664653171E+00;
  COFD[        1455] =   0.1696493513352471E-01;
  COFD[        1456] =  -0.1726613964269955E+02;
  COFD[        1457] =   0.1754487534188974E+01;
  COFD[        1458] =   0.1044739320789194E+00;
  COFD[        1459] =  -0.8933007198223725E-02;
  COFD[        1460] =  -0.2349858517247871E+02;
  COFD[        1461] =   0.5062581035220192E+01;
  COFD[        1462] =  -0.4159619955775484E+00;
  COFD[        1463] =   0.1696464437402649E-01;
  COFD[        1464] =  -0.1871985392716656E+02;
  COFD[        1465] =   0.2406373385171115E+01;
  COFD[        1466] =   0.1073057268211112E-01;
  COFD[        1467] =  -0.4539455264233341E-02;
  COFD[        1468] =  -0.2363569469280461E+02;
  COFD[        1469] =   0.5050661168714743E+01;
  COFD[        1470] =  -0.4076790606763266E+00;
  COFD[        1471] =   0.1634373025964126E-01;
  COFD[        1472] =  -0.2328225161733864E+02;
  COFD[        1473] =   0.5008847120014931E+01;
  COFD[        1474] =  -0.4122620203443967E+00;
  COFD[        1475] =   0.1693220599831681E-01;
  COFD[        1476] =  -0.2177638593119188E+02;
  COFD[        1477] =   0.4563067342827384E+01;
  COFD[        1478] =  -0.3686945303206421E+00;
  COFD[        1479] =   0.1561188821007837E-01;
  COFD[        1480] =  -0.2370839244124438E+02;
  COFD[        1481] =   0.5048027935648303E+01;
  COFD[        1482] =  -0.4040835196333057E+00;
  COFD[        1483] =   0.1606332363762197E-01;
  COFD[        1484] =  -0.1871985825163239E+02;
  COFD[        1485] =   0.2406372669599885E+01;
  COFD[        1486] =   0.1073111423051258E-01;
  COFD[        1487] =  -0.4539506957064400E-02;
  COFD[        1488] =  -0.2363377195653007E+02;
  COFD[        1489] =   0.5047702307335338E+01;
  COFD[        1490] =  -0.4070594187739618E+00;
  COFD[        1491] =   0.1630755566219766E-01;
  COFD[        1492] =  -0.2289906109626653E+02;
  COFD[        1493] =   0.4409354958846299E+01;
  COFD[        1494] =  -0.2878825260320446E+00;
  COFD[        1495] =   0.9793621866625414E-02;
  COFD[        1496] =  -0.2349830619682465E+02;
  COFD[        1497] =   0.5062455579552836E+01;
  COFD[        1498] =  -0.4159432465895745E+00;
  COFD[        1499] =   0.1696373135599558E-01;
  COFD[        1500] =  -0.1726614977220465E+02;
  COFD[        1501] =   0.1754490020116512E+01;
  COFD[        1502] =   0.1044738525084192E+00;
  COFD[        1503] =  -0.8933019725221741E-02;
  COFD[        1504] =  -0.2363542175234041E+02;
  COFD[        1505] =   0.5050538606714960E+01;
  COFD[        1506] =  -0.4076609632033716E+00;
  COFD[        1507] =   0.1634285774662198E-01;
  COFD[        1508] =  -0.2337253695762064E+02;
  COFD[        1509] =   0.4682127905352519E+01;
  COFD[        1510] =  -0.3323278796888292E+00;
  COFD[        1511] =   0.1205004259254705E-01;
  COFD[        1512] =  -0.1462640990483345E+02;
  COFD[        1513] =   0.1200660129297300E+00;
  COFD[        1514] =   0.3363461583567012E+00;
  COFD[        1515] =  -0.1963938399449287E-01;
  COFD[        1516] =  -0.2380856690199094E+02;
  COFD[        1517] =   0.4985767450175743E+01;
  COFD[        1518] =  -0.3856843293132689E+00;
  COFD[        1519] =   0.1488074438089647E-01;
  COFD[        1520] =  -0.1729871483444523E+02;
  COFD[        1521] =   0.4183622454402964E+01;
  COFD[        1522] =  -0.3265686107794125E+00;
  COFD[        1523] =   0.1406921012900881E-01;
  COFD[        1524] =  -0.1725286194864467E+02;
  COFD[        1525] =   0.4167814533234269E+01;
  COFD[        1526] =  -0.3247525848550428E+00;
  COFD[        1527] =   0.1400034963858472E-01;
  COFD[        1528] =  -0.1749395439646448E+02;
  COFD[        1529] =   0.4865192761077812E+01;
  COFD[        1530] =  -0.4135081828421872E+00;
  COFD[        1531] =   0.1775619961599839E-01;
  COFD[        1532] =  -0.1541179728688726E+02;
  COFD[        1533] =   0.3612468014590696E+01;
  COFD[        1534] =  -0.2574709286409149E+00;
  COFD[        1535] =   0.1128548676893656E-01;
  COFD[        1536] =  -0.2025720433242960E+02;
  COFD[        1537] =   0.5159185593402682E+01;
  COFD[        1538] =  -0.4235430393599882E+00;
  COFD[        1539] =   0.1709888299272256E-01;
  COFD[        1540] =  -0.1540309090870326E+02;
  COFD[        1541] =   0.3614696445371578E+01;
  COFD[        1542] =  -0.2577753207878676E+00;
  COFD[        1543] =   0.1129930080443808E-01;
  COFD[        1544] =  -0.2042648432871773E+02;
  COFD[        1545] =   0.5112800688990101E+01;
  COFD[        1546] =  -0.4208144038227885E+00;
  COFD[        1547] =   0.1710457743576905E-01;
  COFD[        1548] =  -0.1645726609264868E+02;
  COFD[        1549] =   0.3866322991209679E+01;
  COFD[        1550] =  -0.2893595968515034E+00;
  COFD[        1551] =   0.1262127247206273E-01;
  COFD[        1552] =  -0.1512111139020410E+02;
  COFD[        1553] =   0.3458819295912228E+01;
  COFD[        1554] =  -0.2375970826594765E+00;
  COFD[        1555] =   0.1042745102315496E-01;
  COFD[        1556] =  -0.1320038983011787E+02;
  COFD[        1557] =   0.3068261259040047E+01;
  COFD[        1558] =  -0.1894129351750649E+00;
  COFD[        1559] =   0.8447150657362242E-02;
  COFD[        1560] =  -0.1672939192904382E+02;
  COFD[        1561] =   0.3975404184226489E+01;
  COFD[        1562] =  -0.3025898020717506E+00;
  COFD[        1563] =   0.1315505412218577E-01;
  COFD[        1564] =  -0.2041541490829319E+02;
  COFD[        1565] =   0.5105115933032721E+01;
  COFD[        1566] =  -0.4196478966186132E+00;
  COFD[        1567] =   0.1704696721515930E-01;
  COFD[        1568] =  -0.1645625614130163E+02;
  COFD[        1569] =   0.3866580965342516E+01;
  COFD[        1570] =  -0.2893531456511056E+00;
  COFD[        1571] =   0.1261944474539078E-01;
  COFD[        1572] =  -0.1921093934794400E+02;
  COFD[        1573] =   0.4770560824446395E+01;
  COFD[        1574] =  -0.3938285475558091E+00;
  COFD[        1575] =   0.1661810137086212E-01;
  COFD[        1576] =  -0.1539426704211137E+02;
  COFD[        1577] =   0.3632170199052923E+01;
  COFD[        1578] =  -0.2601531779068780E+00;
  COFD[        1579] =   0.1140681634760632E-01;
  COFD[        1580] =  -0.2026805290514954E+02;
  COFD[        1581] =   0.5158728558021261E+01;
  COFD[        1582] =  -0.4234667524907835E+00;
  COFD[        1583] =   0.1709481128969107E-01;
  COFD[        1584] =  -0.1641792325775208E+02;
  COFD[        1585] =   0.3854808256627624E+01;
  COFD[        1586] =  -0.2878185797680770E+00;
  COFD[        1587] =   0.1255269968409222E-01;
  COFD[        1588] =  -0.1867449917403642E+02;
  COFD[        1589] =   0.4627416013756130E+01;
  COFD[        1590] =  -0.3794324844557107E+00;
  COFD[        1591] =   0.1616528198979139E-01;
  COFD[        1592] =  -0.2380856690199094E+02;
  COFD[        1593] =   0.4985767450175743E+01;
  COFD[        1594] =  -0.3856843293132689E+00;
  COFD[        1595] =   0.1488074438089647E-01;
  COFD[        1596] =  -0.1728008170220198E+02;
  COFD[        1597] =   0.4182999326187858E+01;
  COFD[        1598] =  -0.3264888667764370E+00;
  COFD[        1599] =   0.1406581894032641E-01;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetKTDIF EGTRANSETKTDIF
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetKTDIF egtransetktdif
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetKTDIF egtransetktdif_
#endif
void egtransetKTDIF(int* KTDIF) {
  KTDIF[           0] =            3;
  KTDIF[           1] =           10;
  KTDIF[           2] =           19;
};
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFTD EGTRANSETCOFTD
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFTD egtransetcoftd
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFTD egtransetcoftd_
#endif
void egtransetCOFTD(double* COFTD) {
  COFTD[           0] =   0.1335613390740672E+00;
  COFTD[           1] =   0.5427093553867389E-03;
  COFTD[           2] =  -0.2560337786112679E-06;
  COFTD[           3] =   0.3734108235989047E-10;
  COFTD[           4] =   0.1388124311125256E+00;
  COFTD[           5] =   0.5423100996522542E-03;
  COFTD[           6] =  -0.2563216336337905E-06;
  COFTD[           7] =   0.3743074786658457E-10;
  COFTD[           8] =   0.0000000000000000E+00;
  COFTD[           9] =   0.0000000000000000E+00;
  COFTD[          10] =   0.0000000000000000E+00;
  COFTD[          11] =   0.0000000000000000E+00;
  COFTD[          12] =   0.2653218947985116E+00;
  COFTD[          13] =   0.3764127896098151E-03;
  COFTD[          14] =  -0.1883983304184878E-06;
  COFTD[          15] =   0.2872283799691947E-10;
  COFTD[          16] =  -0.1446549002209399E+00;
  COFTD[          17] =   0.7715674487897902E-03;
  COFTD[          18] =  -0.3091472931916318E-06;
  COFTD[          19] =   0.4070175242598579E-10;
  COFTD[          20] =   0.2633405620645734E+00;
  COFTD[          21] =   0.3736018682492041E-03;
  COFTD[          22] =  -0.1869914364289777E-06;
  COFTD[          23] =   0.2850834571320456E-10;
  COFTD[          24] =  -0.1281798748010575E+00;
  COFTD[          25] =   0.8029074901776178E-03;
  COFTD[          26] =  -0.3283082526701463E-06;
  COFTD[          27] =   0.4376878388365953E-10;
  COFTD[          28] =   0.2347368576400100E+00;
  COFTD[          29] =   0.4604164702501780E-03;
  COFTD[          30] =  -0.2259286218303810E-06;
  COFTD[          31] =   0.3390102262015666E-10;
  COFTD[          32] =   0.2814197988499757E+00;
  COFTD[          33] =   0.3326523270193756E-03;
  COFTD[          34] =  -0.1683849808056114E-06;
  COFTD[          35] =   0.2592454220086897E-10;
  COFTD[          36] =   0.1504482931689058E+00;
  COFTD[          37] =   0.5840856631801511E-04;
  COFTD[          38] =  -0.3129515477212745E-07;
  COFTD[          39] =   0.5176285460253215E-11;
  COFTD[          40] =   0.2132882440000961E+00;
  COFTD[          41] =   0.4931598944385875E-03;
  COFTD[          42] =  -0.2396666518582280E-06;
  COFTD[          43] =   0.3570013941021640E-10;
  COFTD[          44] =  -0.1284796695276010E+00;
  COFTD[          45] =   0.8047853780428643E-03;
  COFTD[          46] =  -0.3290761195680041E-06;
  COFTD[          47] =   0.4387115292260514E-10;
  COFTD[          48] =   0.2322810318100867E+00;
  COFTD[          49] =   0.4600739597951276E-03;
  COFTD[          50] =  -0.2256287885029010E-06;
  COFTD[          51] =   0.3384082237165661E-10;
  COFTD[          52] =   0.1823923694713997E-01;
  COFTD[          53] =   0.7294390609533908E-03;
  COFTD[          54] =  -0.3255243929397424E-06;
  COFTD[          55] =   0.4579013475655006E-10;
  COFTD[          56] =   0.2558166843155990E+00;
  COFTD[          57] =   0.3629277253771079E-03;
  COFTD[          58] =  -0.1816489221700018E-06;
  COFTD[          59] =   0.2769383652293768E-10;
  COFTD[          60] =  -0.1455172589046899E+00;
  COFTD[          61] =   0.7761671401140823E-03;
  COFTD[          62] =  -0.3109902715659168E-06;
  COFTD[          63] =   0.4094439549991396E-10;
  COFTD[          64] =   0.2336163545989917E+00;
  COFTD[          65] =   0.4582186984122274E-03;
  COFTD[          66] =  -0.2248501644020966E-06;
  COFTD[          67] =   0.3373919801654980E-10;
  COFTD[          68] =   0.6581709312715994E-01;
  COFTD[          69] =   0.6846199844611079E-03;
  COFTD[          70] =  -0.3121490853294717E-06;
  COFTD[          71] =   0.4450966517984661E-10;
  COFTD[          72] =   0.2156069967211590E+00;
  COFTD[          73] =  -0.8350733530287657E-03;
  COFTD[          74] =   0.3124142341339501E-06;
  COFTD[          75] =  -0.3937676136436095E-10;
  COFTD[          76] =   0.1322744734933614E+00;
  COFTD[          77] =   0.5374803423009470E-03;
  COFTD[          78] =  -0.2535668891694801E-06;
  COFTD[          79] =   0.3698130045018236E-10;
  COFTD[          80] =   0.3455103259336228E+00;
  COFTD[          81] =   0.1321913114217400E-03;
  COFTD[          82] =  -0.7086012360787731E-07;
  COFTD[          83] =   0.1173303990698333E-10;
  COFTD[          84] =   0.3538478885885439E+00;
  COFTD[          85] =   0.1302119616495474E-03;
  COFTD[          86] =  -0.6988019757569516E-07;
  COFTD[          87] =   0.1160423794907823E-10;
  COFTD[          88] =  -0.1504482931689058E+00;
  COFTD[          89] =  -0.5840856631801511E-04;
  COFTD[          90] =   0.3129515477212745E-07;
  COFTD[          91] =  -0.5176285460253215E-11;
  COFTD[          92] =   0.4184663328245823E+00;
  COFTD[          93] =   0.1740049929096906E-04;
  COFTD[          94] =  -0.7266155663582315E-08;
  COFTD[          95] =   0.1668860961881446E-11;
  COFTD[          96] =   0.5463813258525711E-01;
  COFTD[          97] =   0.5716146453946012E-03;
  COFTD[          98] =  -0.2605856959213780E-06;
  COFTD[          99] =   0.3715360167358631E-10;
  COFTD[         100] =   0.4121693838336929E+00;
  COFTD[         101] =   0.1713866208243755E-04;
  COFTD[         102] =  -0.7156816851868155E-08;
  COFTD[         103] =   0.1643748470085703E-11;
  COFTD[         104] =   0.1547978351031708E+00;
  COFTD[         105] =   0.5129939480053735E-03;
  COFTD[         106] =  -0.2444497170195778E-06;
  COFTD[         107] =   0.3589955716843627E-10;
  COFTD[         108] =   0.4451226233123079E+00;
  COFTD[         109] =   0.5795895813636342E-04;
  COFTD[         110] =  -0.3071051601810868E-07;
  COFTD[         111] =   0.5599659330462059E-11;
  COFTD[         112] =   0.4050684083869960E+00;
  COFTD[         113] =   0.2269701837029273E-06;
  COFTD[         114] =   0.3185217494944040E-08;
  COFTD[         115] =  -0.1597162014516924E-12;
  COFTD[         116] =   0.0000000000000000E+00;
  COFTD[         117] =   0.0000000000000000E+00;
  COFTD[         118] =   0.0000000000000000E+00;
  COFTD[         119] =   0.0000000000000000E+00;
  COFTD[         120] =   0.4379907438872362E+00;
  COFTD[         121] =   0.7975216050713585E-04;
  COFTD[         122] =  -0.4292321594073912E-07;
  COFTD[         123] =   0.7558155922052741E-11;
  COFTD[         124] =   0.1555253375623714E+00;
  COFTD[         125] =   0.5154048625926477E-03;
  COFTD[         126] =  -0.2455985558916718E-06;
  COFTD[         127] =   0.3606827410240905E-10;
  COFTD[         128] =   0.4402429418906190E+00;
  COFTD[         129] =   0.5862828974396131E-04;
  COFTD[         130] =  -0.3111220911258768E-07;
  COFTD[         131] =   0.5658391439418446E-11;
  COFTD[         132] =   0.3190608345439727E+00;
  COFTD[         133] =   0.3157015231693243E-03;
  COFTD[         134] =  -0.1615480625622576E-06;
  COFTD[         135] =   0.2512460309342888E-10;
  COFTD[         136] =   0.3886123540069732E+00;
  COFTD[         137] =   0.1615912311204916E-04;
  COFTD[         138] =  -0.6747777862907584E-08;
  COFTD[         139] =   0.1549802065388507E-11;
  COFTD[         140] =   0.5529744609569107E-01;
  COFTD[         141] =   0.5785122687327116E-03;
  COFTD[         142] =  -0.2637301604522408E-06;
  COFTD[         143] =   0.3760193089688968E-10;
  COFTD[         144] =   0.4408678331086748E+00;
  COFTD[         145] =   0.5740494628799610E-04;
  COFTD[         146] =  -0.3041696364431554E-07;
  COFTD[         147] =   0.5546133909791147E-11;
  COFTD[         148] =   0.3587095749289513E+00;
  COFTD[         149] =   0.2549270192949080E-03;
  COFTD[         150] =  -0.1329645563167186E-06;
  COFTD[         151] =   0.2110100942106591E-10;
  COFTD[         152] =  -0.4163666475409064E-01;
  COFTD[         153] =  -0.7420443513816751E-03;
  COFTD[         154] =   0.3344525746340108E-06;
  COFTD[         155] =  -0.4734076628212864E-10;
  COFTD[         156] =   0.3387875572167290E+00;
  COFTD[         157] =   0.1296191983867103E-03;
  COFTD[         158] =  -0.6948136243488190E-07;
  COFTD[         159] =   0.1150474423035603E-10;
  COFTD[         160] =  -0.2150948093190753E+00;
  COFTD[         161] =   0.8368367635659997E-03;
  COFTD[         162] =  -0.3135080334769803E-06;
  COFTD[         163] =   0.3955057076088332E-10;
  COFTD[         164] =  -0.2131316243913729E+00;
  COFTD[         165] =   0.8391730916389774E-03;
  COFTD[         166] =  -0.3155200110876533E-06;
  COFTD[         167] =   0.3989852631472832E-10;
  COFTD[         168] =  -0.2156069967211590E+00;
  COFTD[         169] =   0.8350733530287657E-03;
  COFTD[         170] =  -0.3124142341339501E-06;
  COFTD[         171] =   0.3937676136436095E-10;
  COFTD[         172] =  -0.1261195350944164E+00;
  COFTD[         173] =   0.8578171176710015E-03;
  COFTD[         174] =  -0.3536931711268338E-06;
  COFTD[         175] =   0.4739398779251787E-10;
  COFTD[         176] =  -0.1742730942128402E+00;
  COFTD[         177] =   0.3444077011559843E-03;
  COFTD[         178] =  -0.3417522795042020E-07;
  COFTD[         179] =  -0.3052872895702647E-11;
  COFTD[         180] =  -0.1261190225204061E+00;
  COFTD[         181] =   0.8578136313375028E-03;
  COFTD[         182] =  -0.3536917336498642E-06;
  COFTD[         183] =   0.4739379517425008E-10;
  COFTD[         184] =  -0.1967303775760090E+00;
  COFTD[         185] =   0.4203106279105429E-03;
  COFTD[         186] =  -0.7191069880138451E-07;
  COFTD[         187] =   0.2393238182040283E-11;
  COFTD[         188] =  -0.1616001861800044E+00;
  COFTD[         189] =   0.8629977674683832E-03;
  COFTD[         190] =  -0.3458426863256669E-06;
  COFTD[         191] =   0.4553802095640965E-10;
  COFTD[         192] =  -0.1036278124056820E+00;
  COFTD[         193] =   0.8495612640731872E-03;
  COFTD[         194] =  -0.3556336745081299E-06;
  COFTD[         195] =   0.4809664985510671E-10;
  COFTD[         196] =   0.4163666475409064E-01;
  COFTD[         197] =   0.7420443513816751E-03;
  COFTD[         198] =  -0.3344525746340108E-06;
  COFTD[         199] =   0.4734076628212864E-10;
  COFTD[         200] =  -0.1769715465915445E+00;
  COFTD[         201] =   0.8611073025514351E-03;
  COFTD[         202] =  -0.3399354853976130E-06;
  COFTD[         203] =   0.4434361953551270E-10;
  COFTD[         204] =  -0.1967306274104651E+00;
  COFTD[         205] =   0.4203111616770018E-03;
  COFTD[         206] =  -0.7191079012317498E-07;
  COFTD[         207] =   0.2393241221293148E-11;
  COFTD[         208] =  -0.1625657627376637E+00;
  COFTD[         209] =   0.8629643703710789E-03;
  COFTD[         210] =  -0.3455235737063456E-06;
  COFTD[         211] =   0.4547109515793038E-10;
  COFTD[         212] =  -0.2457462032825698E+00;
  COFTD[         213] =   0.7209815038571562E-03;
  COFTD[         214] =  -0.2347526386964109E-06;
  COFTD[         215] =   0.2676902670472271E-10;
  COFTD[         216] =  -0.1261170427389316E+00;
  COFTD[         217] =   0.8578001656167727E-03;
  COFTD[         218] =  -0.3536861814950185E-06;
  COFTD[         219] =   0.4739305119958373E-10;
  COFTD[         220] =  -0.1742736562882112E+00;
  COFTD[         221] =   0.3444088119590448E-03;
  COFTD[         222] =  -0.3417533817428469E-07;
  COFTD[         223] =  -0.3052882741999780E-11;
  COFTD[         224] =  -0.1615997658635892E+00;
  COFTD[         225] =   0.8629955228414649E-03;
  COFTD[         226] =  -0.3458417868009724E-06;
  COFTD[         227] =   0.4553790251361466E-10;
  COFTD[         228] =  -0.2403972795835718E+00;
  COFTD[         229] =   0.7756828293052003E-03;
  COFTD[         230] =  -0.2690332887763941E-06;
  COFTD[         231] =   0.3218675475978092E-10;
  COFTD[         232] =   0.0000000000000000E+00;
  COFTD[         233] =   0.0000000000000000E+00;
  COFTD[         234] =   0.0000000000000000E+00;
  COFTD[         235] =   0.0000000000000000E+00;
  COFTD[         236] =  -0.2150936814154548E+00;
  COFTD[         237] =   0.8368323754023877E-03;
  COFTD[         238] =  -0.3135063895189177E-06;
  COFTD[         239] =   0.3955036336753852E-10;
};




#if 0




\\
\\
\\  This is the mechanism file
\\
\\
ELEMENTS
O  H  C  N   E
END
SPECIES
CH3     CH4     H         OH     H2O     O     HCOp    NO   N
H2      O2      CH2O      CO     CO2     CH    H3Op    N2   NO2
E	CH2	
!CH(S)	C2H3Op	C2H2	CH3p	C3H3p	C
END
!THERMO
! Insert GRI-Mech thermodynamics here or use in default file
!END
REACTIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!31 react.Jing Hu et al. (Experimental Thermal and Fluid Science 21, 200)124-133!!!
!!!!!!!!!M coefficients taken from GRI-Mech 3.0!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CH4<=>CH3+H				4.23E+15     .000        108.69E+3
CH4+OH<=>CH3+H2O	                2.00E+14     .000          8.41E+3
CH4+O<=>CH3+OH          		3.48E+13     .000          8.38E+3
CH4+H<=>CH3+H2                      	4.35E+14     .000         13.74E+3
CH3+O2<=>CH2O+OH                     	5.29E+11     .000          1.70E+3
CH2O+OH<=>CO+H2O+H                   	5.87E+14     .000          4.88E+3
CO+OH=>CO2+H                        	1.30E+12     .000          1.53E+3
CO2+H=>CO+OH                        	1.45E+14     .000         23.76E+3
O2+H=>OH+O                          	2.24E+14     .000         16.80E+3
OH+O=>O2+H                          	1.71E+13     .000          0.87E+3
O+H2=>OH+H                          	1.74E+13     .000          9.45E+3
OH+H=>O+H2                          	7.70E+12     .000          7.58E+3
O+H2O=>2OH                          	5.75E+13     .000         18.10E+3
2OH=>O+H2O                          	5.38E+12     .000          1.05E+3
OH+H2=>H2O+H                        	2.19E+13     .000          5.15E+3
H2O+H=>H2+OH                        	8.41E+13     .000         20.10E+3
H+OH+M<=>H2O+M                    	2.00E+19   -1.000           .00
H2/ .73/ H2O/3.65/ CH4/2.00/ 
O+O+M<=>O2+M                        	8.90E+14   -0.500           .00
H2/ 2.40/ H2O/15.40/ CH4/ 2.00/ CO/ 1.75/ CO2/ 3.60/
H+H+M<=>H2+M				1.00E+18   -1.000           .00	
H2/ .00/ H2O/ .00/ CH4/2.00/ CO2/ .00/ 		
CH3+O<=>CH+H2O				2.80E+08     .000	    .00
CH+O<=>HCOp+E				5.75E+11     .000          6.00E+3
HCOp+H2O<=>CO+H3Op			5.02E+17     .000 	  24.00E+3
H3Op+E<=>H2O+H				1.44E+17     .000           .00
CH+O2<=>CO+OH				6.00E+10     .000           .00
O+N2=>NO+N				1.36E+14     .000         75.40E+3
NO+N=>O+N2				3.10E+13     .000          0.33E+3
N+O2=>NO+O				6.43E+09    1.000          6.25E+3
NO+O=>N+O2				1.55E+09    1.000	  38.64E+3
NO+O+M<=>NO2+M				1.05E+15     .000         -1.87E+3
H2/2.00/ H2O/6.00/ CH4/2.00/ CO/1.50/ CO2/2.00/ 				
NO2+O<=>NO+O2                            2.10E+12     .000           .00
NO2+H<=>NO+OH				3.00E+14     .000	    .00
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!11 react.Pedersen and Brown (Combustion and Flame 94 (1993) 433-448)!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CH+O<=>HCOp+E				2.52E+11     .000       1700.00
!HCOp+H2O<=>CO+H3Op			1.00E+16   -0.100           .00
!H3Op+E<=>H2O+H				2.29E+18   -0.500           .00
!CH(S)+O<=>HCOp+E			5.04E+14     .000       1700.00
!H3Op+C2H2<=>C2H3Op+H2			8.39E+15     .000           .00
!HCOp+CH2<=>CH3p+CO			5.62E+14     .000           .00
!H3Op+CH2<=>CH3p+H2O			6.17E+14     .000           .00 
!CH3p+C2H2<=>C3H3p+H2			7.24E+14     .000	    .00		
!C3H3p+H2O<=>C2H3Op+CH2			7.24E+14     .000           .00
!CH3p+CO2<=>C2H3Op+O			7.24E+14     .000           .00
!CH3p+E<=>CH2+H				2.29E+18   -0.500           .00
END

\\
\\
\\  This is the therm file
\\
\\
THERMO ALL
   200.000  1000.000  5000.000
! GRI-Mech Version 3.0 Thermodynamics released 7/30/99
! NASA Polynomial format for CHEMKIN-II
! see README file for disclaimer
O                 L 1/90O   1               G   200.000  3500.000  1000.000    1
 2.56942078E+00-8.59741137E-05 4.19484589E-08-1.00177799E-11 1.22833691E-15    2
 2.92175791E+04 4.78433864E+00 3.16826710E+00-3.27931884E-03 6.64306396E-06    3
-6.12806624E-09 2.11265971E-12 2.91222592E+04 2.05193346E+00                   4
O2                TPIS89O   2               G   200.000  3500.000  1000.000    1
 3.28253784E+00 1.48308754E-03-7.57966669E-07 2.09470555E-10-2.16717794E-14    2
-1.08845772E+03 5.45323129E+00 3.78245636E+00-2.99673416E-03 9.84730201E-06    3
-9.68129509E-09 3.24372837E-12-1.06394356E+03 3.65767573E+00                   4
H                 L 7/88H   1               G   200.000  3500.000  1000.000    1
 2.50000001E+00-2.30842973E-11 1.61561948E-14-4.73515235E-18 4.98197357E-22    2
 2.54736599E+04-4.46682914E-01 2.50000000E+00 7.05332819E-13-1.99591964E-15    3
 2.30081632E-18-9.27732332E-22 2.54736599E+04-4.46682853E-01                   4
H2                TPIS78H   2               G   200.000  3500.000  1000.000    1
 3.33727920E+00-4.94024731E-05 4.99456778E-07-1.79566394E-10 2.00255376E-14    2
-9.50158922E+02-3.20502331E+00 2.34433112E+00 7.98052075E-03-1.94781510E-05    3
 2.01572094E-08-7.37611761E-12-9.17935173E+02 6.83010238E-01                   4
OH                RUS 78O   1H   1          G   200.000  3500.000  1000.000    1
 3.09288767E+00 5.48429716E-04 1.26505228E-07-8.79461556E-11 1.17412376E-14    2
 3.85865700E+03 4.47669610E+00 3.99201543E+00-2.40131752E-03 4.61793841E-06    3
-3.88113333E-09 1.36411470E-12 3.61508056E+03-1.03925458E-01                   4
H2O               L 8/89H   2O   1          G   200.000  3500.000  1000.000    1
 3.03399249E+00 2.17691804E-03-1.64072518E-07-9.70419870E-11 1.68200992E-14    2
-3.00042971E+04 4.96677010E+00 4.19864056E+00-2.03643410E-03 6.52040211E-06    3
-5.48797062E-09 1.77197817E-12-3.02937267E+04-8.49032208E-01                   4
HO2               L 5/89H   1O   2          G   200.000  3500.000  1000.000    1
 4.01721090E+00 2.23982013E-03-6.33658150E-07 1.14246370E-10-1.07908535E-14    2
 1.11856713E+02 3.78510215E+00 4.30179801E+00-4.74912051E-03 2.11582891E-05    3
-2.42763894E-08 9.29225124E-12 2.94808040E+02 3.71666245E+00                   4
H2O2              L 7/88H   2O   2          G   200.000  3500.000  1000.000    1
 4.16500285E+00 4.90831694E-03-1.90139225E-06 3.71185986E-10-2.87908305E-14    2
-1.78617877E+04 2.91615662E+00 4.27611269E+00-5.42822417E-04 1.67335701E-05    3
-2.15770813E-08 8.62454363E-12-1.77025821E+04 3.43505074E+00                   4
C                 L11/88C   1               G   200.000  3500.000  1000.000    1
 2.49266888E+00 4.79889284E-05-7.24335020E-08 3.74291029E-11-4.87277893E-15    2
 8.54512953E+04 4.80150373E+00 2.55423955E+00-3.21537724E-04 7.33792245E-07    3
-7.32234889E-10 2.66521446E-13 8.54438832E+04 4.53130848E+00                   4
CH                TPIS79C   1H   1          G   200.000  3500.000  1000.000    1
 2.87846473E+00 9.70913681E-04 1.44445655E-07-1.30687849E-10 1.76079383E-14    2
 7.10124364E+04 5.48497999E+00 3.48981665E+00 3.23835541E-04-1.68899065E-06    3
 3.16217327E-09-1.40609067E-12 7.07972934E+04 2.08401108E+00                   4
CH2               L S/93C   1H   2          G   200.000  3500.000  1000.000    1
 2.87410113E+00 3.65639292E-03-1.40894597E-06 2.60179549E-10-1.87727567E-14    2
 4.62636040E+04 6.17119324E+00 3.76267867E+00 9.68872143E-04 2.79489841E-06    3
-3.85091153E-09 1.68741719E-12 4.60040401E+04 1.56253185E+00                   4
CH2(S)            L S/93C   1H   2          G   200.000  3500.000  1000.000    1
 2.29203842E+00 4.65588637E-03-2.01191947E-06 4.17906000E-10-3.39716365E-14    2
 5.09259997E+04 8.62650169E+00 4.19860411E+00-2.36661419E-03 8.23296220E-06    3
-6.68815981E-09 1.94314737E-12 5.04968163E+04-7.69118967E-01                   4
CH3               L11/89C   1H   3          G   200.000  3500.000  1000.000    1
 2.28571772E+00 7.23990037E-03-2.98714348E-06 5.95684644E-10-4.67154394E-14    2
 1.67755843E+04 8.48007179E+00 3.67359040E+00 2.01095175E-03 5.73021856E-06    3
-6.87117425E-09 2.54385734E-12 1.64449988E+04 1.60456433E+00                   4
CH4               L 8/88C   1H   4          G   200.000  3500.000  1000.000    1
 7.48514950E-02 1.33909467E-02-5.73285809E-06 1.22292535E-09-1.01815230E-13    2
-9.46834459E+03 1.84373180E+01 5.14987613E+00-1.36709788E-02 4.91800599E-05    3
-4.84743026E-08 1.66693956E-11-1.02466476E+04-4.64130376E+00                   4
CO                TPIS79C   1O   1          G   200.000  3500.000  1000.000    1
 2.71518561E+00 2.06252743E-03-9.98825771E-07 2.30053008E-10-2.03647716E-14    2
-1.41518724E+04 7.81868772E+00 3.57953347E+00-6.10353680E-04 1.01681433E-06    3
 9.07005884E-10-9.04424499E-13-1.43440860E+04 3.50840928E+00                   4
CO2               L 7/88C   1O   2          G   200.000  3500.000  1000.000    1
 3.85746029E+00 4.41437026E-03-2.21481404E-06 5.23490188E-10-4.72084164E-14    2
-4.87591660E+04 2.27163806E+00 2.35677352E+00 8.98459677E-03-7.12356269E-06    3
 2.45919022E-09-1.43699548E-13-4.83719697E+04 9.90105222E+00                   4
HCO               L12/89H   1C   1O   1     G   200.000  3500.000  1000.000    1
 2.77217438E+00 4.95695526E-03-2.48445613E-06 5.89161778E-10-5.33508711E-14    2
 4.01191815E+03 9.79834492E+00 4.22118584E+00-3.24392532E-03 1.37799446E-05    3
-1.33144093E-08 4.33768865E-12 3.83956496E+03 3.39437243E+00                   4
CH2O              L 8/88H   2C   1O   1     G   200.000  3500.000  1000.000    1
 1.76069008E+00 9.20000082E-03-4.42258813E-06 1.00641212E-09-8.83855640E-14    2
-1.39958323E+04 1.36563230E+01 4.79372315E+00-9.90833369E-03 3.73220008E-05    3
-3.79285261E-08 1.31772652E-11-1.43089567E+04 6.02812900E-01                   4
CH2OH             GUNL93C   1H   3O   1     G   200.000  3500.000  1000.000    1
 3.69266569E+00 8.64576797E-03-3.75101120E-06 7.87234636E-10-6.48554201E-14    2
-3.24250627E+03 5.81043215E+00 3.86388918E+00 5.59672304E-03 5.93271791E-06    3
-1.04532012E-08 4.36967278E-12-3.19391367E+03 5.47302243E+00                   4
CH3O              121686C   1H   3O   1     G   290.00   3000.00   1000.000    1
 0.03770799E+02 0.07871497E-01-0.02656384E-04 0.03944431E-08-0.02112616E-12    2
 0.12783252E+03 0.02929575E+02 0.02106204E+02 0.07216595E-01 0.05338472E-04    3
-0.07377636E-07 0.02075610E-10 0.09786011E+04 0.13152177E+02                   4
CH3OH             L 8/88C   1H   4O   1     G   200.000  3500.000  1000.000    1
 1.78970791E+00 1.40938292E-02-6.36500835E-06 1.38171085E-09-1.17060220E-13    2
-2.53748747E+04 1.45023623E+01 5.71539582E+00-1.52309129E-02 6.52441155E-05    3
-7.10806889E-08 2.61352698E-11-2.56427656E+04-1.50409823E+00                   4
C2H               L 1/91C   2H   1          G   200.000  3500.000  1000.000    1
 3.16780652E+00 4.75221902E-03-1.83787077E-06 3.04190252E-10-1.77232770E-14    2
 6.71210650E+04 6.63589475E+00 2.88965733E+00 1.34099611E-02-2.84769501E-05    3
 2.94791045E-08-1.09331511E-11 6.68393932E+04 6.22296438E+00                   4
C2H2              L 1/91C   2H   2          G   200.000  3500.000  1000.000    1
 4.14756964E+00 5.96166664E-03-2.37294852E-06 4.67412171E-10-3.61235213E-14    2
 2.59359992E+04-1.23028121E+00 8.08681094E-01 2.33615629E-02-3.55171815E-05    3
 2.80152437E-08-8.50072974E-12 2.64289807E+04 1.39397051E+01                   4
C2H3              L 2/92C   2H   3          G   200.000  3500.000  1000.000    1
 3.01672400E+00 1.03302292E-02-4.68082349E-06 1.01763288E-09-8.62607041E-14    2
 3.46128739E+04 7.78732378E+00 3.21246645E+00 1.51479162E-03 2.59209412E-05    3
-3.57657847E-08 1.47150873E-11 3.48598468E+04 8.51054025E+00                   4
C2H4              L 1/91C   2H   4          G   200.000  3500.000  1000.000    1
 2.03611116E+00 1.46454151E-02-6.71077915E-06 1.47222923E-09-1.25706061E-13    2
 4.93988614E+03 1.03053693E+01 3.95920148E+00-7.57052247E-03 5.70990292E-05    3
-6.91588753E-08 2.69884373E-11 5.08977593E+03 4.09733096E+00                   4
C2H5              L12/92C   2H   5          G   200.000  3500.000  1000.000    1
 1.95465642E+00 1.73972722E-02-7.98206668E-06 1.75217689E-09-1.49641576E-13    2
 1.28575200E+04 1.34624343E+01 4.30646568E+00-4.18658892E-03 4.97142807E-05    3
-5.99126606E-08 2.30509004E-11 1.28416265E+04 4.70720924E+00                   4
C2H6              L 8/88C   2H   6          G   200.000  3500.000  1000.000    1
 1.07188150E+00 2.16852677E-02-1.00256067E-05 2.21412001E-09-1.90002890E-13    2
-1.14263932E+04 1.51156107E+01 4.29142492E+00-5.50154270E-03 5.99438288E-05    3
-7.08466285E-08 2.68685771E-11-1.15222055E+04 2.66682316E+00                   4
CH2CO             L 5/90C   2H   2O   1     G   200.000  3500.000  1000.000    1
 4.51129732E+00 9.00359745E-03-4.16939635E-06 9.23345882E-10-7.94838201E-14    2
-7.55105311E+03 6.32247205E-01 2.13583630E+00 1.81188721E-02-1.73947474E-05    3
 9.34397568E-09-2.01457615E-12-7.04291804E+03 1.22156480E+01                   4
HCCO              SRIC91H   1C   2O   1     G   290.00   4000.00   1000.000    1
 0.56282058E+01 0.40853401E-02-0.15934547E-05 0.28626052E-09-0.19407832E-13    2
 0.19327215E+05-0.39302595E+01 0.22517214E+01 0.17655021E-01-0.23729101E-04    3
 0.17275759E-07-0.50664811E-11 0.20059449E+05 0.12490417E+02                   4
HCCOH              SRI91C   2O   1H   2     G   290.000  5000.000  1000.000    1
 0.59238291E+01 0.67923600E-02-0.25658564E-05 0.44987841E-09-0.29940101E-13    2
 0.72646260E+04-0.76017742E+01 0.12423733E+01 0.31072201E-01-0.50866864E-04    3
 0.43137131E-07-0.14014594E-10 0.80316143E+04 0.13874319E+02                   4
H2CN               41687H   2C   1N   1     G   290.00   4000.000  1000.000    1
 0.52097030E+01 0.29692911E-02-0.28555891E-06-0.16355500E-09 0.30432589E-13    2
 0.27677109E+05-0.44444780E+01 0.28516610E+01 0.56952331E-02 0.10711400E-05    3
-0.16226120E-08-0.23511081E-12 0.28637820E+05 0.89927511E+01                   4
HCN               GRI/98H   1C   1N   1     G   200.000  6000.000  1000.000    1
 0.38022392E+01 0.31464228E-02-0.10632185E-05 0.16619757E-09-0.97997570E-14    2
 0.14407292E+05 0.15754601E+01 0.22589886E+01 0.10051170E-01-0.13351763E-04    3
 0.10092349E-07-0.30089028E-11 0.14712633E+05 0.89164419E+01                   4
HNO               And93 H   1N   1O   1     G   200.000  6000.000  1000.000    1
 0.29792509E+01 0.34944059E-02-0.78549778E-06 0.57479594E-10-0.19335916E-15    2
 0.11750582E+05 0.86063728E+01 0.45334916E+01-0.56696171E-02 0.18473207E-04    3
-0.17137094E-07 0.55454573E-11 0.11548297E+05 0.17498417E+01                   4
N                 L 6/88N   1               G   200.000  6000.000  1000.000    1
 0.24159429E+01 0.17489065E-03-0.11902369E-06 0.30226245E-10-0.20360982E-14    2
 0.56133773E+05 0.46496096E+01 0.25000000E+01 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00 0.56104637E+05 0.41939087E+01                   4
NNH               T07/93N   2H   1          G   200.000  6000.000  1000.000    1
 0.37667544E+01 0.28915082E-02-0.10416620E-05 0.16842594E-09-0.10091896E-13    2
 0.28650697E+05 0.44705067E+01 0.43446927E+01-0.48497072E-02 0.20059459E-04    3
-0.21726464E-07 0.79469539E-11 0.28791973E+05 0.29779410E+01                   4
N2O               L 7/88N   2O   1          G   200.000  6000.000  1000.000    1
 0.48230729E+01 0.26270251E-02-0.95850874E-06 0.16000712E-09-0.97752303E-14    2
 0.80734048E+04-0.22017207E+01 0.22571502E+01 0.11304728E-01-0.13671319E-04    3
 0.96819806E-08-0.29307182E-11 0.87417744E+04 0.10757992E+02                   4
NH                And94 N   1H   1          G   200.000  6000.000  1000.000    1
 0.27836928E+01 0.13298430E-02-0.42478047E-06 0.78348501E-10-0.55044470E-14    2
 0.42120848E+05 0.57407799E+01 0.34929085E+01 0.31179198E-03-0.14890484E-05    3
 0.24816442E-08-0.10356967E-11 0.41880629E+05 0.18483278E+01                   4
NH2               And89 N   1H   2          G   200.000  6000.000  1000.000    1
 0.28347421E+01 0.32073082E-02-0.93390804E-06 0.13702953E-09-0.79206144E-14    2
 0.22171957E+05 0.65204163E+01 0.42040029E+01-0.21061385E-02 0.71068348E-05    3
-0.56115197E-08 0.16440717E-11 0.21885910E+05-0.14184248E+00                   4
NH3               J 6/77N   1H   3          G   200.000  6000.000  1000.000    1
 0.26344521E+01 0.56662560E-02-0.17278676E-05 0.23867161E-09-0.12578786E-13    2
-0.65446958E+04 0.65662928E+01 0.42860274E+01-0.46605230E-02 0.21718513E-04    3
-0.22808887E-07 0.82638046E-11-0.67417285E+04-0.62537277E+00                   4
NO                RUS 78N   1O   1          G   200.000  6000.000  1000.000    1
 0.32606056E+01 0.11911043E-02-0.42917048E-06 0.69457669E-10-0.40336099E-14    2
 0.99209746E+04 0.63693027E+01 0.42184763E+01-0.46389760E-02 0.11041022E-04    3
-0.93361354E-08 0.28035770E-11 0.98446230E+04 0.22808464E+01                   4
NO2               L 7/88N   1O   2          G   200.000  6000.000  1000.000    1
 0.48847542E+01 0.21723956E-02-0.82806906E-06 0.15747510E-09-0.10510895E-13    2
 0.23164983E+04-0.11741695E+00 0.39440312E+01-0.15854290E-02 0.16657812E-04    3
-0.20475426E-07 0.78350564E-11 0.28966179E+04 0.63119917E+01                   4
HCNO              BDEA94H   1N   1C   1O   1G   290.000  5000.000  1382.000    1
 6.59860456E+00 3.02778626E-03-1.07704346E-06 1.71666528E-10-1.01439391E-14    2
 1.79661339E+04-1.03306599E+01 2.64727989E+00 1.27505342E-02-1.04794236E-05    3
 4.41432836E-09-7.57521466E-13 1.92990252E+04 1.07332972E+01                   4
HOCN              BDEA94H   1N   1C   1O   1G   290.000  5000.000  1368.000    1
 5.89784885E+00 3.16789393E-03-1.11801064E-06 1.77243144E-10-1.04339177E-14    2
-3.70653331E+03-6.18167825E+00 3.78604952E+00 6.88667922E-03-3.21487864E-06    3
 5.17195767E-10 1.19360788E-14-2.82698400E+03 5.63292162E+00                   4
HNCO              BDEA94H   1N   1C   1O   1G   290.000  5000.000  1478.000    1
 6.22395134E+00 3.17864004E-03-1.09378755E-06 1.70735163E-10-9.95021955E-15    2
-1.66599344E+04-8.38224741E+00 3.63096317E+00 7.30282357E-03-2.28050003E-06    3
-6.61271298E-10 3.62235752E-13-1.55873636E+04 6.19457727E+00                   4
NCO               EA 93 N   1C   1O   1     G   200.000  6000.000  1000.000    1
 0.51521845E+01 0.23051761E-02-0.88033153E-06 0.14789098E-09-0.90977996E-14    2
 0.14004123E+05-0.25442660E+01 0.28269308E+01 0.88051688E-02-0.83866134E-05    3
 0.48016964E-08-0.13313595E-11 0.14682477E+05 0.95504646E+01                   4
CN                HBH92 C   1N   1          G   200.000  6000.000  1000.000    1
 0.37459805E+01 0.43450775E-04 0.29705984E-06-0.68651806E-10 0.44134173E-14    2
 0.51536188E+05 0.27867601E+01 0.36129351E+01-0.95551327E-03 0.21442977E-05    3
-0.31516323E-09-0.46430356E-12 0.51708340E+05 0.39804995E+01                   4
HCNN              SRI/94C   1N   2H   1     G   290.000  5000.000  1000.000    1
 0.58946362E+01 0.39895959E-02-0.15982380E-05 0.29249395E-09-0.20094686E-13    2
 0.53452941E+05-0.51030502E+01 0.25243194E+01 0.15960619E-01-0.18816354E-04    3
 0.12125540E-07-0.32357378E-11 0.54261984E+05 0.11675870E+02                   4
N2                121286N   2               G   290.000  5000.000  1000.000    1
 0.02926640E+02 0.14879768E-02-0.05684760E-05 0.10097038E-09-0.06753351E-13    2
-0.09227977E+04 0.05980528E+02 0.03298677E+02 0.14082404E-02-0.03963222E-04    3
 0.05641515E-07-0.02444854E-10-0.10208999E+04 0.03950372E+02                   4
AR                120186AR  1               G   290.000  5000.000  1000.000    1
 0.02500000E+02 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-0.07453750E+04 0.04366000E+02 0.02500000E+02 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-0.07453750E+04 0.04366000E+02                   4
C3H8              L 4/85C   3H   8          G   290.000  5000.000  1000.000    1
 0.75341368E+01 0.18872239E-01-0.62718491E-05 0.91475649E-09-0.47838069E-13    2
-0.16467516E+05-0.17892349E+02 0.93355381E+00 0.26424579E-01 0.61059727E-05    3
-0.21977499E-07 0.95149253E-11-0.13958520E+05 0.19201691E+02                   4
C3H7              L 9/84C   3H   7          G   290.000  5000.000  1000.000    1
 0.77026987E+01 0.16044203E-01-0.52833220E-05 0.76298590E-09-0.39392284E-13    2
 0.82984336E+04-0.15480180E+02 0.10515518E+01 0.25991980E-01 0.23800540E-05    3
-0.19609569E-07 0.93732470E-11 0.10631863E+05 0.21122559E+02                   4
CH3CHO            L 8/88C   2H   4O   1     G   200.000  6000.000  1000.000    1
 0.54041108E+01 0.11723059E-01-0.42263137E-05 0.68372451E-09-0.40984863E-13    2
-0.22593122E+05-0.34807917E+01 0.47294595E+01-0.31932858E-02 0.47534921E-04    3
-0.57458611E-07 0.21931112E-10-0.21572878E+05 0.41030159E+01                   4
CH2CHO            SAND86O   1H   3C   2     G   290.000  5000.000  1000.000    1
 0.05975670E+02 0.08130591E-01-0.02743624E-04 0.04070304E-08-0.02176017E-12    2
 0.04903218E+04-0.05045251E+02 0.03409062E+02 0.10738574E-01 0.01891492E-04    3
-0.07158583E-07 0.02867385E-10 0.15214766E+04 0.09558290E+02                   4
!
!!!!!!FROM SUNNY
HCOp                    C   1H   1O   1E  -1G   290.00  20000.00  3654.180     1
 0.69157265E+01 0.16631107E-03-0.18476343E-07 0.90856429E-12-0.16454441E-16    2
 0.96129333E+05-0.18089510E+02 0.28762768E+01 0.49486020E-02-0.23085306E-05    3
 0.51121752E-09-0.43520631E-13 0.99130022E+05 0.65020266E+01                   4
OHs               RUS 78O   1H   1          G   200.000  6000.000  1000.000    1
 2.75582920E+00 1.39848756E-03-4.19428493E-07 6.33453282E-11-3.56042218E-15    2
 5.09751756E+04 5.62581429E+00 3.46084428E+00 5.01872172E-04-2.00254474E-06    3
 3.18901984E-09-1.35451838E-12 5.07349466E+04 1.73976415E+00                   4
C3H3                    C   3H   3          G   290.000  5000.00   1000.00     1
 .883104772E+01 .435719410E-02-.410906610E-06-.236872300E-09 .437652000E-13    2
 .384741917E+05-.217791939E+02 .475419900E+01 .110802800E-01 .279332310E-06    3
-.547921220E-08 .194962900E-11 .398888300E+05 .585454700E+00                   4
C2H3Op                  C   2H   3O   1E  -1G   200.00  6000.000  1000.000     1
 0.53137165E+01 0.91737793E-02-0.33220386E-05 0.53947456E-09-0.32452368E-13    2
 0.76901865E+05-0.16757558E+01 0.40358705E+01 0.87729487E-03 0.30710010E-04    3
-0.39247565E-07 0.15296869E-10 0.77864832E+05 0.78617682E+01                   4
!O2n                               O   2E   1G   290.00  6000.000  2008.710     1
! 0.42592867E+01 0.22468072E-03-0.51397955E-07 0.73545978E-11-0.38558652E-15    2
!-0.72426252E+04 0.47599697E+00 0.31021718E+01 0.27980875E-02-0.22651126E-05    3
! 0.86916517E-09-0.12721884E-12-0.68074793E+04 0.67609020E+01                   4
!CHO2n                   H   1C   1O   2E   1G   290.000 5000.000  1000.000     1
! 5.97791811E+00 3.24247847E-03-1.46666291E-06 2.91808902E-10-2.10704956E-14    2
!-5.81813435E+04-7.12854015E+00-3.01936623E+01 2.54607495E-01-6.43484728E-04    3
! 6.92943698E-07-2.65871657E-10-5.36791044E+04 1.47958586E+02                   4
!C2                      C   2               G   300.000 5000.000  1000.000     1
! 0.04135978E+02 0.06531618E-03 0.01837099E-05-0.05295085E-09 0.04712137E-13    2
! 0.09967272E+06 0.07472923E+01 0.06996045E+02-0.07400601E-01 0.03234703E-04    3
! 0.04802535E-07-0.03295917E-10 0.09897487E+06-0.13862268E+02                   4
!CHs                     C   1H   1          G   200.000 6000.000  1000.000     1
! 2.78220752E+00 1.47246754E-03-4.63436227E-07 7.32736021E-11-4.19705404E-15    2
! 1.04547060E+05 5.17421018E+00 3.47250101E+00 4.26443626E-04-1.95181794E-06    3
! 3.51755043E-09-1.60436174E-12 1.04334869E+05 1.44799533E+00                   4
CH3p              A12/04C  1.H  3.E -1.   0.G   200.000  6000.000 1000.        1
 2.41723886E+00 6.40287629E-03-2.21301978E-06 3.46738910E-10-2.02364572E-14    2
 1.31474291E+05 6.78764161E+00 4.73043702E+00-8.66259820E-03 3.12269215E-05    3
-3.13568798E-08 1.09957173E-11 1.31269897E+05-3.03197684E+00 1.32514363E+05    4
!C3H3p  Propargyl  T 7/14C  3.H  3.E -1.   0.G   298.150  6000.000 1000.        1 ! real
C3H3p  Propargyl  T 7/14C  3.H  3.E -1.   0.G   290.000  6000.000 1000.        1
 6.53692307E+00 8.18501375E-03-2.88413266E-06 4.59392185E-10-2.72379381E-14    2
 1.41067534E+05-1.01646822E+01 2.03504329E+00 2.49892888E-02-2.90704802E-05    3
 2.01081492E-08-5.79204271E-12 1.42136831E+05 1.21137159E+01 1.43634441E+05    4
!!!!!!!FROM BURCAT DATABASE
!On                g 1/97O  1.E  1.   0.   0.G   298.150  6000.000 1000.        1 ! real
On                g 1/97O   1E   1          G   290.000  6000.000 1000.        1
 2.54474869E+00-4.66695513E-05 1.84912357E-08-3.18159223E-12 1.98962956E-16    2
 1.15042089E+04 4.52131015E+00 2.90805921E+00-1.69804907E-03 2.98069955E-06    3
-2.43835127E-09 7.61229311E-13 1.14357717E+04 2.80339097E+00 1.22492116E+04    4
!OHn               g 4/02O  1.H  1.E  1.   0.G   298.150  6000.000 1000.        1 ! real
OHn               g 4/02O  1.H  1.E  1.   0.G   290.000  6000.000 1000.        1
 2.80023747E+00 1.13380509E-03-2.99666184E-07 4.01911483E-11-1.78988913E-15    2
-1.82535298E+04 4.69394620E+00 3.43126659E+00 6.31146866E-04-1.92914359E-06    3
 2.40618712E-09-8.66679361E-13-1.85085918E+04 1.07990541E+00-1.74702052E+04    4
!H3Op              ATcT AH  3.O  1.E -1.   0.G   298.150  6000.000 1000.        1 ! real
H3Op              ATcT AH  3.O  1.E -1.   0.G   290.000  6000.000 1000.        1
 2.49647765E+00 5.72844840E-03-1.83953239E-06 2.73577348E-10-1.54093917E-14    2
 7.16244227E+04 7.45850493E+00 3.79295251E+00-9.10852723E-04 1.16363521E-05    3
-1.21364865E-08 4.26159624E-12 7.14027518E+04 1.47156927E+00 7.25739701E+04    4
!!E                 g12/98E  1.   0.   0.   0.G   298.150  6000.000 1000.        1 ! real
E                 g12/98E  1.   0.   0.   0.G   290.000  6000.000 1000.        1
 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-7.45375000E+02-1.17208122E+01 2.50000000E+00 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-7.45375000E+02-1.17208122E+01 0.00000000E+00    4
CH2CO             g 4/02C  2.H  2.O  1.   0.G   200.000  6000.000 1000.        1
 5.75871449E+00 6.35124053E-03-2.25955361E-06 3.62321512E-10-2.15855515E-14    2
-8.08533464E+03-4.96490444E+00 2.13241136E+00 1.81319455E-02-1.74093315E-05    3
 9.35336040E-09-2.01724844E-12-7.14808520E+03 1.33807969E+01-5.84267744E+03    4
CH(S)             EG4/09C  1.H  1.   0.   0.G   200.000  6000.000 1000.        1
 2.78220752E+00 1.47246754E-03-4.63436227E-07 7.32736021E-11-4.19705404E-15    2
 1.04547060E+05 5.17421018E+00 3.47250101E+00 4.26443626E-04-1.95181794E-06    3
 3.51755043E-09-1.60436174E-12 1.04334869E+05 1.44799533E+00 1.05378099E+05    4
!CH excited B2Sig  EG4/09C  1.H  1.   0.   0.G   200.000  6000.000 1000.        1
! 3.10123129E+00 1.48670484E-03-3.88548896E-07 1.04667880E-10-4.71283513E-15    2
! 1.08312371E+05 3.30360492E+00 3.69510692E+00-1.80505435E-03 5.03652887E-06    3
!-3.23608801E-09 6.08848851E-13 1.08268493E+05 7.40672572E-01 1.09328350E+05    4
!
END

\\
\\
\\  This is the tran file
\\
\\
AR                 0   136.500     3.330     0.000     0.000     0.000
C                  0    71.400     3.298     0.000     0.000     0.000 ! *
C2                 1    97.530     3.621     0.000     1.760     4.000
C2O                1   232.400     3.828     0.000     0.000     1.000 ! *
CN2                1   232.400     3.828     0.000     0.000     1.000 ! OIS
C2H                1   209.000     4.100     0.000     0.000     2.500
C2H2               1   209.000     4.100     0.000     0.000     2.500
C2H2OH             2   224.700     4.162     0.000     0.000     1.000 ! *
C2H3               2   209.000     4.100     0.000     0.000     1.000 ! *
C2H4               2   280.800     3.971     0.000     0.000     1.500
C2H5               2   252.300     4.302     0.000     0.000     1.500
C2H6               2   252.300     4.302     0.000     0.000     1.500
C2N                1   232.400     3.828     0.000     0.000     1.000 ! OIS
C2N2               1   349.000     4.361     0.000     0.000     1.000 ! OIS
C3H2               2   209.000     4.100     0.000     0.000     1.000 ! *
C3H4               1   252.000     4.760     0.000     0.000     1.000
C3H6               2   266.800     4.982     0.000     0.000     1.000
C3H7               2   266.800     4.982     0.000     0.000     1.000
C4H6               2   357.000     5.180     0.000     0.000     1.000
I*C3H7             2   266.800     4.982     0.000     0.000     1.000
N*C3H7             2   266.800     4.982     0.000     0.000     1.000
C3H8               2   266.800     4.982     0.000     0.000     1.000
C4H                1   357.000     5.180     0.000     0.000     1.000
C4H2               1   357.000     5.180     0.000     0.000     1.000
C4H2OH             2   224.700     4.162     0.000     0.000     1.000 ! *
C4H8               2   357.000     5.176     0.000     0.000     1.000
C4H9               2   357.000     5.176     0.000     0.000     1.000
I*C4H9             2   357.000     5.176     0.000     0.000     1.000
C5H2               1   357.000     5.180     0.000     0.000     1.000
C5H3               1   357.000     5.180     0.000     0.000     1.000
C6H2               1   357.000     5.180     0.000     0.000     1.000
C6H5               2   412.300     5.349     0.000     0.000     1.000 ! JAM
C6H5O              2   450.000     5.500     0.000     0.000     1.000 ! JAM
C5H5OH             2   450.000     5.500     0.000     0.000     1.000 ! JAM
C6H6               2   412.300     5.349     0.000     0.000     1.000 ! SVE
C6H7               2   412.300     5.349     0.000     0.000     1.000 ! JAM
CH                 1    80.000     2.750     0.000     0.000     0.000
CH2                1   144.000     3.800     0.000     0.000     0.000
CH2(S)             1   144.000     3.800     0.000     0.000     0.000
CH2*               1   144.000     3.800     0.000     0.000     0.000
CH2CHCCH           2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH2CHCCH2          2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH2CHCH2           2   260.000     4.850     0.000     0.000     1.000 ! JAM
CH2CHCHCH          2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH2CHCHCH2         2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH2CO              2   436.000     3.970     0.000     0.000     2.000
CH2O               2   498.000     3.590     0.000     0.000     2.000
CH2OH              2   417.000     3.690     1.700     0.000     2.000
CH3                1   144.000     3.800     0.000     0.000     0.000
CH3CC              2   252.000     4.760     0.000     0.000     1.000 ! JAM
CH3CCCH2           2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH3CCCH3           2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH3CCH2            2   260.000     4.850     0.000     0.000     1.000 ! JAM
CH3CHCH            2   260.000     4.850     0.000     0.000     1.000 ! JAM
CH3CH2CCH          2   357.000     5.180     0.000     0.000     1.000 ! JAM
CH3CHO             2   436.000     3.970     0.000     0.000     2.000
CH2CHO             2   436.000     3.970     0.000     0.000     2.000
CH3CO              2   436.000     3.970     0.000     0.000     2.000
CH3O               2   417.000     3.690     1.700     0.000     2.000
CH3OH              2   481.800     3.626     0.000     0.000     1.000 ! SVE
CH4                2   141.400     3.746     0.000     2.600    13.000
CH4O               2   417.000     3.690     1.700     0.000     2.000
CN                 1    75.000     3.856     0.000     0.000     1.000 ! OIS
CNC                1   232.400     3.828     0.000     0.000     1.000 ! OIS
CNN                1   232.400     3.828     0.000     0.000     1.000 ! OIS
CO                 1    98.100     3.650     0.000     1.950     1.800
CO2                1   244.000     3.763     0.000     2.650     2.100
H                  0   145.000     2.050     0.000     0.000     0.000
H2C4O              2   357.000     5.180     0.000     0.000     1.000 ! JAM
H2                 1    38.000     2.920     0.000     0.790   280.000
H2CCCCH            2   357.000     5.180     0.000     0.000     1.000 ! JAM
H2CCCCH2           2   357.000     5.180     0.000     0.000     1.000 ! JAM
H2CCCH             2   252.000     4.760     0.000     0.000     1.000 ! JAM
H2CN               1   569.000     3.630     0.000     0.000     1.000 ! os/jm
H2NO               2   116.700     3.492     0.000     0.000     1.000 ! JAM
H2O                2   572.400     2.605     1.844     0.000     4.000
H2O2               2   107.400     3.458     0.000     0.000     3.800
HC2N2              1   349.000     4.361     0.000     0.000     1.000 ! OIS
HCCHCCH            2   357.000     5.180     0.000     0.000     1.000 ! JAM
HCCO               2   150.000     2.500     0.000     0.000     1.000 ! *
HCNN               2   150.000     2.500     0.000     0.000     1.000 ! *
HCCOH              2   436.000     3.970     0.000     0.000     2.000
HCN                1   569.000     3.630     0.000     0.000     1.000 ! OIS
HCO                2   498.000     3.590     0.000     0.000     0.000
HE                 0    10.200     2.576     0.000     0.000     0.000 ! *
HCNO               2   232.400     3.828     0.000     0.000     1.000 ! JAM
HOCN               2   232.400     3.828     0.000     0.000     1.000 ! JAM
HNCO               2   232.400     3.828     0.000     0.000     1.000 ! OIS
HNNO               2   232.400     3.828     0.000     0.000     1.000 ! *
HNO                2   116.700     3.492     0.000     0.000     1.000 ! *
HNOH               2   116.700     3.492     0.000     0.000     1.000 ! JAM
HO2                2   107.400     3.458     0.000     0.000     1.000 ! *
N                  0    71.400     3.298     0.000     0.000     0.000 ! *
N2                 1    97.530     3.621     0.000     1.760     4.000
N2H2               2    71.400     3.798     0.000     0.000     1.000 ! *
N2H3               2   200.000     3.900     0.000     0.000     1.000 ! *
N2H4               2   205.000     4.230     0.000     4.260     1.500
N2O                1   232.400     3.828     0.000     0.000     1.000 ! *
NCN                1   232.400     3.828     0.000     0.000     1.000 ! OIS
NCO                1   232.400     3.828     0.000     0.000     1.000 ! OIS
NH                 1    80.000     2.650     0.000     0.000     4.000
NH2                2    80.000     2.650     0.000     2.260     4.000
NH3                2   481.000     2.920     1.470     0.000    10.000
NNH                2    71.400     3.798     0.000     0.000     1.000 ! *
NO                 1    97.530     3.621     0.000     1.760     4.000
NCNO               2   232.400     3.828     0.000     0.000     1.000 ! OIS
NO2                2   200.000     3.500     0.000     0.000     1.000 ! *
O                  0    80.000     2.750     0.000     0.000     0.000
O2                 1   107.400     3.458     0.000     1.600     3.800
OH                 1    80.000     2.750     0.000     0.000     0.000
! Additional terms
H3Op               2   572.400     2.605     1.844     0.000     4.000 ! same as H2O (??)
E                  0   850.        425.      0.000     0.000     1.000 !(singh) ! from CHEMKIN transport data
CHO2n              2   232.400     3.828     0.000     0.000     1.000 ! JAM ! same as HNCO (??)
C2H3Op             2   436.000     3.970     0.000     0.000     2.000 ! same as CH2CHO (??)
HCOp               2   498.000     3.590     0.000     0.000     0.000 ! same as HCO
O2n                1   107.400     3.458     0.000     1.600     3.800 ! same as O2
OHn                1    80.000     2.750     0.000     0.000     0.000 ! same as OH
On                 0    80.000     2.750     0.000     0.000     0.000 ! same as O
CH(S)              1    80.000     2.750     0.000     0.000     0.000 ! same as CH
OHs                1    80.000     2.750     0.000     0.000     0.000 ! same as OH
CH3p               1   144.000     3.800     0.000     0.000     0.000 ! same as CH3
C3H3               1   252.000     4.760     0.000     0.000     1.000 ! from CHEMKIN
C3H3p              1   252.000     4.760     0.000     0.000     1.000 ! same as C3H3

#endif
