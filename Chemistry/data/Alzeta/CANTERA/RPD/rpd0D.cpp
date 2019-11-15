#include "cantera/oneD/Sim1D.h"
#include "cantera/oneD/Inlet1D.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/IdealGasMix.h"
#include "cantera/transport.h"
#include "cantera/kinetics/ReactionPath.h"
#include <fstream>

#include "cantera/zerodim.h"
#include "cantera/IdealGasMix.h"
//#include "example_utils.h"
#include "cantera/base/plots.h"

using namespace Cantera;
using fmt::print;
using std::vector;
using std::string;
using std::map;
using std::endl;
using std::cout;

template<class G, class A>
void saveSoln(double time, const G& gas, A& soln)
{
    soln.resize(static_cast<int>(soln.nRows()),
                static_cast<int>(soln.nColumns()) + 1);
    int back = static_cast<int>(soln.nColumns()) - 1;
    soln(0,back) = time;
    soln(1,back) = gas.temperature();
    soln(2,back) = gas.density();
    soln(3,back) = gas.pressure();
    size_t nsp = gas.nSpecies();
    for (size_t k = 0; k < nsp; k++) {
        soln(4+k,back) = gas.moleFraction(k);
    }
}


template<class G, class A>
void saveSoln(int i, double time, const G& gas, A& soln)
{
    soln(0,i) = time;
    soln(1,i) = gas.temperature();
    soln(2,i) = gas.density();
    soln(3,i) = gas.pressure();
    gas.getMoleFractions(&soln(4,i));
}

string reactionLabel(size_t i, size_t kr, size_t nr,
                     const std::vector<size_t>& slist, const Kinetics& s)
{
  string label = "";
  for (size_t j = 0; j < nr; j++) {
    if (j != kr) {
      label += " + "+ s.kineticsSpeciesName(slist[j]);
    }
  }
  if (s.reactionType(i) == THREE_BODY_RXN) {
    label += " + M ";
  } else if (s.reactionType(i) == FALLOFF_RXN) {
    label += " (+ M)";
  }
  return label;
}

class MyReactionPathBuilder : public ReactionPathBuilder
{
public:
  MyReactionPathBuilder(ReactorBase& r) : ReactionPathBuilder(), m_r(r) {}
  virtual ~MyReactionPathBuilder() {}

  int init(std::ostream& logfile, Kinetics& s) {ReactionPathBuilder::init(logfile,s);}

  int build(Kinetics& s, const std::string& element, std::ostream& output,
	    ReactionPathDiagram& r, const Array2D& soln, bool quiet=false);

  //! Analyze a reaction to determine which reactants lead to which products.
  int findGroups(std::ostream& logfile, Kinetics& s) {ReactionPathBuilder::findGroups(logfile,s);}
protected:
  ReactorBase& m_r;
};

int
MyReactionPathBuilder::build(Kinetics& s, const std::string& element, std::ostream& output,
			     ReactionPathDiagram& r, const Array2D& soln, bool quiet)
{
  int nr = s.nReactions();
  vector_fp ropf(nr), ropfT(nr,0);
  vector_fp ropr(nr), roprT(nr,0);
  int Npts = soln.nColumns();

  size_t nsp = m_r.contents().nSpecies();
  doublereal X[nsp];
  
  for (int i=0; i<Npts; ++i) {
    for (size_t k = 0; k < nsp; k++) {
      X[k] = soln(4+k,i);
    }
    m_r.contents().setState_TPX(soln(1,i),soln(3,i),X);
    s.getFwdRatesOfProgress(ropf.data());
    s.getRevRatesOfProgress(ropr.data());

    doublereal vol = 0;

    if (i<Npts-1) {
      vol += 0.5*(soln(0,i)+soln(0,i+1));
    }
    if (i>0) {
      vol += 0.5*(soln(0,i)+soln(0,i-1));
    }

    for (int j=0; j<nr; ++j) {
      ropfT[j] += ropf[j] * vol;
      roprT[j] += ropr[j] * vol;
    }
  }

  map<size_t, int> warn;
  doublereal threshold = 0.0;
  size_t m = m_enamemap[element]-1;
  r.element = element;
  if (m == npos) {
    return -1;
  }

  // species explicitly included or excluded
  vector<string>& in_nodes = r.included();
  vector<string>& out_nodes = r.excluded();

  vector_int status(s.nTotalSpecies(), 0);
  for (size_t ni = 0; ni < in_nodes.size(); ni++) {
    status[s.kineticsSpeciesIndex(in_nodes[ni])] = 1;
  }
  for (size_t ne = 0; ne < out_nodes.size(); ne++) {
    status[s.kineticsSpeciesIndex(out_nodes[ne])] = -1;
  }

  for (size_t i = 0; i < m_nr; i++) {
    double ropf = ropfT[i];
    double ropr = roprT[i];

    // loop over reactions involving element m
    if (m_elatoms(m, i) > 0) {
      size_t nr = m_reac[i].size();
      size_t np = m_prod[i].size();

      for (size_t kr = 0; kr < nr; kr++) {
	size_t kkr = m_reac[i][kr];
	string fwdlabel = reactionLabel(i, kr, nr, m_reac[i], s);

	for (size_t kp = 0; kp < np; kp++) {
	  size_t kkp = m_prod[i][kp];
	  string revlabel = "";
	  for (size_t j = 0; j < np; j++) {
	    if (j != kp) {
	      revlabel += " + "+ s.kineticsSpeciesName(m_prod[i][j]);
	    }
	  }
	  if (s.reactionType(i) == THREE_BODY_RXN) {
	    revlabel += " + M ";
	  } else if (s.reactionType(i) == FALLOFF_RXN) {
	    revlabel += " (+ M)";
	  }

	  // calculate the flow only for pairs that are not the same
	  // species, both contain atoms of element m, and both are
	  // allowed to appear in the diagram
	  if ((kkr != kkp) && (m_atoms(kkr,m) > 0
			       && m_atoms(kkp,m) > 0)
	      && status[kkr] >= 0 && status[kkp] >= 0) {
	    // if neither species contains the full number of atoms
	    // of element m in the reaction, then we must consider
	    // the type of reaction to determine which reactant
	    // species was the source of a given m-atom in the
	    // product
	    double f;
	    if ((m_atoms(kkp,m) < m_elatoms(m, i)) &&
		(m_atoms(kkr,m) < m_elatoms(m, i))) {
	      map<size_t, map<size_t, Group> >& g = m_transfer[i];
	      if (g.empty()) {
		if (!warn[i] && !quiet) {
		  output << endl;
		  output << "*************** REACTION IGNORED ***************" << endl;
		  output << "Warning: no rule to determine partitioning of " << element
			 << endl << " in reaction " << s.reactionString(i) << "." << endl
			 << "*************** REACTION IGNORED **************" << endl;
		  output << endl;
		  warn[i] = 1;
		}
		f = 0.0;
	      } else {
		if (!g[kkr][kkp]) {
		  f = 0.0;
		} else {
		  f = g[kkr][kkp].nAtoms(m);
		}
	      }
	    } else {
	      // no ambiguity about where the m-atoms come from or
	      // go to. Either all reactant m atoms end up in one
	      // product, or only one reactant contains all the
	      // m-atoms. In either case, the number of atoms
	      // transferred is given by the same expression.
	      f = m_atoms(kkp,m) * m_atoms(kkr,m) / m_elatoms(m, i);
	    }

	    double fwd = ropf*f;
	    double rev = ropr*f;
	    bool force_incl = ((status[kkr] == 1) || (status[kkp] == 1));

	    bool fwd_incl = ((fwd > threshold) ||
			     (fwd > 0.0 && force_incl));
	    bool rev_incl = ((rev > threshold) ||
			     (rev > 0.0 && force_incl));
	    if (fwd_incl || rev_incl) {
	      if (!r.hasNode(kkr)) {
		r.addNode(kkr, s.kineticsSpeciesName(kkr), m_x[kkr]);
	      }
	      if (!r.hasNode(kkp)) {
		r.addNode(kkp, s.kineticsSpeciesName(kkp), m_x[kkp]);
	      }
	    }
	    if (fwd_incl) {
	      r.linkNodes(kkr, kkp, int(i), fwd, fwdlabel);
	    }
	    if (rev_incl) {
	      r.linkNodes(kkp, kkr, -int(i), rev, revlabel);
	    }
	  }
	}
      }
    }
  }
}

int rpd()
{
    try {
        IdealGasMix gas("gri30.cti","gri30_mix");

        size_t nsp = gas.nSpecies();
        vector_fp x(nsp, 0.0);

        gas.setState_TPX(1001.0, OneAtm, "H2:2.0, O2:1.0, N2:4.0");

        // create a reactor
        Reactor r;

        // create a reservoir to represent the environment
        Reservoir env;

        // specify the thermodynamic property and kinetics managers
        r.setThermoMgr(gas);
        r.setKineticsMgr(gas);
        env.setThermoMgr(gas);

        // create a flexible, insulating wall between the reactor and the
        // environment
        Wall w;
        w.install(r,env);

        // set the "Vdot coefficient" to a large value, in order to
        // approach the constant-pressure limit; see the documentation
        // for class Reactor
        w.setExpansionRateCoeff(1.e9);
        w.setArea(1.0);

        // create a container object to run the simulation
        // and add the reactor to it
        ReactorNet sim;
        sim.setVerbose(false);
        sim.addReactor(r);

        double tm;
        double dt = 1.e-5;    // interval at which output is written
        int nsteps = 100;     // number of intervals

        Array2D soln(nsp+4, 1);
        saveSoln(0, 0.0, gas, soln);

        // main loop
        for (int i = 1; i <= nsteps; i++) {
            tm = i*dt;
            sim.advance(tm);
            saveSoln(tm, gas, soln);
        }

        ReactionPathDiagram rpd;
        MyReactionPathBuilder rpb(r);
        std::string rpd_log("rpd0D.log");
        std::ofstream ofs(rpd_log.c_str());
        int init_ok = rpb.init(ofs,gas);
        ofs.close();

        bool build_quiet = false;
        rpb.build(gas,"C",std::cout,rpd,soln,build_quiet);

        rpd.threshold = 1.e-20;
        std::string rpd_dot("rpd0D.dot");
        std::ofstream ofs_dot(rpd_dot.c_str());
        rpd.exportToDot(ofs_dot);
        ofs_dot.close();

    } catch (CanteraError& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << "program terminating." << std::endl;
        return -1;
    }
    return 0;
}

int main()
{
    return rpd();
}
