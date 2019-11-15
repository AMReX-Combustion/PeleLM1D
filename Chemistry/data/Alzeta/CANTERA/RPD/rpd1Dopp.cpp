#include "cantera/oneD/Sim1D.h"
#include "cantera/oneD/Inlet1D.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/IdealGasMix.h"
#include "cantera/transport.h"
#include "cantera/kinetics/ReactionPath.h"
#include <fstream>

using namespace Cantera;
using fmt::print;
using std::vector;
using std::string;
using std::map;
using std::endl;

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
  MyReactionPathBuilder(StFlow& stFlow,Sim1D& flame) : ReactionPathBuilder(), m_stFlow(stFlow), m_flame(flame) {}
  virtual ~MyReactionPathBuilder() {}

  int init(std::ostream& logfile, Kinetics& s) {ReactionPathBuilder::init(logfile,s);}

  int build(Kinetics& s, const std::string& element, std::ostream& output,
	    ReactionPathDiagram& r, bool quiet=false);

  //! Analyze a reaction to determine which reactants lead to which products.
  int findGroups(std::ostream& logfile, Kinetics& s) {ReactionPathBuilder::findGroups(logfile,s);}
protected:
  StFlow& m_stFlow;
  Sim1D& m_flame;
};

int
MyReactionPathBuilder::build(Kinetics& s, const std::string& element, std::ostream& output,
			     ReactionPathDiagram& r, bool quiet)
{
  int nr = s.nReactions();
  vector_fp ropf(nr), ropfT(nr,0);
  vector_fp ropr(nr), roprT(nr,0);
  const vector_fp& grid = m_stFlow.grid();
  int Npts = grid.size();
  
  for (int i=0; i<Npts; ++i) {
    m_stFlow.setGas(m_flame.solution(),i);
    //std::cout << "T at " << i << ": " << m_stFlow.phase().temperature() << std::endl;
    s.getFwdRatesOfProgress(ropf.data());
    s.getRevRatesOfProgress(ropr.data());
    
    doublereal zC = grid[i];
    doublereal zL = (i==0      ? zC : 0.5*(grid[i-1]+zC));
    doublereal zR = (i==Npts-1 ? zC : 0.5*(grid[i+1]+zC));
    doublereal vol = zR - zL;

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

int rpd(double phi)
{
    try {
        IdealGasMix gas("gri30.cti","gri30_mix");

        doublereal temp = 300.0; // K
        doublereal pressure = 1.0*OneAtm; //atm
        doublereal uin = 0.3; //m/sec

        size_t nsp = gas.nSpecies();
        vector_fp x(nsp, 0.0);

        doublereal C_atoms = 1.0;
        doublereal H_atoms = 4.0;
        doublereal ax = C_atoms + H_atoms / 4.0;
        doublereal fa_stoic = 1.0 / (4.76 * ax);
        x[gas.speciesIndex("CH4")] = 1.0;
        x[gas.speciesIndex("O2")] = 0.21 / phi / fa_stoic;
        x[gas.speciesIndex("N2")] = 0.79 / phi/ fa_stoic;

        gas.setState_TPX(temp, pressure, x.data());
        doublereal rho_in = gas.density();

        vector_fp yin(nsp);
        gas.getMassFractions(&yin[0]);

        doublereal phi1 = 0.;
        vector_fp x1(nsp, 0.0);
        x1[gas.speciesIndex("O2")] = 0.21;
        x1[gas.speciesIndex("N2")] = 0.79;
        doublereal temp1 = 1000;
        gas.setState_TPX(temp1, pressure, x1.data());
        vector_fp y1(nsp);
        gas.getMassFractions(&y1[0]);
        doublereal rho1_out = gas.density();

        //=============  build each domain ========================


        //-------- step 1: create the flow -------------

        //FreeFlame flow(&gas);
        AxiStagnFlow flow(&gas);

        // create an initial grid
        int nz = 6;
        doublereal lz = 0.02;
        vector_fp z(nz);
        doublereal dz = lz/((doublereal)(nz-1));
        for (int iz = 0; iz < nz; iz++) {
            z[iz] = ((doublereal)iz)*dz;
        }

        flow.setupGrid(nz, &z[0]);

        // specify the objects to use to compute kinetic rates and
        // transport properties

        std::unique_ptr<Transport> trmix(newTransportMgr("Mix", &gas));
        std::unique_ptr<Transport> trmulti(newTransportMgr("Multi", &gas));

        flow.setTransport(*trmix);
        flow.setKinetics(gas);
        flow.setPressure(pressure);

        //------- step 2: create the inlet  -----------------------

        Inlet1D inlet;

        inlet.setMoleFractions(x.data());
        doublereal mdot=uin*rho_in;
        inlet.setMdot(mdot);
        inlet.setTemperature(temp);



        Inlet1D inlet1;

        inlet1.setMoleFractions(x1.data());
        doublereal mdot1=uin*rho1_out;
        inlet1.setMdot(mdot1);
        inlet1.setTemperature(temp1);


        //------- step 3: create the outlet  ---------------------

        //Outlet1D outlet;

        //=================== create the container and insert the domains =====

        //std::vector<Domain1D*> domains { &inlet, &flow, &outlet };
        std::vector<Domain1D*> domains { &inlet, &flow, &inlet1 };
        Sim1D flame(domains);

        //----------- Supply initial guess----------------------

        vector_fp locs{0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
        vector_fp value;

        double uout = inlet.mdot()/rho1_out;
        //value = {uin, uin, uout, uout};
        value = {uin, uin, uin, uin};
        flame.setInitialGuess("u",locs,value);
        //value = {temp, temp, Tad, Tad};
        value = {temp, temp, temp1, temp1};
        flame.setInitialGuess("T",locs,value);

        for (size_t i=0; i<nsp; i++) {
            value = {yin[i], yin[i], y1[i], y1[i]};
            flame.setInitialGuess(gas.speciesName(i),locs,value);
        }

        inlet.setMoleFractions(x.data());
        inlet.setMdot(mdot);
        inlet.setTemperature(temp);

        // flame.showSolution();

        int flowdomain = 0.02;
        double ratio = 10.0;
        double slope = 0.08;
        double curve = 0.1;

        flame.setRefineCriteria(flowdomain,ratio,slope,curve);

        int loglevel=5;
        flow.solveEnergyEqn();
        bool refine_grid = true;

        //refine_grid=false;

        flame.solve(loglevel,refine_grid);

        ReactionPathDiagram rpd;
        MyReactionPathBuilder rpb(flow,flame);
        std::string rpd_log("rpd.log");
        std::ofstream ofs(rpd_log.c_str());
        int init_ok = rpb.init(ofs,gas);
        ofs.close();

        bool build_quiet = false;
        rpb.build(gas,"C",std::cout,rpd,build_quiet);

        std::string rpd_dot("rpd.dot");
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
    double phi;
    //print("Enter phi: ");
    //std::cin >> phi;
    phi = 0.9;
    return rpd(phi);
}
