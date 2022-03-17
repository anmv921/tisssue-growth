#include <iostream>

using namespace std;

#define NDIM  2

#include "in_mddefs.h"

// TODO move this to a seaparate file
typedef struct {
    VecR ra, fa;
    VecR rb, fb;
    int ncell;
    real rhoa, rhob; // radii
    bool dividing;
    real fna, fnb; // normal force
    real age;
} Cell;


// ************************************************************************************************************************
// Initalization of variables

vector<Cell> cells;
vector<string> inVars;
real dt, D, timeNow, uSum, lCell, gama, kCell, ka, kb, kr, r0, d0, fNcrit, dp0, rad_std, R0fene, \
     kab, rmax, pDivision0, criticalPressure, probDeath0, apoptisisAge;
int stepLimit, Ncells, printStep, stepCount, Nparticles, moreCycles, NcellsMax, Ncells0;
unsigned seed = 1455435;
mt19937 gen(seed);
normal_distribution<real> dist(0.0, 1.0);
ofstream ovito, propsDat;
VecR W, reg;

// ************************************************************************************************************************

int main(int argc, char **argv)
{	
	auto start = chrono::steady_clock::now();
	
	// move this into an init() function
    ReadInput ();
    SetParams ();
    SetupJob ();
    
    moreCycles = 1;
    PrintSummary ();
    while (moreCycles) {
        SingleStep ();
        if (stepCount >= stepLimit) moreCycles = 0;
    }
    
    // TODO move this into a function end()
    ovito.close();
    propsDat.close();
    
    PrintElapsedTime(start);
	
	//return 0;
}

// *****************************************************************************
// functions

void SingleStep ()
{
    ++stepCount;
    timeNow = stepCount * dt;
    computeForcesHarmonic();
    CBD_Step();
    
    if (stepCount > 1000) { // this allows the cells to stop overlapping after initialization
		cellDivision();
		cellDeath();
	}
	
    
    
    if (stepCount % printStep == 0) {
        PrintSummary ();
        writeProps();
    }
}

void SetupJob ()
{	
    AllocArrays ();
    stepCount = 0;
    InitCoordsRand();
    Init_F();
    VZero(W);
    ovito.open("ovito.xyz");
    timeNow = 0.;
    uSum = 0;
    stepCount = 0;
    propsDat.open("props.txt");
    srand(1456359366);
    writePropsHeader();
    
}

void SetParams ()
{
    Nparticles = Ncells * 2;
    d0 = r0*2;
    Ncells0 = Ncells;
    rmax = 3;
}

void AllocArrays ()
{	
	cells.resize(NcellsMax);
}

// harmonic
void pairForceContactHarmonic(VecR& r1, VecR& r2, VecR& f1, VecR& f2, real& rho1, real& rho2, \
				real& fn1, real& fn2) {
	real dx = r1.x - r2.x, dy = r1.y - r2.y, \
	     dr = Dist(dx, dy), \
	     rContact = rho1 + rho2;
	real req = rContact;
	

	// repulsive part of interaction
	if (dr <= req) {
		real F = -kr * (dr - req),
			 Fx = F  * dx/dr, \
			 Fy = F * dy/dr;
			 
		uSum += 0.5 * kr * (dr - req) * (dr - req);  
		AddForce(f1, Fx, Fy);
		AddForce(f2, -Fx, -Fy);
		
		// add normal force
		VecR fij, nVec;
		fij.x = Fx, fij.y = Fy;
		nVec.x = dx/dr; nVec.y = dy/dr;
		fn1 += abs(DotProd(fij, nVec));
		fn2 += abs(DotProd(fij, nVec));
		
	}
	
}

// only force for r>req & r<=rmax
void pairForceHarmonicLongRange(VecR& r1, VecR& r2, VecR& f1, VecR& f2, real& rho1, real& rho2, real k,\
		real& fn1, real& fn2) {
			
	real dx = r1.x - r2.x, dy = r1.y - r2.y, \
	     dr = Dist(dx, dy), \
	     rContact = rho1 + rho2;
	
	real req;
	req = rContact;
	
	if ( dr <= rmax and dr > req) {
		real F = -k * (dr - req),
	         Fx = F  * dx/dr, \
	         Fy = F * dy/dr;
	           
		uSum += 0.5 * k * (dr - req) * (dr - req);  
		AddForce(f1, Fx, Fy);
		AddForce(f2, -Fx, -Fy);
		
		// add normal force
		
		VecR fij, nVec;
		fij.x = Fx, fij.y = Fy;
		nVec.x = dx/dr; nVec.y = dy/dr;
		fn1 += abs(DotProd(fij, nVec));
		fn2 += abs(DotProd(fij, nVec));
	}
}

#define RStep(rinit, ff, dtggi, ww) \
	(rinit).x += (ff).x * dtggi + (ww).x, \
	(rinit).y += (ff).y * dtggi + (ww.y)
	
#define GenerateNoise(ww, s) \
	(ww).x = dist(gen) * s, \
	(ww).y = dist(gen) * s

void CBD_Step ()
{
    int n;
	DO_CELLS {
		// apical pole
		GenerateNoise( W, sqrt(2*dt*D) );
		RStep( cells[n].ra, cells[n].fa, dt/gama, W);
		
		// basal pole
		GenerateNoise( W, sqrt(2*dt*D) );
		RStep( cells[n].rb, cells[n].fb, dt/gama, W);
		
		
	}
}

#define GenRand01 \
	(double)rand() / RAND_MAX
	
void InitCoordsRand() {
	
	int n;
	DO_CELLS {
		cells[n].rhoa = r0, cells[n].rhob = r0, cells[n].ncell = n;
		
		real lo = -reg.x/2, hi = reg.x/2;
		cells[n].ra.x = GenRand01 * (hi - lo) + lo;
		lo = -reg.y/2, hi = reg.y/2;
		cells[n].ra.y = GenRand01 * (hi - lo) + lo;
		
		real angulo = GenRand01 * (2 * 3.1415);
		
		cells[n].rb.x = cells[n].ra.x + lCell * cos(angulo);
		cells[n].rb.y = cells[n].ra.y + lCell * sin(angulo);
		
		cells[n].age = 0;

	}
}

void Init_F ()
{
    int n;
	DO_CELLS VZero (cells[n].fa), VZero (cells[n].fb);
}

void PrintSummary ()
{	
	ovito << Nparticles << endl;
	ovito << "Lattice=\""<<reg.x<<" 0.0 0.0 0.0 "<<reg.y<<" 0.0 0.0 0.0 0.0\" Origin=\""<<-reg.x/2<<" ";
	ovito<<-reg.y/2<<" 0.0\"";
	ovito << " Properties=pos:R:2:Radius:R:1:Type:S:1:Identifier:I:1:Pressure:R:1:Age:R:1";
	ovito << " t=" << timeNow << endl;
	
    int counter=0;
    for( counter=0; counter < Ncells; counter++)
    {	
		Cell& c = cells[counter];
		ovito << c.ra.x << "\t" <<  c.ra.y << "\t" << c.rhoa << "\t" << "A" << "\t" << c.ncell << "\t" << (c.fna + c.fnb) / (NDIM) << "\t" << c.age << endl;
		ovito << c.rb.x << "\t" <<  c.rb.y << "\t" << c.rhob << "\t" << "B" << "\t" << c.ncell << "\t" << (c.fnb + c.fna) / (NDIM) << "\t" << c.age << endl;
    }
}

#define StrToNum(index, var) \
	istringstream ( inVars[index] ) >> var

void ReadInput() {

	ifstream infile("params.in");

	string varName, varString;
	
	while (infile >> varName >> varString) {
		inVars.push_back(varString);
	}
	
	StrToNum(0, Ncells);
	StrToNum(1, dt);
	StrToNum(2, D);
	StrToNum(3, lCell);
	StrToNum(4, gama);
	StrToNum(5, kCell);
	StrToNum(6, ka);
	StrToNum(7, kb);
	StrToNum(8, kr);
	StrToNum(9, r0),
	StrToNum(10, stepLimit);
	StrToNum(11, printStep);
	StrToNum(12, NcellsMax);
	StrToNum(13, R0fene);
	StrToNum(14, pDivision0);
	StrToNum(15, criticalPressure);
	StrToNum(16, reg.x);
	StrToNum(17, reg.y);
	StrToNum(18, kab);
	StrToNum(19, probDeath0);
	StrToNum(20, apoptisisAge);
	infile.close();
}

void writePropsHeader() {
	propsDat << "t" << "\t" << "<U>" << "\t" << "<P>" << endl;
}

void writeProps() {
	propsDat << timeNow << "\t" << uSum / Ncells << "\t" << virialPres() << endl;
}

// não conta para a pressão da célula
void feneForce() {
	int n;
	real sigLJ = lCell;
	real epsLJ = 1;
	
	DO_CELLS {
		
		VecR DR;
		GetDr( DR, cells[n].ra, cells[n].rb );
		real dr = sqrt(DotProd(DR, DR));
		
		real F = - kCell / ( 1 - pow(dr/R0fene, 2) );
		uSum += -0.5 * kCell * R0fene * R0fene * log( 1 - (dr/R0fene)*(dr/R0fene) );
		
		
		//if (dr <= pow(2, 1/6) * sigLJ) {
		// com apenas WCA pot. há uma certa instabilidade
		
		F += (48 * epsLJ / (sigLJ*sigLJ)) * ( pow(sigLJ/dr, 14) - 0.5 * pow(sigLJ/dr, 8) );
			uSum += ( 4 * epsLJ * ( pow(sigLJ/dr, 12) - pow(sigLJ/dr, 6) ) + epsLJ );
		
		//}
		
		AddForce(cells[n].fa, F*DR.x, F*DR.y);
		AddForce(cells[n].fb, -F*DR.x, -F*DR.y);
		
		// normal force for pressure
		//VecR fij, nVec;
		//fij.x = F*DR.x, fij.y = F*DR.y;
		//nVec.x = DR.x/dr; nVec.y = DR.y/dr;
		//cells[n].fna += abs(DotProd(fij, nVec));
		//cells[n].fnb += abs(DotProd(fij, nVec));
	
	}
}

void computeForcesHarmonic() { 
	int j1, j2, n;
    
    DO_CELLS {
		VZero (cells[n].fa);
		VZero (cells[n].fb);
		cells[n].fna = 0;
		cells[n].fnb = 0;
	}
    uSum = 0.;
    
    real rcut = 3;
    
	for (j1 = 0; j1 < Ncells - 1; j1 ++) {
		//cout << j1 << endl;
        for (j2 = j1 + 1; j2 < Ncells; j2 ++) {
			
			// *****************************************************************
			// check if any of the pole distances is lower than rcut
			// aa
			VecR DR_aa;
			GetDr( DR_aa, cells[j1].ra, cells[j2].ra );
			real dr_aa = sqrt(DotProd(DR_aa, DR_aa));
			
			// bb
			VecR DR_bb;
			GetDr( DR_bb, cells[j1].rb, cells[j2].rb );
			real dr_bb = sqrt(DotProd(DR_bb, DR_bb));
			
			// ab
			VecR DR_ab;
			GetDr( DR_ab, cells[j1].ra, cells[j2].rb );
			real dr_ab = sqrt(DotProd(DR_ab, DR_ab));
			
			// ba
			VecR DR_ba;
			GetDr( DR_ba, cells[j1].rb, cells[j2].ra );
			real dr_ba = sqrt(DotProd(DR_ba, DR_ba));
			
			// *****************************************************************
			// 
			
			real r_pole_min = min( min(dr_aa, dr_bb), min(dr_ab, dr_ba) );
			
			if (r_pole_min <= rcut) {
				
				// *************************************************************
				// apical-apical poles (mechanical repulsion)
				pairForceContactHarmonic(cells[j1].ra, cells[j2].ra, cells[j1].fa, cells[j2].fa, cells[j1].rhoa, cells[j2].rhoa, \
				cells[j1].fna, cells[j2].fna);
				
				// a-a atraction
				pairForceHarmonicLongRange(cells[j1].ra, cells[j2].ra, cells[j1].fa, cells[j2].fa, \
				cells[j1].rhoa, cells[j2].rhoa, ka, cells[j1].fna, cells[j2].fna);
				
				
				// *************************************************************
				// basal-basal
				pairForceContactHarmonic(cells[j1].rb, cells[j2].rb, cells[j1].fb, cells[j2].fb, cells[j1].rhob, cells[j2].rhob, \
				cells[j1].fnb, cells[j2].fnb);
				
				// b-b atract
				pairForceHarmonicLongRange(cells[j1].rb, cells[j2].rb, cells[j1].fb, cells[j2].fb, \
				cells[j1].rhob, cells[j2].rhob, kb, cells[j1].fnb, cells[j2].fnb);
				
				
				// *************************************************************
				// apical-basal repulsion
				pairForceContactHarmonic(cells[j1].ra, cells[j2].rb, cells[j1].fa, cells[j2].fb, cells[j1].rhoa, cells[j2].rhob, \
				cells[j1].fna, cells[j2].fnb);
				
				// a-b long-distance repulsion
				pairForceHarmonicLongRange(cells[j1].ra, cells[j2].rb, cells[j1].fa, cells[j2].fb, \
				cells[j1].rhoa, cells[j2].rhob, kab, cells[j1].fna, cells[j2].fnb);
				
				
				// *************************************************************
				// basal-apical
				pairForceContactHarmonic(cells[j1].rb, cells[j2].ra, cells[j1].fb, cells[j2].fa, cells[j1].rhob, cells[j2].rhoa, \
				cells[j1].fnb, cells[j2].fna);
				
				// b-a long-distance rep
				pairForceHarmonicLongRange(cells[j1].rb, cells[j2].ra, cells[j1].fb, cells[j2].fa, \
				cells[j1].rhob, cells[j2].rhoa, kab, cells[j1].fnb, cells[j2].fna);
				
			}
		}
	}
	
	// *************************************************************************
	// internal cell interactiom
	feneForce();
}


real virialPres() {
	 int n;
	 real P = 0.;
	 //real vol = reg.x * reg.y;
	 DO_CELLS {
		P += (cells[n].fna + cells[n].fnb) / (NDIM);
	 }
	 return P/Ncells;
}

void cellDivision() {
	
	int n;
	
	// TODO make pD dependend on pressure
	
	real pDivision = pDivision0;
	
	DO_CELLS {
		
		// TODO vol - use all cells volume or implement periodic boundary conditions
		real cellPressure = (cells[n].fna + cells[n].fnb) / (NDIM);
		
		if (GenRand01 < pDivision and Ncells < NcellsMax and cellPressure < criticalPressure) {
			
			real initialCellLength = lCell; // cells run the risk of exploding if it's low with WCA
			// TODO replace WCA in fene with a softer potencial, maybe harmonic
			
			
			// define new cell
			cells[Ncells].rhoa = r0, cells[Ncells].rhob = r0, cells[Ncells].ncell = Ncells;
			cells[n].rhoa = r0, cells[n].rhob = r0, cells[n].ncell = Ncells+1; // a célula "n" desaparece !!!!!
		
			// ra of 1st new cell is the ra of old cell
			cells[Ncells].ra.x = cells[n].ra.x;
			cells[Ncells].ra.y = cells[n].ra.y;
			// rb of 1st new cell is random
			real angulo = GenRand01 * (2 * 3.1415);
			cells[Ncells].rb.x = cells[Ncells].ra.x + initialCellLength * cos(angulo);
			cells[Ncells].rb.y = cells[Ncells].ra.y + initialCellLength * sin(angulo);
			
			
			// rb of 2nd new cell is the rb of the old cell
			// new ra of 2nd cell
			angulo = GenRand01 * (2 * 3.1415);
			cells[n].ra.x = cells[n].rb.x + initialCellLength * cos(angulo);
			cells[n].ra.y = cells[n].rb.y + initialCellLength * sin(angulo);
			
			VZero (cells[n].fa), VZero (cells[n].fb);
			VZero (cells[Ncells].fa), VZero (cells[Ncells].fb);
			
			Ncells++;
			Nparticles = Ncells*2;
			
			cells[Ncells].age = 0.;
			cells[n].age = 0.;
			
		} else cells[n].age += dt;
	}
}

void cellDeath() {

	int deaths = 0;
	if(!cells.empty()) {
	for(int i = Ncells - 1; i >= 0; i--) {
			if(cells.at(i).age > apoptisisAge and GenRand01 < probDeath0) {
				cells.erase( cells.begin() + i ); 
				deaths ++;
			}
		}
	}
	Ncells -= deaths;
	Nparticles = Ncells*2;
}


void PrintElapsedTime(chrono::steady_clock::time_point start) {
    auto end = chrono::steady_clock::now();
    
    if (chrono::duration_cast<chrono::seconds>(end - start).count() < 60) {
        cout << "Elapsed time : " \
            << chrono::duration_cast<chrono::seconds>(end - start).count() \
            << " sec";
    }
    else if (chrono::duration_cast<chrono::seconds>(end - start).count() > 60) {
        cout << "Elapsed time : " \
            << chrono::duration_cast<chrono::minutes>(end - start).count() << " min";
    }
    else if (chrono::duration_cast<chrono::minutes>(end - start).count() > 60) {
        cout << "Elapsed time : " \
            << chrono::duration_cast<chrono::hours>(end - start).count() << " hours";
    }
}
