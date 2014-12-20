/**********************************************************************
obrotate = rotate a tortional bond matched by a SMART pattern
Copyright (C) 2003 Fabien Fontaine
Some portions Copyright (C) 2004-2005 Geoffrey R. Hutchison
Some portions Copyright (C) 2008 Tim Vandermeersch
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

/*
  Require a SMART pattern, a file containing molecule coordinates
  4 atoms of the SMART pattern to define the tortional, an angle value
  The angle value must be in degree
  the 2 atoms of the rotating bond must be bonded but the 2 others not
  the part of the molecule on the side of the second atom is kept fixed
  whereas the part on the side of the third atom is rotated.
  example of command line:
  obrotate "[nH]ccccc[O,C][C,O]" test.sdf 1 6 7 8 180.0
*/


// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <pthread.h>
#include <fstream>
#include <sstream>

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/parsmart.h>
#include <openbabel/rotamer.h>
//#include <unistd.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>

#include "TemplateHDVector.hpp"

using namespace std;
using namespace OpenBabel;

static int sum=0;
char *FileIn = NULL;
float angle = 0;      // tortional angle value to set in degree
int angleSum = 0;
long double totalSum=1;
bool isEnergyCalcing = false;
OBConversion conv; //NF...
OBForceField *forceField=NULL;
MultiVector<double> *energySheet=NULL;
int* indexSheet=NULL;
int degrees=0;
pthread_mutex_t sum_lock = PTHREAD_MUTEX_INITIALIZER;

bool isContained(vector<int*>& mols,int* indexes)
{
	for(vector<int*>::iterator ind=mols.begin();ind!=mols.end();++ind)
	{
		bool set=true;
		for (int i=0;i<degrees;++i)
			if (indexes[i]!=(*ind)[i])
				set=false;
		if (set)
			return true;
	}
	return false;
}

void outputMol(vector<int*>& mols,OBMol& mol,vector< vector <int> >& maplist,int k)
{
	if (k<0)
		return;
	if (k>=maplist.size())
		return;

	OBAtom *a1, *a2, *a3, *a4;
	double energy;

	string fileName;
	ostringstream oss(fileName);
	for (int i=0;i<angleSum;++i) {
		fileName = "";
		a2 = mol.GetAtom(maplist[k][0]);
		a3 = mol.GetAtom(maplist[k][1]);
		std::vector<OBEdgeBase*>::iterator temp = a2->BeginBonds();
		a1 = a2->BeginNbrAtom(temp);
		if (a1==a3)
			a1 = a2->NextNbrAtom(temp);
		temp = a3->BeginBonds();
		a4 = a3->BeginNbrAtom(temp);
		if (a4==a2)
			a4 = a3->NextNbrAtom(temp);
		
		mol.SetTorsion(a1, a2, a3, a4, i * angle * DEG_TO_RAD);

		indexSheet[k]=i;
		if (isContained(mols,indexSheet))
		{
			oss.str("");
			pthread_mutex_lock(&sum_lock);
			oss << sum++ << "," << FileIn;
			pthread_mutex_unlock(&sum_lock);
			ofstream ofs(oss.str().c_str());
			cerr << "Outputing low energy file no." << sum << "/" << totalSum << endl;
			conv.Write(&mol,&ofs); //NF
		}
		outputMol(mols,mol,maplist,k-1);
	}
}

void turnMol(OBMol& mol,vector< vector <int> >& maplist,int k)
{
	if (k<0)
		return;
	if (k>=maplist.size())
		return;

	OBAtom *a1, *a2, *a3, *a4;
	double energy;

	string fileName;
	ostringstream oss(fileName);
	for (int i=0;i<angleSum;++i) {
		fileName = "";
		a2 = mol.GetAtom(maplist[k][0]);
		a3 = mol.GetAtom(maplist[k][1]);
		std::vector<OBEdgeBase*>::iterator temp = a2->BeginBonds();
		a1 = a2->BeginNbrAtom(temp);
		if (a1==a3)
			a1 = a2->NextNbrAtom(temp);
		temp = a3->BeginBonds();
		a4 = a3->BeginNbrAtom(temp);
		if (a4==a2)
			a4 = a3->NextNbrAtom(temp);
		if ( !a2->IsConnected(a3) ) {
			cerr << "obrotate: The atoms of the rotating bond must be bonded." << endl;
			exit(-1);
		}

		mol.SetTorsion(a1, a2, a3, a4, i * angle * DEG_TO_RAD);

		if (isEnergyCalcing)
		{
			indexSheet[k]=i;
			forceField->Setup(mol);
			energy = forceField->Energy(false);
			energySheet->getVectorValue(indexSheet)=energy;
		}
		else
		{
			oss.str("");
			pthread_mutex_lock(&sum_lock);
			oss << sum++ << "," << FileIn ;
			pthread_mutex_unlock(&sum_lock);
			ofstream ofs(oss.str().c_str());
			cerr << "Outputing file no." << sum << "/" << totalSum << endl;
			conv.Write(&mol,&ofs); //NF
		}
		turnMol(mol,maplist,k-1);
	}
}

class MolData
{
public:
	OBMol mol;
	vector< vector <int> >& maplist;
	int k;
	MolData(OBMol& molin,vector< vector <int> >& maplistin,int kin)
		: mol(molin)
		, maplist(maplistin)
		, k(kin)
	{
	}
};

void* turnMolSingle(void* MolDataIn)
{
	MolData &data=*((MolData*)MolDataIn);
	OBMol& mol=data.mol;
	vector< vector <int> >& maplist=data.maplist;
	int &k=data.k;
	
	turnMol(mol,maplist,k);

	delete &data;	
	return NULL;
}

void turnMolMT(OBMol& mol,vector< vector <int> >& maplist,int k)
{
	if (k<0)
		return;
	if (k>=maplist.size())
		return;

	OBAtom *a1, *a2, *a3, *a4;
	pthread_t threadList[angleSum];

	for (int i=0;i<angleSum;++i) {
		pthread_create(threadList+i,NULL,turnMolSingle,new MolData(mol,maplist,k-1));

		a2 = mol.GetAtom(maplist[k][0]);
		a3 = mol.GetAtom(maplist[k][1]);
		std::vector<OBEdgeBase*>::iterator temp = a2->BeginBonds();
		a1 = a2->BeginNbrAtom(temp);
		if (a1==a3)
			a1 = a2->NextNbrAtom(temp);
		temp = a3->BeginBonds();
		a4 = a3->BeginNbrAtom(temp);
		if (a4==a2)
			a4 = a3->NextNbrAtom(temp);
		mol.SetTorsion(a1, a2, a3, a4, i * angle * DEG_TO_RAD);
	}
	for (int i=0;i<angleSum;++i) {
		pthread_join(threadList[i],NULL);
	}
}

///////////////////////////////////////////////////////////////////////////////
//! \brief Set a tortional bond to a given angle
int main(int argc,char **argv)
{
  const char *Pattern=NULL;
  unsigned int i, t, errflg = 0;
  int c;
  char flags[255];
  string err;
  bool graphOutput=false;

  // parse the command line -- optional -a flag to change all matching torsions
  if (argc < 3 || argc > 4) {
    errflg++;
  } else {
    FileIn = argv[1];
    Pattern = "[!$(*#*)&!D1]-!@[!$(*#*)&!D1]";
    // Read the atom position
    c = sscanf(argv[2], "%d", &angleSum);
	angle = 360./angleSum;
    if (argc == 4)
	{
    		c = sscanf(argv[3], "%s", flags);
		int flagid=1;
    		while (flags[flagid]!=0)
			switch (flags[flagid++])
			{
			case 'g':
				graphOutput=true;
			case 'e':
				forceField=OBForceField::FindForceField("MMFF94");
				isEnergyCalcing=true;
				break;
    			}
 	}
  }
  if (errflg) {
    cerr << "Usage: rkrotate <filename> <angle> [options]" << endl;
    exit(-1);
  }

  // create pattern
  OBSmartsPattern sp;
  sp.Init(Pattern);

  OBFormat* format = conv.FormatFromExt(FileIn);
  if(!(format && conv.SetInAndOutFormats(format, format))) { //in and out formats same
    cerr << "obrotate: cannot read and/or write this file format!" << endl;
    exit (-1);
  } //...NF

  //Open the molecule file
  ifstream ifs;

  // Read the file
  ifs.open(FileIn);
  if (!ifs) {
    cerr << "obrotate: cannot read input file!" << endl;
    exit (-1);
  }

  OBMol mol;
  vector< vector <int> > maplist;      // list of matched atoms
//  vector< vector <int> >::iterator m;  // and its iterators
  //   int tindex;
  
  // Set the angles
  for (;;) {
    mol.Clear();
    //NF      ifs >> mol;                   // Read molecule
    conv.Read(&mol,&ifs); //NF
    if (mol.Empty())
      break;

    if (sp.Match(mol)) {          
      // if match perform rotation
      maplist = sp.GetUMapList(); // get unique matches
      
      if (maplist.size() > 1)
        cerr << "obrotate: Found " << maplist.size() << " matches." << endl;

	energySheet=new MultiVector<double>(degrees=maplist.size(),angleSum);
	indexSheet=new int[maplist.size()];

      for (int EXO=0;EXO<maplist.size();++EXO)
	totalSum*=angleSum+EXO;
      // look at all the mapping atom but save only the last one.
	turnMol(mol,maplist,maplist.size()-1);
	
      if (graphOutput)
      {
	ofstream ofs("energyGraph.mlog");
	int ind[degrees];
	for (int i=0;i<degrees;++i)
		ind[i]=0;
	do
	{
		for (int i=0;i<degrees;++i)
			ofs<<ind[i]<<'\t';
		ofs<<energySheet->getVectorValue(ind)<<endl;
	}
	while(energySheet->incressIndex(ind));
      }

	if (isEnergyCalcing)
	{
		std::vector<int*> lowEnergySheet;
		totalSum=energySheet->getMinValues(lowEnergySheet);
		if (totalSum)
			outputMol(lowEnergySheet,mol,maplist,maplist.size()-1);
		else
			cerr << "rkrotate: No low energy conformation found." << endl;
	}

	cout << sum;
    } else {
      cerr << "obrotate: Found 0 matches for the SMARTS pattern." << endl;
      exit(-1);
    }
    //NF      cout << mol;
  }

  return(0);
}


/* obrotate man page*/
/** \page obrotate batch-rotate dihedral angles matching SMARTS patterns
*
* \n
* \par SYNOPSIS
*
* \b obrotate '<SMARTS-pattern>' \<filename\> \<atom1\> \<atom2\> \<atom3\> \<atom4\> \<angle\>
*
* \par DESCRIPTION
*
* The obrotate program rotates the torsional (dihedral) angle of a specified 
* bond in molecules to that defined by the user. In other words, it does the
* same as a user setting an angle in a molecular modelling package, but much 
* faster and in batch mode.
* \n\n
* The four atom IDs required are indexes into the SMARTS pattern, which starts
* at atom 1. The angle supplied is in degrees. The two atoms used to set
* the dihedral angle \<atom1\> and \<atom4\> do not need to be connected 
* to the atoms of the bond \<atom2\> and \<atom3\> in any way.
*\n\n
* The order of the atoms matters -- the portion of the molecule attached to
* \<atom1\> and \<atom2\> remain fixed, but the portion bonded to \<atom3\> and
& \<atom4\> moves.
* 
* \par EXAMPLES
*  - Let's say that you want to define the conformation of a large number of
*  molecules with a pyridyl scaffold and substituted with an aliphatic chain
*  at the 3-position, for example for docking or 3D-QSAR purposes.
* \n\n
*    To set the value of the first dihedral angle to 90 degrees:\n
*   obrotate "c1ccncc1CCC" pyridines.sdf 5 6 7 8 90
* \n
* Here 6 and 7 define the bond to rotate in the SMARTS patter, i.e., c1-C and 
* atoms 5 and 8 define the particular dihedral angle to rotate.
*  - Since the atoms to define the dihedral do not need to be directly
*  connected, the nitrogen in the pyridine can be used:\n
*   obrotate "c1ccncc1CCC" pyridines.sdf 4 6 7 8 90
*
*  - Keep the pyridyl ring fixed and moves the aliphatic chain:\n
*   obrotate "c1ccncc1CCC" pyridines.sdf 5 6 7 8 90
*  - Keep the aliphatic chain fixed and move the pyridyl ring:\n
*   obrotate "c1ccncc1CCC" pyridines.sdf 8 7 6 5 90
*
* \par AUTHORS
*
* The obrotate program was contributed by \b Fabien \b Fontaine.
*
* Open Babel is currently maintained by \b Geoff \b Hutchison, \b Chris \b Morley and \b Michael \b Banck.
*
* For more contributors to Open Babel, see http://openbabel.org/THANKS.shtml
*
* \par COPYRIGHT
*  Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
*  Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison \n \n
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation version 2 of the License.\n \n
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
* \par SEE ALSO
*   The web pages for Open Babel can be found at: http://openbabel.org/ \n
*   A guide for constructing SMARTS patterns can be found at: http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
**/
