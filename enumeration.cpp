#include "enumeration.h"
#include <cassert>
#include "printer.h"
#include "log.h"


EnumerationFilePrinter::EnumerationFilePrinter():
  initialisedFile(0),
  file(0),
  filename("")
{
}


EnumerationFilePrinter::~EnumerationFilePrinter()
{
  assert(file==0);
}


void EnumerationFilePrinter::open(string filename)
{
  this->filename=filename;
  string name=filename+extension();
  initialisedFile=fopen(name.c_str(),"w");
  file=initialisedFile;

  assert(file);

  fprintf(Stderr,"Output file opened: \"%s\"\n",name.c_str());

  onOpened();
}


void EnumerationFilePrinter::open(FILE *file)
{
  initialisedFile=0;
  this->file=file;

  assert(file);

  onOpened();
}


void EnumerationFilePrinter::close()
{
  onClose();

  if(initialisedFile)
    {
      fclose(initialisedFile);
      onClosed();
      initialisedFile=0;
    }
  file=0;
}


string EnumerationFilePrinter::extension()
{
  return "";
}


void EnumerationAlgorithm::printProgress(int step)
{
  while(step>0)
    {
      progressCounter++;
      //      if(!(progressCounter&4095))
	//      if(!(progressCounter&255))
      if(!(progressCounter&15))
	log2 fprintf(Stderr,"Number of Gr\"obner Bases found %i\n",progressCounter);
      fflush(Stderr);
      step--;
    }
}


//--------------------------------------
// EnumerationTargetCollector
//--------------------------------------

void EnumerationTargetCollector::beginEnumeration(PolynomialSet const &g)
{
  theList=PolynomialSetList();
}


void EnumerationTargetCollector::endEnumeration()
{
}


bool EnumerationTargetCollector::basis(const PolynomialSet &groebnerBasis)
{
  theList.push_back(groebnerBasis);
  return true;
}


PolynomialSetList EnumerationTargetCollector::getList()
{
  return theList;
}

#include "traverser_groebnerfan.h"

TargetGlue::TargetGlue(EnumerationTarget &target_):
	target(target_)
	{
	}

bool TargetGlue::process(ConeTraverser &traverser)
{
	GroebnerFanTraverser &r=dynamic_cast<GroebnerFanTraverser&>(traverser);
	return target.basis(r.refToGroebnerBasisRepresentation());
}
