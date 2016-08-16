#include "macaulay2.h"

#include <assert.h>
#include <unistd.h>
#include "printer.h"
#include "parser.h"
#include "programs.h"

//--------------------------------------------------
// Pipe
//--------------------------------------------------

Pipe::Pipe()
{
  int r;
  r=pipe(fdM2Input);
  assert(r==0);
  r=pipe(fdM2Output);
  assert(r==0);
  
  int pid1=fork();
  
  if(pid1==0)//child
    {
      fprintf(stderr,"Child\n");
      int r;
      r=dup2(fdM2Input[0],fileno(stdin));
      assert(r>=0);
      r=dup2(fdM2Output[1],fileno(stdout));
      assert(r>=0);
      
      static char M2Name[1024];
      sprintf(M2Name,"%s",programNameM2);//"M2";
      char *argv[2];
      argv[0]=M2Name;
      argv[1]=0;
      execvp(argv[0], argv);
      
      assert(0);
    }
  else //parent
    {
      fprintf(stderr,"Parent\n");
      
      pipeInput=fdopen(fdM2Input[1],"w");
      assert(pipeInput);
      pipeOutput=fdopen(fdM2Output[0],"r");
      assert(pipeOutput);
    }
}

Pipe::~Pipe()
{
  //fclose?
  fclose(pipeInput);
  fclose(pipeOutput);
}


//--------------------------------------------------
// Macaulay2Pipe
//--------------------------------------------------

Macaulay2Pipe::Macaulay2Pipe():
  printer(pipeInput)
{
  readLine();
}


Macaulay2Pipe::~Macaulay2Pipe()
{
  fprintf(pipeInput,"exit\n");

  /*  while(!feof(pipeOutput))
    {
      fprintf(stdout,"%c ",(getc(pipeOutput)));
    }
  */
}


void Macaulay2Pipe::skipStartOfLine()
{
  int c;
  do
    {
      c=fgetc(pipeOutput);
      if(c==EOF)
	{
	  fprintf(stderr,"EOF\n");
	}
      //      fprintf(stderr,"skipStartOfLine c=\'%c\',%i\n",c,c);
    }
  while(c!='=');
}


char *Macaulay2Pipe::readLine()
{
  static char line[2048];
  fgets(line,2048,pipeOutput);
  assert(line[strlen(line)-1]=='\n');
  //  fprintf(stderr,"readLine:%s",line);
  return line;
}


int Macaulay2Pipe::readInt()
{
  char *line=readLine();
  while(*line!=0 && *line!='=')line++;
  assert(*line!=0);
  int ret;
  sscanf(line+1,"%i",&ret);
  return ret;
}


bool Macaulay2Pipe::readBool()
{
  char *line=readLine();
  while(*line!=0 && *line!='=')line++;
  assert(*line!=0);
  int ret;
  return 'f'!=line[2];
}


int Macaulay2Pipe::getPdimCokerGensMonomial(gbasis const &monomialIdeal)
{
  fprintf(pipeInput,"pdim coker gens monomialIdeal ");
  printer.printMonomialIdeal(monomialIdeal);
  fprintf(pipeInput,"\n");
  fflush(pipeInput);
  
  readLine();
  readLine();
  int pdim=readInt();
  readLine();
  return pdim;
}


bool Macaulay2Pipe::isHomogeneousIdeal(gbasis &ideal)
{
  fprintf(pipeInput,"isHomogeneous ideal ");
  printer.printGroebnerBasis(ideal);
  fprintf(pipeInput,"\n");
  fflush(pipeInput);
  
  readLine();
  readLine();
  int pdim=readBool();
  readLine();
  return pdim;
}


void Macaulay2Pipe::setPolynomialRing(int numberOfVariables)
{
  fprintf(pipeInput,"R = ZZ/32003[");
  if(numberOfVariables>26)
    {
      fprintf(stderr,"Variable index out of range!\n");
      assert(0);
    }
  for(int i=0;i<numberOfVariables;i++)fprintf(pipeInput,(i==0)? "%c":",%c",i+'a');
  fprintf(pipeInput,"]\n");
  fflush(pipeInput);
  
  readLine();
  readLine();
  readLine();
  readLine();
  readLine();
  readLine();
}

void Macaulay2Pipe::setPolynomialRing(gbasis const &g)
{
  assert(g.begin()!=g.end());
  setPolynomialRing(g.begin()->size());
}


StandardPairList Macaulay2Pipe::standardPairs(gbasis &monomialIdeal)
{
  assert(monomialIdeal.begin()!=monomialIdeal.end());
  StandardPairList s;
  setPolynomialRing(monomialIdeal);

  fprintf(pipeInput,"toString standardPairs monomialIdeal ");
  printer.printMonomialIdeal(monomialIdeal);
  fprintf(pipeInput,"\n");
  fflush(pipeInput);
  
  skipStartOfLine();
  //  skipStartOfLine();
  return FileParser(pipeOutput).parseStandardPairList(monomialIdeal.begin()->size());
}
