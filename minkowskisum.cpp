#include "minkowskisum.h"
#include "printer.h"
#include "parser.h"

#define MINKOWSKIFILEINPUT "minkowski.input"
#define MINKOWSKIFILEOUTPUT "minkowski.output"
#define MINKOWSKIPROGRAM "~/math/software/minkowski/TOPCOM-0.13.2/weibel/essai"

static IntegerVector const& index(IntegerVectorList const &l, int index)
{
  for(IntegerVectorList::const_iterator i=l.begin();i!=l.end();i++)
    {
      if(index==0)return *i;
      index--;
    }
  assert(0);
  return *l.end();
}

IntegerVectorList minkowski(IntegerVectorListList const &polytopes, IntegerVectorListList *sums)
{
  if(sums)*sums=IntegerVectorListList();

  IntegerVectorList ret;
  {
    FILE *f=fopen(MINKOWSKIFILEINPUT,"w");

    assert(f);
    {
      fprintf(f,"%i",polytopes.size());
      TopcomPrinter p(f);
      //AsciiPrinter p(f);

      for(IntegerVectorListList::const_iterator i=polytopes.begin();i!=polytopes.end();i++)
	p.printVectorList(*i);
    }
    fclose(f);
  }
  system(MINKOWSKIPROGRAM" -v <" MINKOWSKIFILEINPUT " >" MINKOWSKIFILEOUTPUT);

  {
    FILE *f=fopen(MINKOWSKIFILEOUTPUT,"r");
    assert(f);


    {
      FileParser p(f);
      while(p.nextNonBlankDoNotGet())
	{
	  IntegerVector indices=p.parseIntegerVector();
	  if(sums)
	    {
	      IntegerVectorList sum;
	      int i=0;
	      for(IntegerVectorListList::const_iterator I=polytopes.begin();I!=polytopes.end();I++)
		{
		  sum.push_back(index(*I,indices[i]-1));
		  i++;
		}

	      sums->push_back(sum);
	    }
	  int c=p.nextNonBlank();
	  assert(c==':');
	  ret.push_back(p.parseIntegerVector());
	  c=p.nextNonBlank();
	  while(c!=']')c=p.nextNonBlank();
	}
    }
    fclose(f);
  }
  return ret;
}
