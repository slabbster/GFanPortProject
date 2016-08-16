#include "application.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <list>

#include "field.h"  //REMOVE THIS INCLUDE

using namespace std;

#define APPLICATIONNAME "gfan"

Application *Application::applicationList;
Application::Option *Application::Option::constructionList;


static char *findName(char *p)
{
  // extracts filename from argv[0]
  int l=strlen(p);
  assert(l>0);
  do
    {
      l--;
    }
  while(p[l]!=0 && p[l]!='/');

  return p+l+1;
}

static char* tail(char *p)
{
  const char *n=APPLICATIONNAME;

  const char *m=n;
  assert(p);
  while(*m)
    {
      assert(*p==*m);
      m++;
      p++;
    }
  return p;
}

//--------------------------------------------------
// Main - starts an application
//--------------------------------------------------

int main(int argc, char *argv[])
{
  bool helpOption=false;
  if(argc==2)if(strcmp(argv[1],"--help")==0)helpOption=true;
  if(argc==3)if(strcmp(argv[2],"--help")==0)helpOption=true;
  /*  if(argc==2)
    {
      if(strcmp(argv[1],"installlinks")==0)
	{
	  Application::makeSymbolicLinks(APPLICATIONNAME,false,"");
	  return 0;
	}
      if(strcmp(argv[1],"installlinksall")==0)
	{
	  Application::makeSymbolicLinks(APPLICATIONNAME,true,"");
	  return 0;
	}
      if(strcmp(argv[1],"documentation")==0)
	{
	  Application::produceLatexDocumentation(false);
	  return 0;
	}
      if(strcmp(argv[1],"--help")==0)helpOption=true;
    }
  if(argc==3)
    {
      if(strcmp(argv[1],"installlinks")==0)
	{
	  Application::makeSymbolicLinks(APPLICATIONNAME,false,argv[2]);
	  return 0;
	}
      if(strcmp(argv[1],"installlinksall")==0)
	{
	  Application::makeSymbolicLinks(APPLICATIONNAME,true,argv[2]);
	  return 0;
	}
      if(strcmp(argv[2],"--help")==0)helpOption=true;
    }
  */
  //  fprintf(Stderr,"Extracted name=%s",);
  //  Application *app=Application::applicationList;
  int argumentsToSkip=0;
  Application *app=0;
  Application *app2=0;
  if(argc>1)app2=Application::findApplication(argv[1]);
  Application *app3=Application::findApplication(tail(findName(argv[0])));

  if(app2)
    {
      argumentsToSkip=1;
      app=app2;
    }
  else
    app=app3;

  if(app==0)
    {
      fprintf(stderr,"Application not found!\n");
      assert(0);

      return 0;
    }

  /*  if((!app) || app->next)
    {
      fprintf(Stderr,app?"Error: multiple applications defined.\n"
	                :"Error: no applications defined.\n");
      assert(0);
      return 0;
      }*/

  if(helpOption)
    {
      app->printHelp();
      return 0;
    }

  if(app->parseOptions(argc,argv,argumentsToSkip))
    {
      int ret=app->main();
      //      fprintf(Stderr,"Number of rationals living:%i\n",FieldElementRationalsLiving);
      app->onExit();
      return ret;
    }


  return 0;
}


//--------------------------------------------------
// A few internal applications for installation
// and generating documentation.
//--------------------------------------------------

class DocumentationApplication : public Application
{
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "Internal application for generating the LaTeX documentation.\n";
  }
  const char *name()
  {
    return "documentation";
  }
  int main()
  {
    Application::produceLatexDocumentation(false);
    return 0;
  }
};

DocumentationApplication theDocumentationApplication;


class InstallationApplication : public Application
{
  SimpleOption optionAll;
  StringOption optionPath;
  bool includeInDefaultInstallation()
  {
    return false;
  }
public:
  InstallationApplication():
    optionAll("--all","Install all links to all programs. Also the ones only known by the author.\n"),
    optionPath("--path","Specify the installation path.\n","")
  {
    registerOptions();
  }
  const char *helpText()
  {
    return "Internal application for installing symboilic links.\n";
  }
  const char *name()
  {
    return "installlinks";
  }
  int main()
  {
    Application::makeSymbolicLinks(APPLICATIONNAME,optionAll.getValue(),optionPath.getValue());
    return 0;
  }
};

InstallationApplication theInstallationApplication;


//--------------------------------------------------
// Application::Option
//--------------------------------------------------

Application::Option::Option()
{
  bool hidden=false;
  next=constructionList;constructionList=this;
}


Application::Option *Application::Option::getOptionList()
{
  Option *ret=constructionList;constructionList=0;return ret;
}


void Application::Option::onOptionsParsed()
{
}


bool Application::Option::isHidden()const
{
  return hidden;
}


void Application::Option::hide(bool b)
{
  hidden=b;
}


//--------------------------------------------------
// Application::StringMatchingOption
//--------------------------------------------------

bool Application::StringMatchingOption::isExactMatch(const char *s)
{
  return std::string(s)==std::string(matchString);
}


Application::StringMatchingOption::StringMatchingOption(const char *s, const char *description_):
  matchString(s),
  description(description_)
{
}


bool Application::StringMatchingOption::matches(const char *s)
{
  for(int i=0;matchString[i];i++)
    if(matchString[i]!=s[i])return false;
  return true;
}


std::string Application::StringMatchingOption::documentationMatchString()
{
  return std::string(matchString);
}


std::string Application::StringMatchingOption::documentationDescription()
{
  return description;
}

//--------------------------------------------------
// Application::SimpleOption
//--------------------------------------------------

Application::SimpleOption::SimpleOption(const char *s, const char *description):
  StringMatchingOption(s, description),
  value(false)
{
}


void Application::SimpleOption::parseValue(int t, char **argv, bool *ok, int *numberOfArgumentsParsed)
{
  *numberOfArgumentsParsed=0;
  if(isExactMatch(argv[t]))
    {
      value=true;
      *ok=true;
    }
  else
    *ok=false;
}


bool Application::SimpleOption::getValue()
{
  return value;
}


//--------------------------------------------------
// Application::ValueOption
//--------------------------------------------------

Application::ValueOption::ValueOption(const char *s, const char *description):
  StringMatchingOption(s, description)
{
}


std::string Application::ValueOption::documentationMatchString()
{
  return std::string(matchString)+" value";
}


void Application::ValueOption::parseValue(int t, char **argv, bool *ok, int *numberOfArgumentsParsed)
{
  if(isExactMatch(argv[t]))
    {
      if(argv[t+1]==0 || argv[t+1][0]=='-')
        {
          assignValue("");
          *ok=true;
          *numberOfArgumentsParsed=0;
        }
      else
        {
          assignValue(argv[t+1]);
          if(argv[t+1])
            {
              *ok=true;
              *numberOfArgumentsParsed=1;
            }
          else
            {
              *ok=false;
              *numberOfArgumentsParsed=0;
            }
	}
    }
  else
    {
      assignValue(argv[t]+std::string(matchString).length());
      *ok=true;
      *numberOfArgumentsParsed=0;
    }
}


//--------------------------------------------------
// Application::StringOption
//--------------------------------------------------

Application::StringOption::StringOption(const char *s, const char *description, const char *initialValue):
  ValueOption(s, description),
  value(initialValue)
{
}


void Application::StringOption::assignValue(const char *s)
{
  value=s;
}


const char *Application::StringOption::getValue()
{
  return value;
}


//--------------------------------------------------
// Application::IntegerOption
//--------------------------------------------------

Application::IntegerOption::IntegerOption(const char *s, const char *description, int initialValue):
  ValueOption(s, description),
  value(initialValue),
  hasRange(false)
{
}


Application::IntegerOption::IntegerOption(const char *s, const char *description, int initialValue, int lower_, int upper_):
  ValueOption(s, description),
  value(initialValue),
  hasRange(true),
  lower(lower_),
  upper(upper_)
{
}


void Application::IntegerOption::assignValue(const char *s)
{
  bool onlyDigits=true;
  value=0;
  for(int t=0;s[t];t++)
    {
      if(s[t]<'0' || s[t]>'9')onlyDigits=false;
      value*=10;
      value+=s[t]-'0';
    }
  assert(onlyDigits);

  if(hasRange)
    if(lower>value || upper<value)
      {
	fprintf(stderr,"Integer option %s has out of range value. Range={%i,...,%i}\n",matchString,lower,upper);
	assert(0);
      }
}


int Application::IntegerOption::getValue()
{
  return value;
}


//--------------------------------------------------
// Application::ZeroOneOption
//--------------------------------------------------

Application::ZeroOneOption::ZeroOneOption(const char *s, const char *description, int initialValue):
  IntegerOption(s, description,initialValue,0,1)
{
}


//--------------------------------------------------
// Application
//--------------------------------------------------

Application::Application()
{
  next=applicationList;applicationList=this;optionList=Option::getOptionList();
}


bool Application::includeInDefaultInstallation()
{
  return true;
}


bool Application::parseOptions(int argc, char **argv, int argumentsToSkip)
{
  for(int t=1+argumentsToSkip;t<argc;t++)
    {
      int numberOfMatchesFound=0;

      for(Option *i=optionList;i;i=i->next)
	{
	  if(i->matches(argv[t]))numberOfMatchesFound++;
	}
      //      fprintf(Stderr,"NumberOfMatches %i\n",numberOfMatchesFound);
      assert(numberOfMatchesFound<2);
      if(numberOfMatchesFound==0)
	{
	  fprintf(stderr,"UNKNOWN OPTION: %s.\n",argv[t]);
	  fprintf(stderr,"USE --help AS A SINGLE OPTION TO VIEW THE HELP TEXT.\n");
	  return false;
	}

      for(Option *i=optionList;i;i=i->next)
	if(i->matches(argv[t]))
	  {
	    bool ok=false;
	    int argumentsParsed=1;
	    i->parseValue(t,argv,&ok,&argumentsParsed);
	    if(!ok)
	      {
		fprintf(stderr,"PARSE ERROR PARSING ARGUMENTS IN OPTION %s.\n",argv[t]);
		fprintf(stderr,"USE --help AS A SINGLE OPTION TO VIEW THE HELP TEXT.\n");
		//		assert(0);
		//		exit(0);
		return false;
	      }
	    t+=argumentsParsed;
	    break;
	  }
    }
  for(Option *i=optionList;i;i=i->next)i->onOptionsParsed();

  return true;
}


void Application::registerOptions()
{//merge current list of options with new list
  Option *l=Option::getOptionList();
  while(l)
    {
      Option *o=l;
      l=o->next;
      o->next=optionList;
      optionList=o;
    }
}


void Application::printHelp()
{
  fprintf(stderr,"%s",helpText());

  Application *p=applicationList;

  FILE *f=stderr;
  Option *l=optionList;
  if(l)
    {
      fprintf(f,"Options:\n");

      while(l)
	{
	  if(!l->isHidden())fprintf(f,"%s:\n %s\n",l->documentationMatchString().c_str(),l->documentationDescription().c_str());

	  l=l->next;
	}
    }
};


class Application *Application::findApplication(char *name)
{
  Application *p=applicationList;
  while(p)
    {
      if(strcmp(name,p->name())==0)return p;
      p=p->next;
    }
  return 0;
}

void Application::makeSymbolicLinks(const char *name, bool all, const char *path)
{
  Application *p=applicationList;
  while(p)
    {
      if(all || p->includeInDefaultInstallation())
	if(strlen(p->name())>0)
	  {
	    char c[1024];
	    sprintf(c,"ln -s %s%s %s%s%s\n",path,name,path,name,p->name());
	    fprintf(stderr,"%s",c);
	    system(c);
	  }
      p=p->next;
    }
}


static int substituteSingleString(FILE *f,const char* s, const char* pattern, const char* substitute)
{
  int n=0;
  while(*pattern)
    {
      if(*s!=*pattern)return 0;
      pattern++;
      s++;
      n++;
    }
  fprintf(f,"%s",substitute);
  return n;
}

static void quoteLatexPrint(FILE *f, const char *s)
{
  while(*s)
    {
      s+=substituteSingleString(f,s,"\\omega","\\omega");
      if(s[0]=='-' && s[1]=='-')
	{
	  fprintf(f,"-\\hspace{0.013cm}-");
	  s++;
	}
      else if(s[0]=='G' && s[1]=='r' && s[2]=='o' && s[3]=='e')
	{
	  fprintf(f,"Gr\\\"o");
	  s+=3;
	}
      else if(s[0]=='~')
	fprintf(f,"\\~{}");
      else if(s[0]=='{')
	fprintf(f,"\\{");
      else if(s[0]=='}')
	fprintf(f,"\\}");
      else if(s[0]=='_')
	fprintf(f,"\\_");
      else if(s[0]=='\\')
	fprintf(f,"\\backslash");
      else if(s[0]=='<')
	fprintf(f,"\\symbol{60}");
      else
	fprintf(f,"%c",s[0]);
      s++;
    }
}

static bool compare_appname(Application *a, Application *b)
{
  return string(a->name())<string(b->name());
}

std::list<Application*> Application::getSortedApplicationList()
{
  Application *p=applicationList;
  list<Application*> alist2;
  while(p)
    {
      alist2.push_back(p);
      p=p->next;
    }

  alist2.sort(compare_appname);
  return alist2;
}

void Application::produceLatexDocumentation(bool all)
{
  FILE *f=stdout;

  std::list<Application*> alist2=getSortedApplicationList();

  //  while(p)
  for(list<Application*>::const_iterator i=alist2.begin();i!=alist2.end();i++)
    {
      Application *p=*i;
      if(all || p->includeInDefaultInstallation())
	{
	  fprintf(f,"{\\subsection{%s",APPLICATIONNAME);
	  quoteLatexPrint(f,p->name());
	  fprintf(f,"}");

	  fprintf(f,"\\label{applist:%s}\n",p->name());
	  quoteLatexPrint(f,p->helpText());

	  Option *l=p->optionList;

	  bool containsNonHidden=false;
	  Option *l2=l;
	  while(l2){if(!l2->isHidden())containsNonHidden=true;l2=l2->next;}
	  if(containsNonHidden)
	    {
	      fprintf(f,"\\newline\n");
	      fprintf(f,"{\\bf Options:}\n");

	      fprintf(f,"\\begin{description}\n");
	      while(l)
		{
		  if(!l->isHidden())
		    {
		      fprintf(f,"\\item[");
		      quoteLatexPrint(f,l->documentationMatchString().c_str());
		      fprintf(f,"]");
		      quoteLatexPrint(f,l->documentationDescription().c_str());
		    }
		  l=l->next;
		}
	      fprintf(f,"\\end{description}\n");
	    }
          fprintf(f,"\n\n");
	}
      //      p=p->next;
    }
}


void Application::onExit()
{
}
