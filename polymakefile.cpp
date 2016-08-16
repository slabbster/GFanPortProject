#include "polymakefile.h"
#include "printer.h"
#include "log.h"

#include <assert.h>
#include <sstream>

bool PolymakeFile::xmlForced=false;

static void eatComment2(int &c, stringstream &s)
{
  if(c=='#')
    {
      do
	c=s.get();
      while(c!='\n' && !s.eof());
    }
}

static void eatComment(stringstream &s)
{
  int c=s.get();
  while(c==' '||c=='\t')c=s.get();
  eatComment2(c,s);
  s.unget();
}

PolymakeProperty::PolymakeProperty(const string &name_, const string &value_, bool embedded_):
  name(name_),
  value(value_),
  embedded(embedded_)
{
}


static string readUntil(FILE *f, int c)
{
  stringstream ret;
  int c2;
  c2=fgetc(f);
  while(c2!=c && c2!=EOF)
    {
      ret<<char(c2);
      c2=fgetc(f);
    }
  return ret.str();
}


const char *PolymakeFile::mapToPolymakeNames(const char *s)
{
  if(!isXml)return s;
  if(strcmp(s,"CONES_ORBITS")==0)return "CONES_REPS";
  if(strcmp(s,"MAXIMAL_CONES_ORBITS")==0)return "MAXIMAL_CONES_REPS";

  return s;
}


void PolymakeFile::open(const char *fileName_)
{
  isXml=false;
  fileName=string(fileName_);

  FILE *f=fopen(fileName.c_str(),"r");
  if(!f)fprintf(Stderr,"Could not open file:\"%s\"\n",fileName_);
  assert(f);

  int c=fgetc(f);
  while(c!=EOF)
    {
      if(c=='_')
	{
	  readUntil(f,'\n');
	}
      else if(c!='\n')
	{
	  ungetc(c,f);
	  string name=readUntil(f,'\n');

log1	  fprintf(Stderr,"Reading:\"%s\"\n",name.c_str());
	  stringstream value;
	  while(1)
	    {
	      string l=readUntil(f,'\n');
	      if(l.size()==0)break;
	      value << l <<endl;
	    }
	  properties.push_back(PolymakeProperty(name.c_str(),value.str().c_str()));
	}
      c=fgetc(f);
    }
}


void PolymakeFile::create(const char *fileName_, const char *application_, const char *type_, bool isXml_)
{
  fileName=string(fileName_);
  application=string(application_);
  type=string(type_);
  isXml=isXml_;
  if(xmlForced){isXml=true;}
}



void PolymakeFile::close()
{
  FILE *f=fopen(fileName.c_str(),"w");
  assert(f);
  if(isXml)
    {
      fprintf(f,"<properties>\n");

      for(list<PolymakeProperty>::const_iterator i=properties.begin();i!=properties.end();i++)
	{
	  if(i->embedded)
	    {
              fprintf(f,"<property name=\"%s\" value=\"%s\" />\n",mapToPolymakeNames(i->name.c_str()),i->value.c_str());
	    }
	  else
	    {
	      fprintf(f,"<property name=\"%s\">\n",mapToPolymakeNames(i->name.c_str()));
	      fprintf(f,"%s",i->value.c_str());
	      fprintf(f,"</property>\n");
	    }
	 }
      fprintf(f,"</properties>\n");
    }
  else
    {
      fprintf(f,"_application %s\n",application.c_str());
      fprintf(f,"_version 2.2\n");
      fprintf(f,"_type %s\n",type.c_str());

      for(list<PolymakeProperty>::const_iterator i=properties.begin();i!=properties.end();i++)
	{
	  fprintf(f,"\n%s\n",i->name.c_str());
	  fprintf(f,"%s",i->value.c_str());
	}
    }
  fclose(f);
}


void PolymakeFile::writeStream(ostream &file)
{
  if(isXml)
    {
      file << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
      file << "<object name=\"1\" type=\""<<application<<"::"<<type<<"&lt;Rational&gt;\" version=\"2.9.9\" xmlns=\"http://www.math.tu-berlin.de/polymake/#3\">\n";

      for(list<PolymakeProperty>::const_iterator i=properties.begin();i!=properties.end();i++)
	{
	  if(i->embedded)
	    {
              file << "<property name=\""<<mapToPolymakeNames(i->name.c_str())<<  "\" value=\""<<i->value.c_str()<<"\" />\n";
	    }
	  else
	    {
	      file << "<property name=\"" << mapToPolymakeNames(i->name.c_str()) << "\">\n";
	      file << i->value.c_str();
	      file << "</property>\n";
	    }
	}
      file << "</object>\n";
    }
  else
    {
      file << "_application " << application << endl;
      file << "_version 2.2\n";
      file << "_type " << type << endl;

      for(list<PolymakeProperty>::const_iterator i=properties.begin();i!=properties.end();i++)
	{
	  file << endl << i->name.c_str() << endl;
	  file << i->value;
	}
    }
}


list<PolymakeProperty>::iterator PolymakeFile::findProperty(const char *p)
{
  string s(p);
  for(list<PolymakeProperty>::iterator i=properties.begin();i!=properties.end();i++)
    {
      if(s==i->name)return i;
    }

  return properties.end();
}


void PolymakeFile::writeProperty(const char *p, const string &data, bool embedded)
{
  if(hasProperty(p))
    {
      assert(0);
    }
  properties.push_back(PolymakeProperty(string(p),data,embedded));
}


bool PolymakeFile::hasProperty(const char *p, bool doAssert)
{
  if(doAssert)
    if(findProperty(p)==properties.end())
      {
	fprintf(stderr,"Property: \"%s\" not found in file.\n",p);
	assert(0);
      }

  return findProperty(p)!=properties.end();
}


int PolymakeFile::readCardinalProperty(const char *p)
{
  assert(hasProperty(p,true));

  list<PolymakeProperty>::iterator prop=findProperty(p);
  stringstream s(prop->value);

  int ret;
  s>>ret;

  return ret;
}


void PolymakeFile::writeCardinalProperty(const char *p, int n)
{
  stringstream t;

  if(isXml)
    {
      t<<n;
      writeProperty(p,t.str(),true);
    }
  else
    {
      t<<n<<endl;
      writeProperty(p,t.str());
    }
}


bool PolymakeFile::readBooleanProperty(const char *p)
{
  return false;
}


void PolymakeFile::writeBooleanProperty(const char *p, bool n)
{
  stringstream t;

if(isXml)
    {
      t<<(n?"true":"false");
      writeProperty(p,t.str(),true);
    }
  else
    {
      t<<n<<endl;
      writeProperty(p,t.str());
    }
}


IntegerMatrix PolymakeFile::readMatrixProperty(const char *p, int height, int width)
{
  IntegerMatrix ret(height,width);

  assert(hasProperty(p,true));
  list<PolymakeProperty>::iterator prop=findProperty(p);
  stringstream s(prop->value);
  for(int i=0;i<height;i++)
    for(int j=0;j<width;j++)
      {
	int v;
	eatComment(s);
	s>>v;
	ret[i][j]=v;
      }

  return ret;
}


void PolymakeFile::writeMatrixProperty(const char *p, const IntegerMatrix &m, bool indexed, const vector<string> *comments)
{
  stringstream t;

  if(comments)assert(comments->size()>=m.getHeight());
  if(isXml)
    {
      t << "<m>\n";
      for(int i=0;i<m.getHeight();i++)
	{
	  t << "<v>";
	  for(int j=0;j<m.getWidth();j++)
	    {
	      if(j!=0)t<<" ";
	      t<<m[i][j];
	    }
	  t << "</v>\n";
	}
      t << "</m>\n";
    }
  else
    {
      for(int i=0;i<m.getHeight();i++)
	{
	  for(int j=0;j<m.getWidth();j++)
	    {
	      if(j!=0)t<<" ";
	      t<<m[i][j];
	    }
	  if(indexed)t<<"\t# "<<i;
	  if(comments)t<<"\t# "<<(*comments)[i];
	  t<<endl;
	}
    }
  writeProperty(p,t.str());
}

IntegerMatrix PolymakeFile::readArrayArrayIntProperty(const char *p, int width)
{
  assert(0);//Not implemented yet.
}


void PolymakeFile::writeArrayArrayIntProperty(const char *p, const IntegerMatrix &m)
{
  stringstream t;

  if(isXml)
    {
      t << "<m>\n";
      for(int i=0;i<m.getHeight();i++)
        {
          t << "<v>";
          for(int j=0;j<m.getWidth();j++)
            {
              if(j!=0)t<<" ";
              t<<m[i][j];
            }
          t << "</v>\n";
        }
      t << "</m>\n";
    }
  else
    {
      for(int i=0;i<m.getHeight();i++)
        {
          for(int j=0;j<m.getWidth();j++)
            {
              if(j!=0)t<<" ";
              t<<m[i][j];
            }
          t<<endl;
        }
    }
  writeProperty(p,t.str());
}


static list<int> readIntList(istream &s)
{
  list<int> ret;
  int c=s.peek();
  while((c>='0') && (c<='9')|| (c==' '))
    {
      //      fprintf(Stderr,"?\n");
      int r;
      s >> r;
      ret.push_back(r);
      c=s.peek();
    }
  return ret;
}

vector<list<int> > PolymakeFile::readMatrixIncidenceProperty(const char *p)
{
  vector<list<int> > ret;
  assert(hasProperty(p,true));
  list<PolymakeProperty>::iterator prop=findProperty(p);
  stringstream s(prop->value);

  while((s.peek()!=-1)&&(s.peek()!='\n')&&(s.peek()!=0))
    {
      //      fprintf(Stderr,"!\n");
      int c=s.get();
      //fprintf(Stderr,"%i",c);
      assert(c=='{');
      ret.push_back(readIntList(s));
      c=s.get();
      assert(c=='}');
      c=s.get();
      while(c==' ' || c=='\t')c=s.get();
      eatComment2(c,s);
      assert(c=='\n');
    }
  return ret;
}


void PolymakeFile::writeIncidenceMatrixProperty(const char *p, const vector<list<int> > &m, int baseSetSize)
{
  stringstream t;

  if(isXml)
    {
      t<<"<m cols=\""<<baseSetSize<<"\">";
      for(int i=0;i<m.size();i++)
	{
	  t<<"<v>";
	  list<int> temp=m[i];
	  temp.sort();
	  for(list<int>::const_iterator j=temp.begin();j!=temp.end();j++)
	    {
	      if(j!=temp.begin())t<<' ';
	      t<< *j;
	    }
	  t<<"<v>\n"<<endl;
	}
      t<<"</m>\n";
    }
  else
    {
      for(int i=0;i<m.size();i++)
	{
	  t<<'{';
	  list<int> temp=m[i];
	  temp.sort();
	  for(list<int>::const_iterator j=temp.begin();j!=temp.end();j++)
	    {
	      if(j!=temp.begin())t<<' ';
	      t<< *j;
	    }
	  t<<'}'<<endl;
	}
    }
  writeProperty(p,t.str());
}


IntegerVector PolymakeFile::readCardinalVectorProperty(const char *p)
{
  assert(hasProperty(p,true));
  list<PolymakeProperty>::iterator prop=findProperty(p);
  stringstream s(prop->value);

  list<int> temp=readIntList(s);

  IntegerVector ret(temp.size());
  int I=0;
  for(list<int>::const_iterator i=temp.begin();i!=temp.end();i++,I++)ret[I]=*i;

  return ret;
}


void PolymakeFile::writeCardinalVectorProperty(const char *p, IntegerVector const &v)
{
  stringstream t;

  if(isXml)
    {
      t<<"<v>";
      for(int i=0;i<v.size();i++)
	{
	  if(i!=0)t<<" ";
	  t<<v[i];
	}
      t<<"</v>\n";
    }
  else
    {
      for(int i=0;i<v.size();i++)
	{
	  if(i!=0)t<<" ";
	  t<<v[i];
	}
      t<<endl;
    }
  writeProperty(p,t.str());
}


void PolymakeFile::writeStringProperty(const char *p, const string &s)
{
  if(isXml)
    writeProperty(p,s);
  else
    writeProperty(p,s);
}
