#ifndef APPLICATION_H_INCLUDED
#define APPLICATION_H_INCLUDED

#include <string>
#include <list>

// A general application.
// This could solve several problems:
// - the help text could be automatically generated
// - duplicate main code could be removed
// - it would be easy to parse options
// - all applications could be linked to a single file thereby saving disk-space


/**
   @brief A command line application.

   The class Application describes the general command line
   application. The compiled Gfan executable consists of several
   command line applications. When the static main() is run it is
   decided which application to invoke by investigating argv[0] and
   other command line options. A command line application will
   register itself during construction. Typically the Appliction is a
   static object. Therefore the Application has been registered when
   main() is executed. The static main() is a friend of Application
   which makes it possible for main() to search through the linked
   list of Applications. The Application::next member is used as a
   next pointer in the linked list of registered Application s.

   An Application can have a set of associated Option s. When an
   Option is constructed it is put into a static temporary linked
   list. When an Application has constructed all its Options it should
   call registerOptions() to move the linked list of Options into the
   Application s linked lists of Option s.
 */
class Application{
  // Nested Option Classes
  /**
     The general commmand line option.
   */
 protected:
  /**
     This abstract superclass describes the general Option. An Option object is
     intended to be part of an Application object. Before the
     Application::main() is run the command line is parsed by matching
     argv[1], argv[2]... with the Option s of the Application in
     question. The matches() method is called with argv[i] and by the
     return value the Option tells if it matches the argv[i]
     string. If this is the case then the parseValue() method will be
     called to let the Option take values from argv[i+1], argv[i+2]... as
     it likes.

     All this has already been taken care of in the code of the subclasses
     of this class. These subclasses can be used without worrying
     about the implementation details. See the documentation for the
     respective subclasses.
   */
  class Option{
    bool hidden;
  protected:
    /**
       The base pointer for a temporary linked list of Option s. When
       an Option is created it is added to the list by its
       constructor. The content of this list is later moved to a
       linked list for the Application of which the Option is a part.
     */
    static Option *constructionList;
  public:
    Option();
    /**
       An Option can be hidden which means that it does not show up in
       the help text of the Application. This function tells whether
       the Option is hidden or not.
     */
    bool isHidden()const;
    /**
       This method hides the Option. See isHidden().
     */
    void hide(bool b=true);
    /**
       This method tells whether the Option matches a given string or
       not.  Notice that a Option may match more than a single string
       allowing command line parameters such as -n2 and -n7 to match
       the same Option.  @param s The string to be matched.
     */
    virtual bool matches(const char *s)=0;
    /**
       @todo document this method
    */
    virtual void onOptionsParsed();
    /**
       This function returns a string that will be used in the
       documentation to denote the Option.
     */
    virtual std::string documentationMatchString()=0;
    /**
       This function returns the a string describing the Option to be
       used for documentation.
     */
    virtual std::string documentationDescription()=0;
    /**
       When this method is called the Option has the possibility to
       read off its value from the command line through argv[]. This
       only happens if the Option already accepted the argv[t] string
       as a match through matches(). How many additional strings the
       Option reads its value from is written to *numberOfArguments.

       @param t the index of the matching command line string.

       @param argv the command line strings.

       @param ok is considered to be a return value telling whether
       the commandline strings were parsed correctly. The pointer is
       not allowed to be null.

       @param numberOfArguments is considered to be return value
       telling how many additional (besides argv[t]) command line
       strings the method parsed. This allows each option to "eat" a
       different number of command line strings. The pointer is not
       allowed to be null.
    */
    virtual void parseValue(int t, char **argv, bool *ok, int *numberOfArgumentsParsed)=0;
    static Option *getOptionList();
    Option *next;
  };
  class StringMatchingOption: public Option
    {
    protected:
      const char *matchString;
      const char *description;
      bool isExactMatch(const char *s);
    public:
      virtual std::string documentationMatchString();
      virtual std::string documentationDescription();
      StringMatchingOption(const char *s, const char *description_="");
      bool matches(const char *s);
    };

  class SimpleOption: public StringMatchingOption{
    bool value;
  public:
    SimpleOption(const char *s, const char *description);
    void parseValue(int t, char **argv, bool *ok, int *numberOfArgumentsParsed);
    bool getValue();
  };

  class ValueOption: public StringMatchingOption{
  public:
    virtual std::string documentationMatchString();
    virtual void assignValue(const char *s)=0;
    ValueOption(const char *s, const char *description);
    void parseValue(int t, char **argv, bool *ok, int *numberOfArgumentsParsed);
  };

  class StringOption: public ValueOption{
    const char *value;
  public:
    StringOption(const char *s, const char *description, const char *initialValue=0);
    void assignValue(const char *s);
    const char *getValue();
  };

  class IntegerOption: public ValueOption{
    int value;
    bool hasRange;
    int lower,upper;
  public:
    IntegerOption(const char *s, const char *description, int initialValue=0);
    IntegerOption(const char *s, const char *description, int initialValue, int lower, int upper);
    void assignValue(const char *s);
    int getValue();
  };

  class ZeroOneOption: public IntegerOption{
  public:
    ZeroOneOption(const char *s, const char *description, int initialValue=0);
  };

  // Application members
  private:
    static class Application *applicationList;
  static class Application *findApplication(char *name);
 protected:
   static std::list<Application*> getSortedApplicationList();
  /**
     This static procedure makes a symbolic link (on the file system) for registered
     Application s to the Gfan executable in the specified path. This
     procedure is supposed to be called during the installation of Gfan.
     @param name The name of the Gfan executable.
     @param all If false an Application only gets a symbolic link if its includeInDefaultInstallation() returns true.
     @param path The path to the Gfan executable and directory of the symbolic links. Must be terminated by a '/'.
   */
  static void makeSymbolicLinks(const char *name, bool all, const char *path);
  /**
     This procedure produces the list of Gfan
     Application s in the appendix of the Gfan user's manual. The
     contents is written as LaTeX to stdout. The output contains one
     subsection for each Application. An application is not included
     in the list if its includeInDefaultInstallation() returns false.
     @param all Forces all Application s to be documented.
   */
  static void produceLatexDocumentation(bool all);
  /**
     The base pointer for a the linked list of the Option s of the Application.
   */
  class Application::Option* optionList;
 public:
  /**
     Superconstructor. Adds the Application to the static linked list of Application s.
   */
  Application();
  /**
     The next pointer for the linked list of existing Application s.
   */
  class Application *next;
  /**
     This procedure parses the arguments for the static main() and
     assigns values of the appropriate Option s of the Application.
     @param argc The number of arguments on the command line
     (including the name of the command). See K&R:"The C Programming
     Language".  @param argv The arguments. See K&R:"The C Programming
     Language".  @param argumentsToSkip The number of arguments to
     skip (excluding the name of the executable). Usually no options
     are skipped, but if the program is not invoked using a symbolic
     link, the first (index 0) argument is the executable name and the
     second (index 1) argument is the Application name which should be
     skipped when parsing Option s.
   */
  bool parseOptions(int argc, char **argv, int argumentsToSkip);
  /**
     @return true if the Application should appear in the documentation / be installed during default installation of Gfan.
   */
  virtual bool includeInDefaultInstallation();
  /**
     After construction of its Application::Option s an Application should call this procedure to collect the Option s in the optionList.
   */
  void registerOptions();
  /**
     This procedure writes the help text of the Application to stderr and lists the Option s.
   */
  virtual void printHelp();
  /**
     This virtual method contains the code to be executed when the Application is run.
     @return The value to be passed to the shell when the program finishes execution.
  */
  virtual int main()=0;
  virtual void onExit();
  /**
     @return The help text for the documentation. The format is usual ASCII.
   */
  virtual const char *helpText()=0;
  friend int main(int argc, char *argv[]);
  /**
     This function returns the name of the Application. This name is
     used for matching with arg[0] and other options to decide which
     application to run.
  */
  virtual const char *name()=0;
};


#endif
