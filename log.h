#ifndef LOG_H_INCLUDED
#define LOG_H_INCLUDED

/* The main purpose of this file is to define notation for identifying the
   parts of the Gfan code that prints out debug information.
   Later this notation will probably change.
   A second purpose is to be able to disable logging.*/

extern int logLevel;
void setLogLevel(int l);

/* Use log0 for temporary output while debugging. This makes it easy
   to remove debugging commands after the bug has been fixed. */
#define log0 if(logLevel>=0)

#define log1 if(logLevel>=1)
#define log2 if(logLevel>=2)
#define log3 if(logLevel>=3)

#endif
