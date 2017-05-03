/*****************************************************************************/
/*****************************************************************************/

#include <stdio.h>
#include <string.h>

/*****************************************************************************/

void position(),
     checkentry();
int
setup(fptr)
FILE *fptr;
{
  int read;

  position(fptr);
  checkentry(fptr, &read);

  return read;
}

/*****************************************************************************/

void
position(fptr)
FILE *fptr;
{
  int c;

  while ((c=getc(fptr))!=':');

}

/*****************************************************************************/

void
checkentry(fptr, read)
FILE *fptr;
int *read;
{
  int c, found;

  found = 0;
  while (!found) {
    c=getc(fptr);
    if ((c==' ') || (c=='\t')) { found = 0; }
    else { found = 1; }
  }

  if ((c=='\n')||(c==',')) {
    *read = 0;
  }
  else {
    *read = 1;
    fseek(fptr, ftell(fptr)-1, 0);
  }

}

/*****************************************************************************/

/*
int
readentries(variables, fptr)
int *variables;
FILE *fptr;
{
  int i, c, found, true;

  i = 0;
  true = 1;

  while (true) {

    found = 0;
    while (!found) {
      c=getc(fptr);
      if ((c==' ') || (c=='\t')) { found = 0; }
      else { found = 1; }
    }

    if ((c=='\n')||(c==',')) {
      true = 0;
      variables[i] = -1;
    }
    else {
      fseek(fptr, ftell(fptr)-1, 0);
      fscanf(fptr, "%d", &variables[i++]);
    }

  }

}
*/

/*****************************************************************************/

void
myfscanf(string, fptr)
char *string;
FILE *fptr;
{
  int i, c;

  i = 0;
  while ((c=getc(fptr))!=',') { ++i; }
  fseek(fptr, ftell(fptr)-i-1, 0);
  fgets(string,i+2,fptr);
  strncpy(strrchr(string,','),"\0",1);

}

/*****************************************************************************/
/*****************************************************************************/
