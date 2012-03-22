/* Procedures to return parameters from the input file
 * Integer parameters should be accessed with get_iparam
 * double-precision parameters should be accessed with get_dparam
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int get_iparam(char paramname[], char filename[])
{
  FILE *inputfile;
  char inputstring[256];
  char searchstring[256];
  char scanstring[256];
  int returnvalue = 0;
  
  strcpy(searchstring, paramname);
  strcat(searchstring, "=");

  strcpy(scanstring, paramname);
  strcat(scanstring, "=%i");

  inputfile = fopen(filename, "r");
  
  if (inputfile == NULL) {
    fprintf(stderr,"Input file failed to open.\n");
  }

  while(fgets(inputstring, 256, inputfile) != NULL) {
    if((strstr(inputstring, searchstring) != NULL) &&
       (strstr(inputstring, searchstring) - inputstring) == 0) {
      sscanf(inputstring, scanstring, &returnvalue);
      break;
    }
  }

  fclose(inputfile);
  return returnvalue;
}

double get_dparam(char paramname[], char filename[])
{
  FILE *inputfile;
  char inputstring[256];
  char searchstring[256];
  char scanstring[256];
  double returnvalue;

  strcpy(searchstring, paramname);
  strcat(searchstring, "=");

  strcpy(scanstring, paramname);
  strcat(scanstring, "=%le");

  inputfile = fopen(filename, "r");
  
  if (inputfile == NULL) {
    fprintf(stderr,"Input file failed to open.\n");
  }
  
  while(fgets(inputstring, 256, inputfile) != NULL) {
    if((strstr(inputstring, searchstring) != NULL) &&
       (strstr(inputstring, searchstring) - inputstring) == 0) {
      sscanf(inputstring, scanstring, &returnvalue);
      break;
    }
  }

  fclose(inputfile);
  return returnvalue;
}

void get_sparam(char paramname[], char filename[], char destination[])
{
  FILE *inputfile;
  char inputstring[256];
  char searchstring[256];
  char scanstring[256];
  char returnstring[256];

  strcpy(searchstring, paramname);
  strcat(searchstring, "=");

  strcpy(scanstring, paramname);
  strcat(scanstring, "=%s");

  inputfile = fopen(filename, "r");
  
  if (inputfile == NULL) {
    fprintf(stderr,"Input file failed to open.\n");
  }
  
  while(fgets(inputstring, 256, inputfile) != NULL) {
    if((strstr(inputstring, searchstring) != NULL) &&
       (strstr(inputstring, searchstring) - inputstring) == 0) {
      sscanf(inputstring, scanstring, returnstring);
      break;
    }
  }

  fclose(inputfile);
  strcpy(destination, returnstring);
}
