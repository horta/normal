#include "normal.h"
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>

double* read_file(char *filename)
{
  FILE   *fp = fopen(filename, "r");
  char    line[27];
  char   *end;
  int     nelems  = 50000;
  double *numbers = malloc(sizeof(double) * nelems);

  for (int i = 0; i < nelems; i++)
  {
    fread(line, sizeof(char), 26, fp);
    nelems[i] = strtod(line, &end);
  }

  fclose(fp);
  return nelems;
}

int main()
{
  double *x           = read_file("x.txt");
  double *logsf_vals  = read_file("logsf.txt");
  double *logcdf_vals = read_file("logcdf.txt");

  return 0;
}
