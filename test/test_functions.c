#include "normal.h"
#include <stdio.h>
#include <stdlib.h>

double* read_file(char *filename)
{
  FILE *fp = fopen(filename, "r");

  char line[40];
  char   *end;
  int     nelems  = 50000;
  double *numbers = malloc(sizeof(double) * nelems);


  int i = 0;

  while (fgets(line, sizeof(line), fp)) {
    numbers[i] = strtod(line, &end);
    i++;
  }

  fclose(fp);

  return numbers;
}

int main()
{
  double *x           = read_file("test/x.txt");
  double *logcdf_vals = read_file("test/logcdf.txt");

  int nelems = 50000;

  for (int i = 0; i < nelems; i++)
  {
    if (fabs(logcdf(x[i]) - logcdf_vals[i]) > 1e-13) return 1;
  }


  free(x);
  free(logcdf_vals);

  return 0;
}
