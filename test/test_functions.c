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
  double *logsf_vals  = read_file("test/logsf.txt");
  double *logcdf_vals = read_file("test/logcdf.txt");

  // printf("%.70f\n",   x[0]);
  // printf("%.70f\n\n", x[50000 - 1]);
  //
  // printf("%.70f\n",   logsf_vals[0]);
  // printf("%.70f\n\n", logsf_vals[50000 - 1]);
  //
  // printf("%.70f\n",   logcdf_vals[0]);
  // printf("%.70f\n\n", logcdf_vals[50000 - 1]);

  int nelems = 50000;

  for (int i = 0; i < nelems; i++)
  {
    if (fabs(logcdf(x[i]) - logcdf_vals[i]) > 1e-13) return 1;

    // {
    //   printf("x: %.30f\n",              x[i]);
    //   printf("logcdf(x[i]): %.30f\n",   logcdf(x[i]));
    //   printf("logcdf_vals[i]: %.30f\n", logcdf_vals[i]);
    //   return 1;
    // }
  }


  free(x);
  free(logsf_vals);
  free(logcdf_vals);


  return 0;
}
