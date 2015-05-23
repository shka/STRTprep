#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUFSIZE 256

main() {
  char header[BUFSIZE];
  char sequence[BUFSIZE];
  char dummy[BUFSIZE];
  char quality[BUFSIZE];
  char *acc;
    
  while(fgets(header, BUFSIZE, stdin) != NULL) {
    if(fgets(sequence, BUFSIZE, stdin) == NULL
       || fgets(dummy, BUFSIZE, stdin) == NULL
       || fgets(quality, BUFSIZE, stdin) == NULL) {
      exit(1);
    }
    acc = strtok(header+1, " \t");
    sequence[strcspn(sequence, "\r\n")] = 0;
    quality[strcspn(quality, "\r\n")] = 0;
    fprintf(stdout, "%s\t%s\t%s\n", acc, sequence, quality);
  }
}
