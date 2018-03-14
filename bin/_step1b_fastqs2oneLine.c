#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUFSIZE 512

int main(int argc, char *argv[]) {
  char header[BUFSIZE];
  char sequence[BUFSIZE];
  char dummy[BUFSIZE];
  char quality[BUFSIZE];
  char *acc;
  char *tmp;
  char sequence_umi[BUFSIZE];
  char *sequence_gapcdna;
  char *sequence_barcode;
  char quality_umi[BUFSIZE];
  char *quality_gapcdna;
  char *quality_barcode;
  int umi;

  umi = atoi(argv[1]);
  memset(sequence_umi, '\0', sizeof(sequence_umi));
  memset(quality_umi, '\0', sizeof(quality_umi));
  
  while(fgets(header, BUFSIZE, stdin) != NULL) {
    if(fgets(sequence, BUFSIZE, stdin) == NULL
       || fgets(dummy, BUFSIZE, stdin) == NULL
       || fgets(quality, BUFSIZE, stdin) == NULL) {
      exit(1);
    }
    header[strcspn(header, "\r\n")] = 0;
    acc = strtok(header+1, " \t\r\n");
    sequence[strcspn(sequence, "\r\n")] = 0;
    tmp = strtok(sequence, "\t");
    strncpy(sequence_umi, tmp, umi);
    sequence_gapcdna = tmp+umi;
    sequence_barcode = strtok(NULL, "\t");
    quality[strcspn(quality, "\r\n")] = 0;
    tmp = strtok(quality, "\t");
    strncpy(quality_umi, tmp, umi);
    quality_gapcdna = tmp+umi;
    quality_barcode = strtok(NULL, "\t");
    fprintf(stdout, "%s\t%s%s%s\t%s%s%s\n", acc, sequence_umi, sequence_barcode, sequence_gapcdna, quality_umi, quality_barcode, quality_gapcdna);
  }
}
