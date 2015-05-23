#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUFSIZE 256*3

main(int argc, char *argv[]) {
  char line[BUFSIZE];
  FILE *fp;
  char *accession;
  char *sequence;
  char *quality;
  int qualified;
  
  char *libid = argv[1];
  int len = atoi(argv[2]);
  char *stat = argv[3];

  sprintf(line, "gsort -n --parallel=`gnproc` | guniq -c | sed -e 's/^ *//' > %s", stat);
  fp = popen(line, "w");
  while(fgets(line, BUFSIZE, stdin) != NULL) {
    accession = strtok(line, "\t");
    sequence = strtok(NULL, "\t");
    quality = strtok(NULL, "\t");
    qualified = strcspn(quality, "!\"#") ;
    fprintf(fp, "%d\n", qualified);
    if(qualified >= len) {
      quality[len] = 0;
      sequence[len] = 0;
      fprintf(stdout, "%s\t%s:1-%i\t%s\t%s\n", libid, accession, len, quality, sequence);
    }
  }
  pclose(fp);
}
