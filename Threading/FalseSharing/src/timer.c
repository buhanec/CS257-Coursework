#include <sys/time.h>
#include <time.h>
#include "example.h"

/* Return gettimeofday value in seconds */

double timer()
{
	struct timeval t;
	gettimeofday(&t,NULL);
	double microsec = t.tv_usec * 1E-6;
	double sec = (double) t.tv_sec;
	return sec+microsec;
	
}
