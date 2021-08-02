#include <stdio.h>
#include <stdlib.h>


int main ()
{
	int year,x,y;

	int pp[10][7];


	int day,yy,ly,month;
	ly=0;

	FILE *fp;
	char ds[100];

	month=9;
	day=12;

	int dayn,starty,count;

        for (starty=1000;starty<3000;starty++)
        {

		for (x=0;x<10;x++) { for (y=0;y<7;y++) { pp[x][y]=0; } }

		count=0;
		year=starty;
		while (count < 70)
		{
  			char path[1035];
			sprintf(ds,"date -d \"%02d/%02d %04d\" +%%w",month,day,year);
		  	/* Open the command for reading. */
  			fp = popen(ds, "r");
  			if (fp == NULL) {
    				printf("Failed to run command\n" );
    			exit(1); }

  			/* Read the output a line at a time - output it. */
  			while (fgets(path, sizeof(path), fp) != NULL) {
    				//printf("%s", path); 
			}
			int ret;
			ret=WEXITSTATUS(pclose(fp));
			if (ret !=0 )
			{
				printf("BAD DATE\n");
				exit (1);
			}
			dayn=atoi(path);
			//printf("year %d day %d\n",year,dayn);
			pp[year%10][dayn]=1;

			count=0;
		        for (x=0;x<10;x++)
                	{
                        	for (y=0;y<7;y++)
                        	{
                                	if (pp[x][y]==1){count ++;}
                        	}
                	}
			year++;
		}
               	printf ("%02d/%02d/%04d Age %d \n",day,month,starty,year-starty-1);
  	}

  /* close */

}
