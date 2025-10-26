/* Generate brown trout data with detection heterogeneity with psi=.6, p=2. or .6 and 4 sample occasions */
/* About 1/3 of the sites had low detectability */

data history;;
   file 'browntrout.inp' noprint;
   length history $4;
   do site=1 to 500;
      history = '0000';
	  occupy = ranuni(23432)<0.5; /* is the site occupied */
	  group = ranuni(24234) < .333; /* is this group 0 or group 1 */
	  p = rand('BETA',2,2);
	  if occupy = 0 then p = 0;
	  do i=1 to 4;
	     if ranuni(23423) < p then substr(history,i,1)='1';
      end;
	  put history " 1;";
	  output history;
   end;
run;

proc tabulate data=history;
  class history occupy group;
  table history ALL, occupy*group*n*f=5.0;
run; 

