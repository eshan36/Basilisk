#include "iso3D.h"
event figure_velocityyy (t = 0.1; t <= MAXTIME; t += 0.1)
{
  view (camera = "top", fov = 3, quat = {0,0,0,1}, tx = 0.0, ty = -0.03, bg = {1,1,1}, width = 1080, height = 1920, samples = 4);
\\setting the position of camera
  clear();

  scalar omegay[],mag_u[];
  double yb = 0., sb = 0.;
	 
  \\position of center of bubble in y axis
   foreach(reduction(+:yb) 
	  reduction(+:sb)) {
    double dv = (1. - f[])*dv();
    yb += y*dv;
    sb += dv;
  }
 
  foreach()
    mag_u[] = sqrt(sq(u.x[])+sq(u.y[])+sq(u.z[]));
  
  boundary ({mag_u});

  squares("mag_u" ,n = {0,1,0} ,alpha= yb/sb, min = 0, max =2.8, map = jet, linear = true);
  isoline2 ("f", val = 0.5,
	      np = {0,1,0}, alpha = yb/sb, lc = {1,1,0}, lw = 3);
  char name[100];
  sprintf(name,"velocityyy-%.1f.png",t);

  save (name);
  }
