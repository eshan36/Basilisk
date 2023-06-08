event figure_vorticityyy (t = 0.1; t <= MAXTIME; t += 0.1)
{
  view (camera = "top", fov = 3, quat = {0,0,0,1}, tx = 0.01, ty = -0.01, bg = {1,1,1}, width = 1080, height = 1920, samples = 4);

  clear();
 
  double yb = 0., sb = 0.;
  scalar omegay[];
   foreach(reduction(+:yb) 
	  reduction(+:sb)) {
    double dv = (1. - f[])*dv();
    yb += y*dv;
    sb += dv;
  }

   foreach()
      omegay[] = (u.x[0,0,1] - u.x[0,0,-1] - u.z[1,0,0] + u.z[-1,0,0])/(2.*Delta);
  boundary ({omegay});

  squares("omegay" ,n = {0,1,0} ,alpha= (yb/sb)-2, min=-1, max=1, map = jet, linear = true);
 
  char name[100];
  sprintf(name,"vorticityyy-%.1f.png",t);

  save (name);
}
