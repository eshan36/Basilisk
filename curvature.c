event figure_curveiso (t = 0.1; t <= MAXTIME; t += 0.1)
{
  view (camera = "iso", fov = 3, quat = {0,0,0,1}, tx = -0.015, ty = -0.34, tz= -0.110, bg = {1,1,1}, width = 1080, height = 1920, samples = 4);
  
 
  scalar kappa[], kappa1[];
  curvature(f,kappa);
  boundary({kappa});
  foreach()
    kappa1[] = fabs(kappa[]);
  draw_vof ("f",color = "kappa1", min=0.0, max= 10.0,  map = jet);


  char name[100];
  sprintf(name,"curveiso-%.1f.png",t);

  save (name);
  }
