#include "grid/octree.h"
#include "output_vtu_foreach.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"

#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))

#include "two-phase.h"

#include "tension.h"

#include "lambda2.h"

#include "view.h"
#include "iso3D.h"


#include "maxruntime.h"


#define RHOR 1000.
#define MUR 100.

# define Ga 201.8
# define Bo 16.
# define MAXTIME 100


#define WIDTH 100
#define Z0 4.
int LEVEL = 13;


/*u.t[left] = dirichlet(0.);
u.r[left] = dirichlet(0.);
u.t[right] = dirichlet(0.);
u.r[right] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.r[bottom] = dirichlet(0.);
u.t[back] = dirichlet(0.);
u.r[back] = dirichlet(0.);
u.t[front] = dirichlet(0.);
u.r[front] = dirichlet(0.);*/

int main (int argc, char * argv[]) {
  maxruntime (&argc, argv);
  if (argc > 1)
    LEVEL = atoi (argv[1]);
  
  
  size (WIDTH);
  origin (-0.50*L0, 0, -0.50*L0);
  init_grid (128);

  rho1 = 1.;
  rho2 = 1./RHOR;
  mu1 = 1./Ga;
  mu2 = 1./(MUR*Ga);
  f.sigma = 1./Bo;

  TOLERANCE = 1e-4;
  run();
}


event init (t = 0) {
  if (!restore (file = "restart")) {
    refine (sq(x) + sq(y - Z0) + sq(z) - sq(0.75) < 0 && level < LEVEL);
    fraction (f, sq(x) + sq(y - Z0) + sq(z) - sq(.5));
  }
}


event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] -= 1.;
}


event adapt (i++) {
  double uemax = 1e-2;
  #if 1
  scalar f1[];
  foreach()
    f1[] = f[];
  boundary ({f1});
  adapt_wavelet ({f1,u}, (double[]){0.01,uemax,uemax,uemax}, LEVEL, 5);
  #else
  adapt_wavelet ({f, u}, (double[]){0.01,uemax,uemax,uemax}, LEVEL, 5);
  #endif
}


event logfile (i += 10) {

  static FILE * fp_bubble = fopen("log", "w");
  static FILE * fp_vorticity = fopen("vorticity", "w");


  double xb = 0., yb = 0., zb = 0., sb = 0.;
  double vbx = 0., vby = 0., vbz = 0., pb=0;

  double rmax = -HUGE, rmin = HUGE, oy_max = -HUGE, o_max = -HUGE;
  scalar omegax[], omegay[], omegaz[];
  double sum_yp = 0., sum_yn = 0., sum_o = 0., area = 0.; //sum_yp is the positive Omega_y, sum_yn is the negative Omega_y, sum_o is the total Omega, area is the interface area

  foreach(reduction(+:xb) reduction(+:yb) reduction(+:zb)
	  reduction(+:vbx) reduction(+:vby) reduction(+:vbz)
	  reduction(+:pb) reduction(+:sb)) {
    double dv = (1-f[])*dv();
    xb += x*dv;
    yb += y*dv;
    zb += z*dv;
    vbx += u.x[]*dv;
    vby += u.y[]*dv;
    vbz += u.z[]*dv;
    pb += p[]*dv;
    sb += dv;
    omegax[] = (u.z[0,1,0] - u.z[0,-1,0] - u.y[0,0,1] + u.y[0,0,-1])/(2.*Delta);
    omegay[] = (u.x[0,0,1] - u.x[0,0,-1] - u.z[1,0,0] + u.z[-1,0,0])/(2.*Delta);
    omegaz[] = (u.y[1,0,0] - u.y[-1,0,0] - u.x[0,1,0] + u.x[0,-1,0])/(2.*Delta);
  }

    boundary ({omegax, omegay,omegaz});

   foreach(reduction(max:rmax) reduction(min:rmin) reduction(+:area) reduction(+:sum_yp) reduction(+:sum_yn)  reduction(max:oy_max) reduction(+:sum_o) reduction(max:o_max)) {
    if (f[] > 0 && f[] < 1) {
      coord p;
      coord n = mycs (point, f);
      double alpha = plane_alpha (f[], n);
      double s = plane_area_center (n, alpha, &p);

    double rad  = sqrt(sq(x + Delta*p.x-xb/sb) + sq(y + Delta*p.y-yb/sb) + sq(z + Delta*p.z-zb/sb));
      if (rad > rmax)
	     rmax = rad;
      if (rad < rmin)
	     rmin = rad;

      area += pow(Delta,2)*s;

      double inter_ox = interpolate(omegax,(x+Delta*p.x),(y+Delta*p.y),(z+Delta*p.z));
      double inter_oy = interpolate(omegay,(x+Delta*p.x),(y+Delta*p.y),(z+Delta*p.z));
      double inter_oz = interpolate(omegaz,(x+Delta*p.x),(y+Delta*p.y),(z+Delta*p.z));
      if (inter_ox > 1000. || inter_ox < -1000.)
         inter_ox = omegax[];
      if (inter_oy > 1000. || inter_oy < -1000.)
         inter_oy = omegay[];
      if (inter_oz > 1000. || inter_oz < -1000.)
         inter_oz = omegaz[];

      if (inter_oy > 0)
        sum_yp += pow(Delta,2)*s*inter_oy;
      else if (inter_oy < 0)
        sum_yn += pow(Delta,2)*s*inter_oy;

      if (fabs(inter_oy) > oy_max)
	      oy_max = fabs(inter_oy);

      double inter_o = sqrt(pow(inter_ox,2)+pow(inter_oy,2)+pow(inter_oz,2));
      sum_o += pow(Delta,2)*s*inter_o;

      if (inter_o > o_max)
        o_max = inter_o;
    }
 }

  // 1-time 2-volume 3-x 4-y 5-z 6-vx 7-vy 8-vz 9-dmax/2.0 10-dmin/2.0
  fprintf (fp_bubble,
	   "%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n",
	   t, sb, xb/sb, yb/sb, zb/sb, vbx/sb, vby/sb, vbz/sb, rmax, rmin, pb/sb);
  fflush (fp_bubble);

  // 1-time 2-area 3-omega-y-positive 4-omega-y-negative 5-omegay-max 6-omega 7-omega-max
  fprintf (fp_vorticity,
	   "%.8f %.8f %.8f %.8f %.8f %.8f %.8f\n",
	   t, area, sum_yp, sum_yn, oy_max, sum_o, o_max);
  fflush (fp_vorticity);

}


  
event snapshot (t = 0; t <= MAXTIME; t += 1)
{
  scalar l2[], omegay[], kappa[];
  lambda2 (u, l2);
  foreach()
    omegay[] = (u.z[1] - u.z[-1] - u.x[0,0,1] + u.x[0,0,-1])/(2.*Delta);
  boundary ({omegay});
  curvature(f,kappa);
  char name[80];
  sprintf (name, "dump-%03d", (int) t);
  dump (file = name);

}

void save_data (scalar f, scalar p, scalar omegay, vector u, int t) {
    
  char name[80], subname[80];
  
  FILE * fp ;
  
  t > 0 ? sprintf(name, "profiles_%3.3d_n%3.3d.vtu", t,pid()) : sprintf(name, "profiles_n%3.3d.vtu",pid());
    fp = fopen(name, "w");
    output_vtu_ascii_foreach ((scalar *) {f, p, omegay}, (vector *) {u}, N, fp, false);
    fclose (fp);


    if (pid()==0){
      t > 0 ? sprintf(name, "profiles_%3.3d.pvtu", t) : sprintf(name, "profiles.pvtu");
      t > 0 ? sprintf(subname, "profiles_%3.3d",  t) : sprintf(subname, "profiles");
        fp = fopen(name, "w");
        output_pvtu_ascii ((scalar *) {f, p, omegay}, (vector *) {u}, N, fp, subname);
        fclose (fp);

    }
 
MPI_Barrier(MPI_COMM_WORLD);
  }

event loggfilek (t = 1; t <= MAXTIME; t += 1){
  
  scalar kappa[],  omegay[];

  foreach()
    omegay[] = (u.x[0,0,1] - u.x[0,0,-1] - u.z[1,0,0] + u.z[-1,0,0])/(2.*Delta);
  
 

  boundary((scalar *){omegay});
  
  save_data(f, p, omegay, u, (int) t);
   
  }


event figure_lambda2 (t = 0.1; t <= MAXTIME; t += 0.1)
{
  view (camera = "front", fov = 6, quat = {0,0,0,1}, tx = 0.025, ty = -0.35, bg = {1,1,1}, width = 1080, height = 1920, samples = 4);

  
  
  clear();
 
  draw_vof ("f", fc = {0.410156,0.410156,0.410156});

  scalar l2[];
  lambda2 (u, l2);
  isosurface ("l2", -0.0002);

  char name[100];
  sprintf(name,"lambda2-%.1f.png",t);

  save (name);
}

event figure_lambdaa2 (t = 0.1; t <= MAXTIME; t += 0.1)
{
  view (camera = "right", fov = 6, quat = {0,0,0,1}, tx = -0.0275, ty = -0.10, bg = {1,1,1}, width = 1080, height = 1920, samples = 4);

  travelling (23, 82, ty = -0.6);
    
  clear();
 
  draw_vof ("f", fc = {0.410156,0.410156,0.410156});

  scalar l2[];
  lambda2 (u, l2);
  isosurface ("l2", -0.0002);

  char name[100];
  sprintf(name,"lambdaa2-%.1f.png",t);

  save (name);
}

event figure_omega (t = 0.1; t <= MAXTIME; t += 0.1)
{
  view (camera = "front", fov = 6, quat = {0,0,0,1}, tx = 0.0275, ty = -0.35, bg = {1,1,1}, width = 1080, height = 1920, samples = 4);

  

  clear();
 
  draw_vof ("f", fc = {0.410156,0.410156,0.410156});

  scalar omegay[];
  foreach()
    omegay[] = -(u.z[1] - u.z[-1] - u.x[0,0,1] + u.x[0,0,-1])/(2.*Delta);
  boundary ({omegay});
  isosurface ("omegay", 0.15, fc = {1,0,0.0});
  isosurface ("omegay", -0.15, fc = {0,1.0,0.0});

  char name[100];
  sprintf(name,"omegay-%.1f.png",t);

  save (name);
}

event figure_omegaa (t = 0.1; t <= MAXTIME; t += 0.1)
{
  view (camera = "right", fov = 6, quat = {0,0,0,1}, tx = -0.025, ty = -0.35, bg = {1,1,1}, width = 1080, height = 1920, samples = 4);
  
  clear();
 
  draw_vof ("f", fc = {0.410156,0.410156,0.410156});

  scalar omegay[];
  foreach()
    omegay[] = -(u.z[1] - u.z[-1] - u.x[0,0,1] + u.x[0,0,-1])/(2.*Delta);
  boundary ({omegay});
  isosurface ("omegay", 0.15, fc = {1,0,0.0});
  isosurface ("omegay", -0.15, fc = {0,1.0,0.0});

  char name[100];
  sprintf(name,"omegaay-%.1f.png",t);


  save (name);
}

event figure_omegaaa (t = 0.1; t <= MAXTIME; t += 0.1)
{
  view (camera ="iso", fov = 6, quat = {0,0,0,1}, tx = -0.025, ty = -0.30, bg = {1,1,1}, width = 1080, height = 1920, samples = 4);

  
  clear();
 
  draw_vof ("f", fc = {0.410156,0.410156,0.410156});

  scalar omegay[];
  foreach()
    omegay[] = -(u.z[1] - u.z[-1] - u.x[0,0,1] + u.x[0,0,-1])/(2.*Delta);
  boundary ({omegay});
  isosurface ("omegay", 0.15, fc = {1,0,0.0});
  isosurface ("omegay", -0.15, fc = {0,1.0,0.0});

  char name[100];
  sprintf(name,"omegaaay-%.1f.png",t);


  save (name);
}
/*
event figure_bubble (t = 0.1; t <= MAXTIME; t += 0.0001)
{
  view (camera = "front", fov = 2, quat = {0,0,0,1}, tx = 0.0, ty = -0.38, bg = {1,1,1}, width = 1080, height = 1920, samples = 4);

 
  draw_vof ("f", fc = {0.410156,0.410156,0.410156});


  char name[100];
  sprintf(name,"bubble-%.1f.png",t);

  save (name);
  }
event figure_bubblee (t = 0.1; t <= MAXTIME; t += 0.0001)
{
  view (camera = "right", fov = 2, quat = {0,0,0,1}, tx = 0.015, ty = -0.38, bg = {1,1,1}, width = 1080, height = 1920, samples = 4);

 
  draw_vof ("f", fc = {0.410156,0.410156,0.410156});


  char name[100];
  sprintf(name,"bubblee-%.1f.png",t);

  save (name);
  }
event figure_bubbleee (t = 0.1; t <= MAXTIME; t += 0.1)
{
  view (camera = "top", fov = 3, quat = {0,0,0,1}, tx = 0.0, ty = -0.03, bg = {1,1,1}, width = 1080, height = 1920, samples = 4);

 
  draw_vof ("f", fc = {0.410156,0.410156,0.410156});


  char name[100];
  sprintf(name,"bubbleee-%.1f.png",t);

  save (name);
  }


event figure_vorticityyy (t = 0.1; t <= MAXTIME; t += 0.1)
{
  view (camera = "top", fov = 3, quat = {0,0,0,1}, tx = 0.01, ty = -0.01, bg = {1,1,1}, width = 1080, height = 1920, samples = 4);

  clear();
 
  // draw_vof ("u.x", fc = {0.410156,0.410156,0.410156});
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
event figure_vbubblebot (t = 0.1; t <= MAXTIME; t += 0.1)
{
  view (camera = "bottom", fov = 3, quat = {0,0,0,1}, tx = 0.0, ty = -0.03, bg = {1,1,1}, width = 1080, height = 1920, samples = 4);

 
   scalar omegay[];
   foreach()
      omegay[] = -(u.z[1] - u.z[-1] - u.x[0,0,1] + u.x[0,0,-1])/(2.*Delta);
  boundary ({omegay});
  
  draw_vof ("f",color = "omegay", min=-5.0, max=5.0,  map = jet, linear = true);


  char name[100];
  sprintf(name,"vbubblebot-%.1f.png",t);

  save (name);
  }
event figure_vbubbletop (t = 0.1; t <= MAXTIME; t += 0.1)
{
  view (camera = "top", fov = 3, quat = {0,0,0,1}, tx = 0.0, ty = -0.03, bg = {1,1,1}, width = 1080, height = 1920, samples = 4);

 
   scalar omegay[];
   foreach()
      omegay[] = -(u.z[1] - u.z[-1] - u.x[0,0,1] + u.x[0,0,-1])/(2.*Delta);
  boundary ({omegay});
  
  draw_vof ("f",color = "omegay", min=-5.0, max=5.0,  map = jet, linear = true);


  char name[100];
  sprintf(name,"vbubbletop-%.1f.png",t);

  save (name);
  }

event figure_vbubble (t = 0.1; t <= MAXTIME; t += 0.1)
{
  view (camera = "iso", fov = 3, quat = {0,0,0,1}, tx = -0.015, ty = -0.18, tz= -0.110, bg = {1,1,1}, width = 1080, height = 1920, samples = 4);

 
   scalar omegay[];
   foreach()
      omegay[] = -(u.z[1] - u.z[-1] - u.x[0,0,1] + u.x[0,0,-1])/(2.*Delta);
  boundary ({omegay});
  
  draw_vof ("f",color = "omegay", min=-5.0, max=5.0,  map = jet, linear = true);


  char name[100];
  sprintf(name,"vbubble-%.1f.png",t);

  save (name);
  }
event figure_curvebot (t = 0.1; t <= MAXTIME; t += 0.1)
{
   view (camera = "bottom", fov = 3, quat = {0,0,0,1}, tx = 0.0, ty = -0.03, bg = {1,1,1}, width = 1080, height = 1920, samples = 4);

 
  scalar kappa[], kappa1[];
  curvature(f,kappa);
  boundary({kappa});
  foreach()
    kappa1[] = fabs(kappa[]);
  draw_vof ("f",color = "kappa1", min=0.0, max= 10.0,  map = jet);

  char name[100];
  sprintf(name,"curvebot-%.1f.png",t);

  save (name);
  }
event figure_curvetop (t = 0.1; t <= MAXTIME; t += 0.1)
{
   view (camera = "top", fov = 3, quat = {0,0,0,1}, tx = 0.0, ty = -0.03, bg = {1,1,1}, width = 1080, height = 1920, samples = 4);

 
  scalar kappa[], kappa1[];
  curvature(f,kappa);
  boundary({kappa});
  foreach()
    kappa1[] = fabs(kappa[]);
  draw_vof ("f",color = "kappa1", min=0.0, max= 10.0,  map = jet);


  char name[100];
  sprintf(name,"curvetop-%.1f.png",t);

  save (name);
  }


event figure_curveiso (t = 0.1; t <= MAXTIME; t += 0.1)
{
  view (camera = "iso", fov = 3, quat = {0,0,0,1}, tx = -0.015, ty = -0.18, tz= -0.110, bg = {1,1,1}, width = 1080, height = 1920, samples = 4);

 
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
