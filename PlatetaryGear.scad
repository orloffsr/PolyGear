use <PolyGear.scad>
use <shortcuts.scad>
use <star.scad>

//module planetary_gear(modul, sun_teeth, planet_teeth, number_planets, width, rim_width, bore, pressure_angle=20, helix_angle=0, together_built=true, optimized=true)

module planetary_gear(
//basic options
  n_planets = 3,
  s = 11, //sun teeth
  p = 11, //planet teeth
  z = 1, //thickness
  m = 1, //module
  pressure_angle = 20,
  helix_angle    = 0,   // the sign gives the handiness, can be a list
  backlash       = 0.1, // in module units
//advanced options
  chamfer       = 0, // degrees, should be in [0:90[
  chamfer_shift = 0, // from the pitch radius in module units
  add = 0, // add to addendum
  ded = 0, // subtract to the dedendum
  x   = 0, // profile shift
  together_built=true, //if false, will build exploded model suitable for making individual parts
  od_padding = 1, //amount to add to radius of ring gear
//finesse options
  $fn=5,     // tooth profile subdivisions
){
  
  
  
  //dimensions
  sun_d = m*s;
  planet_d = m*p;
  r = s+2*p; //ring teeth
  ring_d = m*r;
  carrier_radius = together_built ? (sun_d + planet_d)/2 
                                  : ring_d + planet_d/2 + od_padding + 2*m;
  
  
  //sun gear
  spur_gear(n=s, m=m, z=z, 
    pressure_angle=pressure_angle,
    helix_angle=helix_angle,   
    backlash=backlash,
    chamfer = chamfer,
    chamfer_shift = chamfer_shift,
    add=add,
    ded=ded,
    x=x,
    $fn=$fn
    );
    
  
  //planets
  for (n=[0:1:n_planets-1]){
    thetaN = n*360/n_planets; //rotation on carrier axis
    theta = thetaN*s/p ; //rotation on gear axis
    theta_tooth = 360/p/2; //rotate 1/2 tooth
    
      Rz(thetaN)
      Tx(carrier_radius)
      Rz(theta + theta_tooth)//teeth  collide, need extra 1/2-tooth rotation to correct
      Mx()// this switches handedness of helix (I think).
        spur_gear(
          n=p, m=m, z=z, 
          pressure_angle=pressure_angle,
          helix_angle=helix_angle, 
          backlash=backlash,
          chamfer = chamfer,
          chamfer_shift=chamfer_shift,
          add=add,
          ded=ded,
          x=x,
          $fn=$fn
          );
  }//end of planets
    
  //ring
  Rz(r%2<0.0001 ? 0 : 360/r/2)  //rotate 1/2 tooth if r is odd
  //Mx()//switches handedness
  difference(){ //switches handedness
    //Tz(0.005)
      Cy(d= ring_d + 2*(m + od_padding), h=z, $fn=50);
    
    spur_gear(
      n=r, m=m, z=z+1,
      pressure_angle=pressure_angle,
      helix_angle=-helix_angle,
      backlash=-backlash,//negative backlash 
      chamfer = chamfer,
      chamfer_shift=chamfer_shift,
      add= add +0.1,
      ded= ded -0.1,
      x=x,
      $fn=$fn
      );
    }//end of difference
    
}


//planetary_gear(s=23, p=16, n_planets=6, z=10,
//              helix_angle=herringbone(helix=30, $fn=10), 
//              together_built=true,
//              $fn=10);


    
module compound_planetary_gear(
//basic options
  n_planets = 3,
  s = [], //sun teeth in each gearset
  p = [], //planet teeth in each gearset
  z = [], //thickness of each gearset
  r = undef, //optional. if not specified, r = s+2*p.
  //m = 1, //module of first gearset -- other modules will be calculated
  od = 100,//if defined, it m will be overridden based on meeting od and padding
  pressure_angle = 20,
  helix_angle    = 0,   // the sign gives the handiness, can be a list
  backlash       = 0.1, // in module units
  gap = 3, //gap between sections. Should be > (m1*p1 - m2*p2)/2 for printability
  connect_suns = true,//will connect suns 1 and 3 (note s[1,3] must equal p[1,3].
  split_ring = true, //will put a split in the ring for assembly
//advanced options
  chamfer       = 0, // degrees, should be in [0:90[
  chamfer_shift = 0, // from the pitch radius in module units
  add = 0, // add to addendum
  ded = 0, // subtract to the dedendum
  x   = 0, // profile shift
  together_built=true, //if false, will build exploded model suitable for making individual parts
  od_padding = 3, //amount to add to radius of ring gear
  od_bumps = 16,//number of protrusions on the OD
  od_bumps_w = 3,
  od_bumps_h = 3,// width and height of the protrusions
//finesse options
  $fn=5,     // tooth profile subdivisions
){
  assert(len(s)==len(p) && len(s)==len(z));
  
  rings_ideal = [ for (i = [0 : 1 : len(s)-1]) (s[i]+2*p[i]) ];
  rings = is_undef(r) ? rings_ideal : r;
  //rings = [ for (i = [0 : 1 : len(s)-1]) (s[i]+2*p[i]) ];
  
  d_list = [ for (i= [0 : 1 : len(s)-1]) rings_ideal[i]/(s[i]+p[i]) ]; //m[i]*rings[i]
  m_idx = index_max(d_list);
  m = (od-od_padding*2)/(rings_ideal[m_idx]+2);
  m_vec = [ for (i = [0 : 1 : len(s)-1]) m*(s[m_idx]+p[m_idx])/(s[i]+p[i]) ];
  //echo(m_idx);
  //echo(m);
  //echo(d_list);
  //echo(m_vec);
  //echo(rings_ideal);
  
  //r1 = s[0] + 2*p[0];
  //m = is_undef(od) ? m
  //    : od/(r1+1+od_padding); //m*(r1+1)+od_padding = od 
  //od =  m*(r1+1)+od_padding;
  extra_radius = together_built ? 0 : (od/2+2*m);
  carrier_radius = m*(s[m_idx]+p[m_idx])/2 + extra_radius;

  

hex_shaft_hole(size=0.25*25.4 + .1) //1/4" hex, in mm plus clearance
etch_star(n_planets, r1 = od+carrier_radius, r2=m, h=m, z=add(z)+gap*(len(z)-1)-m/2)
{
  //planets
  {
    etch_star(n_planets, r1 = od+carrier_radius, r2=m, h=m/2)
    for (n=[0:1:n_planets-1]){
      thetaN = n*360/n_planets;
      
      Rz(thetaN)
      Tx(carrier_radius)
      label(text=str(n+1), size=10*m, depth=m/2, theta=-thetaN)
      for (i = [0:1:len(p)-1]){
        pi = p[i];
        si = s[i];
        zi = z_shift(z, i, gap);
        mi = m_vec[i];
        theta = thetaN*si/pi;
        theta_tooth = 360/pi/2; //rotate 1/2 tooth
        bevel = mi+ded+backlash;
  
        Tz(zi)      
        Rz(theta + theta_tooth)//teeth  collide, need extra 1/2-tooth rotation to correct
        Mx()// this switches handedness of helix (I think).
        union(){
        cylindrical_bevel(r=mi*(pi)/2 , h=z[i], bevel=bevel)
          spur_gear(
            n=pi, m=mi, z=z[i], 
            pressure_angle=pressure_angle,
            helix_angle=helix_angle, 
            backlash=backlash,
            chamfer = chamfer,
            chamfer_shift=chamfer_shift,
            add=add,
            ded=ded,
            x=x
            );
          if(i>0){ //add connecting bits in gaps
              mim1 = m*(s[m_idx] + p[m_idx])/(s[i-1]+p[i-1]);
              bevelim1 = mim1+ded+backlash;
              r1 = mi*(pi)/2-bevel;
              r2 = mim1*(p[i-1])/2 - bevelim1;
              Tz(z[i]/2+gap/2+0.001) Cy(r1=r1, r2=r2, h=gap+0.002, $fn=0);
            }
        }//end union
        
        //label
        
      }//end this planet
      
    }}//end of planets
  
  ///////////////////sun(s)
  
  if(connect_suns){
    assert(len(s)>=3 && len(p) >= 3);
    assert(s[0]==s[2] && p[0] == p[2]);
    
    r0 = m_vec[0]*(s[0]-2)/2; //r and m same for sun 1 and 3
    r1 = m_vec[1]*(s[1] - 2)/2 - add - backlash - 2.1*m_vec[1]; //min radius, with clearance of 2m subtracted 
    
    //etch_star(n_planets, r1 = od+carrier_radius, r2=m, h=m/2)
    union(){
    Tz(z_shift(z, 0, gap))
    cylindrical_bevel(r=r0+m_vec[0], h=z[0], bevel=m_vec[0]+ded+backlash)
    spur_gear(
            n=s[0], m=m_vec[0], z=z[0], 
            pressure_angle=pressure_angle,
            helix_angle=helix_angle,  backlash=backlash,
            chamfer = chamfer, chamfer_shift=chamfer_shift,
            add=add, ded=ded, x=x
            );
    
    //connecting shaft
    Tz(z_shift(z, 1, gap))
    union(){
      r_shaft = min(r0,r1);
      Cy(r=r_shaft , h=z[1], $fn=0);
      Tz((z[1]+gap)/2 + 0.001)
        Cy(r1=r_shaft, r2=r0, h=gap+0.002 ,$fn=0);
      Tz(-(z[1]+gap)/2+ 0.001)
        Cy(r2=r_shaft, r1=r0, h=gap+ 0.002, $fn=0);
      }
      
    Tz(z_shift(z, 2, gap))
    cylindrical_bevel(r=r0+m_vec[2], h=z[2], bevel=m_vec[2]+ded+backlash)
    spur_gear(
            n=s[2], m=m_vec[2], z=z[2], 
            pressure_angle=pressure_angle,
            helix_angle=helix_angle,  backlash=backlash,
            chamfer = chamfer, chamfer_shift=chamfer_shift,
            add=add, ded=ded, x=x
            );
      }//end union
  }//end connected suns
  else{
    for (i = [0:1:len(s)-1]){
      r0 = m_vec[i]*s[i]/2;
      last = i==len(s)-1;
      //etch_depth = last ? m/2 : 0;
      x_offset = last? 0 : max(od,carrier_radius);
      z_offset = together_built ? z_shift(z, i, gap) : z[i]/2; 
      
      Tx(together_built ? 0 : (len(s)-1-i)* max(m_vec)*max(s) + x_offset )
      //label(text=str(len(s)- i), size=10*m, depth=m/2)
      //etch_star(n_planets, r1 = od+carrier_radius, r2=m, h=etch_depth)
      Tz( z_offset )
      cylindrical_bevel(r=r0, h=z[i], bevel=m_vec[i]+ded+backlash)
      spur_gear(
              n=s[i], m=m_vec[i], z=z[i], 
              pressure_angle=pressure_angle,
              helix_angle=helix_angle,  backlash=backlash,
              chamfer = chamfer, chamfer_shift=chamfer_shift,
              add=add, ded=ded, x=x
              );
    
    }//end loop over suns
  }
//end else //ends of suns

//rings

  for (i = [0:1:len(rings)-1]){
    half_tooth = p[i]%2<0.0001 ? 360/rings[i]/2 :0 ; //rotate 1/2 tooth if p is even
    last = i==len(rings)-1;
    y_offset = last ? 0 : max(od, carrier_radius) ;
    z_offset = together_built ?  z_shift(z, i, gap) : z[i]/2;
    m_scale = rings_ideal[i]/rings[i];
    
    Tz(z_offset)
    Ty(together_built ? 0 : (len(rings)-1-i)* (od+2*m+6) + y_offset )
    split_ring(gap=together_built? 0:0.2)
    Rz(half_tooth)
    difference(){ //switches handedness
      bumpyCy(d=od, h=z[i], n_bumps=od_bumps, h_bumps=od_bumps_h, w_bumps=od_bumps_w, $fn=0, $fa=9);//, h_bumps=3, w_bumps=3, $fn=0, $fa=6)//bumpyCy(d= od, h=z[i], $fn=0, $fa=6);
      spur_gear(
        n=rings[i], 
        m=m_vec[i] * m_scale,
        z=z[i]+0.1*m_vec[i],
        pressure_angle=pressure_angle,
        helix_angle=-helix_angle,
        backlash=-backlash,//negative backlash 
        chamfer = chamfer,
        chamfer_shift=chamfer_shift,
        add= add +0.15*m_vec[i], //add clearance after taking difference
        ded= ded -0.15*m_vec[i], //add clearance after taking difference
        x=x
        );
      }//end of difference
  }//end rings
  

}// end of etch



}//end diff planet

//find index of max value in list:
function index_max(l) = search(max(l), l)[0];

//add elements in vector, per openscad manual:
function add(v, i = 0, r = 0) = i < len(v) ? add(v, i + 1, r + v[i]) : r;

function z_shift(z, i, gap) = add(z, i+1)+ (len(z)-i-1)*gap + z[i]/2;

module cylindrical_bevel(r=undef, h=undef, bevel=undef){
        rmin = r- bevel;
        intersection(){
        Tz(-h/4)union()
          {
          Cy(h = h/2, r2=rmin+h/2, r1=rmin, $fn=0);
          Tz(h/2)Mz()
            Cy(h = h/2, r2=rmin+h/2, r1=rmin, $fn=0);
          }
        children();  
        }//end intersection
  }//end bevel_cylinder

module label(text="", size=1, depth=1,theta=0){
  difference(){
    children();
    //Tz(depth)
    Rz(theta)
    //Mz()
    Mx()
    linear_extrude(height = depth)
      text(text, size=size, halign="center", valign="center");   
    }
}//end label

module etch_star(points=5, r1=0, r2=1, h=1, z=0) {
  /*if(h!=0){
    difference(){
      children();
      Tz(z)
      Rz(360/points/2) Tz(-h/10) star(points, r1, r2, h=h*1.1);
    }
  }
  else {children();}
  */
  children(); //need to fix
  
}

module hex_shaft_hole(size=1){ //size is distance flat-to-flat
  difference()
  {
    children();
    Cy(r=size/sqrt(3), h=10000, $fn=6);
  }
}

module split_ring(gap=0.2){
  if(gap>0){
   difference(){
     children();
     Tz(-10000)
     cube([10000,gap,20000]);
  }
  }else {children();}
}

module bumpyCy(d=1, h=1, n_bumps=2, h_bumps=3, w_bumps=3, $fn=0, $fa=6){
  delta_w = w_bumps/(d/2)*(180/PI)/2;
  cylindrical_bevel(r=d/2+h_bumps, h=h, bevel=h_bumps/2)
  union(){
    for(i = [0.5:1:n_bumps-.5])
    {
    w=i*360/n_bumps;
    //echo(delta_w, w);
    CyS(r=d/2 + h_bumps, h=h, w1 =w-delta_w, w2 = w+delta_w);
    } 
    Cy(r=d/2,h=h);
  }
}
//intersection()
{
//cube(1000);

//compound_planetary_gear(s=[11,15,11,12], p=[13,11,13,14], n_planets=4, z=[20,40,20,20], gap=2,
//compound_planetary_gear(s=[25,31,25,34], p=[17,17,17,23], n_planets=6, z=[10,20,10,5], gap=2,
//compound_planetary_gear(s=[12], p=[13],r=[36], n_planets=4, z=[5], gap=2, 
compound_planetary_gear(s=[17, 20], p=[18, 17],r=[53,50], n_planets=5, z=[20,20], gap=2, 
              helix_angle=herringbone(helix=30), 
              together_built=false,
              connect_suns=false, 
              od=55, od_padding=1.5,
              $fn=5);
}
