// Bolidozor short-circuited QHA antenna

include <shortcuts.scad>  // http://www.thingiverse.com/thing:644830
$fn = 100;

wire_r = 19/2;
base_height = 50;
base_side = 2000;
base_margin = 10; 

height = 0.6*1000;
diameter = 0.3899*1000;

helix_circumference = (diameter*PI);

helix_tube_lenght = sqrt(pow(helix_circumference/2,2) + pow(height,2));    // pythogarian calculation of helix_tube_lenght. 

echo(helix_tube_lenght);

helix_angle = atan(height/(helix_circumference/2));
echo(helix_angle);


module elliptic_ring(r1 = 10, r2 = 5, r = 2, slices = 100, h = 0, w = 360)
{
  dz = h/slices; 
   dwx = atan(h/(r1+r2)/PI); 
  for (i = [0:slices-1])
  hull()
  {
    T(r1*cos((i+1)*w/slices), r2*sin((i+1)*w/slices), (i+1)*dz)
     R(90+dwx, 0, (i+1)*w/slices)
    Cy(r = r, h = .01); 
  
    T(r1*cos((i)*w/slices), r2*sin((i)*w/slices), i*dz)
     R(90+dwx, 0, i*w/slices)
    Cy(r = r, h = .01); 
  }
}


module loop() {
  union() {
    rotate([0,90,0])
      cylinder(r=wire_r, h=diameter + 2 * wire_r, center=true);
    translate([0,0,height]) 
      rotate([0,90,0]) 
        cylinder(r=wire_r, h=diameter + 2 * wire_r, center=true);
    
    elliptic_ring(h = height, w = 180, r1 = diameter/2, r2 = diameter/2, r = wire_r);  // partial ring
    rotate([0,0,180])
      elliptic_ring(h = height, w = 180, r1 = diameter/2, r2 = diameter/2, r = wire_r);  // partial ring
  }
}

module quad_helix() {   
  union() {
    //z = (big_height - 0.238*lambda)/2;
    loop();
    rotate([0,0,90]) loop();
  }
}



module base() {
translate([base_margin/2, base_margin/2,0]) cube([base_side/2, base_side/2, base_height], center = false);
translate([-base_margin/2 - base_side/2, base_margin/2,0]) cube([base_side/2, base_side/2, base_height], center = false);
translate([-base_margin/2 - base_side/2,-base_margin/2 - base_side/2, 0]) cube([base_side/2, base_side/2, base_height], center = false);
translate([base_margin/2,-base_margin/2 - base_side/2, 0]) cube([base_side/2, base_side/2, base_height], center = false);
}

module antenna(){
  union() {
      translate([0,0,2*base_height]) 
        quad_helix();
      translate([0,0,base_height + 40])
        color("blue") cube([150, 150, 80], center = true);
      base();
  }
}  

module manufacturing_helper(part = 1){
x_size = diameter/5;
y_size = diameter/5;

  difference() {
    translate([0, diameter/2, height/2])
      rotate([0, helix_angle, 0])
      {
        cube([x_size, y_size, wire_r*5], center = true);

      }
    loop();
    translate([0, diameter/2, height/2])
      rotate([0, helix_angle, 0])
      {
        if (part == 1)
          translate([0, 0, (wire_r*5)/2])        
            cube([diameter, diameter, wire_r*5], center = true);  // half cut
        else
          translate([0, 0, -(wire_r*5)/2])        
            cube([diameter, diameter, wire_r*5], center = true);  // half cut


        // mountig holes
        translate([x_size/3, y_size/3, 0])        
          cylinder(r=3, h=diameter, center=true);

        translate([x_size/3, -y_size/3, 0])        
          cylinder(r=3, h=diameter, center=true);

        translate([-x_size/3, -y_size/3, 0])        
          cylinder(r=3, h=diameter, center=true);

        translate([-x_size/3, y_size/3, 0])        
          cylinder(r=3, h=diameter, center=true);
      }
  }  
}

//manufacturing_helper(part = 1);
antenna();
