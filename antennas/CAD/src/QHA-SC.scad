
// http://www.rish.kyoto-u.ac.jp/digitalbeacon/information/Building_QFH_Antenna_Guide.pdf
freq = 143050000;
c = 300000000000.0;  // speed of light in mm per second. 
lambda =  c / freq;  
big_height = 0.26*lambda;
small_height =  0.238*lambda;
wire_r = (0.0088 * lambda)/2;
base_height = 50;
base_side = 2000;
base_margin = 10; 

echo(lambda);


module loop() {
  height = big_height;
  diameter = 0.173*lambda;
  leg_size = 0.560*lambda;

  union() {
    translate([0,0,height]) rotate([0,90,0]) cylinder(r=wire_r, h=diameter + 2 * wire_r, center=true);
    linear_extrude(height = height, twist = 180)
      union () {
        translate ([diameter/2,0]) circle (r=wire_r);
        translate ([-diameter/2,0]) circle (r=wire_r);
      }
  }
}

module quad_helix() {   
  union() {
    z = (big_height - 0.238*lambda)/2;
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

union() {
    translate([0,0,base_height]) 
        quad_helix();
    base();
}
  
