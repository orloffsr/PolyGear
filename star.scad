// With recursive functions very complex results can be generated,
// e.g. calculating the outline of a star shaped polygon. Using
// default parameters can help with simplifying the usage of functions
// that have multiple parameters.
function point(angle) = [ sin(angle), cos(angle) ];
function radius(i, r1, r2) = (i % 2) == 0 ? r1 : r2;
function star_f(count, r1, r2, i = 0, result = []) = i < count
    ? star_f(count, r1, r2, i + 1, concat(result, [ radius(i, r1, r2) * point(360 / count * i) ]))
    : result;


module star(points, r1, r2, h=1){
    count = points*2;
        linear_extrude(height = h)
            polygon(star_f(count, r1, r2));
}

//color("yellow")
//    //translate([0, 50, 0])
//        linear_extrude(height = 1)
//            polygon(star_f(30, 40, 10));

star(6, 40, 1, h=1);
