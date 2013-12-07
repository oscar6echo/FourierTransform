
// ********************************************************************************global variables

var Fourier_Series = {
	'square' : {
		coef: [0.75, 0, 0.25, 0, 0.15, 0, 0.107, 0],
		phi: [-0.5, 0, -0.5, 0, -0.5, 0, -0.5, 0]},
	'saw' : {
		coef: [-0.75, -0.375, -0.25, -0.188, -0.15, -0.125, -0.107, -0.094],
		phi: [-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5]},
	'triangle' : {
		coef: [0.75, 0, -0.083, 0, 0.03, 0., -0.015, 0],
		phi: [-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5]},
	'pulse' : {
		coef: [1, 1, 1, 1, 1, 1, 1, 1],
		phi: [0, 0, 0, 0, 0, 0, 0, 0]},
	'manual' : {
		coef: [0.414, 0.721, 0.378, 0.033, 0.688, 0.867, 0.702, 0.788],
		phi: [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]}
}

var circle_area_dim = {width: 400, height: 250},
	curve_area_dim = {width: 400, height: 250},
	margin = {top: 15, right: 10, bottom: 15, left: 10},
	width_curve_area = curve_area_dim.width-margin.left-margin.right,
	width_circle_area = circle_area_dim.width-margin.left-margin.right
	height = curve_area_dim.height-margin.top-margin.bottom,
	width_indiv_area = width_circle_area+width_curve_area;
	height_indiv_area = 100;
	height_coef_area = 50;


var x_curve = d3.scale.linear().range([0, width_curve_area]),
	x_circle_area = d3.scale.linear().range([0, width_circle_area]),
	r_circle = d3.scale.linear(),
	y = d3.scale.linear().range([height, 0]);
	x_indiv = d3.scale.linear(),
	y_indiv = d3.scale.linear().range([height_indiv_area, 0]),
	r_indiv = d3.scale.linear();


var x_Axis = d3.svg.axis().scale(x_curve).orient("bottom"),
	y_Axis = d3.svg.axis().scale(y).orient("left");

var svg = d3.select("#vis")
		.append("svg:svg")
			.attr("width", margin.left + width_circle_area + width_curve_area + margin.right)
			.attr("height", margin.top + height + margin.bottom + height_indiv_area + height_coef_area);

var circle_area = svg.append("svg:g")
		.attr("id", "circle_area")
		.attr("transform", "translate(" + margin.left + "," + margin.top + ")");

var curve_area = svg.append("svg:g")
		.attr("id", "curve_area")
		.attr("transform", "translate("+ (margin.left + circle_area_dim.width) + "," + margin.top + ")");

var indiv_circle_area = svg.append("svg:g")
		.attr("id", "indiv_circle_area")
		.attr("transform", "translate("+ margin.left  + "," + circle_area_dim.height + ")");

var coef_area = svg.append("svg:g")
		.attr("id", "coef_area")
		.attr("transform", "translate("+ margin.left  + "," + (circle_area_dim.height + height_indiv_area) +")");

var selected_coef_amp = null;


curve_area.append("svg:g")
		.attr("class", "x_Axis")
		.attr("transform", "translate("+ 0 +", "+ height/2 + ")");

curve_area.append("svg:g")
		.attr("class", "y_Axis")
		.attr("transform", "translate(" + 0+ "," + 0 + ")");

var curve_name = 'square',
	fourier_series = Fourier_Series[curve_name],
	dataset=[],
	time_series,
	fourier_vector_series,
	fourier_sum_series,
	r_max,
	Nb_Fourier_Coef = 8,
	timestep = 1/80,
	timestep_0 = timestep,
	nb_timestep = Nb_Fourier_Coef / timestep,
	period = 6000,
	time_series = d3.range(0, Nb_Fourier_Coef, timestep),
	x_Axis_time_span = period*1.8,
	nb_index_displayed = round(nb_timestep * x_Axis_time_span /period),
	t = 0.0,
	dir = +1,
	trail_fraction = 0.3;


var curve_line = d3.svg.line()
	.x(function(d, i) { return x_curve(i); })
	.y(function(d) { return y(dataset[d].full_sum.re); });

var trail_line = d3.svg.line()
	.x(function(d, i) { return x_circle_area(-dataset[d].full_sum.i); })
	.y(function(d) { return y(dataset[d].full_sum.re); });


// ********************************************************************************functions

function mod(a, b) { return ((a%b)+b)%b; }
function cos(d){ return Math.cos(d); }
function sin(d){ return Math.sin(d); }
function exp(d){ return Math.exp(d); }
function abs(d){ return Math.abs(d); }
function round(d){ return Math.round(d); }
function max(a, b){ return Math.max(a, b); }
function min(a, b){ return Math.min(a, b); }


bisect_time = d3.bisector(function(d) { return d; }).left;
var pi = Math.PI;

var	stroke = d3.scale.category20b();
function color(d, i) { return i==0  ? "black ": stroke(i); }

function transpose_fourier_series(d) {
	// return transposed fourier_series
	var transposed = [];
	for (var k=0; k<d.coef.length;  k++) {
		var obj_ = {};
		obj_.coef = d.coef[k];
		obj_.phi = d.phi[k];
		transposed.push(obj_);
	}
	return transposed;
}

function get_fourier_vector(fourier_series, t){
	// returns a dict of coef, cumul_sum; the length of the dict is that of the Fourier_Coef input
	var coef = fourier_series.coef;
	phi = fourier_series.phi;
	var N = coef.length;
	var coef_ = [];
	for (var k=0; k<N;  k++) {
		coef_[k] = Complex(coef[k], 0).mult(Complex.Polar(1, 2*pi*t*(1+k)/N+phi[k]*pi));
	}
	var cumul_sum_ = [];
	for (var k=0; k<N;  k++) {
		if (k>0){
			cumul_sum_[k] = coef_[k].add(cumul_sum_[k-1]);
		}
		else{
			cumul_sum_[k] = coef_[k];
		}
	}

	var obj_ = {};
	obj_.coef = coef_;
	obj_.cumul_sum = cumul_sum_;

	return obj_;
}

function get_fourier_data(curve_name){
	// return a dict of t, vector (it self a dict of coef, cumul_sum)
	var fourier_series = Fourier_Series[curve_name];
	var coef = fourier_series.coef;
	function get_fourier_vector_given_coef(t) {
		return get_fourier_vector(fourier_series, t);
	}

	var fourier_vector_series = time_series.map(get_fourier_vector_given_coef);
	var dataset = [];
	N = fourier_series.coef.length;

	for (var k=0; k<time_series.length;  k++) {
		var obj_ = {};
		var cumul_sum_ = fourier_vector_series[k].cumul_sum;
		var coef_ = fourier_vector_series[k].coef;
		full_sum_ = cumul_sum_.pop();
		cumul_sum_.unshift(Complex(0, 0));

		var vector_ = [];
		for (var q=0; q<N;  q++) {
			var obj2_ = {};
			// obj2_ = {};
			obj2_.coef = coef_[q];
			obj2_.cumul_sum = cumul_sum_[q];
			obj2_.initcoef = Complex.Polar(fourier_series.coef[q], pi*fourier_series.phi[q]);
			vector_.push(obj2_);
		}

		obj_.t = time_series[k];
		obj_.full_sum = full_sum_;
		obj_.vector = vector_;
		dataset.push(obj_);
	}
	return dataset;
}


function update_square() {
	console.log("square");
	curve_name = "square";
	dataset = [];
	// t=0;
	update();
}

function update_saw() {
	console.log("saw");
	curve_name = "saw";
	dataset = [];
	// t=0;
	update();
}

function update_triangle() {
	console.log("triangle");
	curve_name = "triangle";
	dataset = [];
	// t=0;
	update();
}

function update_pulse() {
	console.log("pulse");
	curve_name = "pulse";
	dataset = [];
	// t=0;
	update();
}

function update() {

	if (dataset.length==0) {
		dataset = get_fourier_data(curve_name);
	}

	// update displayed curve points
	var start_index = bisect_time(time_series, t);
	var end_index = start_index + nb_index_displayed;
	var rng_index_displayed = d3.range(start_index, end_index).map(function (d){ return mod(d, nb_timestep); })

	var time_series_displayed = rng_index_displayed.map(function (d){ return dataset[d].t; })
	var fourier_vector_displayed = dataset[start_index%nb_timestep];

	// update displayed trail points
	end_index = bisect_time(time_series, t);
	start_index = end_index - round(nb_index_displayed*trail_fraction);
	var rng_index_trail = d3.range(start_index, end_index).map(function (d){ return mod(d, nb_timestep); });


	// update scales
	var buffer = 1.15

	x_curve.domain([0, time_series_displayed.length-1]);

	var max_re =d3.max(dataset.map(function (d){ return abs(d.full_sum.re); }))
	var max_i =d3.max(dataset.map(function (d){ return abs(d.full_sum.re); }))
	if (max_re/height > max_i/width_circle_area) {
		y.domain([-max_re * buffer, max_re * buffer]);
		x_circle_area.domain([-max_re * buffer * width_circle_area/height, max_re * buffer * width_circle_area/height]);
	}
	else {
		x_circle_area.domain([-max_i * buffer, max_i * buffer]);
		y.domain([-max_i * buffer * height/width_circle_area, max_i * buffer * height/width_circle_area]);
	}

	r_circle.domain([0, d3.max(y.domain())-d3.min(y.domain())]);
	r_circle.range([0, d3.max(y.range())-d3.min(y.range())]);

	width_unit_indiv_area = width_indiv_area/fourier_vector_displayed.vector.length;
	x_indiv.range([-width_unit_indiv_area/2, width_unit_indiv_area/2]);
	max_r =d3.max(fourier_vector_displayed.vector.map(function (d){ return d.coef.r; }))
	if (width_unit_indiv_area<height_indiv_area) {
		x_indiv.domain([-max_r * buffer, max_r * buffer]);
		y_indiv.domain([-max_r * buffer * height_indiv_area/width_unit_indiv_area, max_r * buffer * height_indiv_area/width_unit_indiv_area]);
	}
	else {
		y_indiv.domain([-max_r * buffer, max_r * buffer]);
		x_indiv.domain([-max_r * buffer * width_unit_indiv_area/height_indiv_area, max_r * buffer * width_unit_indiv_area/height_indiv_area]);
	}

	r_indiv.domain([0, d3.max(x_indiv.domain())-d3.min(x_indiv.domain())]);
	r_indiv.range([0, d3.max(x_indiv.range())-d3.min(x_indiv.range())]);

	// draw curve
	var curve = curve_area.selectAll(".curve")
			.data([rng_index_displayed]);

	curve.enter()
		.append("path")
			.attr("class", "curve")
			.attr("d", curve_line);

	curve.transition()
			.duration(0)
			.attr("d", curve_line);

	// draw trail
	var trail = circle_area.selectAll(".trail")
			.data([rng_index_trail]);

	trail.enter()
		.append("path")
			.attr("class", "trail")
			.attr("d", trail_line);

	trail.transition()
			.duration(0)
			.attr("d", trail_line);


	// draw coef circles in circle_area
	var circles = circle_area.selectAll(".coef_circle")
			.data(fourier_vector_displayed.vector);

	circles.enter()
			.append("svg:circle")
			.attr("class", "coef_circle")
			.attr("r", function(d) { return r_circle(d.coef.r); })
			.attr("cx", function(d) { return x_circle_area(-d.cumul_sum.i); })
			.attr("cy", function(d) { return y(d.cumul_sum.re); });

	circles.transition()
			.duration(0)
			.attr("r", function(d) { return r_circle(d.coef.r); })
			.attr("cx", function(d) { return x_circle_area(-d.cumul_sum.i); })
			.attr("cy", function(d) { return y(d.cumul_sum.re); });

	// draw coef value points in circle_area
	var points = circle_area.selectAll(".coef_point")
			.data(fourier_vector_displayed.vector);

	points.enter()
			.append("svg:circle")
			.attr("class", "coef_point")
			.attr("r", function (d, i) { return 1.5; })
			.attr("cx", function(d) { return x_circle_area( - d.cumul_sum.i - d.coef.i ); })
			.attr("cy", function(d) { return y( d.cumul_sum.re + d.coef.re ); });

	points.transition()
			.duration(0)
			.attr("r", function (d, i) { return 1.5; })
			.attr("cx", function(d) { return x_circle_area( - d.cumul_sum.i - d.coef.i ); })
			.attr("cy", function(d) { return y( d.cumul_sum.re + d.coef.re ); });

	// draw lines from cumul_sum to coef in circle_area
	var coef_lines = circle_area.selectAll(".coef_line")
			.data(fourier_vector_displayed.vector);

	coef_lines.enter()
			.append("svg:line")
			.attr("class", "coef_line")
			.attr("x1", function (d){ return x_circle_area( - d.cumul_sum.i ); })
			.attr("y1", function (d){ return y( d.cumul_sum.re ); })
			.attr("x2", function (d){ return x_circle_area( - d.cumul_sum.i - d.coef.i ); })
			.attr("y2", function (d){ return y( d.cumul_sum.re + d.coef.re ); });

	coef_lines.transition()
			.duration(0)
			.attr("x1", function (d){ return x_circle_area( - d.cumul_sum.i ); })
			.attr("y1", function (d){ return y( d.cumul_sum.re ); })
			.attr("x2", function (d){ return x_circle_area( - d.cumul_sum.i - d.coef.i ); })
			.attr("y2", function (d){ return y( d.cumul_sum.re + d.coef.re ); });

	// draw one link line between circle and curve areas
	var link_line = circle_area.selectAll(".link_line")
			.data([fourier_vector_displayed]);

	link_line.enter()
			.append("svg:line")
			.attr("class", "link_line")
			.attr("x1", function (d){ return x_circle_area( - d.full_sum.i ); })
			.attr("y1", function (d){ return y( d.full_sum.re ); })
			.attr("x2", function (d){ return circle_area_dim.width; })
			.attr("y2", function (d){ return y( d.full_sum.re ); });

	link_line.transition()
			.duration(0)
			.attr("x1", function (d){ return x_circle_area( - d.full_sum.i ); })
			.attr("y1", function (d){ return y( d.full_sum.re ); })
			.attr("x2", function (d){ return circle_area_dim.width; })
			.attr("y2", function (d){ return y( d.full_sum.re ); });

	// draw link point at the end of link line on y axis
	var link_point = circle_area.selectAll(".link_point")
			.data([fourier_vector_displayed]);

	link_point.enter()
			.append("svg:circle")
			.attr("class", "link_point")
			.attr("r", function (d, i) { return 4; })
			.attr("cx", function(d) { return circle_area_dim.width; })
			.attr("cy", function(d) { return y( d.full_sum.re ); });

	link_point.transition()
			.duration(0)
			.attr("r", function (d, i) { return 4; })
			.attr("cx", function(d) { return circle_area_dim.width; })
			.attr("cy", function(d) { return y( d.full_sum.re ); });

	// draw x axis
	curve_area.selectAll(".x_Axis")
		.transition()
			.duration(0)
			.call(x_Axis.ticks(0))
			.selectAll("text")
				.attr("class", "font_axis")
				.attr("class", "invisible")
				.style("text-anchor", "center")
				.attr("dx", "0.0em")
				.attr("dy", "0.8em");

	// draw y axis
	curve_area.selectAll(".y_Axis")
		.transition()
			.duration(0)
			.call(y_Axis.ticks(0))
			.selectAll("text")
				.attr("class", "font_axis")
				.attr("class", "invisible")
				.style("text-anchor", "end")
				.attr("dx", "-0.2em")
				.attr("dy", "0em");

	// draw individual circles
	var indiv_circle_obj = indiv_circle_area.selectAll("g")
			.data(fourier_vector_displayed.vector);

	indiv_circle_obj.enter()
		.append("g")
			.attr("transform", function (d, i) { return "translate(" + (width_indiv_area/fourier_series.coef.length*(i+0.5)) + "," + 0 + ")"});

	var indiv_circle = indiv_circle_obj.selectAll(".indiv_circle")
			.data(function (d, i) { return [d] });

	indiv_circle.enter()
		.append("svg:circle")
			.attr("class", "indiv_circle")
			.attr("r", function(d, i) { return r_indiv(d.coef.r); })
			.attr("cx", function(d, i) { return x_indiv(0); })
			.attr("cy", function(d, i) { return y_indiv(0); });

	indiv_circle.transition()
			.duration(0)
			.attr("r", function(d, i) { return r_indiv(d.coef.r); })
			.attr("cx", function(d, i) { return x_indiv(0); })
			.attr("cy", function(d, i) { return y_indiv(0); });


	var indiv_initphi_line = indiv_circle_obj.selectAll(".indiv_initphi_line")
			.data(function (d, i) { return [d] });

	indiv_initphi_line.enter()
		.append("svg:line")
			.attr("class", "indiv_initphi_line")
			.attr("x1", function (d){ return x_indiv(0); })
			.attr("y1", function (d){ return y_indiv(0); })
			.attr("x2", function (d, i){ return x_indiv(-d.initcoef.i); })
			.attr("y2", function (d, i){ return y_indiv(d.initcoef.re); });

	indiv_initphi_line.transition()
			.duration(0)
			.attr("x1", function (d){ return x_indiv(0); })
			.attr("y1", function (d){ return y_indiv(0); })
			.attr("x2", function (d, i){ return x_indiv(-d.initcoef.i); })
			.attr("y2", function (d, i){ return y_indiv(d.initcoef.re); });


	var indiv_circle_point = indiv_circle_obj.selectAll(".indiv_circle_point")
			.data(function (d, i) { return [d] });

	indiv_circle_point.enter()
		.append("svg:circle")
			.attr("class", "indiv_circle_point")
			.attr("r", 2)
			.attr("cx", function(d, i) { return x_indiv(-d.coef.i); })
			.attr("cy", function(d, i) { return y_indiv(d.coef.re); });

	indiv_circle_point.transition()
			.duration(0)
			.attr("r", 1.5)
			.attr("cx", function(d, i) { return x_indiv(-d.coef.i); })
			.attr("cy", function(d, i) { return y_indiv(d.coef.re); });

	var indiv_circle_line = indiv_circle_obj.selectAll(".indiv_circle_line")
			.data(function (d, i) { return [d] });

	indiv_circle_line.enter()
		.append("svg:line")
			.attr("class", "indiv_circle_line")
			.attr("x1", function (d){ return x_indiv(0); })
			.attr("y1", function (d){ return y_indiv(0); })
			.attr("x2", function (d){ return x_indiv(-d.coef.i); })
			.attr("y2", function (d){ return y_indiv(d.coef.re); });

	indiv_circle_line.transition()
			.duration(0)
			.attr("x1", function (d){ return x_indiv(0); })
			.attr("y1", function (d){ return y_indiv(0); })
			.attr("x2", function (d){ return x_indiv(-d.coef.i); })
			.attr("y2", function (d){ return y_indiv(d.coef.re); });

	// draw buttons
	var coef_obj = coef_area.selectAll("g")
			.data(transpose_fourier_series(fourier_series));

	coef_obj.enter()
		.append("g")
			.attr("transform", function (d, i) { return "translate(" + (width_indiv_area/fourier_series.coef.length*(i+0.5)) + "," + 0 + ")"})
		.append("rect")
			.attr("class", "coef_box")
			.attr("x", -0.9*width_indiv_area/fourier_series.coef.length/2)
			.attr("y", 0)
			.attr("width", 0.9*width_indiv_area/fourier_series.coef.length)
			.attr("height", height_indiv_area/2.5);

	var coef_display = coef_obj.selectAll(".coef")
			.data(function(d) { return [d]; });

	coef_display.enter()
		.append("svg:text")
			.attr("class", "coef" )
			.attr("x", 0)
			.attr("y", height_coef_area/2)
			.text(function(d) { return d.coef+", "+d.phi; });

	coef_display.transition()
			.duration(0)
			.text(function(d) { return d.coef+", "+d.phi; });

}


// key function definitions and attach key function to body
function keydown_function() {
	console.log("d3.event.keyCode = ", d3.event.keyCode);
	if (!selected_coef_amp) {
		// right
		if (d3.event.keyCode==39) {
			timestep *= 1.1;
		}
		// left
		if (d3.event.keyCode==37) {
			timestep /= 1.1;
		}
		// down
		if (d3.event.keyCode==40) {
			timestep = timestep_0 / 5;
		}
		// up
		if (d3.event.keyCode==38) {
			timestep = timestep_0;
		}
		// L
		if (d3.event.keyCode==76) {
			trail_fraction = min(trail_fraction+0.02, 1.0);
		}
		// S
		if (d3.event.keyCode==83) {
			trail_fraction = max(0, trail_fraction-0.02);
		}
	}
}
d3.select("body").on("keydown", keydown_function);


// timing
function next_step() {
	t += dir*timestep;
	if (t<0) {
		t = Nb_Fourier_Coef;
	}
	if (t>Nb_Fourier_Coef) {
		t = 0;
	}
	update();
}
var myTimer = setInterval(function() { next_step() }, period*timestep_0/Nb_Fourier_Coef);

