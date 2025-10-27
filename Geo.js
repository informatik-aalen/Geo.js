
/*
 Stand 20241127
 
 ToDo:
 - Klasse Algo:
	- Kruskal
	- Diverse Grafik-Algos:
		- Orientierung Polygon
		- Finde nächstes Objekt
 Klasse Datenstrukturen:
	- KdTree
 
 Nomenklatur für Winkel:
    phi, lambda, alpha: math., Bogenmaß
	lat, lon: WGS85, Gradmaß
 
 
 */
const cos = Math.cos, sin = Math.sin, tan = Math.tan;
const acos = x => Math.acos(Math.min(Math.max(x, -1), 1));
const asin = x => Math.asin(Math.min(Math.max(x, -1), 1));
const atan = Math.atan, atan2 = Math.atan2;
const abs = Math.abs, signum = Math.sign, log = Math.log;
const M_PI = Math.PI,
	M_PI_2 = Math.PI / 2,
	M_PI_4 = Math.PI / 4,
	M_PI_3_4 = Math.PI * 3 / 4;

class Coordinates {
	/**
	 */
	constructor() {}

	/**
	 * International airports
	 * @type {object}
	 */
	get airports() {
		return {
			EDDS: [9.221944,48.69], // Stuttgart
			KJFK: [-73.77888888888889, 40.63972222222222], // NY JFK
			YSSY: [151.177222, -33.946111] // Sydney
		};
	}
	/**/
	/**
	 * Buildings
	 * @type {object}
	 */
	get buildings () {
		return {
			HS_AA_G2: [[10.067109,48.841427], // Nord
					   [10.0664974,48.840523], // West
					   [10.0667579,48.8404466], // // Sued
					   [10.0673695,48.8413506] // Ost
			],
			TV_Tower_Stuttgart: [[9.190166666666666,48.75575]]
		};
	}
	get nature() {
		return {
		lake_of_constance: {west:  [9.0307915, 47.8071598], east: [9.7537997, 47.5151325]},
		};
	}

	/**
	 * Cities
	 * @type {object}
	 */
	get cities()  {
		return {
			berlin: [13.40,52.517], // https://de.wikipedia.org/wiki/Orthodrome
			tokio: [139.767, 35.70]// https://de.wikipedia.org/wiki/Orthodrome
		}
	}

	
	/**
	 * @param {array} c - Coordinates (each coordinate is array with length 2)
	 * @return {array} Bounding Box [left, bottom, right, ceiling]
	 */
bbox (c)  {
	let res = [c[0][0],c[0][1],c[0][0],c[0][1]]
	for (let i = 1; i < c.length; i++) {
		if (c[i][0] < res[0])
			res[0] = c[i][0];
		else if (c[i][0] > res[2])
			res[2] = c[i][0];
		if (c[i][1] < res[1])
			res[1] = c[i][1];
		else if (c[i][1] > res[3])
			res[3] = c[i][1];
	}
	return res;
}
}


class c_Linesegment{
	#P1; #P2;
	constructor(P1, P2) {
		this.#P1 = new c_Point(P1), this.#P2 = new c_Point(P2);
	}
	dist(o) {
		if (o instanceof c_Point) {
			const v_p1_p2 = new c_Vector(this.#P1, this.#P2);
			let v = new c_Vector(this.#P1, o);
			return (v_p1_p2.dot(v) <= 0) ? v.abs : (
				(v_p1_p2.neg.dot(v = new c_Vector(this.#P2, o)) <= 0) ? v.abs
				: v_p1_p2.neg.cross(v).abs / v_p1_p2.abs);
		}
	}
	intersection(o) {
		if (o instanceof c_Linesegment)
			;
	}
}

class c_Point {
	#c;
	constructor(c) { // Bisher nur: Vektor, Punkt, Array
		if (c instanceof c_Vector || c instanceof c_Point)
			this.#c = c.c.slice();
		else
			this.#c = c.slice();
	}
	get c() {
		return this.#c;
	}
};

class c_Line {
	/** https://de.wikipedia.org/wiki/Gerade
	 
	 */
	#a;#b;#c;
	/*
	 
	 */
	constructor(x, y, z) {
		if (x instanceof c_Point && y instanceof c_Point) {// Zweipunktform
			this.#a = x.c[1] - y.c[1], this.#b = y.c[0] - x.c[0], this.#c = y.c[0] * x.c[1] - x.c[0] * y.c[1];
		}
		else if (x instanceof Array && y instanceof Array) {// Zweipunktform
			this.#a = x[1] - y[1], this.#b = y[0] - x[0], this.#c = y[0] * x[1] - x[0] * y[1];
		}
		else if  (x instanceof c_Line) // Gerade
			this.#a = x.#a, this.#b = x.#b, this.#c = x.#c;
		else if (x instanceof c_Point && y instanceof c_Vector) {// Punkt-Vektor
			const p = new c_Point(new c_Vector(x).add(y));
			const l = new c_Line(x, p);
			this.#a = l.#a, this.#b = l.#b, this.#c = l.#c;
		}
		else if (x instanceof c_Point && typeof(y) == 'number') {// Punkt-Winkel
			y = new c_Vector([cos(y), sin(y)]);
			const p = new c_Point(new c_Vector(x).add(y));
			const l = new c_Line(x, p);
			this.#a = l.#a, this.#b = l.#b, this.#c = l.#c;
		}
		else if (typeof x == 'number' && typeof y == 'number')
		{
			if (typeof z == 'number') // Koordinatenform
				this.#a = x, this.#b = y, this.#c = z;
			else //
				this.#a = x, this.#b = -1, this.#c = -y; // tbd
		}
		else
			console.log("error");
//		console.log(this.#a, this.#b, this.#c);
	}
	get c() {
		return [this.#a, this.#b, this.#c];
	}
	y(x) {
		return (this.#b) ? (this.#c - this.#a * x) / this.#b : undefined;
	}
	x(y) {
		return (this.#a) ? (this.#c - this.#b * y) / this.#a : undefined;
	}
	intersection(line) {
		const denominator = this.#a * line.#b - line.#a * this.#b;
		return (denominator) ? new c_Point([
			(this.#c * line.#b - line.#c * this.#b) / denominator,
			(this.#a * line.#c - line.#a * this.#c) / denominator])
		: undefined;
	}
	get s() {
		return -this.#a / this.#b;
	}
	get alpha() {
		return atan2(-this.#a, this.#b);
	}
}

class c_Circle{
	#C; #r;
	constructor(C, r) {
		if (C instanceof c_Point)
			console.log(C.c), this.#C = C.c.slice(), this.#r = r;
		else
		this.#C = C, this.#r = r;
		console.log("===>", this.#C);
	}
	
	get C() {
		return this.#C;
	}
	get r() {
		return this.#r;
	}
	powerline(re) {
		return new c_Line(2 * (re.#C[0] - this.#C[0]), 2 * (re.#C[1] - this.#C[1]),
						  this.#r * this.#r
						  - this.#C[0] * this.#C[0]
						  - this.#C[1] * this.#C[1]
						  - re.#r * re.#r
						  + re.#C[0] * re.#C[0]
						  + re.#C[1] * re.#C[1]);
	}
};

class c_Vector {
	#c;
	constructor(c, d) {
		if ((c instanceof c_Vector || c instanceof c_Point) && d == undefined)
			this.#c = c.c.slice();
		else if (c.length) // Array
			this.#c = c.slice();
		else if (c instanceof c_Point && d instanceof c_Point)
			this.#c = [d.c[0] - c.c[0], d.c[1] - c.c[1]];
		else
			console.log("error constructor c_Vector")
/*		else if (typeof y == "number") // polar, r, alpha
			this.c = [d * cos(]*/
	}
	get c() {
		return this.#c;
	}

	add(b) {
		let res = [];
		for (let i = 0; i < this.#c.length; i++)
			res.push(this.c[i] + b.c[i]);
		return new c_Vector(res);
	}

	sub(b) {
		let res = [];
		for (let i = 0; i < this.#c.length; i++)
			res.push(this.c[i] - b.c[i]);
		return new c_Vector(res);
	}

	mul(b) { // b ist skalar
		let res = [];
		for (let i = 0; i < this.#c.length; i++)
			res.push(this.c[i] * b);
		return new c_Vector(res);
	}
	
	dot(right) { // Skalarprodukt zweier Vektoren
		let res = 0;
		for (let i = 0; i < this.#c.length; i++)
			res += this.#c[i] * right.#c[i];
		return res;
	}
	get neg() { // Antivector
		return new c_Vector([-this.#c[0], -this.#c[1]]);
	}
	
	angle(right) { // Skalarprodukt zweier Vektoren tbd
		//console.log(this.abs, right.abs, 23123123);
		return acos(this.dot(right) / (this.abs * right.abs));
	}
	
	cross(right) { // Kreuzproduk // Skalarprodukt zweier Vektoren
		const a = this.c.slice(), b = right.c.slice();
		if (a.length == 2)
			a.push(0), b.push(0);
		return new c_Vector([a[1] * b[2] - a[2] * b[1],
					 a[2] * b[0] - a[0] * b[2],
					 a[0] * b[1] - a[1] * b[0]]);
	}
	
	get abs() {
		return Math.sqrt(this.dot(this));
	}
	
	absqr(x) {
		return this.dot(this);
	}
	N() {
		if (o.length == 2) { // Polarkoordinaten Radiant
			if (o[1] == M_PI_2 || o[1] == -M_PI_2)
				return undefined;
			if (o[0] >= 0)
				return [M_PI_2 - o[0], o[1] + (o[1] > 0) ? - M_PI : M_PI]
				return [o[0] + M_PI_2, o[1]];
		}
		else if (o.length == 3) { // Kartesische Koordinaten
			if (o[0] == 0 && o[1] == 0)
				return undefined;
			let r = Math.sqrt(o[0] * o[0] + o[1] * o[1]), c = -o[2] / r;
			return [o[0] * c , o[1] * c, r];
		}
		return undefined;
	}
	/** // Ostkursvektor
	 * @param {array} c
	 */
	E() {
		if (o.length == 2) { // Polarkoordinaten Radiant
			if (o[1] == M_PI_2 || o[1] == -M_PI_2)
				return undefined;
			return [0, o[1] + (o[1] > M_PI_4) ? - M_PI_3_4 : M_PI_4]
		}
		else if (o.length == 3) { // Kartesische Koordinaten
			if (o[0] == 0 && o[1] == 0)
				return undefined;
			return geo.vector.norm([-o[1], o[0], 0]);
		}
		return undefined;
	}
	rotation_matrix(v, alpha) {
		v = this.norm(v);
		let s = Math.sin(alpha), c = Math.cos(alpha), c_ = 1 - c;
		let v_0_c_ = v[0] * c_, v_1_c_ = v[1] * c_, v_2_c_ = v[2] * c_;
		return [
			[v[0]*v_0_c_ + c, 		v[0]*v_1_c_ - v[2]*s, 	v[0]*v_2_c_ + v[1]*s],
			[v[1]*v_0_c_ + v[2]*s, 	v[1]*v_1_c_ + c, 		v[1]*v_2_c_ - v[0]*s],
			[v[2]*v_0_c_ - v[1]*s, 	v[2]*v_1_c_ + v[0]*s, 	v[2]*v_2_c_ + c]];
	}
	norm(x) {
		return this.div(x, this.abs(x));
	}
	get alpha() {
//		console.log(this.#b, this.#a);
//		return atan2(this.#b, this.#a) / M_PI * 180;
	}
	get rwk() {
		return 90 - this.alpha;
	}

}

class c_Matrix  {
	#c;
	constructor(o) {
		if (o instanceof c_Matrix) {
			this.#c = [];
			o.c.forEach(v => this.#c.push(v.slice()));
		}
		else if (o instanceof c_Vector) {
			this.#c = [];
			o.c.forEach(v => this.#c.push([v]));
		}
		else { // Array
			this.#c = [];
			o.forEach(v => this.#c.push(v.slice()));
		}
	}
	get c() {
		return this.#c;
	}
	row(r) {
		return new c_Vector(this.#c[r]);
	}
	col(c) {
		let dummy = new Array(this.#c.length);
		for (let r = 0; r < this.#c.length; r++)
			dummy[r] = this.#c[r][c];
		return new c_Vector(dummy);
	}

	mul(y) {
		if (y instanceof c_Matrix) { // result c_Matrix
			//console.log("==>", this.c, y.c);
			let dummy = [];
			for (let row = 0; row < this.#c.length; row++) {
				let dummydummy = [];
				for (let col = 0; col < y.#c[0].length; col++)
					dummydummy.push(this.row(row).dot(y.col(col)));
				dummy.push(dummydummy);
			}
			return new c_Matrix(dummy);
		}
		else if (y instanceof c_Vector) { // result c_Vector
			let dummy = [];
			for (let row = 0; row < this.#c.length; row++) {
				dummy.push(this.row(row).dot(y));
			}
			return new c_Vector(dummy);

		}
		else if (typeof y == "number"){ // result c_Matrix
			let res = new c_Matrix(this);
			for (let row = 0; row < this.#c.length; row++)
				for (let col = 0; col < this.#c[0].length; col++)
					res.#c[row][col] *= y;
			return res;
		}
		else
			return undefined;
	}
	/**
	 */
	/**
	 */
/*	rot(v,r, a) {
		return this.mul(this.rotation_matrix(r, a), v);
	}*/
	/**
	 */
	inv() {
		let det_M = this.det();
		if (!det_M)
			return undefined;
		if (this.#c[0].length == 3) {
			// https://www.wikihow.com/Find-the-Inverse-of-a-3x3-Matrix
			
			let M_T = [
				[this.#c[0][0], this.#c[1][0], this.#c[2][0]],
				[this.#c[0][1], this.#c[1][1], this.#c[2][1]],
				[this.#c[0][2], this.#c[1][2], this.#c[2][2]]
			];
			var M_Adj = [
				[+ M_T[1][1] * M_T[2][2] - M_T[1][2] * M_T[2][1],
				   - M_T[1][0] * M_T[2][2] + M_T[1][2] * M_T[2][0],
				   + M_T[1][0] * M_T[2][1] - M_T[1][1] * M_T[2][0]
				],
				[- M_T[0][1] * M_T[2][2] + M_T[0][2] * M_T[2][1],
				   + M_T[0][0] * M_T[2][2] - M_T[0][2] * M_T[2][0],
				   - M_T[0][0] * M_T[2][1] + M_T[0][1] * M_T[2][0]
				],
				[+ M_T[0][1] * M_T[1][2] - M_T[0][2] * M_T[1][1],
				   - M_T[0][0] * M_T[1][2] + M_T[0][2] * M_T[1][0],
				   + M_T[0][0] * M_T[1][1] - M_T[0][1] * M_T[1][0]
				]
			];
		}
		else { // 2*2
			var M_Adj = new c_Matrix([
				[ this.#c[1][1], -this.#c[0][1]],
				[-this.#c[1][0],  this.#c[0][0]]]);
				//console.log("==>", M_Adj.c);
		}
		return M_Adj.mul(1 / det_M);
	}
	/**
	 */
	det() {
		if (this.#c[0].length == 3)
			return this.#c[0][0] * this.#c[1][1] * this.#c[2][2]
			+ this.#c[0][1] * this.#c[1][2] * this.#c[2][0]
			+ this.#c[0][2] * this.#c[1][0] * this.#c[2][1]
			- this.#c[0][0] * this.#c[1][2] * this.#c[2][1]
			- this.#c[0][1] * this.#c[1][0] * this.#c[2][2]
			-  this.#c[0][2] * this.#c[1][1] * this.#c[2][0];
		return this.#c[0][0] * this.#c[1][1] - this.#c[0][1] * this.#c[1][0];
	}
	/**
	 */
	intersect(P1, v1, P2, v2) {
		let A_inv = this.inv([
			[v1[0], -v2[0]],
			[v1[1], -v2[1]]
		]);
		//console.log(A_inv);
		let a_b = this.mul(A_inv, this.sub(P2, P1));
		return this.add(P1, this.mul(v1, a_b[0]));
	}
}

class c_Transform {
	/**
	 */
	constructor() {}
	
	get grad2rad_f() {return M_PI / 180;}

	/**
	 */
linear_transform_2d (M, p) {
	var res;
	if (typeof p[0] == "number")
		res =  [M[0][0]*p[0] + M[0][1]*p[1] + M[0][2],
				M[1][0]*p[0] + M[1][1]*p[1] + M[1][2]];
	else {
		res = [];
		p.forEach((p) => res.push([M[0][0]*p[0] + M[0][1]*p[1] + M[0][2],
								   M[1][0]*p[0] + M[1][1]*p[1] + M[1][2]]));
	}
	return res;
}
	/**
	 A, B; Bounding Boxen
	 */

linear_transform_2d_fit (A, B) {
	let s_x=((B[0] - B[2]) / (A[0]-A[2])), d_x;
	let s_y=((B[1] - B[3]) / (A[1]-A[3])), d_y;
	if (Math.abs(s_x) < Math.abs(s_y)) {
		s_y = Math.sign(s_x) * s_x;
		d_x = B[2] - (A[2] * s_x);
		d_y = (B[1] + B[3] - (A[1] + A[3]) * s_y) / 2;
	}
	else {
		s_x = abs(s_y);
		d_y = B[3] - (A[3] * s_y);
		d_x = (B[2] + B[0] - (A[2] + A[0]) * s_x) / 2;
	}
	return [new c_Matrix([[s_x, 0], [0, s_y]]), new c_Vector([d_x, d_y])];
}
	/**
	 */
grad2rad(o) {
	let e = undefined;
	if (typeof(o.length) != "undefined") {
		e = [];
		for (let p of o)
			e.push(this.grad2rad(p));
	}
	else if (typeof o == "object") {
		e = {}
		for (let p in o)
			e[p] = this.grad2rad(o[p]);
	}
	else if (typeof o == "number")
		e = o * this.grad2rad_f;
	return e;
}
	/**
	 */
rad2grad (o) {
	let e = undefined;
	if (typeof o == "object") {
		e = [];
		for (let p of o)
			e.push(this.rad2grad(p));
	}
	else if (typeof o == "number")
		e = o / this.grad2rad_f;
	return e;
}
	/**
	 * @param {array} polar coordinates (lon, lat) in radiant
	 * @return {array} cartesian ccordinates
	 */
polar2kart_3d(k) {
	let res;
	if (typeof k[0] == "number") {
		let c = Math.cos(k[1]);
		res = [c * Math.cos(k[0]),
				c * Math.sin(k[0]),
				Math.sin(k[1])];
	}
	else {
		res = [];
		for (let kk of k)
			res.push(this.polar2kart_3d(kk));
	}
	return res;
}
	/**
	 */
kart2polar_3d(k) {
		return [Math.atan2(k[1], k[0]), Math.asin(k[2])];
	}
	/**
	 */
gms2grad(g,m,s,hemi) {
		return ("WwSs".indexOf(hemi) >= 0)
			? -(g + ((m) ? m / 60 : 0) + (s ? s / 3600 : 0))
			: +(g + ((m) ? m / 60 : 0) + (s ? s / 3600 : 0));
	}
	/**
	 */
grad2gms(t) {
		return [t >> 0, t * 60 % 60 >> 0, t * 3600 % 60];
	}
	/**
	 */
EPSG_4326_2_EPSG_3857(k) {
		return [
			geo.transform.grad2rad(k[0]) * geo.nav.r,
			Math.atanh(Math.sin(geo.transform.grad2rad(k[1]))) * geo.nav.r
		]
	}
	/**
	 */
EPSG_3857_2_EPSG_4326(k) {
		return [
			geo.transform.rad2grad(k[0] / geo.nav.r),
			geo.transform.rad2grad(Math.atan(Math.sinh(k[1] / geo.nav.r)))
		];
	}
}


class c_Nav {
	/**
	 Blubber
	 */
	constructor() {}
	get r() {return 6378137;}
	get M_EPSG_3857_MAX() {return this.r * M_PI;}

	
	delta_lambda_wgs84_if_else = (l_S, l_Z) => {
		let delta = l_Z - l_S;
		if (delta > 180)
			delta -= 360;
		else if (delta < -180)
			delta += 360;
		return delta;
	}
	
	delta_lambda_wgs84   = (l_S, l_Z) => (l_Z - l_S + 540) % 360 - 180;
	delta_lambda_wgs84_2 = (l_S, l_Z) => (180 - (540 - (l_Z - l_S)) % 360);
	lambda_plus_delta_wgs84   = (l, d_l) => this.lambda_norm_wgs84(l + d_l);
	lambda_plus_delta_wgs84_2 = (l, d_l) => this.lambda_norm_wgs84_2(l + d_l);
	lambda_norm_wgs84   = (l) => 180 - (540 - (l)) % 360;
	lambda_norm_wgs84_2 = (l) => (l + 540) % 360 - 180;
	lambda_avg_wgs84 = (l1, l2) => this.lambda_plus_delta(l1, this.delta_lambda_wgs84(l1, l2) / 2);
	
	delta_lambda   = (l_S, l_Z) => (l_Z - l_S + M_PI * 3) % (2 * M_PI) - M_PI;
	delta_lambda_2 = (l_S, l_Z) => (M_PI - (M_PI * 3 - (l_Z - l_S)) % (2 * M_PI));
	lambda_plus_delta   = (l, d_l) => this.lambda_norm(l + d_l);
	lambda_plus_delta_2 = (l, d_l) => this.lambda_norm_2(l + d_l);
	lambda_norm   = (l) => M_PI - (M_PI * 3 - (l)) % (2 * M_PI);
	lambda_norm_2 = (l) => (l + M_PI * 3) % (2 * M_PI) - M_PI;
	lambda_avg = (l1, l2) => this.lambda_plus_delta(l1, this.delta_lambda(l1, l2) / 2);
	
	delta_lambda_epsg3827  = (l_S, l_Z) => this.r * this.delta_lambda(l_S / this.r, l_Z / this.r);

	*longitude_loop_wgs84(l_S, l_Z, n, start = 0, end = n) {
		let step = this.delta_lambda_wgs84(l_S, l_Z) / n;
		for (let i = start; i <= end; i++)
			yield this.lambda_plus_delta_wgs84(l_S, i * step);
	}
	*longitude_loop(l_S, l_Z, n, start = 0, end = n) {
		let step = this.delta_lambda(l_S, l_Z) / n;
		for (let i = start; i <= end; i++)
			yield this.lambda_plus_delta(l_S, i * step);
	}
	/**
		Rückgabe in Grad
	 */
	rho = (x, y) => (((x instanceof c_Point) ? atan2(x.c[0], x.c[1]) :
	(y == undefined) ? x : atan2(x, y))  / geo.transform.grad2rad_f + 360) % 360
	//rwk = (alpha) => (alpha + 2 * M_PI) / geo.transform.grad2rad_f  % 360;

/*	*longitude_interpolation_wgs84(l_S, l_Z, f, far_away) {
		for (let p of geo.interpolation.midway(l_S, l_Z, f, geo.nav.lambda_avg, far_away)
			 yield(p);
	/*	for (let p of geo.interpolation.midway(l_S, l_Z, f
			(a, b) => geo.nav.lambda_plus_delta(a, geo.nav.delta_lambda_wgs84(a, b) / 2.),
			(a, b) => Math.abs(geo.nav.delta_lambda_wgs84(a, b)) > delta)
			 )
			yield p;
	}*/

	/**
	 Berechnet die Distanz zwischen den Punkten A und B
	 * @param {array} A - Coordinates point A WGS84
	 * @param {array} B - Coordinates point B WGS84
	 * @return {number} Distance in km between A and B
	 */
	dist(A, B) {
		let a = geo.transform.grad2rad(A), b = geo.transform.grad2rad(B);
		return this.r *
		acos(sin(b[1]) * sin(a[1]) + cos(b[1]) * cos(a[1]) * cos(b[0] - a[0]));
	}
	/**
	 Berechnet die Distanz zwischen den Punkten A und B auf Einheitskugel dadiant!
	* @param {array} A - Coordinates point A rad
	* @param {array} B - Coordinates point B rad
	* @return {number} Distance in km between A and B
	*/
	dist_norm(a, b) {
		return acos(sin(b[1]) * sin(a[1]) + cos(b[1]) * cos(a[1]) * cos(b[0] - a[0]));
	}
	/**
	 * @param {array} B - Coordinates points B WGS84
	 @return {object} sin and cos of latitude B
	 */
	dist_d(B) {
		let b = geo.transform.grad2rad(B);
		return {sin_lon: sin(b[1]), cos_lon: cos(b[1]), b: b};
	}
	/**
	 * @param {array} A - Coordinates point A WGS84
	 * @param {dest} B - precalculates values from destinarion (return value from dist_d)
	 * @return {number} Distance in km between A and dest
	 */
	dist_2(A, dest) { // A in WGS84; dest aus dist_d!
		let a = geo.transform.grad2rad(A);
		let dist = geo.nav.r *
			acos(dest.sin_lon * sin(a[1]) + dest.cos_lon * cos(a[1]) * cos(dest.b[0] - a[0]));
		return dist;
	}

	dist_haversin(A, B) {
		const a = geo.transform.grad2rad(A), b = geo.transform.grad2rad(B);
		const d_lat = a[1] - b[1], d_lon = a[0] - b[0];
		const diskr = sin((a[1] - b[1]) / 2) ** 2 + sin((a[0] - b[0]) / 2) ** 2 * cos(a[1]) * cos(b[1]);
		return 2 * geo.nav.r * asin(Math.sqrt(diskr));
	}

	dist_haversin2(A, B) { // tbd: Genau gleiche Ergebnise wie dist !?!?!?
		const a = geo.transform.grad2rad(A), b = geo.transform.grad2rad(B);
		return 2 * this.r * asin(
		 Math.sqrt(
		 (1 - cos(a[1] - b[1]) + cos(a[1]) * cos(b[1]) * (1 - cos(a[0] - b[0])))
		 / 2
		 )
		 );
	}

	dist_pythagoras(A, B) {
		let a = geo.transform.grad2rad(A), b = geo.transform.grad2rad(B);
									 

		return geo.nav.r * Math.sqrt(
			((a[0] - b[0]) * cos((a[0] + a[1])/2)) ** 2 + (a[1] - b[1]) ** 2
		);
	}
	
	/**
	 * @param {array} v - vector for calculation
	 * @param {array} rot_M - Rotation Matrix for Vector
	 * @param {int} n - number of waypoints
	 * @return {array} ToDo
	 */
	rwK(v, rot_M, rot) {
		let rho =  geo.transform.rad2grad(
			Math.acos(
				Math.min(
					geo.vector.dot(geo.vector.mul(rot_M, v), geo.vector.N(v))
				, 1)
			)
		);
		return (rot[2] > 0) ? rho : 360 - rho;
	}
	/**
	 * @param {array} A - WGS64 starting point
	 * @param {array} B - WGS64 destination point
	 * @param {int} n - number of waypoints
	 * @return {array} ToDo
	 */
	orthodrome_kart (A, B, n = 0) {
		let a = geo.transform.polar2kart_3d(geo.transform.grad2rad(A));
		let b = geo.transform.polar2kart_3d(geo.transform.grad2rad(B));
		let dist = Math.acos(geo.vector.dot(a, b)), delta = dist / (n + 1);
		let rot = geo.vector.cross(a, b);
		let ortho = [A];
		let rot_M = geo.vector.rotation_matrix(rot, M_PI_2);
		let rho = [geo.nav.rwK(a, rot_M, rot)];
		for (let i = 1; i <= n; i++) {
			let p = geo.vector.mul(geo.vector.rotation_matrix(rot, delta * i), a);
			ortho.push(geo.transform.rad2grad(geo.transform.kart2polar_3d(p)));
			rho.push(geo.nav.rwK(p, rot_M, rot));
		}
		ortho.push(B);
		rho.push(geo.nav.rwK(b, rot_M, rot));
		return {track: ortho, rwK: rho, dist: geo.r * dist};
	}
	


} //
												
												
/* Beispiel für gleitenden Mittelwert
 q = new G.Queue((that, z) => {that.sum = z + ((typeof(that.sum) == "undefined") ? 0 : that.sum); if (that.length == 5) that.pop();},
	that => that.sum -= that.start.val)
 */

/*
class Queue { // Altlast
				constructor(push_fn, pop_fn) {
					this.push_fn = push_fn, this.pop_fn = pop_fn;
					this.start = this.end = null;
					this.length = 0;
				}
				push(o) {
					if (this.push_fn)
						this.push_fn(this, o);
					this.length++;
					if (!this.start)
						this.start = this.end = {val: o, next: null};
					else
						this.end = this.end.next = {val: o, next: null};
				}
				pop() {
					let rc = undefined;
					if (this.start) {
						if (this.pop_fn)
							this.pop_fn(this);
						rc = this.start.val, this.start = this.start.next, this.length--;
						if (!this.start)
							this.end = null;
					}
					return rc;
				}
}
*/
class c_Stack {
	constructor() {
		this.array = [];
	}
	push(data) {
		this.array.push(o);
	}
	pop() {
		return this.pop();
	}
	get length () {
		return this.array.length;
	}
}
class c_Queue {
	#length;
	constructor() {
		this.start = this.end = null;
		this.#length = 0;
	}
	push(data) {
		if (this.start == null)
			this.end = this.start = {data: data, next: null};
		else
			this.end = this.end.next = {data: data, next: null};
		this.#length++;
	}
	pop(data) {
		var res = undefined;
		if (this.start) {
			res = this.start.data;
			if (this.start != this.end)
				this.start = this.start.next;
			else
				this.start = this.end = null;
			this.#length--;
		}
		return res;
	}
	get length() {return this.#length;}
};

class c_Binheap {
	constructor(wrong) {
		this.array = [];
		this.wrong = wrong; //1. Parameter: Parent; 2. Parameter: Child; return true if wrong
	}
	push(o) {
		this.array.push(o);
		this.swim_up(this.array.length - 1);
	}
	swim_up(ind) {
		var p, tmp = this.array[ind];
		while (ind > 0 && this.wrong(this.array[p = (ind - 1) >> 1], tmp))
			this.array[ind] = this.array[p], ind = p;
		this.array[ind] = tmp;
		// tbd: return ind;
	}
	sink_down(i) {
		let tmp = this.array[i];
		while (true) {
			let child = 2 * i + 1, cmp_obj = tmp, ind = undefined;
			if (child < this.array.length && this.wrong(cmp_obj, this.array[child]))
				cmp_obj = this.array[ind = child];
			if (++child < this.array.length && this.wrong(cmp_obj, this.array[child]))
				ind = child;
			if (ind == undefined)
				break;
			this.array[i] = this.array[ind], i = ind;
		}
		this.array[i] = tmp;
		// tbd: return i;
	}
	set(ind, o) {
		console.log("set", this);
		if (this.wrong(this.array[ind], o)) // o is less
			this.array[ind] = o, this.swim_up(ind);
		else
			this.array[ind] = o, this.sink_down(ind);
		console.log("set", this);
	}
	pop() {
		return this.array.pop();
	}
	get length() {
			return this.array.length;
	}
};

class c_Extrememem {
	constructor(n, wrong) {
		this.n = n, this.heap = new c_Binheap(wrong);
	}
	push(o) {
		if (this.heap.length < this.n)
			this.heap.push(o);
		else if (this.heap.wrong(o, this.heap.array[0]))
			this.heap.set(0, o);
	}
}

class c_Priorityqueue {
	constructor(wrong) {
		this.heap = new c_Binheap(wrong);
	}
	push(o) {
		this.heap.push(o);
	}
	pop() {
		var r = undefined;
		if (this.heap.length > 1)
			r = this.heap.array[0], this.heap.set(0, this.heap.pop());
		else if (this.heap.length == 1) // Letztes / einziges Element
			r = this.heap.pop();
		return r;
	}
	get length() {
		return this.heap.length
	}
};

class c_Dijkstra_Queue {
	constructor() {
		this.array = [];
		//this.map = new Map();
		this.map = {};
	}
	push(o) {
		const key = o.key, data = o.data;
		var ind;
		if (data in this.map)
			ind = this.map[data], this.array[ind] = {key: key, data: data};
		else
			ind = this.array.length, this.array.push({key: key, data: data});
		let tmp = this.array[ind], p;
		while (ind > 0 && this.array[p = (ind - 1) >> 1].key > tmp.key)
			this.array[ind] = this.array[p], this.map[this.array[ind].data] = ind, ind = p;
		this.array[ind] = tmp, this.map[tmp.data] = ind;
	}
	pop() {
		var r = undefined;
			if (this.array.length > 1) {
				r = this.array[0], delete this.map[r.data];
				let tmp = this.array.pop(), i = 0;
				while (true) {
					let left = 2 * i + 1, right = left + 1, cmp_obj = tmp, ind = undefined;
					if (left < this.array.length && this.array[left].key < cmp_obj.key)
						ind = left, cmp_obj = this.array[ind];
					if (right < this.array.length && this.array[right].key < cmp_obj.key)
						ind = right;
					if (ind == undefined)
						break;
					this.array[i] = this.array[ind], this.map[this.array[i].data] = i, i = ind;
				}
				this.array[i] = tmp, this.map[tmp.data] = i;
			}
			else if (this.array.length == 1) // Letztes / einziges Element
				r = this.array.pop(), delete this.map[r.data];
		 return r;
	}
	get length() {
		return this.array.length
	}
};



class c_Graph {
	/**
	 *
	 */
	constructor(nodes = {}, edges = {}) {
		this.edges = edges;
		this.nodes = nodes;
	}
	/**
	 s: Startknoten
	 z: Zielknoten
	 q: Warteschlange
	 f_w: Gewichtsfunktion
	 f_q: Heuristikfunkktion Distanz zum Ziel
	 *
	 */
	route (s, z, q, f_w, f_q) {
		console.time("route")
		const visited = {}, path = [];
		q.push({key: 0, data: s}), visited[s] = {pred: null, dist: 0};
		for (let nd; (nd = q.pop()) != undefined && (nd = nd.data) != z;) {
			for (let a in this.edges[nd]) {
				const d2 = visited[nd].dist + f_w(nd, a);
				if (!(a in visited) || visited[a].dist > d2) {
					visited[a] = {pred: nd, dist: d2};
					q.push({key: d2 + f_q(a), data: a})
				}
			}
		}
		if (visited[z])
			for (let nd = z; nd!= null; nd = visited[nd].pred)
				path.push(nd);
		path.reverse();
		console.timeEnd("route")
		return path;
	};
	/**
	s: Startknoten
	z: Zielknoten
	*/
	route_breadth_search(s, z) {
		return this.route(s, z, new c_Queue,
			() => 1, // Distanz zwischen Knoten
			() => 0);// Distanz zum Ziel
	}
	/**
	s: Startknoten
	z: Zielknoten
	*/
	route_dijkstra(s, z) {
		return this.route(s, z, new c_Dijkstra_Queue,
			(a, b) => this.edges[a][b].w,// Distanz zwischen Knoten
			() => 0);// Distanz zum Ziel
	}
	/**
	s: Startknoten
	z: Zielknoten
	*/
	route_a_star(s, z) {
		const dist_d = geo.nav.dist_d( this.nodes[z]);
		return this.route(s, z, new c_Dijkstra_Queue,
			(a, b) => this.edges[a][b].w,// Distanz zwischen Knoten
			(a) => geo.nav.dist_2(this.nodes[a], dist_d));// Distanz zum Ziel
	}
/*	route(id_from, id_to, f_w, f_q) {
		var u, res = {path: [], w: 0}, h = new geo.binheap(), nodes = {};
		if (!this.n[id_from] || !this.n[id_to])
			return res;
		nodes[id_from] =  {pred: null, w: 0};
		h.push({id: id_from, w: 0});
		while (!h.empty() && (u = h.pop().id) != id_to) {
			for (const v in this.n[u].a) {
				let w_v = f_w(nodes[u].w, u, v);
				if (nodes[v] == undefined || nodes[v].w > w_v) {
					nodes[v] = {pred: u, w: w_v};
					h.push({id: v, w: f_q(w_v, v, id_to)});
				}
			}
		}
		if (u == id_to) {
			res.path.push({u: u, e: null});
			while (nodes[u].pred != null) {
				res.path.push({u: nodes[u].pred, e: this.n[nodes[u].pred].a[u]});
				res.w += this.e[this.n[nodes[u].pred].a[u]].w
				u = nodes[u].pred;
			}
		}
		res.path.reverse()
		return res;
	}*/
	

	fromOpenStreetMap(data) {
		this.edges = {};
		this.nodes = {};

		data.elements.filter(e => e.type == "node").forEach(node =>{ // tbd: for ... of
			this.nodes[node.id] = [node.lon, node.lat];
			this.edges[node.id] = {}
		});
		data.elements.filter(e => e.type == "way").forEach(way => {
			for (let i = 0; i < way.nodes.length - 1; i++) {
				const s = way.nodes[i], z = way.nodes[i + 1];
				if (this.nodes[s] && this.nodes[z]) {
					this.edges[s][z] = this.edges[z][s] =
					{w: geo.nav.dist_haversin(this.nodes[s], this.nodes[z])};
				}
			}
		});
	}

	to_OL_VectorLayer(vec, style_fn) {
		for (let s in this.edges) {
			for (let z in this.edges[s]) {
				this.edges[s][z].disp = this.edges[z][s].disp = new ol.Feature({
					geometry: new ol.geom.LineString([this.nodes[s], this.nodes[z]])});
				this.edges[s][z].disp.getGeometry().transform('EPSG:4326', 'EPSG:3857');
				if (style_fn)
					this.edges[s][z].disp.setStyle(style_fn(s, z));
				vec.addFeature(this.edges[s][z].disp);
			}
		}
	}
};



class kd_tree {
	/**
	 data: Zu indizierendes Feld
	 fn_xy: Funktion die zu einem Datenelemnt aus data das Array [x, y] liefert
	 fn_dist: Distanz zwischen zwei Punken
		Default: Pythagoras ** 2
	 */
	constructor (data, fn_xy = e => e,
		fn_dist = (a, b) => (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) {
		this.fn_xy = fn_xy, this.fn_dist = fn_dist;
		this.keys = (data instanceof Array) ? data :  Object.values(data);
		let stack = [{d: 0, start: 0, l: this.keys.length}], o;
		while (o = stack.pop()) {
			if (o.l <= 1)
				continue;
			let m = o.l >> 1;
			this.select(o.d % 2, o.start + m, o.start, o.start + o.l - 1);
			stack.push({d: o.d + 1, start: o.start + m + 1, l: o.l - m - 1});
			stack.push({d: o.d + 1, start: o.start, l: m});
		}
	}
	mem_nn() {
		class c_Kdtree_mem_nn {
			#dist;
			constructor() {
				this.#dist = +Infinity, this.ind = undefined;
			}
			push(ind, dist) {
				if (dist < this.#dist)
					this.#dist = dist, this.ind = ind;
			}
			get dist () {
				return this.#dist;
			}
		};
		return new c_Kdtree_mem_nn;
	}
	/**
	 p: Punkt zu dem NN gesucht wird
	 arr: Startindex und Länge als Array
	 Alle anderen per Defaultwert!
	 */
	#NN_recursive(p,  mem, s = 0, l = this.keys.length, dim = 0,
			bbox = [-Infinity, -Infinity, +Infinity, +Infinity]) {
		if (l < 1 || this.dist_bbox_p(bbox, p) > mem.dist)
			return null;
		const m = s + (l >> 1), pos = this.fn_xy(this.keys[m]);
		mem.push(m, this.fn_dist(pos, p));
		var bbox_1 = bbox.slice(), bbox_2 = bbox.slice();
		if (dim % 2 == 0)// aktuell x
			bbox_1[2] = bbox_2[0] = pos[0]; // right half , left half
		else // aktuell y
			bbox_1[3] = bbox_2[1] = pos[1]; // upper half, left half
		if (p[dim % 2] < pos[dim % 2])
			this.#NN_recursive(p, mem, s,	m - s, dim + 1, bbox_1),
			this.#NN_recursive(p, mem, m + 1, s + l - 1 - m, dim + 1, bbox_2);
		else
			this.#NN_recursive(p, mem, m + 1, s + l - 1 - m, dim + 1, bbox_2),
			this.#NN_recursive(p, mem, s,	m - s, dim + 1, bbox_1);
		return mem;
	}
	#NN(p) {
		const stack = [{start: 0, length: this.keys.length, dim: 0,
			bbox: [-Infinity, -Infinity, +Infinity, +Infinity]}];
		let merk = {dist: Infinity, ind: undefined}, e;
		while(e = stack.pop()) {
			if (e.length < 1 || this.dist_bbox_p(e.bbox, p) > merk.dist)
				continue;
			let m = e.start + (e.length >> 1), pos = this.fn_xy(this.keys[m]), d2;
			if ((d2 = this.fn_dist(pos, p)) < merk.dist)
				merk = {dist: d2, ind: m};
			var bbox_1 = e.bbox.slice(), bbox_2 = e.bbox.slice();
			if (e.dim % 2 == 0)// aktuell x
				bbox_1[2] = bbox_2[0] = pos[0]; // right half , left half
			else // aktuell y
				bbox_1[3] = bbox_2[1] = pos[1]; // upper half, left half
			stack.push({start: m + 1, length: e.start + e.length - 1 - m, dim: e.dim + 1, bbox: bbox_2});
			stack.push({start: e.start,	length: m - e.start, dim: e.dim + 1, bbox: bbox_1});
		}
		return merk;
	}

	NN_k(p) {
		return this.fn_xy(this.NN_e(p));
	}
	NN_e(p) {
		return this.keys[this.#NN_recursive(p, this.mem_nn()).ind];
	}

	dist_bbox_p(bb, p) {  // Only 2D !!!
		return (bb[0] <= p[0] && bb[1] <= p[1] && bb[2] >= p[0] && bb[3] >= p[1]) ? 0 :// p innerhalb bb
		this.fn_dist([
			(bb[2] < p[0]) ?  bb[2] : ((bb[0] > p[0]) ? bb[0] : p[0]),
			(bb[3] < p[1]) ?  bb[3] : ((bb[1] > p[1]) ? bb[1] : p[1])
		], p);
	}

	
	select(dim, k, li, re) { // Analog "Numerical Recipes in C"
		const val = (idx) => this.fn_xy(this.keys[idx])[dim];
		for (;;) {
			if (re <= li + 1) { // Active partition contains 1 or 2elements.
				if (re == li + 1 && val(re) < val(li)) // Case of 2elements.
					this.swap(li, re);
				break;
			} else {
				let mid=(li + re) >> 1;
				this.swap(mid, li + 1);
				if (val(li) > val(re))
					this.swap(li, re);
				if (val(li + 1) > val(re))
					this.swap(li + 1, re);
				if (val(li) > val(li + 1))
					this.swap(li, li + 1);
				let i=li + 1, j=re;
				for (;;) {
					while (val(++i) < val(li + 1));
					while (val(--j) > val(li + 1));
					if (j < i)
						break;
					this.swap(i,j);
				}
				this.swap(li + 1, j);
				if (j >= k)
					re = j-1;
				if (j <= k)
					li = i;
			}
		}
	}

	swap(a, b) {
		[this.keys[a], this.keys[b]] = [this.keys[b], this.keys[a]];
	}
	print(start = 0, l = this.keys.length, level = 1) {
		if (l) {
			const m = start + (l >> 1) ;
			console.log("".padEnd(level," "), `${start}...${l} m=${m}:`, this.fn_xy(this.keys[m]));
			this.print(start, m - start, level + 1);//
			this.print(m + 1, l - (m - start) -1 ,level + 1);
		}
	}
	
	tree(start = 0, l = this.keys.length) {
		const m = start + (l >> 1) ;
		return (!l) ? null : {
			val: this.fn_xy(this.keys[m]),
			l: this.tree(start, m - start),
			r: this.tree(m + 1, l - (m - start) -1)
		}
	}

	/*	select(arr, k, li, re) { // Analog "Numerical Recipes in C"
	 const f_xy = (idx) => this.fn_xy(this.keys[idx])[k];
	 var i, j, mid;
	 for (;;) {
	 if (re <= li + 1) { // Active partition contains 1 or 2elements.
	 //if (re == li + 1 && arr[re] < arr[li]) { // Case of 2elements.
	 if (re == li + 1 && arr[re] < arr[li]) { // Case of 2elements.
	 this.swap(li, re);
	 }
	 return;
	 } else {
	 mid=(li + re) >> 1;
	 this.swap(mid, li + 1);
	 if (arr[li] > arr[re])
	 this.swap(li, re);
	 if (arr[li + 1] > arr[re])
	 this.swap(li + 1, re);
	 if (arr[li] > arr[li + 1]) {
	 this.swap(li, li + 1);
	 }
	 i=li + 1, j=re;
	 for (;;) { // tbd: wile (j >= i) und if--break raus!
	 do i++; while (arr[i] < arr[li + 1]); // tbd: while (arr[++i] < arr[li + 1])
	 do j--; while (arr[j] > arr[li + 1]); // tbd: while (arr[--j] > arr[li + 1])
	 if (j < i) break;
	 this.swap(i,j);
	 }
	 this.swap(li + 1, j);
	 if (j >= k) re = j-1;
	 if (j <= k) li = i;
	 }
	 }
	 }
	


	
	tikz_tree(d, start, length) {
		var s = "";
		if (typeof d !== 'number' || d === 0) {
			d = 0, start = 0, length = this.keys.length;
			s += "\\documentclass{standalone}\n\\usepackage{tikz}\\usepackage{tikz-qtree}\n\\begin{document}\n\\begin{tikzpicture}[every node/.style = {anchor=north,align=center,shape=rectangle, rounded corners,	draw, align=center,top color=white, bottom color=blue!20},blank/.style={draw=none, top color=white!0, bottom color=white!0}]\n\\Tree";
		}
		let half = length >> 1, m = start + half;
		if (length > 1) {
			s += "     ".substring(0,d) + `[.\\node{${this.indK[0][m]}/${this.indK[1][m]}};\n`;
			s += this.tikz_tree(d + 1, start, half);
			if (length > 2)
				s += this.tikz_tree(d + 1, m + 1, length - half - 1);
			else
				s += "     ".substring(0,d+1) + "\\edge[blank];\\node[blank]{};\n";
			s +="     ".substring(0,d) + "]\n"
		}
		else if (length == 1)
			s += "     ".substring(0,d) + `\\node{${this.indK[0][start]}/${this.indK[1][start]}};\n`;
		if (!d)
			s += "\\end{tikzpicture}\n\\end{document}";
		return s;
	}
	
	tikz_plane(d, start, length, bbox) {
		var s = "";
		if (typeof d !== 'number' || d === 0) {
			d = 0, start = 0, length = this.keys.length, bbox = [this.bbox[0] - 1,this.bbox[1] - 1,this.bbox[2] + 1,this.bbox[3] + 1];
			s = "\\documentclass{standalone}\n";
			s += "\\usepackage{tikz}\n";
			s += "\\begin{document}\n\\begin{tikzpicture}\n";
			s += `\\draw[very thick] (${bbox[0]},${bbox[1]})--(${bbox[2]},${bbox[1]})--(${bbox[2]},${bbox[3]}) -- (${bbox[0]},${bbox[3]})--cycle;\n`;
		}
		if (!length)
			return "";
		let half = length >> 1, m = start + half, p = [this.indK[0][m], this.indK[1][m]];
		console.log("p=", p);
		if (d % 2 == 0)
			s += `\\draw (${p[0]},${bbox[1]})--(${p[0]},${bbox[3]});\n`;
		else
			s += `\\draw (${bbox[0]},${p[1]})--(${bbox[2]},${p[1]});\n`;
		s += `\\filldraw[red] (${p[0]},${p[1]}) circle (2pt) node[black,above right ] {\\tiny ${p[0]}/${p[1]}};\n`;
		if (d % 2 == 0)
			s += this.tikz_plane(d + 1, start, half, [bbox[0],bbox[1],p[0],bbox[3]]),
			s += this.tikz_plane(d + 1, m + 1, length - half -1, [p[0],bbox[1],bbox[2],bbox[3]]);
		else
			s += this.tikz_plane(d + 1, start, half, [bbox[0], bbox[1],bbox[2], p[1]]),
			s += this.tikz_plane(d + 1, m + 1, length - half -1, [bbox[0],p[1], bbox[2], bbox[3]]);
		
		if (!d)
			s += "\\end{tikzpicture}\n\\end{document}";
		return s;
	}
	 
	 
	 /*	minK(dim, d, start, l) {
	 if (typeof d === 'undefined')
	 d = 0, start = 0, l = this.n;
	 if (!l)
	 return undefined;
	 // console.log("\t\t\t\t\t\t\t\t\t\t".substring(0, d), d, "==> ", start, l);
	 let m = l >> 1, erg = this.indK[dim][start + m];
	 var erg1;
	 if ((erg1 = this.minK(dim, d + 1, start, m)) < erg)
	 erg = erg1;
	 if (d % this.dim != dim && (erg1 = this.minK(dim, d + 1, start + m + 1, l - m -1)) < erg)
	 erg = erg1;
	 return erg;
	 }
	 
	 maxK(dim, d, start, l) {
	 if (typeof d === 'undefined')
	 d = 0, start = 0, l = this.n;
	 if (!l)
	 return undefined;
	 let m = l >> 1, erg = this.indK[dim][start + m];
	 //console.log("\t\t\t\t\t\t\t\t\t\t".substring(0, d), d, "==> ", start, l);
	 var erg1;
	 if ((erg1 = this.maxK(dim, d + 1, start + m + 1, l - m - 1)) > erg)
	 erg = erg1;
	 if (d % this.dim != dim) {
	 if ((erg1 = this.maxK(dim, d + 1, start, m)) > erg)
	 erg = erg1;
	 }
	 return erg;
	 }*/
	 
};



class c_Orthodrome{
	#alpha_A; #A; #B;#cos_d;#d;#NN;#tan_phi_NN;#sin_phi_NN;#sin_phi_A;#cos_phi_A;#cos_alpha_A;#sin_alpha_A;
	constructor(A, B) {// ToDo: B als Objekt {alpha: ..., d: ...} übergeben und B über P(d) berechnen
		this.#A = geo.transform.grad2rad(A);
		this.#sin_phi_A = sin(this.#A[1]);
		this.#cos_phi_A = cos(this.#A[1]);
		if (!B.alpha) { // B ist Punkt
			this.#B = geo.transform.grad2rad(B);
			console.log(this.A, this.B);
			this.#cos_d = sin(this.#B[1]) * sin(this.#A[1]) + cos(this.#B[1]) * cos(this.#A[1]) * cos(this.#B[0] - this.#A[0]);
			this.#d = acos(this.#cos_d);
			//this.#alpha_A = asin(cos(this.#B[1]) * sin(geo.nav.delta_lambda(this.#A[0], this.#B[0])) / Math.sqrt(1 - this.#cos_d * this.#cos_d)); // Sinussatz, händisch Seite 3
			//console.log(this.#alpha_A);
			this.#alpha_A = atan2( // Bronshtein 3.208 p. 176
			 cos(this.#A[1]) * cos(this.#B[1]) * sin(geo.nav.delta_lambda(this.#A[0], this.#B[0])),
			 sin(this.#B[1]) - sin(this.#A[1]) * this.#cos_d
			 );
		}
		else { // B ist Kurs / Dist
			this.#d = B.d / geo.nav.r;
			this.#cos_d = cos(this.#d);
			this.#alpha_A = B.alpha;
			this.#B = geo.transform.grad2rad(this.P(b.d));
		}
		this.#NN = [0.0, acos(sin(abs(this.#alpha_A))* cos(this.#A[1]))]; // 3.204a ergibt sich aus Sinussatz, Händisch Seite 4
		this.#NN[0] = this.#A[0] + signum(this.#alpha_A) * abs(acos(tan(this.#A[1]) / tan(this.#NN[1]))); // 3.204b aus 3.202e gepimpt
		this.#tan_phi_NN = tan(this.#NN[1]);
		this.#sin_phi_NN = sin(this.#NN[1]);
		this.#sin_alpha_A = sin(this.#alpha_A);
		this.#cos_alpha_A = cos(this.#alpha_A);
	}
	get alpha_A() {
		return geo.nav.rwk(this.#alpha_A);
	}
	get A() {
		return geo.transform.rad2grad(this.#A);
	}
	get B() {
		return geo.transform.rad2grad(this.#B);
	}
	get d() {
		return this.#d * geo.nav.r;
	}
	/**
	d: integer: Anzahl Schritte; float: Max. Distanz
	*/
	track(d) {
		if (!Number.isInteger(d))
			d = Math.ceil(this.d / d);
		const t = [this.A], step = this.d / d;
		for (let i = 1; i < d; i++)
			t.push(this.P(i * step));
		t.push(this.B);
		return t;
	}
	/**
	lambda: rad
	*/
	#lat(lambda) {
		return atan(this.#tan_phi_NN * cos(geo.nav.delta_lambda(lambda,this.#NN[0])));
	}
	lat(lon) {
		return this.#lat(lon  * geo.transform.grad2rad_f) / geo.transform.grad2rad_f;
	}
	/**
	lon: Länge an der das Bearing berechnet wird
	 tbd: Was ist bei Westkursen???
	*/
	rwk(lon) {//  ρ=arccos(sin(λP −λNN)sinφNN)
		const rho = acos(sin(geo.nav.delta_lambda(lon *  geo.transform.grad2rad_f, this.#NN[0])) * this.#sin_phi_NN);
		//return M_PI - rho;
		console.log(rho/geo.transform.grad2rad_f);
		return rho/geo.transform.grad2rad_f;//geo.nav.rwk(M_PI - rho);
	}
	/**
	x: Distanz des return-Punktes zu A
	 return: Punkt auf
	*/
	P(x, r = geo.nav.r) {
		const a = acos(this.#sin_phi_A * cos(x / r) + this.#cos_phi_A * sin(x / r) * this.#cos_alpha_A);
		return geo.transform.rad2grad([geo.nav.lambda_plus_delta(this.#A[0], asin(this.#sin_alpha_A * sin(x / r) / sin(a))), M_PI_2 - a]);
	}
};
			

class c_Loxodrome {
			#A; #line; #alpha; #cos_alpha; #phi_A; #phi_B;
			/**
			 ToDo: Statt Zielpunkt B
			 */
			constructor(A, B, C) {
				var B;
				this.#A = geo.transform.EPSG_4326_2_EPSG_3857(A);
				this.#phi_A = A[1], this.#phi_B = B[1];
				if (arguments.length == 2) {
					B = geo.transform.EPSG_4326_2_EPSG_3857(B);
					const v = new c_Vector([
						geo.nav.delta_lambda_epsg3827(this.#A[0], B[0]),
						B[1] - this.#A[1]
					]);
					this.#line = new c_Line(new c_Point(this.#A), v);
					this.#alpha = Math.PI / 2 - this.#line.alpha;
					this.#cos_alpha = cos(this.#alpha);
				}
				else { // A ist Array, B ist Winkel rho, C ist dist d
	
				}
			}
			get rho() {
					return geo.nav.rho(this.#alpha);
			}
			get d() {
				return abs((this.#phi_A - this.#phi_B) / this.#cos_alpha) * geo.transform.grad2rad_f  * geo.nav.r;
			}
			lat(lon) {
				lon = geo.transform.EPSG_4326_2_EPSG_3857([lon, 0])[0];
				//const d_x = geo.nav.delta_lambda_epsg3827(this.#A[0], lon)
				const d_x = lon - this.#A[0];
				const y = this.#line.y(this.#A[0] + d_x);
				return geo.transform.EPSG_3857_2_EPSG_4326([0, y])[1];
			}
			lon(lat) {
				const y = geo.transform.EPSG_4326_2_EPSG_3857([0, lat])[1];
				const x = this.#line.x(y);
				return geo.transform.EPSG_3857_2_EPSG_4326([x, 0])[0]; // ggf normieren
			}
			P(s) {
				const phi_B = s * this.#cos_alpha / (geo.transform.grad2rad_f  * geo.nav.r) + this.#phi_A;
				return [this.lon(phi_B), phi_B];
			}

};
class c_N_S_course{
	constructor(A, B) {
		
	}
	/**
	 d: integer: Anzahl Schritte; float: Max. Distanz
	 */
	/**
	 
	 */
	lat(lon) {
		
	}
	track(d) {
		
	}
	/**
	 lon: Länge an der das Bearing berechnet wird
	 */
	bearing(lon) {
		
	}
	/**
	 d: Disanz des return-Punktes zu A
	 */
	P(d) {
		
	}
};
												

class c_Raster_tile  {
	/**
	 * calculates earth-pixel from EPSG:3857
	 *
	 * @param {array} c - web-mercator-coordinates
	 * @param {int} z - zoom level
	 * @return {array} - Pixel over all tiles (over world)
	 */
	px_world(c, z) {
		const n_2 = 1 << (z + 7); // Halbe Pixelzahl
		return [
			((c[0] / geo.nav.M_EPSG_3857_MAX + 1) * n_2) >> 0,
			((1 - c[1] / geo.nav.M_EPSG_3857_MAX) * n_2) >> 0,
		];
	}

	px_tile(c, z) {
		return [
				[c[0] >> 8, c[0] % 256],
				[c[1] >> 8, c[1] % 256]
		];
	}
	/** constructor
	 */

	constructor() {
	}
				
	/**
	* calculatesEPSG:3857 for Pixel depending of a zoomlevel
	*
	* @param {double[]} {px}  Tilepixels if type is [[x_nr, x_px], [[y_nr, y_px]]] else Worldpixel [x, y]
	* @param {int} {z}  Zoomlevel
	* @return {array} [x, y]
	*/
	epsg_3857(px, z) {
		const n_2 = 1 << (z+7);
		const p_welt = (typeof(px[0]) == 'number') ? px :
			[256 * px[0][0] + px[0][1], 256 * px[1][0] + px[1][1]];
		return [
			(p_welt[0] / n_2 - 1) * geo.nav.M_EPSG_3857_MAX,
			(1 - p_welt[1] / n_2) * geo.nav.M_EPSG_3857_MAX
		];
	}
	

	xyz_valid (x, y, z) {
		return (z < 0 || z >= 30 || x < 0 || x >= (1 << z)
			|| y < 0 || y >= (1 << z)) ? false : true;
	}


	/**
	 * calculates coordinate-delta of one tile depending of zoomleve.
	 *
	 * @param {int} {z}  Zoomlevel
	 * @return {double} delta
	 */
	delta (z) {
		return 2 * geo.M_EPSG_3857_MAX / (1 << z);
	}
	
	/**
	 * calculates EPSG:3857-coordinate of top-left corner of a tile
	 *
	 * @param {int} {x}  tile-x
	 * @param {int} {x}  tile-y
	 * @param {int} {z}  Zoomlevel
	 * @return {double[]} top-left coordinate
	 */

}
		
		
class c_Interpolation {
	*midway(x1, x2, f, mid, far_away, l = 0) {
		//console.log("===>", x1, x2, l);
		const f_x_1 = f(x1), f_x_2 = f(x2);
		if (l == 0)
			yield f_x_1;
		if (far_away(f_x_1, f_x_2)) {
			let m = mid(x1, x2);
			yield* this.midway(x1, m, f, mid, far_away, l + 1);
			yield f(m);
			yield* this.midway(m, x2, f, mid, far_away, l + 1);
		}
		if (l == 0)
			yield f_x_2;
	}
	
			dim_1(p1, p2) {
				class dim_1 {
				constructor(p1, p2) {
					this.a = (p2[1] - p1[1]) / (p2[0] - p1[0]);
					this.b = p1[1] - this.a * p1[0];
				}
					y(x) {
						return this.a * x + this.b;
					}
			};
			return new dim_1(p1, p2)}

			/*	*midway(x1, x2, f, mid, far_away, l = 0) {
					//console.log(x1, x2,  Math.abs(x1, x2), far_away(x1, x2))
					if (l == 0)
						yield f(x1);// way.push(f(x1));
					
					if (far_away(x1, x2)) {
						let m = mid(x1, x2);
						yield* this.midway(x1, m, f, mid, far_away, l + 1);
						yield(f(m)); //way.push(f(m));
						yield* this.midway(m, x2, f, mid, far_away, l + 1);
					}
					if (l == 0) {
						yield(f(x2));//way.push(f(x2));
						
						//			return way;
					}
				}*/
};
			
												
												



class Orthodrome_trig {
	constructor(A, B) { // A, b in WGS84
		this.geo = new Geo();
		this.A = this.geo.transform.grad2rad(A);
		this.B = this.geo.transform.grad2rad(B);
		this.a = this.geo.transform.grad2rad(90 - B[1]);
		this.b = this.geo.transform.grad2rad(90 - A[1]);
		this.delta_ = Math.abs(this.A[0] - this.B[0]);
		this.d_ = Math.acos(Math.cos(this.a) * Math.cos(this.b) +
			Math.sin(this.a) * Math.sin(this.b) * Math.cos(this.delta_));
		let dummy = this.delta_ / Math.sin(this.d_);
		this.alpha_= dummy * Math.cos(this.B[1]);
		this.beta_ = dummy * Math.cos(this.A[1]);
	}
	get d() {return this.d_;}
	get alpha() {
		return this.alpha_ / this.geo.grad2rad_f;
	}
	get beta() {return this.beta_ / this.geo.grad2rad_f;}
	P(dist) {}

};

class Geo {
	constructor() {
	}
	get raster_tile() {return new c_Raster_tile();}
	get coordinates() {return new Coordinates();}
	get transform() {return new c_Transform();}
	get nav() {return new c_Nav()}
	get interpolation() {return new c_Interpolation;}
	get version() {return 0.1;}
};
var geo = new Geo();
	if (typeof module != "undefined") {// node.js
			//	module.exports = geo = new Geo;
			geo = new Geo(); // tbd zeiel löschen
			module.exports = geo;
			}
//else // browser
//	geo = new Geo;



 
 
 
 /**
 *
 *
 * @param various x1 Startpunkt
 * @param various x2 Startpunkt
 * @param function f Berechnung eines Punktes
 * @param function mid Berechnung der Mitte
 * @param function far_away Berechnung des Max Distanz
 * @param {array} way Punktfolge
 */
/*midway(x1, x2, f, mid, far_away, way = [], l = 0) {
				if (l == 0)
					way.push(f(x1));
				if (far_away(x1, x2)) {
					let m = mid(x1, x2);
					this.midway(x1, m, f, mid, far_away, way, l + 1);
					way.push(f(m));
					this.midway(m, x2, f, mid, far_away, way, l + 1);
				}
				if (l == 0) {
					way.push(f(x2));
					return way;
				}
			}
												circle(P1, P2, P3) {
				let c = geo.vector.intersect(
											 geo.vector.mul(geo.vector.add(P1, P2), .5),
											 [P1[1] - P2[1], P2[0] - P1[0]],
											 geo.vector.mul(geo.vector.add(P2, P3), .5),
											 [P2[1] - P3[1], P3[0] - P2[0]]
											 );
				return {c: c, r: geo.vector.abs(geo.vector.sub(P1, c))};
			}
												/*catmull_rom : {
												 A(tau) {
												 return [
												 [0, 1, 0, 0],
												 [-tau, 0, tau, 0],
												 [2 * tau, tau - 3, 3 - 2 * tau, -tau],
												 [-tau,2 - tau, tau - 2, tau]
												 ];
												 },
												 P(t, A) {
												 return geo.vector.mul([1, t, t ** 2, t ** 3], A)[0];
												 },
												 extend(x) {
												 x.unshift(geo.vector.sub(geo.vector.mul(x[0], 2), x[1]));
												 x.push(geo.vector.sub(geo.vector.mul(x[x.length - 1], 2), x[x.length - 2]));
												 },
												 unextend(x) {
												 x.shift();
												 },
												 circle(x) {
												 x.push(x[0], x[1], x[2]);
												 },
												 uncircle(x) {
												 x.pop(), x.pop(), x.pop();
												 },
												 spline(x, tau, delta) {
												 const spline = [];
												 var last;
												 for (let i = 1; i < x.length - 2; i++) {
												 const A = geo.vector.mul(this.A(tau), [x[i-1], x[i], x[i+1], x[i+2]]);
												 const dummy = geo.interpolation.midway(0, 1,
												 t => this.P(t, A),
												 (a, b) => (a + b) / 2,
												 (a, b) => delta(this.P(a, A), this.P(b, A))
												 );
												 last = dummy.pop();
												 spline.push(...dummy);
												 }
												 spline.push(last);
												 return spline;
												 }
												 
												 }
												 }
												 binheap = BinHeap
												 graph = {
												 Graph() { return Graph;},
												 Edge(a) {return new Edge(a);},
												 Node(a) {return new Node(a);},
												 ortsplan () {
												
												}
 
 
 */
