//import DollarRecognizer from './oneDoller.js';
// jsClipper uses X/Y instead of x/y...
function toClipperCoordinates(polygon){
	var clone = [];
	for(var i=0; i<polygon.length; i++){
		clone.push({
			X: polygon[i].x,
			Y: polygon[i].y
		});
	}
	
	return clone;
};

function toNestCoordinates(polygon, scale){
	var clone = [];
	for(var i=0; i<polygon.length; i++){
		clone.push({
			x: polygon[i].X/scale,
			y: polygon[i].Y/scale
		});
	}
	
	return clone;
};

function rotatePolygon(polygon, degrees){
	var rotated = [];
	var angle = degrees * Math.PI / 180;
	for(var i=0; i<polygon.length; i++){
		var x = polygon[i].x;
		var y = polygon[i].y;
		var x1 = x*Math.cos(angle)-y*Math.sin(angle);
		var y1 = x*Math.sin(angle)+y*Math.cos(angle);
						
		rotated.push({x:x1, y:y1});
	}
	
	if(polygon.children && polygon.children.length > 0){
		rotated.children = [];
		for(var j=0; j<polygon.children.length; j++){
			rotated.children.push(rotatePolygon(polygon.children[j], degrees));
		}
	}
	
	return rotated;
};

///////  recognizer

function Point(x, y) // constructor
{
	this.X = x;
	this.Y = y;
}

function Rectangle(x, y, width, height) // constructor
{
	this.X = x;
	this.Y = y;
	this.Width = width;
	this.Height = height;
}

function Compare(ObjectX, ObjectY)
{
	
	const NumPoints = 64;
	const SquareSize = 250.0;
	const Origin = new Point(0,0);
	const Diagonal = Math.sqrt(SquareSize * SquareSize + SquareSize * SquareSize);
	const HalfDiagonal = 0.5 * Diagonal;
	const AngleRange = Deg2Rad(45.0);
	const AnglePrecision = Deg2Rad(2.0);
	const Phi = 0.5 * (-1.0 + Math.sqrt(5.0)); // Golden Ratio
	
	

	function Unistroke(name, points) // constructor
	{
		this.Name = name;
		this.Points = Resample(points, NumPoints);
		var radians = IndicativeAngle(this.Points);
		this.Points = RotateBy(this.Points, -radians);
		this.Points = ScaleTo(this.Points, SquareSize);
		this.Points = TranslateTo(this.Points, Origin);
		this.Vector = Vectorize(this.Points); // for Protractor
	}

	let points=[]
	for(let k = 0; k<ObjectX[0].length ; k++)
	{
		points.push(new Point(ObjectX[0][k].x,ObjectX[0][k].y));
	}
	let UniX = new Unistroke("X", points); 


	points=[]
	for(let k = 0; k<ObjectY[0].length ; k++)
	{
		points.push(new Point(ObjectY[0][K].x,ObjectY[0][k].y));
	}
	let UniY = new Unistroke("Y", points); 



	let d = OptimalCosineDistance(UniX.Vector, UniY.Vector);
	return (1.0 - d)
}

function PlacementWorker(binPolygon, paths, ids, rotations, config, nfpCache){
	this.binPolygon = binPolygon;
	this.paths = paths;
	this.ids = ids;
	this.rotations = rotations;
	this.config = config;
	this.nfpCache = nfpCache || {};

	// return a placement for the paths/rotations given
	// happens inside a webworker
	this.placePaths = function(paths){

		//document.write("hello");
		//console.log("heellloo");
		var self = global.env.self;

		if(!self.binPolygon){
			return null;
		}		
		
		var i, j, k, m, n, path;
		
		// rotate paths by given rotation
		var rotated = [];
		for(i=0; i<paths.length; i++){
			var r = rotatePolygon(paths[i], paths[i].rotation);
			r.rotation = paths[i].rotation;
			r.source = paths[i].source;
			r.id = paths[i].id;
			rotated.push(r);
		}
		
		paths = rotated;
		
		var allplacements = [];
		var fitness = 0;
		var binarea = Math.abs(GeometryUtil.polygonArea(self.binPolygon));
		var key, nfp;
		var SimCost = 0
		while(paths.length > 0){
			
			var placed = [];
			var placements = [];
			fitness += 1; // add 1 for each new bin opened (lower fitness is better)

			for(i=0; i<paths.length; i++){
				path = paths[i];
				
				// inner NFP
				key = JSON.stringify({A:-1,B:path.id,inside:true,Arotation:0,Brotation:path.rotation});
				var binNfp = self.nfpCache[key];
				
				// part unplaceable, skip
				if(!binNfp || binNfp.length == 0){
					continue;
				}
				
				// ensure all necessary NFPs exist
				var error = false;
				for(j=0; j<placed.length; j++){			
					key = JSON.stringify({A:placed[j].id,B:path.id,inside:false,Arotation:placed[j].rotation,Brotation:path.rotation});
					nfp = self.nfpCache[key];
										
					if(!nfp){
						error = true;
						break;
					}	
				}
				
				// part unplaceable, skip
				if(error){
					continue;
				}
				
				var position = null;
				if(placed.length == 0){
					// first placement, put it on the left
					for(j=0; j<binNfp.length; j++){
						for(k=0; k<binNfp[j].length; k++){
							if(position === null || binNfp[j][k].x-path[0].x < position.x ){
								position = {
									x: binNfp[j][k].x-path[0].x,
									y: binNfp[j][k].y-path[0].y,
									id: path.id,
									rotation: path.rotation
								}
							}
						}
					}
					
					placements.push(position);
					placed.push(path);
					
					continue;
				}
				
				var clipperBinNfp = [];
				for(j=0; j<binNfp.length; j++){
					clipperBinNfp.push(toClipperCoordinates(binNfp[j]));
				}
				
				ClipperLib.JS.ScaleUpPaths(clipperBinNfp, self.config.clipperScale);
				
				var clipper = new ClipperLib.Clipper();
				var combinedNfp = new ClipperLib.Paths();
				
				
				for(j=0; j<placed.length; j++){			
					key = JSON.stringify({A:placed[j].id,B:path.id,inside:false,Arotation:placed[j].rotation,Brotation:path.rotation});
					nfp = self.nfpCache[key];
										
					if(!nfp){
						continue;
					}
					
					for(k=0; k<nfp.length; k++){
						var clone = toClipperCoordinates(nfp[k]);
						for(m=0; m<clone.length; m++){
							clone[m].X += placements[j].x;
							clone[m].Y += placements[j].y;
						}
						
						ClipperLib.JS.ScaleUpPath(clone, self.config.clipperScale);
						clone = ClipperLib.Clipper.CleanPolygon(clone, 0.0001*self.config.clipperScale);
						var area = Math.abs(ClipperLib.Clipper.Area(clone));
						if(clone.length > 2 && area > 0.1*self.config.clipperScale*self.config.clipperScale){
							clipper.AddPath(clone, ClipperLib.PolyType.ptSubject, true);
						}
					}		
				}
				
				if(!clipper.Execute(ClipperLib.ClipType.ctUnion, combinedNfp, ClipperLib.PolyFillType.pftNonZero, ClipperLib.PolyFillType.pftNonZero)){
					continue;
				}
				
				// difference with bin polygon
				var finalNfp = new ClipperLib.Paths();
				clipper = new ClipperLib.Clipper();
				
				clipper.AddPaths(combinedNfp, ClipperLib.PolyType.ptClip, true);
				clipper.AddPaths(clipperBinNfp, ClipperLib.PolyType.ptSubject, true);
				if(!clipper.Execute(ClipperLib.ClipType.ctDifference, finalNfp, ClipperLib.PolyFillType.pftNonZero, ClipperLib.PolyFillType.pftNonZero)){
					continue;
				}
				
				finalNfp = ClipperLib.Clipper.CleanPolygons(finalNfp, 0.0001*self.config.clipperScale);
				
				for(j=0; j<finalNfp.length; j++){
					var area = Math.abs(ClipperLib.Clipper.Area(finalNfp[j]));
					if(finalNfp[j].length < 3 || area < 0.1*self.config.clipperScale*self.config.clipperScale){
						finalNfp.splice(j,1);
						j--;
					}
				}
				
				if(!finalNfp || finalNfp.length == 0){
					continue;
				}
				
				var f = [];
				for(j=0; j<finalNfp.length; j++){
					// back to normal scale
					f.push(toNestCoordinates(finalNfp[j], self.config.clipperScale));
				}
				finalNfp = f;
				
				// choose placement that results in the smallest bounding box
				// could use convex hull instead, but it can create oddly shaped nests (triangles or long slivers) which are not optimal for real-world use
				// todo: generalize gravity direction
				var minwidth = null;
				var minarea = null;
				var minx = null;
				var nf, area, shiftvector;

				for(j=0; j<finalNfp.length; j++){
					nf = finalNfp[j];
					if(Math.abs(GeometryUtil.polygonArea(nf)) < 2){
						continue;
					}
					
					for(k=0; k<nf.length; k++){
						var allpoints = [];
						for(m=0; m<placed.length; m++){
							for(n=0; n<placed[m].length; n++){
								allpoints.push({x:placed[m][n].x+placements[m].x, y: placed[m][n].y+placements[m].y});
							}
						}
						
						shiftvector = {
							x: nf[k].x-path[0].x,
							y: nf[k].y-path[0].y,
							id: path.id,
							rotation: path.rotation,
							nfp: combinedNfp
						};
						
						for(m=0; m<path.length; m++){
							allpoints.push({x: path[m].x+shiftvector.x, y:path[m].y+shiftvector.y});
						}
						
						var rectbounds = GeometryUtil.getPolygonBounds(allpoints);
						
						// weigh width more, to help compress in direction of gravity
						area = rectbounds.width*2 + rectbounds.height;
						
						if(minarea === null || area < minarea || (GeometryUtil.almostEqual(minarea, area) && (minx === null || shiftvector.x < minx))){
							minarea = area;
							minwidth = rectbounds.width;
							position = shiftvector;
							minx = shiftvector.x;
							var contained_area = rectbounds.width*rectbounds.height;
						}
					}
				}
				if(position){
					placed.push(path);
					placements.push(position);
				}

			}
			
			if(minwidth){
				fitness += minwidth/binarea;
			}
			



	
			
			
			
			for(i=0; i<placed.length; i++){
				var index = paths.indexOf(placed[i]);
				if(index >= 0){
					paths.splice(index,1);
				}
			}
			
			if(placements && placements.length > 0){
				allplacements.push(placements);
			}
			else{
				break; // something went wrong
			}
		}
		
		// there were parts that couldn't be placed
		fitness += 2*paths.length;
		

		
		// decreaseing simCost 
		for (i =0 ; i<allplacements.length -1 ;i++ )
		{
			SimCost+= Compare(allplacements[i],allplacements[i+1]) // sum all the similarties
		}
		SimCost = SimCost / allplacements.length // ratio overall simalarty 

		SimCost = 1 - SimCost // ratio overall diffrence

		fitness += SimCost // because less fitness is better

		return {placements: allplacements, fitness: fitness, paths: paths, area: binarea, contained_area:contained_area };
	};

}
(typeof window !== 'undefined' ? window : self).PlacementWorker = PlacementWorker;

// clipperjs uses alerts for warnings
function alert(message) { 
    console.log('alert: ', message);
}


function Resample(points, n)
{
	var I = PathLength(points) / (n - 1); // interval length
	var D = 0.0;
	var newpoints = new Array(points[0]);
	for (var i = 1; i < points.length; i++)
	{
		var d = Distance(points[i-1], points[i]);
		if ((D + d) >= I)
		{
			var qx = points[i-1].X + ((I - D) / d) * (points[i].X - points[i-1].X);
			var qy = points[i-1].Y + ((I - D) / d) * (points[i].Y - points[i-1].Y);
			var q = new Point(qx, qy);
			newpoints[newpoints.length] = q; // append new point 'q'
			points.splice(i, 0, q); // insert 'q' at position i in points s.t. 'q' will be the next i
			D = 0.0;
		}
		else D += d;
	}
	if (newpoints.length == n - 1) // somtimes we fall a rounding-error short of adding the last point, so add it if so
		newpoints[newpoints.length] = new Point(points[points.length - 1].X, points[points.length - 1].Y);
	return newpoints;
}
function IndicativeAngle(points)
{
	var c = Centroid(points);
	return Math.atan2(c.Y - points[0].Y, c.X - points[0].X);
}
function RotateBy(points, radians) // rotates points around centroid
{
	var c = Centroid(points);
	var cos = Math.cos(radians);
	var sin = Math.sin(radians);
	var newpoints = new Array();
	for (var i = 0; i < points.length; i++) {
		var qx = (points[i].X - c.X) * cos - (points[i].Y - c.Y) * sin + c.X
		var qy = (points[i].X - c.X) * sin + (points[i].Y - c.Y) * cos + c.Y;
		newpoints[newpoints.length] = new Point(qx, qy);
	}
	return newpoints;
}
function ScaleTo(points, size) // non-uniform scale; assumes 2D gestures (i.e., no lines)
{
	var B = BoundingBox(points);
	var newpoints = new Array();
	for (var i = 0; i < points.length; i++) {
		var qx = points[i].X * (size / B.Width);
		var qy = points[i].Y * (size / B.Height);
		newpoints[newpoints.length] = new Point(qx, qy);
	}
	return newpoints;
}
function TranslateTo(points, pt) // translates points' centroid
{
	var c = Centroid(points);
	var newpoints = new Array();
	for (var i = 0; i < points.length; i++) {
		var qx = points[i].X + pt.X - c.X;
		var qy = points[i].Y + pt.Y - c.Y;
		newpoints[newpoints.length] = new Point(qx, qy);
	}
	return newpoints;
}
function Vectorize(points) // for Protractor
{
	var sum = 0.0;
	var vector = new Array();
	for (var i = 0; i < points.length; i++) {
		vector[vector.length] = points[i].X;
		vector[vector.length] = points[i].Y;
		sum += points[i].X * points[i].X + points[i].Y * points[i].Y;
	}
	var magnitude = Math.sqrt(sum);
	for (var i = 0; i < vector.length; i++)
		vector[i] /= magnitude;
	return vector;
}
function OptimalCosineDistance(v1, v2) // for Protractor
{
	var a = 0.0;
	var b = 0.0;
	// 1  2  3  4  v1
	// 5  6  7  8  v2
	for (var i = 0; i < v1.length; i += 2) {
		a += v1[i] * v2[i] + v1[i+1] * v2[i+1]; // dot product
		b += v1[i] * v2[i+1] - v1[i+1] * v2[i]; // cross product
	}
	var angle = Math.atan(b / a);
	// tan = cross/dot
	return Math.acos(a * Math.cos(angle) + b * Math.sin(angle)); // this is the dashed line
}
function DistanceAtBestAngle(points, T, a, b, threshold)
{
	var x1 = Phi * a + (1.0 - Phi) * b;
	var f1 = DistanceAtAngle(points, T, x1);
	var x2 = (1.0 - Phi) * a + Phi * b;
	var f2 = DistanceAtAngle(points, T, x2);
	while (Math.abs(b - a) > threshold)
	{
		if (f1 < f2) {
			b = x2;
			x2 = x1;
			f2 = f1;
			x1 = Phi * a + (1.0 - Phi) * b;
			f1 = DistanceAtAngle(points, T, x1);
		} else {
			a = x1;
			x1 = x2;
			f1 = f2;
			x2 = (1.0 - Phi) * a + Phi * b;
			f2 = DistanceAtAngle(points, T, x2);
		}
	}
	return Math.min(f1, f2);
}
function DistanceAtAngle(points, T, radians)
{
	var newpoints = RotateBy(points, radians);
	return PathDistance(newpoints, T.Points);
}
function Centroid(points)
{
	var x = 0.0, y = 0.0;
	for (var i = 0; i < points.length; i++) {
		x += points[i].X;
		y += points[i].Y;
	}
	x /= points.length;
	y /= points.length;
	return new Point(x, y);
}
function BoundingBox(points)
{
	var minX = +Infinity, maxX = -Infinity, minY = +Infinity, maxY = -Infinity;
	for (var i = 0; i < points.length; i++) {
		minX = Math.min(minX, points[i].X);
		minY = Math.min(minY, points[i].Y);
		maxX = Math.max(maxX, points[i].X);
		maxY = Math.max(maxY, points[i].Y);
	}
	return new Rectangle(minX, minY, maxX - minX, maxY - minY);
}
function PathDistance(pts1, pts2)
{
	var d = 0.0;
	for (var i = 0; i < pts1.length; i++) // assumes pts1.length == pts2.length
		d += Distance(pts1[i], pts2[i]);
	return d / pts1.length;
}
function PathLength(points)
{
	var d = 0.0;
	for (var i = 1; i < points.length; i++)
		d += Distance(points[i - 1], points[i]);
	return d;
}
function Distance(p1, p2)
{
	var dx = p2.X - p1.X;
	var dy = p2.Y - p1.Y;
	return Math.sqrt(dx * dx + dy * dy);
}
function Deg2Rad(d) { return (d * Math.PI / 180.0); }