var canvas = document.getElementById('canvas');
var canvas_width = document.documentElement.scrollWidth;
var canvas_height = document.documentElement.scrollHeight;
canvas.width = canvas_height;
canvas.height = canvas_height;
var context = canvas.getContext('2d');
var buffer = context.getImageData(0, 0, canvas.width, canvas.height);
var step = buffer.width * 4;

var viewport_size = 1;
var projection_plane_z = 1;

var Point = function(x,y,h,color) {
    if (!(this instanceof Point)) {return new Point(x,y,h,color);}
    this.x = x;
    this.y = y;
    this.h = h;
    this.color = color;
}

var Vertex = function(x,y,z) {
    if (!(this instanceof Vertex)) {return new Vertex(x,y,z);}
    this.x = x;
    this.y = y;
    this.z = z;
}


var CANVAS = {
    /**
    * Inserts a pixel of the desired color into the canvas.
    * @param {number} x - coordinate x on the canvas.
    * @param {number} y - coordinate y on the canvas.
    * @param {number[]} color - [R, G, B, A] or [R, G, B] for full opacity. Default - black.
    */
    pixel: function(x,y,color=[0,0,0,255]) {
        x = (canvas.width/2 + x)|0;
        y = (canvas.height/2 - y - 1)|0;

        if (x < 0 ||  y < 0 || x >= canvas.width || y >= canvas.height) {
            return;
        }

        var offset = 4*x + step*y;

        buffer.data[offset++] = color[0];
        buffer.data[offset++] = color[1];
        buffer.data[offset++] = color[2];
        buffer.data[offset++] = color[3] || 255;
    },

    viewport_to_canvas: function(point_2d) {
        return Point(
            point_2d.x*canvas.width/viewport_size,
            point_2d.y*canvas.height/viewport_size);
    },

    project_vertex: function(vertex) {
        return CANVAS.viewport_to_canvas(
            Point(
                vertex.x*projection_plane_z/vertex.z,
                vertex.y*projection_plane_z/vertex.z));
    },

    update: function() {
        context.putImageData(buffer, 0, 0);
    }
}

var DRAW = {
    interpolate: function(i0,d0,i1,d1) {
        if (i0 == i1) {
            return [d0];
        }

        var values = [];
        var a = (d1-d0)/(i1-i0);
        var d = d0;

        for (var i = i0; i <= i1; i++) {
            values.push(d);
            d += a;
        }

        return values;
    },

    line: function(p0,p1,color=[0,0,0,255]) {
        if (Math.abs(p1.x-p0.x) > Math.abs(p1.y-p0.y)) {
            if (p0.x > p1.x) {[p0,p1] = [p1,p0];}

            var yi = DRAW.interpolate(p0.x,p0.y,p1.x,p1.y);
            
            for (var x = p0.x; x <= p1.x; x++) {
                CANVAS.pixel(x, yi[(x-p0.x)|0], color);
            }
        } else {
            if (p0.y > p1.y) {[p0,p1] = [p1,p0];}

            var xi = DRAW.interpolate(p0.y,p0.x,p1.y,p1.x);
            
            for (var y = p0.y; y <= p1.y; y++) {
                CANVAS.pixel(xi[(y-p0.y)|0], y, color);
            }
        }
    },

    triangle: {
        wireframe: function(p0,p1,p2, color=[0,0,0,255]) {
            DRAW.line(p0, p1, color);
            DRAW.line(p1, p2, color);
            DRAW.line(p2, p0, color);
        },

        filled: function(p0,p1,p2, color = [0,0,0,255]) {
            if (p1.y < p0.y) {[p0,p1] = [p1,p0];}
            if (p2.y < p0.y) {[p0,p2] = [p2,p0];}
            if (p2.y < p1.y) {[p1,p2] = [p2,p1];}
            
            var x01 = DRAW.interpolate(p0.y,p0.x,p1.y,p1.x);
            var x12 = DRAW.interpolate(p1.y,p1.x,p2.y,p2.x);
            var x02 = DRAW.interpolate(p0.y,p0.x,p2.y,p2.x);
            
            x01.pop();
            var x012 = x01.concat(x12);
            var x_middle = (x012.length/2)|0;
            
            if (x02[x_middle] < x012[x_middle]) {
                var x_left = x02;
                var x_right = x012;
            } else {
                var x_left = x012;
                var x_right = x02;
            }

            for (var y = p0.y; y <= p2.y; y++) {
                var xl = x_left[y-p0.y]|0;
                var x_right_current = x_right[y-p0.y]|0;
                for (var x = xl; x <= x_right_current; x++) {
                    CANVAS.pixel(x,y,color);
                }
            }
        },

        shaded: function(p0,p1,p2,color = [0,0,0,255]) {
            if (p1.y < p0.y) {[p0,p1] = [p1,p0];}
            if (p2.y < p0.y) {[p0,p2] = [p2,p0];}
            if (p2.y < p1.y) {[p1,p2] = [p2,p1];}

            var x01 = DRAW.interpolate(p0.y,p0.x,p1.y,p1.x);
            var h01 = DRAW.interpolate(p0.y,p0.h,p1.y,p1.h);

            var x12 = DRAW.interpolate(p1.y,p1.x,p2.y,p2.x);
            var h12 = DRAW.interpolate(p1.y,p1.h,p2.y,p2.h);

            var x02 = DRAW.interpolate(p0.y,p0.x,p2.y,p2.x);
            var h02 = DRAW.interpolate(p0.y,p0.h,p2.y,p2.h);

            x01.pop();
            h01.pop();

            var x012 = x01.concat(x12);
            var h012 = h01.concat(h12);

            var x_middle = (x012.length/2)|0;

            if (x02[x_middle] < x012[x_middle]) {
                var x_left = x02; var x_right = x012;
                var h_left = h02; var h_right = h012;
            } else {
                var x_left = x012; var x_right = x02;
                var h_left = h012; var h_right = h02;
            }

            for (var y = p0.y; y <= p2.y; y++) {
                var x_left_current = x_left[y-p0.y]|0;
                var x_right_current = x_right[y-p0.y]|0;
                var h_segment = DRAW.interpolate(x_left_current,h_left[y-p0.y],x_right_current,h_right[y-p0.y]);
                for (var x = x_left_current; x <= x_right_current; x++) {
                    var hue = h_segment[x-x_left_current];
                    var shaded_color = [
                        color[0]*hue,
                        color[1]*hue,
                        color[2]*hue,
                        255
                    ]; 
                    CANVAS.pixel(x,y,shaded_color);
                }
            }
        },

        gradient: function(p0,p1,p2) {
            if (p1.y < p0.y) {[p0,p1] = [p1,p0];}
            if (p2.y < p0.y) {[p0,p2] = [p2,p0];}
            if (p2.y < p1.y) {[p1,p2] = [p2,p1];}

            var x01 = DRAW.interpolate(p0.y,p0.x,p1.y,p1.x);
            var h01 = DRAW.interpolate(p0.y,p0.h,p1.y,p1.h);
            var r01 = DRAW.interpolate(p0.y,p0.color[0],p1.y,p1.color[0]);
            var g01 = DRAW.interpolate(p0.y,p0.color[1],p1.y,p1.color[1]);
            var b01 = DRAW.interpolate(p0.y,p0.color[2],p1.y,p1.color[2]);
            var a01 = DRAW.interpolate(p0.y,p0.color[3],p1.y,p1.color[3]);

            var x12 = DRAW.interpolate(p1.y,p1.x,p2.y,p2.x);
            var h12 = DRAW.interpolate(p1.y,p1.h,p2.y,p2.h);
            var r12 = DRAW.interpolate(p1.y,p1.color[0],p2.y,p2.color[0]);
            var g12 = DRAW.interpolate(p1.y,p1.color[1],p2.y,p2.color[1]);
            var b12 = DRAW.interpolate(p1.y,p1.color[2],p2.y,p2.color[2]);
            var a12 = DRAW.interpolate(p1.y,p1.color[3],p2.y,p2.color[3]);

            var x02 = DRAW.interpolate(p0.y,p0.x,p2.y,p2.x);
            var h02 = DRAW.interpolate(p0.y,p0.h,p2.y,p2.h);
            var r02 = DRAW.interpolate(p0.y,p0.color[0],p2.y,p2.color[0]);
            var g02 = DRAW.interpolate(p0.y,p0.color[1],p2.y,p2.color[1]);
            var b02 = DRAW.interpolate(p0.y,p0.color[2],p2.y,p2.color[2]);
            var a02 = DRAW.interpolate(p0.y,p0.color[3],p2.y,p2.color[3]);

            x01.pop();
            h01.pop();
            r01.pop();
            g01.pop();
            b01.pop();
            a01.pop();

            var x012 = x01.concat(x12);
            var h012 = h01.concat(h12);
            var r012 = r01.concat(r12);
            var g012 = g01.concat(g12);
            var b012 = b01.concat(b12);
            var a012 = a01.concat(a12);

            var x_middle = (x012.length/2)|0;

            if (x02[x_middle] < x012[x_middle]) {
                var x_left = x02; var x_right = x012;
                var h_left = h02; var h_right = h012;
                var r_left = r02; var r_right = r012;
                var g_left = g02; var g_right = g012;
                var b_left = b02; var b_right = b012;
                var a_left = a02; var a_right = a012;
            } else {
                var x_left = x012; var x_right = x02;
                var h_left = h012; var h_right = h02;
                var r_left = r012; var r_right = r02;
                var g_left = g012; var g_right = g02;
                var b_left = b012; var b_right = b02;
                var a_left = a012; var a_right = a02;
            }

            for (var y = p0.y; y <= p2.y; y++) {
                var x_left_current = x_left[y-p0.y]|0;
                var x_right_current = x_right[y-p0.y]|0;
                var h_segment = DRAW.interpolate(x_left_current,h_left[y-p0.y],x_right_current,h_right[y-p0.y]);
                var r_segment = DRAW.interpolate(x_left_current,r_left[y-p0.y],x_right_current,r_right[y-p0.y]);
                var g_segment = DRAW.interpolate(x_left_current,g_left[y-p0.y],x_right_current,g_right[y-p0.y]);
                var b_segment = DRAW.interpolate(x_left_current,b_left[y-p0.y],x_right_current,b_right[y-p0.y]);
                var a_segment = DRAW.interpolate(x_left_current,a_left[y-p0.y],x_right_current,a_right[y-p0.y]);
                for (var x = x_left_current; x <= x_right_current; x++) {
                    var hue = h_segment[x-x_left_current];
                    var r = r_segment[x-x_left_current];
                    var g = g_segment[x-x_left_current];
                    var b = b_segment[x-x_left_current];
                    var a = a_segment[x-x_left_current];
                    var shaded_color = [
                        r*hue,
                        g*hue,
                        b*hue,
                        a
                    ]; 
                    CANVAS.pixel(x,y,shaded_color);
                }
            }
        },
    }
}

function triangle_coordinates(r) {
    var xa = r*Math.cos(210*Math.PI/180);
    var ya = r*Math.sin(210*Math.PI/180);

    var xb = 0;
    var yb = r;

    var xc = r*Math.cos(30*Math.PI/180);
    var yc = -r*Math.sin(30*Math.PI/180);

    var shift = (2*r-(r-yc))/2

    return [[xa,ya-shift],[xb,yb-shift],[xc,yc-shift]];
}

var points = triangle_coordinates(300);
var P0 = Point(points[0][0],points[0][1], 1.0, [255,0,0,255]);
var P1 = Point(points[1][0],points[1][1], 1.0, [0,255,0,255]);
var P2 = Point(points[2][0],points[2][1], 1.0, [0,0,255,255]);

//DRAW.triangle.gradient(P0,P1,P2,[255,0,0,255]);

var vAf = Vertex(-2, -0.5, 5);
var vBf = Vertex(-2, 0.5, 5);
var vCf = Vertex(-1, 0.5, 5);
var vDf = Vertex(-1, -0.5, 5);

var vAb = Vertex(-2, -0.5, 6);
var vBb = Vertex(-2, 0.5, 6);
var vCb = Vertex(-1, 0.5, 6);
var vDb = Vertex(-1, -0.5, 6);

//front side
DRAW.line(CANVAS.project_vertex(vAf),CANVAS.project_vertex(vBf), [0,0,255,255]);
DRAW.line(CANVAS.project_vertex(vBf),CANVAS.project_vertex(vCf), [0,0,255,255]);
DRAW.line(CANVAS.project_vertex(vCf),CANVAS.project_vertex(vDf), [0,0,255,255]);
DRAW.line(CANVAS.project_vertex(vDf),CANVAS.project_vertex(vAf), [0,0,255,255]);

//back side
DRAW.line(CANVAS.project_vertex(vAb),CANVAS.project_vertex(vBb), [255,0,0,255]);
DRAW.line(CANVAS.project_vertex(vBb),CANVAS.project_vertex(vCb), [255,0,0,255]);
DRAW.line(CANVAS.project_vertex(vCb),CANVAS.project_vertex(vDb), [255,0,0,255]);
DRAW.line(CANVAS.project_vertex(vDb),CANVAS.project_vertex(vAb), [255,0,0,255]);

//lateral side
DRAW.line(CANVAS.project_vertex(vAf),CANVAS.project_vertex(vAb), [0,255,0,255]);
DRAW.line(CANVAS.project_vertex(vBf),CANVAS.project_vertex(vBb), [0,255,0,255]);
DRAW.line(CANVAS.project_vertex(vCf),CANVAS.project_vertex(vCb), [0,255,0,255]);
DRAW.line(CANVAS.project_vertex(vDf),CANVAS.project_vertex(vDb), [0,255,0,255]);


CANVAS.update();
