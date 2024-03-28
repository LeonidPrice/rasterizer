var canvas = document.getElementById('canvas');
var canvas_width = document.documentElement.scrollWidth;
var canvas_height = document.documentElement.scrollHeight;
canvas.width = canvas_width;
canvas.height = canvas_height;
var context = canvas.getContext('2d');
var buffer = context.getImageData(0, 0, canvas.width, canvas.height);
var step = buffer.width * 4;




var CANVAS = {
    /**
    * Inserts a pixel of the desired color into the canvas.
    * @param {number} x - coordinate x on the canvas.
    * @param {number} y - coordinate y on the canvas.
    * @param {number[]} color - [R, G, B, A] or [R, G, B] for full opacity. Default - black.
    */
    pixel: function(x, y, color=[0,0,0,255]) {
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

    line: function(p0, p1, color=[0,0,0,255]) {
        if (Math.abs(p1[0]-p0[0]) > Math.abs(p1[1]-p0[1])) {
            if (p0[0] > p1[0]) {
                [p0,p1] = [p1,p0];
            }

            var yi = DRAW.interpolate(p0[0],p0[1],p1[0],p1[1]);
            
            for (var x = p0[0]; x <= p1[0]; x++) {
                CANVAS.pixel(x, yi[(x-p0[0])|0], color);
            }
        } else {
            if (p0[1] > p1[1]) {
                [p0,p1] = [p1,p0];
            }

            var xi = DRAW.interpolate(p0[1],p0[0],p1[1],p1[0]);
            
            for (var y = p0[1]; y <= p1[1]; y++) {
                CANVAS.pixel(xi[(y-p0[1])|0], y, color);
            }
        }
    },

    triangle: {
        wireframe: function(p0, p1, p2, color=[0,0,0,255]) {
            DRAW.line(p0, p1, color);
            DRAW.line(p1, p2, color);
            DRAW.line(p2, p0, color);
        },

        filled: function(p0, p1, p2, color = [0,0,0,255]) {
            if (p1[1] < p0[1]) {[p0,p1] = [p1,p0]}
            if (p2[1] < p0[1]) {[p0,p2] = [p2,p0]}
            if (p2[1] < p1[1]) {[p1,p2] = [p2,p1]}
            
            var x01 = DRAW.interpolate(p0[1],p0[0],p1[1],p1[0]);
            var x12 = DRAW.interpolate(p1[1],p1[0],p2[1],p2[0]);
            var x02 = DRAW.interpolate(p0[1],p0[0],p2[1],p2[0]);
            
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

            for (var y = p0[1]; y <= p2[1]; y++) {
                var xl = x_left[y-p0[1]]|0;
                var xr = x_right[y-p0[1]]|0;
                for (var x = xl; x <= xr; x++) {
                    CANVAS.pixel(x,y,color);
                }
            }
        },

        shaded: function(p0,p1,p2,color = [0,0,0,255]) {
            if (p1[1] < p0[1]) {[p0,p1] = [p1,p0]}
            if (p2[1] < p0[1]) {[p0,p2] = [p2,p0]}
            if (p2[1] < p1[1]) {[p1,p2] = [p2,p1]}

            var x01 = DRAW.interpolate(p0[1],p0[0],p1[1],p1[0]);
            var h01 = DRAW.interpolate(p0[1],p0[2],p1[1],p1[2]);

            var x12 = DRAW.interpolate(p1[1],p1[0],p2[1],p2[0]);
            var h12 = DRAW.interpolate(p1[1],p1[2],p2[1],p2[2]);

            var x02 = DRAW.interpolate(p0[1],p0[0],p2[1],p2[0]);
            var h02 = DRAW.interpolate(p0[1],p0[2],p2[1],p2[2]);

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

            for (var y = p0[1]; y <= p2[1]; y++) {
                var x_left_current = x_left[y-p0[1]]|0;
                var x_right_current = x_right[y-p0[1]]|0;
                var h_segment = DRAW.interpolate(x_left_current,h_left[y-p0[1]],x_right_current,h_right[y-p0[1]]);
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
    }
}

DRAW.triangle.shaded([-200,-250,0.3],[200,50,0.1],[20,250,1.0],[0,255,0,255]);
CANVAS.update();
