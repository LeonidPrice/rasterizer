var canvas = document.getElementById('canvas');
var canvas_width = document.documentElement.scrollWidth;
var canvas_height = document.documentElement.scrollHeight;
// var canvas_width = 400;
// var canvas_height = 400;
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
    pixel: function(x, y, color = [0,0,0,255]) {
        var x = Math.round(canvas.width / 2 + x);
        var y = Math.round(canvas.height / 2 - y - 1);

        if (x < 0 ||  y < 0 || x >= canvas.width || y >= canvas.height) {
            return;
        }

        var offset = 4 * x + step * y;

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

    line: function(p0, p1, color = [0,0,0,255]) {
        if (Math.abs(p1[0]-p0[0]) > Math.abs(p1[1]-p0[1])) {
            if (p0[0] > p1[0]) {
                [p0,p1] = [p1,p0];
            }

            var yi = this.interpolate(p0[0],p0[1],p1[0],p1[1]);
            
            for (var x = p0[0]; x <= p1[0]; x++) {
                CANVAS.pixel(x, yi[x-p0[0]], color);
            }
        } else {
            if (p0[1] > p1[1]) {
                [p0,p1] = [p1,p0];
            }

            var xi = this.interpolate(p0[1],p0[0],p1[1],p1[0]);
            
            for (var y = p0[1]; y <= p1[1]; y++) {
                CANVAS.pixel(xi[y-p0[1]], y, color);
            }
        }
    },

    triangle: {
        wireframe: function(p0, p1, p2, color = [0,0,0,255]) {
            DRAW.line(p0, p1, color);
            DRAW.line(p1, p2, color);
            DRAW.line(p2, p0, color);
        },

        fill: function(p0, p1, p2, color = [0,0,0,255]) {
            if (p1[1] < p0[1]) {
                [p0,p1] = [p1,p0];
            }
            if (p2[1] < p0[1]) {
                [p0,p2] = [p2,p0];
            }
            if (p2[1] < p1[1]) {
                [p1,p2] = [p2,p1];
            }
            
            var x01 = DRAW.interpolate(p0[1],p0[0],p1[1],p1[0]);
            var x12 = DRAW.interpolate(p1[1],p1[0],p2[1],p2[0]);
            var x02 = DRAW.interpolate(p0[1],p0[0],p2[1],p2[0]);
            
            x01.pop();
            var x012 = x01.concat(x12);
            var x_middle = Math.floor(x012.length/2);
            
            if (x02[x_middle] < x012[x_middle]) {
                var x_left = x02;
                var x_right = x012;
            } else {
                var x_left = x012;
                var x_right = x02;
            }

            for (var y = p0[1]; y <= p2[1]; y++) {
                for (var x = x_left[y-p0[1]]; x <= x_right[y-p0[1]]; x++) {
                    CANVAS.pixel(x,y,color);
                }
            }
        },
    }
}



CANVAS.pixel(0,0,[255,0,0,255]);
DRAW.triangle.fill([-200,-250],[200,50],[20,250],[0,255,0,255]);
DRAW.triangle.wireframe([-200,-250],[200,50],[20,250],[0,0,0,255]);


CANVAS.update();

