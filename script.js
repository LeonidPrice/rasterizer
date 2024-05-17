const projectionPlaneZ = 1;
const epsilon = 1e-8;


/**
 * Represents a three-dimensional vector or point in space.
 * @constructor
 * @param {Number} x - The x-coordinate of the vertex.
 * @param {Number} y - The y-coordinate of the vertex.
 * @param {Number} z - The z-coordinate of the vertex.
 * @param {Number} w - Optional w-coordinate of the vertex (defaults to null).
 */
class Vertex {
    /**
     * Creates a new Vertex instance.
     * @param {Number} x - The x-coordinate of the vertex.
     * @param {Number} y - The y-coordinate of the vertex.
     * @param {Number} z - The z-coordinate of the vertex.
     * @param {Number} w - Optional w-coordinate of the vertex (defaults to null).
     * @param {Number[]} color - An array of the RGBA color channels. Default - black.
     */
    constructor (x, y, z, w = 1, color = [0,0,0,255]) {
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;
        this.color = color;
    }

    /**
     * Computes the length (magnitude) of the vector.
     * @description Algorithmic complexity: O(1).
     * @returns {Number} The length of the vector.
     */
    Length() {
        return Math.sqrt(Math.pow(this.x,2)+Math.pow(this.y,2)+Math.pow(this.z,2));
    }

    /**
     * Computes the dot product of this vector with another vector.
     * @description Algorithmic complexity: O(1).
     * @param {Vertex} vertex - The other vector.
     * @returns {Vertex} The dot product of the two vectors.
     */
    DotProduct(vertex) {
        return this.x*vertex.x+this.y*vertex.y+this.z*vertex.z;
    }

    /**
     * Adds another vector to this vector.
     * @description Algorithmic complexity: O(1).
     * @param {Vertex} vertex - The vector to be added.
     * @returns {Vertex} The result of adding the two vectors.
     */
    Addition(vertex) {
        return new Vertex(
            this.x+vertex.x,
            this.y+vertex.y,
            this.z+vertex.z);
    }

    /**
     * Subtracts another vector from this vector.
     * @description Algorithmic complexity: O(1).
     * @param {Vertex} vertex - The vector to be subtracted.
     * @returns {Vertex} The result of subtracting the two vectors.
     */
    Subtraction(vertex) {
        return new Vertex(
            this.x-vertex.x,
            this.y-vertex.y,
            this.z-vertex.z);
    }

    /**
     * Multiplies the vector by a scalar.
     * @description Algorithmic complexity: O(1).
     * @param {Number} scalar - The scalar value.
     * @returns {Vertex} The result of multiplying the vector by the scalar.
     */
    MultiplyScalar(scalar) {
        return new Vertex(
            this.x*scalar,
            this.y*scalar,
            this.z*scalar);
    }

    /**
     * Multiplies this vector component-wise with another vector.
     * @description Algorithmic complexity: O(1).
     * @param {Vertex} vertex - The other vector.
     * @returns {Vertex} The result of component-wise multiplication.
     */
    MultiplyVector(vertex) {
        return new Vertex(
            this.x*vertex.x,
            this.y*vertex.y,
            this.z*vertex.z);
    }

    /**
     * Computes the cross product of this vector with another vector.
     * @description Algorithmic complexity: O(1).
     * @param {Vertex} vertex - The other vector.
     * @returns {Vertex} The cross product of the two vectors.
     */
    CrossProduct(vertex) {
        return new Vertex(
            this.y*vertex.z-this.z*vertex.y,
            this.z*vertex.x-this.x*vertex.z,
            this.x*vertex.y-this.y*vertex.x
        )
    }

    /**
     * Finds the midpoint of two vectors.
     * @description Algorithmic complexity: O(1).
     * @param {Vertex} vertex - The other vector.
     * @returns {Vertex} The cross product of the two vectors.
     */
    MiddlePoint(vertex) {
        return new Vertex((vertex.x-this.x)/2, (vertex.y-this.y)/2, (vertex.z-this.z)/2, this.color);
    }
}

/**
 * Presents a mathematical matrix and methods of its processing.
 * @constructor
 * @param {Number[][]} data - Two-dimensional array representing the matrix data.
 */
class Matrix {
    /**
     * Creates a new Matrix instance.
     * @param {Number[][]} data - Two-dimensional array representing the matrix data.
     */
    constructor(data) {
        this.data = data;
    }

    static errors = {
        degree() {console.log('Matrix error: degree is negative or equal to zero.\nThe degree of the matrix must be positive and greater than one.')},
        compatible() {console.log("Matrix error: the matrices aren't compatible.\nThe number of rows and columns of first and second matrices must be equal.")},
        multiply() {console.log("Matrix error: the matrices aren't compatible.\nThe number of columns of first matrix must be equal to the number of rows of the second matrix.")},
        determinant() {console.log("Matrix error: the matrices aren't compatible.\nThe number of rows and columns mast be equal.")},
        iterations() {console.log("Matrix error: the number of iterations is not enough to obtain a solution.")}
    }

    /**
     * Performs linear interpolation between two values.
     * @description Algorithmic complexity: O(i1 - i0), where i0 and i1 are the indices of the first and second data points, respectively.
     * @param {Number} i0 - The index of the first data point.
     * @param {Number} d0 - The value of the first data point.
     * @param {Number} i1 - The index of the second data point.
     * @param {Number} d1 - The value of the second data point.
     * @returns {Number[]} An array of interpolated values between the two data points.
     */
    static Interpolation(i0, d0, i1, d1) {
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
    }

    /**
     * Sorts an array of numbers.
     * @description Algorithmic complexity: O(n^2), where n is the length of an array.
     * @param {Matrix} matrix - Array as a one-dimensional matrix.
     * @returns {[Matrix|Number[][], number]} An array where the first element is sorted array, and the second element represents the index of the largest element in the source array.
     */
    static SelectionSort(matrix) {
        var data = matrix.data.flat();
        var sortedArray = [];
        var index = [];
    
        while (data.length > 0) {
            var largest = data[0];
            var largestIndex = 0;
            for (var i = 1; i < data.length; i++) {
                if (data[i] > largest) {
                    largest = data[i];
                    largestIndex = i;
                }
            }
            index.push(largestIndex);
            sortedArray.push(data.splice(largestIndex, 1)[0]);
        }

        return [new Matrix([sortedArray]), index[0]];
    }
    
    /**
     * Creates a unit matrix.
     * @description Algorithmic complexity: O(n^2), where n is the dimensionality of the square matrix.
     * @param {Number} order - Dimension of the matrix.
     * @returns {Matrix|Number[][]}
     */
    static Identity(order) {
        var result = [];
        for (var i = 0; i < order; i++) {
            result[i] = [];
            for (var j = 0; j < order; j++) {
                result[i][j] = i === j ? 1 : 0;
            }
        }
        return new Matrix(result);
    }

    /**
     * Creates a rotation matrix around the x-axis.
     * @description Algorithmic complexity: O(1).
     * @param {Number} angle - Angle in degrees.
     * @returns {Matrix|Number[][]}
     */
    static RotationX(angle) {
        angle *= Math.PI/180;
        var c = Math.cos(angle);
        var s = Math.sin(angle);

        return new Matrix([
            [1, 0, 0, 0],
            [0, c,-s, 0],
            [0, s, c, 0],
            [0, 0, 0, 1]
        ]);
    }

    /**
     * Creates a rotation matrix around the y-axis.
     * @description Algorithmic complexity: O(1).
     * @param {Number} angle - Angle in degrees.
     * @returns {Matrix|Number[][]}
     */
    static RotationY(angle) {
        angle *= Math.PI/180;
        var c = Math.cos(angle);
        var s = Math.sin(angle);

        return new Matrix([
            [c, 0, s, 0],
            [0, 1, 0, 0],
            [-s,0, c, 0],
            [0, 0, 0, 1]
        ]);
    }

    /**
     * Creates a rotation matrix around the z-axis.
     * @description Algorithmic complexity: O(1).
     * @param {Number} angle - Angle in degrees.
     * @returns {Matrix|Number[][]}
     */
    static RotationZ(angle) {
        angle *= Math.PI/180;
        var c = Math.cos(angle);
        var s = Math.sin(angle);

        return new Matrix([
            [c,-s, 0, 0],
            [s,c, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
        ]);
    }

    /**
     * Creates a complex rotation matrix around three axes.
     * @description Algorithmic complexity: O(1).
     * @param {Number} rx - Angle x in degrees. Default - 0.
     * @param {Number} ry - Angle y in degrees. Default - 0.
     * @param {Number} rz - Angle z in degrees. Default - 0.
     * @returns {Matrix|Number[][]}
     */
    static Rotation(rx = 0, ry = 0, rz = 0) {
        rx *= Math.PI/180;
        ry *= Math.PI/180;
        rz *= Math.PI/180;
        var cAlpha = Math.cos(rx);
        var sAlpha = Math.sin(rx);
        var cBetha = Math.cos(ry);
        var sBetha = Math.sin(ry);
        var cGamma = Math.cos(rz);
        var sGamma = Math.sin(rz);

        return new Matrix([
            [cBetha*cGamma,                      -sGamma*cBetha,                      sBetha,         0],
            [sAlpha*sBetha*cGamma+sGamma*cAlpha, -sAlpha*sBetha*sGamma+cAlpha*cGamma, -sAlpha*cBetha, 0],
            [sAlpha*sGamma-sBetha*cAlpha*cGamma, sAlpha*cGamma+sBetha*sGamma*cAlpha,  cAlpha*cBetha,  0],
            [0,                                  0,                                   0,              1]
        ]);
    }

    /**
     * Creates a scaling matrix on three axes.
     * @description Algorithmic complexity: O(1).
     * @param {Number} sx - Scaling factor on the x-axis. Default - 1.
     * @param {Number} sy - Scaling factor on the y-axis. Default - 1.
     * @param {Number} sz - Scaling factor on the z-axis. Default - 1.
     * @returns {Matrix|Number[][]}
     */
    static Scale(sx = 0, sy = 0, sz = 0) {
        return new Matrix([
            [sx,0, 0, 0],
            [0,sy, 0, 0],
            [0, 0,sz, 0],
            [0, 0, 0, 1]
        ]);
    }

    /**
     * Creates a translation matrix on three axes.
     * @description Algorithmic complexity: O(1).
     * @param {Vertex} vertex - Translation factor on the x,y,z-axes.
     * @returns {Matrix|Number[][]}
     */
    static Translation(vertex) {
        var tx = vertex.x;
        var ty = vertex.y;
        var tz = vertex.z;
        return new Matrix([
            [1, 0, 0,tx],
            [0, 1, 0,ty],
            [0, 0, 1,tz],
            [0, 0, 0, 1]
        ]);
    }

    /**
     * Adds compatible matrices.
     * @description Algorithmic complexity: O(n^2), where n is the dimensionality of the square matrix.
     * @param {Matrix} matrix - Second matrix.
     * @returns {Matrix|Number[][]}
     */
    Addition(matrix) {
        if (!this || !matrix || this.data.length === 0 || matrix.data.length === 0 || this.data.length !== matrix.data.length || this.data[0].length !== matrix.data[0].length) {
            Matrix.errors.compatible();
            return;
        }

        var result = [];
        for (var i = 0; i < this.data.length; i++) {
            result[i] = [];
            for (var j = 0; j < this.data.length; j++) {
                result[i][j] = this.data[i][j] + matrix.data[i][j];
            }
        }
        return new Matrix(result);
    }

    /**
     * Subtracts compatible matrices.
     * @description Algorithmic complexity: O(n^2), where n is the dimensionality of the square matrix.
     * @param {Matrix} matrix - Second matrix.
     * @returns {Matrix|Number[][]}
     */
    Subtraction(matrix) {
        if (!this || !matrix || this.data.length === 0 || matrix.data.length === 0 || this.data.length !== matrix.data.length || this.data[0].length !== matrix.data[0].length) {
            Matrix.errors.compatible();
            return;
        }

        var result = [];
        for (var i = 0; i < this.data.length; i++) {
            result[i] = [];
            for (var j = 0; j < this.data.length; j++) {
                result[i][j] = this.data[i][j] - matrix.data[i][j];
            }
        }
        return new Matrix(result);
    }

    /**
     * Transposes a matrix of size i by j into a matrix of size j by i.
     * @description Algorithmic complexity: O(i*j), where i - amount of rows, j - amount of columns.
     * @returns {Matrix|Number[][]}
     */
    Transpose() {
        var result = [];
        for (var i = 0; i < this.data[0].length; i++) {
            result[i] = [];
            for (var j = 0; j < this.data.length; j++) {
                result[i][j] = this.data[j][i];
            }
        }
        return new Matrix(result);
    }

    /**
     * Multiplies a matrix by a number.
     * @description Algorithmic complexity: O(n^2), where n is the dimensionality of the square matrix.
     * @param {Number} scalar - Any number.
     * @returns @returns {Matrix|Number[][]}
     */
    MultiplyScalar(scalar) {
        var result = [];
        for (var i = 0; i < this.data.length; i++) {
            result[i] = [];
            for (var j = 0; j < this.data[0].length; j++) {
                result[i][j] = this.data[i][j] * scalar;
            }
        }
        return new Matrix(result);
    }

    /**
     * Multiplies compatible matrices.
     * @description Algorithmic complexity: O(i*j*k), where i - amount of rows in first matrix, j - amount of columns in second matrix, k - amount of rows in second matrix.
     * @param {Matrix} matrix - Second matrix.
     * @returns {Matrix|Number[][]}
     */
    MultiplyMatrix(matrix) {
        if (!this || !matrix || this.data.length === 0 || matrix.data.length === 0 || this.data[0].length !== matrix.data.length) {
            Matrix.errors.multiply();
            return;
        }

        var result = [];
        for (var i = 0; i < this.data.length; i++) {
            result[i] = [];
            for (var j = 0; j < matrix.data[0].length; j++) {
                result[i][j] = 0;
            }
        }
        for (var i = 0; i < this.data.length; i++) {
            for (var j = 0; j < matrix.data[0].length; j++) {
                var s = 0;
                for (var k = 0; k < matrix.data.length; k++) {
                    s += this.data[i][k] * matrix.data[k][j];
                }
                result[i][j] = s;
            }
        }
        return new Matrix(result);
    }

    /**
     * Elevates the matrix to positive degree.
     * @description Algorithmic complexity: O(n * i * j * k), where n - value of the degree,i - amount of rows in this matrix, j - amount of columns in second matrix, k - amount of rows in second matrix.
     * @param {Number} degree - Degree of the matrix.
     * @returns {Matrix|Number[][]}
     */
    Degree(degree) {
        if (degree <= 0) {
            Matrix.errors.degree();
            return;
        }
        if (degree == 1) {
            return new Matrix(this.data);
        }

        var previousMatrix = new Matrix(this.data).MultiplyMatrix(new Matrix(this.data));

        if (degree == 2) {
            return previousMatrix;
        } else {
            var nextMatrix;
            for (var i = 0; i <= degree-3; i++) {
                nextMatrix = new Matrix(this.data).MultiplyMatrix(previousMatrix);
                previousMatrix = nextMatrix;
            }
            return nextMatrix;
        }
    }

    /**
     * Calculates the rank of the matrix.
     * @description lgorithmic complexity: O(min(m,n) * m * n), where n - count of rows, m - count of columns.
     * @param {Number} epsilon - Error magnitude. Default: 1e-6.
     * @returns {Number} Rank of the matrix.
     */
    Rank(epsilon = 1e-6) {
        if (!this || this.data.length === 0 || this.data[0].length === 0) {
            return 0;
        }
    
        var lineUsed = [];
        var rank = 0;
    
        for (var i = 0; i < this.data.length; i++) {
            lineUsed[i] = false;
        }
    
        for (var i = 0; i < this.data[0].length; i++) {
            var j;
            for (j = 0; j < this.data.length; j++) {
                if (!lineUsed[j] && Math.abs(this.data[j][i]) > epsilon) {
                    break;
                }
            } 
            if (j !== this.data.length) {
                rank++;
                lineUsed[j] = true;
                for (var p = i + 1; p < this.data[0].length; p++) {
                    this.data[j][p] /= this.data[j][i];
                }
                for (var k = 0; k < this.data.length; k++) {
                    if (k !== j && Math.abs(this.data[k][i]) > epsilon) {
                        for (var p = i + 1; p < this.data[0].length; p++) {
                            this.data[k][p] -= this.data[j][p] * this.data[k][i];
                        }
                    }
                }
            }
        }
        return rank;
    }

    /**
     * Calculates the first norma of the matrix.
     * @description Algorithmic complexity: O(n^2), where n is the dimensionality of the square matrix. 
     * @returns {Number} First norm of the matrix.
     */
    FirstNorma() {
        if (!this || this.data.length === 0 || this.data[0].length === 0) {
            return;
        }

        var firstNorma = [];
        var sum = 0;
        for (var i = 0; i < this.data[0].length; i++) {
            for (var j = 0; j < this.data.length; j++) {
                sum += Math.abs(this.data[j][i]);
            }
            firstNorma.push(sum);
            sum = 0;
        }
        return Matrix.SelectionSort(new Matrix([firstNorma]))[0].data[0][0];
    }  
    
    /**
     * Calculates the infinity norma of the matrix.
     * @description Algorithmic complexity: O(n^2), where n is the dimensionality of the square matrix. 
     * @returns {Number} Infinity norm of the matrix.
    */
    InfinityNorma() {
        if (!this || this.data.length === 0 || this.data[0].length === 0) {
            return;
        }

        var infinityNorma = [];
        var sum = 0;
        for (var i = 0; i < this.data.length; i++) { // Проходим по строкам
            for (var j = 0; j < this.data[0].length; j++) { // Проходим по столбцам
                sum += Math.abs(this.data[i][j]); // Суммируем абсолютные значения элементов в строке
            }
            infinityNorma.push(sum);
            sum = 0;
        }
        return Matrix.SelectionSort(new Matrix([infinityNorma]))[0].data[0][0];
    }

    /**
     * Solves the determinant of the source matrix by Gauss method with selection of the main element.
     * @description Algorithmic complexity: O(n^3), where n is the dimensionality of the square matrix. 
     * @param {Number} epsilon - Error magnitude. Default: 1e-6.
     * @returns {Number} Determinant value.
     */
    Determinant(epsilon = 1e-6) {
        if (!this || this.data.length === 0 || this.data[0].length !== this.data.length) {
            Matrix.errors.determinant();
            return;
        }

        var determinant = 1;
        for (var i = 0; i < this.data.length; i++) {
            var k = i;
            for (var j = i+1; j < this.data[0].length; ++j) {
                if (Math.abs(this.data[j][i]) > Math.abs(this.data[k][i])) {
                    k = j;
                }
            }

            if (Math.abs(this.data[k][i]) < epsilon) {
                determinant = 0;
                break;
            }

            [this.data[i], this.data[k]] = [this.data[k], this.data[i]];

            if (i != k) {
                determinant = -determinant;
            }

            determinant *= this.data[i][i];

            for (var j = i+1; j < this.data[0].length; ++j) {
                this.data[i][j] /= this.data[i][i];
            }

            for (var j = 0; j < this.data[0].length; ++j) {
                if (j != i && Math.abs(this.data[j][i]) > epsilon) {
                    for (var k = i+1; k < this.data[0].length; ++k) {
                        this.data[j][k] -= this.data[i][k] * this.data[j][i];
                    }
                }
            }
        }
        return determinant;
    }

    /**
    * Solves the system of equations AX = B with respect to X. Successive Over-Relaxation (SOR) method.
    * @description Algorithmic complexity: O(k * n^2), where k - amount of iterations, n - size of the matrix A.
    * @param {Matrix|Number[][]} matrix - Right-hand side (RHS).
    * @param {Number} epsilon - Error magnitude. Default: 1e-6.
    * @param {Number} omega - Relaxation parameter in the range 0..1. Default: 0.5.
    * @param {Number} max_iterations - Maximum number of iterations. Default: 1000.
    * @returns {Matrix|Number[][]} Vector of roots X.
    */
    Roots(matrix, epsilon = 1e-6, omega = 0.5, maxIterations = 1000) {
        if (!this || !matrix || this.data.length === 0 || matrix.data.length === 0 || this.data[0].length !== matrix.data.length) {
            Matrix.errors.multiply();
            return;
        }

        var x0 = [];
        var xTau = [];
        var x = [];
        var error = [];
        var k = 0;
        while (k <= maxIterations) {
            for (var i = 0; i < matrix.data.length; i++) {
                x0[i] = x[i] || [0];
                xTau[i] = x[i] || [0];
                x[i] = [0];
            }

            for (var i = 0; i < matrix.data.length; i++) {
                var sum = 0;
                for (var j = 0; j < this.data.length; j++) {
                    if (i !== j) {
                        sum += this.data[i][j] * xTau[j];
                    }
                }

                xTau[i] = 1 / this.data[i][i] * (matrix.data[i][0] - sum);

                for (var j = 0; j < this.data.length; j++) {
                    x[j][0] = omega * xTau[j] + (1 - omega) * x0[j][0];
                }
                error.push(Math.abs(x[i][0] - x0[i][0]));
            }

            error = Matrix.SelectionSort(new Matrix([error]))[0].data[0][matrix.data.length - 1];

            if (error < epsilon) {
                return new Matrix(x);
            }

            error = [];
            k++
        }
        Matrix.errors.iterations();
        return;
    }

    /**
     * Computes the inverse matrix using Schultz's iterative method.
     * @description Algorithmic complexity: O(k * n^3), where k - amount of iterations, n - size of the matrix A.
     * @param {Number} epsilon - Error magnitude. Default: 1e-6.
     * @param {Number} m - Order of convergence of the method. Default 3.
     * @param {Number} max_iterations - Maximum number of iterations. Default: 1000.
     * @returns {Matrix|Number[][]} Inversed matrix.
     */
    Inverse(epsilon = 1e-6, m = 3, max_iterations = 1000) {
        var e = Matrix.Identity(this.data.length);
        var aNorms = this.FirstNorma() * this.InfinityNorma();
        var aTransposed = this.Transpose();
        var u0 = [0];
        
        for (var i = 0; i < this.data.length; i++) {
            u0[i] = [];
            for (var j = 0; j < this.data.length; j++) {
                u0[i][j] = aTransposed.data[i][j] / aNorms;
            }
        }

        u0 = new Matrix(u0);
        
        var k = 0;
        while (k <= max_iterations) {
            var psi = e.Subtraction(this.MultiplyMatrix(u0));
            var psiNorma = 0;
            for (var i = 0; i < psi.data.length; i++) {
                for (var j = 0; j < psi.data.length; j++) {
                    psiNorma += Math.abs(psi.data[i][j])**m;
                }
            }

            psiNorma = Math.sqrt(psiNorma);

            if (psiNorma <= epsilon) {
                return u0;
            } else {
                var psiPrevious = psi;
                var psiM;
                for (var i = 1; i <= m-1; i++) {
                    psiM = psiPrevious.Degree(i);
                    psiM = psiPrevious.Addition(psiM);
                }
                var u1 = u0.MultiplyMatrix(e.Addition(psiPrevious));
                u0 = u1;
                k++
            }
        }
        Matrix.errors.iterations();
        return;
    }

    /**
     * Multiply matrix by vector.
     * @description Algorithmic complexity: O(1).
     * @returns {Vertex} Vector
     */
    MultiplyVector(vertex) {
        var result = [];
        var vector = [vertex.x, vertex.y, vertex.z, vertex.w];

        for (var i = 0; i < this.data.length; i++) {
            result[i] = 0;
            for (var j = 0; j < this.data.length; j++) {
                result[i] += this.data[i][j]*vector[j];
            }
        }
        return new Vertex(result[0], result[1], result[2], result[3], vertex.color);
    }
}

/**
 * Represents a canvas element with methods to manipulate pixels.
 * @constructor
 * @param {HTMLCanvasElement} canvas - The canvas element to use. If not provided, a new canvas will be created.    
 */
class Canvas {
    /**
     * Creates an instance of Canvas.
     * If no canvas element is provided, a new canvas element is created and appended to the body.
     * @constructor
     * @param {HTMLCanvasElement} canvas - The canvas element to use. If not provided, a new canvas will be created.
     */
    constructor(canvas) {
        if (!canvas) {
            this.canvas = document.createElement('canvas');
            this.canvas.id = "Canvas";
            document.body.appendChild(this.canvas);
        } else {
            this.canvas = canvas;
        }
        this.context = this.canvas.getContext('2d');
        this.canvas.width = document.documentElement.clientWidth;
        this.canvas.height = document.documentElement.clientHeight;
        this.buffer = this.context.getImageData(0, 0, this.canvas.width, this.canvas.height);
        this.step = this.buffer.width * 4;
        this.context.imageSmoothingEnabled = false;
        this.dpi = 1/2.54*this.context.miterLimit;
        this.ratio = {
            x: this.canvas.width > this.canvas.height ? this.canvas.width/this.canvas.height : this.canvas.height/this.canvas.width,
            y: this.canvas.width > this.canvas.height ? this.canvas.height/this.canvas.width : this.canvas.width/this.canvas.height
        }
        this.multiplyer = {
            mm: this.dpi/this.canvas.height,
        }
    }

    /**
     * Clear the canvas.
     * @description Algorithmic complexity: O(1).
     */
    Clear() {
        var imageData = this.context.createImageData(this.canvas.width, this.canvas.height);
        var data = new Uint32Array(imageData.data.buffer);
        data.fill(0);
        this.buffer = imageData;
        this.reserveBuffer = imageData;
    }

    /**
     * Insert the ImageData object to the canvas.
     * @description Algorithmic complexity: O(1).
     */
    Update() { 
        if (this.reserveBuffer) {
            this.context.putImageData(this.reserveBuffer, 0, 0);
        } else {
            console.error("Reserve buffer is not set.");
        }
    }

    /**
    * Inserts a pixel of the desired color into the canvas.
    * @description Algorithmic complexity: O(1).
    * @param {Number} x - coordinate x on the canvas.
    * @param {Number} y - coordinate y on the canvas.
    * @param {Number[]} Color - [R, G, B, A] or [R, G, B] for full opacity. Default - black.
    */
    PutPixel(x, y, color = [0,0,0,255]) {
        x = (this.canvas.width/2 + x)|0;
        y = (this.canvas.height/2 - y - 1)|0;

        if (x < 0 ||  y < 0 || x >= this.canvas.width || y >= this.canvas.height) {
            return;
        }

        var offset = 4*x + this.step*y;

        this.buffer.data[offset++] = color[0];
        this.buffer.data[offset++] = color[1];
        this.buffer.data[offset++] = color[2];
        this.buffer.data[offset++] = color[3] || 255;
    }
}

/**
 * The Draw class provides methods for drawing on a specified canvas.
 */
class Draw {
    /**
     * Creates a new instance of Draw.
     * @constructor
     * @param {HTMLCanvasElement} canvas - The canvas element to use. If not provided, a new canvas will be created.
     */
    constructor(canvas) {
        this.canvas = canvas;
    }

    /**
     * A static class providing methods for drawing lines.
     */
    static Line = class {
        /**
         * Creates a new instance of Line.
         * @param {Canvas} canvas - The canvas to draw on.
         */
        constructor(canvas) {
            this.canvas = canvas;
        }

        /**
         * Draws a line between two vertices with solid color filling.
         * @param {Vertex} vertex0 - The starting vertex.
         * @param {Vertex} vertex1 - The ending vertex.
         */
        static Filled(vertex0, vertex1, color = vertex0.color) {
            if (Math.abs(vertex1.x - vertex0.x) > Math.abs(vertex1.y - vertex0.y)) {
                if (vertex0.x > vertex1.x) {
                    [vertex0, vertex1] = [vertex1, vertex0];
                }
    
                var yi = Matrix.Interpolation(vertex0.x, vertex0.y, vertex1.x, vertex1.y);
                
                for (var x = vertex0.x; x <= vertex1.x; x++) {
                    canvas.PutPixel(x, yi[(x - vertex0.x) | 0], color);
                }
            } else {
                if (vertex0.y > vertex1.y) {
                    [vertex0, vertex1] = [vertex1, vertex0];
                }
    
                var xi = Matrix.Interpolation(vertex0.y, vertex0.x, vertex1.y, vertex1.x);
                
                for (var y = vertex0.y; y <= vertex1.y; y++) {
                    canvas.PutPixel(xi[(y - vertex0.y) | 0], y, color);
                }
            }
        }

        /**
         * Draws a line between two vertices with color gradient filling.
         * @param {Vertex} vertex0 - The starting vertex.
         * @param {Vertex} vertex1 - The ending vertex.
         */
        static Gradient(vertex0, vertex1) {
            if (Math.abs(vertex1.x - vertex0.x) > Math.abs(vertex1.y - vertex0.y)) {
                if (vertex0.x > vertex1.x) {
                    [vertex0, vertex1] = [vertex1, vertex0];
                }
    
                var yi = Matrix.Interpolation(vertex0.x, vertex0.y, vertex1.x, vertex1.y);
                var hi = Matrix.Interpolation(vertex0.x, vertex0.w, vertex1.x, vertex1.w);
                var ri = Matrix.Interpolation(vertex0.x, vertex0.color[0], vertex1.x, vertex1.color[0]);
                var gi = Matrix.Interpolation(vertex0.x, vertex0.color[1], vertex1.x, vertex1.color[1]);
                var bi = Matrix.Interpolation(vertex0.x, vertex0.color[2], vertex1.x, vertex1.color[2]);
                var ai = Matrix.Interpolation(vertex0.x, vertex0.color[3], vertex1.x, vertex1.color[3]);

                for (var x = vertex0.x; x <= vertex1.x; x++) {
                    var hue = hi[(x - vertex0.x) | 0];
                    var color = [
                        ri[(x - vertex0.x) | 0]*hue,
                        gi[(x - vertex0.x) | 0]*hue,
                        bi[(x - vertex0.x) | 0]*hue,
                        ai[(x - vertex0.x) | 0]
                    ]
                    canvas.PutPixel(x, yi[(x - vertex0.x) | 0], color);
                }
            } else {
                if (vertex0.y > vertex1.y) {
                    [vertex0, vertex1] = [vertex1, vertex0];
                }
    
                var xi = Matrix.Interpolation(vertex0.y, vertex0.x, vertex1.y, vertex1.x);
                var hi = Matrix.Interpolation(vertex0.y, vertex0.w, vertex1.y, vertex1.w);
                var ri = Matrix.Interpolation(vertex0.y, vertex0.color[0], vertex1.y, vertex1.color[0]);
                var gi = Matrix.Interpolation(vertex0.y, vertex0.color[1], vertex1.y, vertex1.color[1]);
                var bi = Matrix.Interpolation(vertex0.y, vertex0.color[2], vertex1.y, vertex1.color[2]);
                var ai = Matrix.Interpolation(vertex0.y, vertex0.color[3], vertex1.y, vertex1.color[3]);

                
                for (var y = vertex0.y; y <= vertex1.y; y++) {
                    var hue = hi[(y - vertex0.y) | 0];
                    var color = [
                        ri[(y - vertex0.y) | 0]*hue,
                        gi[(y - vertex0.y) | 0]*hue,
                        bi[(y - vertex0.y) | 0]*hue,
                        ai[(y - vertex0.y) | 0]
                    ]
                    canvas.PutPixel(xi[(y - vertex0.y) | 0], y, color);
                }
            }
        }
    }
    
    /**
     * A static class providing methods for drawing triangles.
     */
    static Triangle = class {
        /**
         * Draws a line between two vertices with solid color filling.
         * @param {Vertex} vertex0 - The starting vertex.
         * @param {Vertex} vertex1 - The ending vertex.
         */
        constructor(canvas) {
            this.canvas = canvas;
        }

        /**
         * Draws the wireframe of a triangle by connecting its vertices with lines.
         * @param {Vertex} vertex0 - The first vertex.
         * @param {Vertex} vertex1 - The second vertex.
         * @param {Vertex} vertex2 - The third vertex.
         */
        static Wireframe(vertex0, vertex1, vertex2, color = vertex0.color) {
            Draw.Line.Filled(vertex0, vertex1, color);
            Draw.Line.Filled(vertex1, vertex2, color);
            Draw.Line.Filled(vertex2, vertex0, color);
        }

        /**
         * Draws a filled triangle.
         * @param {Vertex} vertex0 - The first vertex.
         * @param {Vertex} vertex1 - The second vertex.
         * @param {Vertex} vertex2 - The third vertex.
         */
        static Filled(vertex0, vertex1, vertex2, color = vertex0.color) {
            if (vertex1.y < vertex0.y) {[vertex0,vertex1] = [vertex1,vertex0];}
            if (vertex2.y < vertex0.y) {[vertex0,vertex2] = [vertex2,vertex0];}
            if (vertex2.y < vertex1.y) {[vertex1,vertex2] = [vertex2,vertex1];}
            
            var x01 = Matrix.Interpolation(vertex0.y, vertex0.x, vertex1.y, vertex1.x);
            var x12 = Matrix.Interpolation(vertex1.y, vertex1.x, vertex2.y, vertex2.x);
            var x02 = Matrix.Interpolation(vertex0.y, vertex0.x, vertex2.y, vertex2.x);
            
            x01.pop();
            var x012 = x01.concat(x12);
            var xMiddle = (x012.length/2)|0;
            
            if (x02[xMiddle] < x012[xMiddle]) {
                var xLeft = x02;
                var xRight = x012;
            } else {
                var xLeft = x012;
                var xRight = x02;
            }

            for (var y = vertex0.y; y <= vertex2.y; y++) {
                var xl = xLeft[y-vertex0.y]|0;
                var xRightCurrent = xRight[y-vertex0.y]|0;
                for (var x = xl; x <= xRightCurrent; x++) {
                    canvas.PutPixel(x,y,color);
                }
            }
        }

        /**
         * Draws a shaded triangle with color interpolation.
         * @param {Vertex} vertex0 - The first vertex.
         * @param {Vertex} vertex1 - The second vertex.
         * @param {Vertex} vertex2 - The third vertex.
         */
        static Shaded(vertex0, vertex1, vertex2, color = vertex0.color) {
            if (vertex1.y < vertex0.y) {[vertex0,vertex1] = [vertex1,vertex0];}
            if (vertex2.y < vertex0.y) {[vertex0,vertex2] = [vertex2,vertex0];}
            if (vertex2.y < vertex1.y) {[vertex1,vertex2] = [vertex2,vertex1];}

            var x01 = Matrix.Interpolation(vertex0.y, vertex0.x, vertex1.y, vertex1.x);
            var h01 = Matrix.Interpolation(vertex0.y, vertex0.w, vertex1.y, vertex1.w);

            var x12 = Matrix.Interpolation(vertex1.y, vertex1.x, vertex2.y, vertex2.x);
            var h12 = Matrix.Interpolation(vertex1.y, vertex1.w, vertex2.y, vertex2.w);

            var x02 = Matrix.Interpolation(vertex0.y, vertex0.x, vertex2.y, vertex2.x);
            var h02 = Matrix.Interpolation(vertex0.y, vertex0.w, vertex2.y, vertex2.w);

            x01.pop();
            h01.pop();

            var x012 = x01.concat(x12);
            var h012 = h01.concat(h12);

            var xMiddle = (x012.length/2)|0;

            if (x02[xMiddle] < x012[xMiddle]) {
                var xLeft = x02; var xRight = x012;
                var hLeft = h02; var hRight = h012;
            } else {
                var xLeft = x012; var xRight = x02;
                var hLeft = h012; var hRight = h02;
            }

            for (var y = vertex0.y; y <= vertex2.y; y++) {
                var xLeftCurrent = xLeft[y-vertex0.y]|0;
                var xRightCurrent = xRight[y-vertex0.y]|0;
                var h_segment = Matrix.Interpolation(xLeftCurrent, hLeft[y-vertex0.y], xRightCurrent, hRight[y-vertex0.y]);
                for (var x = xLeftCurrent; x <= xRightCurrent; x++) {
                    var hue = h_segment[x-xLeftCurrent];
                    var shaded_color = [
                        color[0]*hue,
                        color[1]*hue,
                        color[2]*hue,
                        255
                    ]; 
                    canvas.PutPixel(x,y,shaded_color);
                }
            }
        }

        /**
         * Draws a triangle with color gradient filling.
         * @param {Vertex} vertex0 - The first vertex.
         * @param {Vertex} vertex1 - The second vertex.
         * @param {Vertex} vertex2 - The third vertex.
         */
        static Gradient(vertex0, vertex1, vertex2) {
            if (vertex1.y < vertex0.y) {[vertex0,vertex1] = [vertex1,vertex0];}
            if (vertex2.y < vertex0.y) {[vertex0,vertex2] = [vertex2,vertex0];}
            if (vertex2.y < vertex1.y) {[vertex1,vertex2] = [vertex2,vertex1];}

            var x01 = Matrix.Interpolation(vertex0.y, vertex0.x, vertex1.y, vertex1.x);
            var h01 = Matrix.Interpolation(vertex0.y, vertex0.w, vertex1.y, vertex1.w);
            var r01 = Matrix.Interpolation(vertex0.y, vertex0.color[0], vertex1.y, vertex1.color[0]);
            var g01 = Matrix.Interpolation(vertex0.y, vertex0.color[1], vertex1.y, vertex1.color[1]);
            var b01 = Matrix.Interpolation(vertex0.y, vertex0.color[2], vertex1.y, vertex1.color[2]);
            var a01 = Matrix.Interpolation(vertex0.y, vertex0.color[3], vertex1.y, vertex1.color[3]);

            var x12 = Matrix.Interpolation(vertex1.y, vertex1.x, vertex2.y, vertex2.x);
            var h12 = Matrix.Interpolation(vertex1.y, vertex1.w, vertex2.y, vertex2.w);
            var r12 = Matrix.Interpolation(vertex1.y, vertex1.color[0], vertex2.y, vertex2.color[0]);
            var g12 = Matrix.Interpolation(vertex1.y, vertex1.color[1], vertex2.y, vertex2.color[1]);
            var b12 = Matrix.Interpolation(vertex1.y, vertex1.color[2], vertex2.y, vertex2.color[2]);
            var a12 = Matrix.Interpolation(vertex1.y, vertex1.color[3], vertex2.y, vertex2.color[3]);

            var x02 = Matrix.Interpolation(vertex0.y, vertex0.x, vertex2.y, vertex2.x);
            var h02 = Matrix.Interpolation(vertex0.y, vertex0.w, vertex2.y, vertex2.w);
            var r02 = Matrix.Interpolation(vertex0.y, vertex0.color[0], vertex2.y, vertex2.color[0]);
            var g02 = Matrix.Interpolation(vertex0.y, vertex0.color[1], vertex2.y, vertex2.color[1]);
            var b02 = Matrix.Interpolation(vertex0.y, vertex0.color[2], vertex2.y, vertex2.color[2]);
            var a02 = Matrix.Interpolation(vertex0.y, vertex0.color[3], vertex2.y, vertex2.color[3]);

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

            var xMiddle = (x012.length/2)|0;

            if (x02[xMiddle] < x012[xMiddle]) {
                var xLeft = x02; var xRight = x012;
                var hLeft = h02; var hRight = h012;
                var rLeft = r02; var rRight = r012;
                var gLeft = g02; var gRight = g012;
                var bLeft = b02; var bRight = b012;
                var aLeft = a02; var aRight = a012;
            } else {
                var xLeft = x012; var xRight = x02;
                var hLeft = h012; var hRight = h02;
                var rLeft = r012; var rRight = r02;
                var gLeft = g012; var gRight = g02;
                var bLeft = b012; var bRight = b02;
                var aLeft = a012; var aRight = a02;
            }

            for (var y = vertex0.y; y <= vertex2.y; y++) {
                var xLeftCurrent = xLeft[y-vertex0.y]|0;
                var xRightCurrent = xRight[y-vertex0.y]|0;
                var hSegment = Matrix.Interpolation(xLeftCurrent, hLeft[y-vertex0.y], xRightCurrent, hRight[y-vertex0.y]);
                var rSegment = Matrix.Interpolation(xLeftCurrent, rLeft[y-vertex0.y], xRightCurrent, rRight[y-vertex0.y]);
                var gSegment = Matrix.Interpolation(xLeftCurrent, gLeft[y-vertex0.y], xRightCurrent, gRight[y-vertex0.y]);
                var bSegment = Matrix.Interpolation(xLeftCurrent, bLeft[y-vertex0.y], xRightCurrent, bRight[y-vertex0.y]);
                var aSegment = Matrix.Interpolation(xLeftCurrent, aLeft[y-vertex0.y], xRightCurrent, aRight[y-vertex0.y]);
                for (var x = xLeftCurrent; x <= xRightCurrent; x++) {
                    var hue = hSegment[x-xLeftCurrent];
                    var r = rSegment[x-xLeftCurrent];
                    var g = gSegment[x-xLeftCurrent];
                    var b = bSegment[x-xLeftCurrent];
                    var a = aSegment[x-xLeftCurrent];
                    var shaded_color = [
                        r*hue,
                        g*hue,
                        b*hue,
                        a
                    ]; 
                    canvas.PutPixel(x,y,shaded_color);
                }
            }
        }
    }
}

class Triangle {
    constructor(vertex0, vertex1, vertex2, color = vertex0.color) {
        this.vertex0 = vertex0;
        this.vertex1 = vertex1;
        this.vertex2 = vertex2;
        this.color = color;
    }
}

class Polygon {
    constructor(data) {
        this.data = data;
    }
}

/**
 * Model representation
 * @param {Vertex[]} vertices - array of the vertices.
 * @param {Triangle[]} triangles - array of the triangles.
 */
class Model {
    constructor(vertices, polygons) {
        this.vertices = vertices;
        this.triangles = polygons;
    }
}

/**
 * @param {Model} model - vertices and triangles.
 * @param {Vertex} position - model center coordinates.
 * @param {Matrix} orientation - rotation of the model.
 * @param {Matrix} scale - scale matrix.
 */
class Instance {
    constructor(model, position, orientation, scale) {
        this.model = model;
        this.position = position || Matrix.Identity(4);
        this.orientation = orientation || Matrix.Identity(4);
        this.scale = scale || Matrix.Scale();
        this.transform = Matrix.Translation(this.position).MultiplyMatrix(this.orientation.MultiplyMatrix(this.scale));
    }
}

/**
 * @param {Vertex} position - model center coordinates.
 * @param {Matrix} orientation - rotation of the model.
 */
class Camera {
    constructor(position, orientation) {
        this.position = position;
        this.orientation = orientation;
    }
}

class Render {
    constructor(canvas) {
        this.canvas = canvas;
        this.canvas.width = document.documentElement.clientWidth;
        this.canvas.height = document.documentElement.clientHeight;
        this.dpi = 1/2.54*this.canvas.context.miterLimit;
        this.ratio = {
            x: this.canvas.width > this.canvas.height ? this.canvas.width/this.canvas.height : this.canvas.height/this.canvas.width,
            y: this.canvas.width > this.canvas.height ? this.canvas.height/this.canvas.width : this.canvas.width/this.canvas.height
        }
        this.multiplyer = {
            mm: this.dpi/this.canvas.height,
        }
    }
    
    ViewportToCanvas(vertex) {
        if (this.canvas.width > this.canvas.height) {
            return new Vertex(
                vertex.x*this.canvas.width/this.ratio.x|0,
                vertex.y*this.canvas.height|0,
                vertex.z,
                vertex.w,
                vertex.color);
        } else {
            return new Vertex(
                vertex.x*this.canvas.width*this.ratio.x|0,
                vertex.y*this.canvas.height|0,
                vertex.z,
                vertex.w,
                vertex.color);
        }
    }

    ProjectVertex(vertex) {
        return this.ViewportToCanvas(
            new Vertex(
                vertex.x*projectionPlaneZ/vertex.z, 
                vertex.y*projectionPlaneZ/vertex.z, 
                vertex.z*this.canvas.multiplyer.mm,
                vertex.w,
                vertex.color));
    }

    WireframeTriangle(polygon, projected) {
        Draw.Triangle.Wireframe(
            projected[polygon.vertex0],
            projected[polygon.vertex1],
            projected[polygon.vertex2],
            polygon.color
        )
    }

    Model(model, transform) {
        var projected = [];
        for (var i = 0; i < model.vertices.length; i++) {
            var vertex = model.vertices[i];
            var vertexH = new Vertex(vertex.x, vertex.y, vertex.z);
            projected.push(this.ProjectVertex(transform.MultiplyVector(vertexH)));
        }

        for (var i = 0; i < model.triangles.length; i++) {
            this.WireframeTriangle(model.triangles[i], projected)
        }
    }

    Object(vertices, triangles) {
        var projected = [];
        for (var i = 0; i < vertices.length; i++) {
            projected.push(this.ProjectVertex(vertices[i]))
        }
        for (var i = 0; i < triangles.length; i++) {
            this.WireframeTriangle(model.triangles[i], projected);
        }
    }

    Instance(instanse) {
        var projected = [];
        var model = instanse.model;

        for (var i = 0; i < model.vertices.length; i++) {
            var vertex = model.vertices[i].Addition(instanse.position);
            projected.push(this.ProjectVertex(vertex))
        }

        for (var i = 0; i < model.triangles.length; i++) {
            this.WireframeTriangle(model.triangles[i], projected);
        }
    }

    Scene(camera, instances) {
        var cameraMatrix = camera.orientation.Transpose().MultiplyMatrix(Matrix.Translation(camera.position.MultiplyScalar(-1)));
        for (var i = 0; i < instances.length; i++) {
            var Transform = cameraMatrix.MultiplyMatrix(instances[i].transform);
            this.Model(instances[i].model, Transform);
        } 
    }

    Cloner(model, offsets, counts) {
        var instances = [];
        for (var x = 0; x < counts.x; x++) {
            for (var y = 0; y < counts.y; y++) {
                for (var z = 0; z < counts.z; z++) {
                    instances.push(new Instance(model, new Vertex(x*offsets.x, y*offsets.y, z*offsets.z)));
                }
            }
        }
        return instances;
    }
}

class Mesh {
    static Triangle(r) {
        var xa = r*Math.cos(210*Math.PI/180);
        var ya = r*Math.sin(210*Math.PI/180);
    
        var xb = 0;
        var yb = r;
    
        var xc = r*Math.cos(30*Math.PI/180);
        var yc = -r*Math.sin(30*Math.PI/180);
    
        var shift = (2*r-(r-yc))/2

        var vertices = [
            new Vertex(xa,ya-shift,0),
            new Vertex(xb,yb-shift,0),
            new Vertex(xc,yc-shift,0)
        ]

        var triangles = [
            new Triangle(0,1,2)
        ]
    
        //return [[xa,ya-shift],[xb,yb-shift],[xc,yc-shift]];
        return [vertices, triangles];
    }

    static Cube(x,y,z, color = [0,0,0,255]) {
        x /= 2;
        y /= 2;
        z /= 2;
        return [
            [
                new Vertex(-x,y,-z),
                new Vertex(x,y,-z),
                new Vertex(x,-y,-z),
                new Vertex(-x,-y,-z),
                new Vertex(-x,y,z),
                new Vertex(x,y,z),
                new Vertex(x,-y,z),
                new Vertex(-x,-y,z),
            ],
            [
                new Triangle(0,1,2,color),
                new Triangle(0,2,3,color),
                new Triangle(0,1,5,color),
                new Triangle(0,4,5,color),
                new Triangle(1,2,5,color),
                new Triangle(2,5,6,color),
                new Triangle(3,2,6,color),
                new Triangle(3,7,6,color),
                new Triangle(0,3,4,color),
                // new Triangle(3,5,7,color),
                new Triangle(4,5,6,color),
                new Triangle(4,6,7,color),
            ]
        ]
    }

    static Icosahedron(subdivisions, color = [0,0,0,255]) {
        var phi = (1 + Math.sqrt(5))/2;
        var triangles = 4**subdivisions;
        var sidePoints = 2**subdivisions+1
        var points = 0;

        var vertices = [
            new Vertex(-1, phi, 0),
            new Vertex(1, phi, 0),
            new Vertex(-1,-phi, 0),
            new Vertex(1, -phi, 0),

            new Vertex(0, -1, phi),
            new Vertex(0, 1, phi),
            new Vertex(0, -1, -phi),
            new Vertex(0, 1, -phi),

            new Vertex(phi, 0, -1),
            new Vertex(phi, 0, 1),
            new Vertex(-phi, 0, -1),
            new Vertex(-phi, 0, 1),
        ];

        var triangles = [
            new Triangle(0,11,5,color),
            new Triangle(0,5,1,color),
            new Triangle(0,1,7,color),
            new Triangle(0,7,10,color),
            new Triangle(0,10,11,color),

            new Triangle(1,5,9,color),
            new Triangle(5,11,4,color),
            new Triangle(11,10,2,color),
            new Triangle(10,7,6,color),
            new Triangle(7,1,8,color),

            new Triangle(3,9,4,color),
            new Triangle(3,4,2,color),
            new Triangle(3,2,6,color),
            new Triangle(3,6,8,color),
            new Triangle(3,8,9,color),

            new Triangle(4,9,5,color),
            new Triangle(2,4,11,color),
            new Triangle(6,2,10,color),
            new Triangle(8,6,7,color),
            new Triangle(9,8,1,color),
        ];

        return [vertices, triangles];
    }

    static Star(n, innerDiameter, outerDiameter, frontThickness=innerDiameter/8, backThickness=innerDiameter/4) {
        var outerRadius = outerDiameter/2;
        var innerRadius = innerDiameter/2;
        if (n%2 === 0) {
            var startAngle = (360/n-90)+(360/n/2);
        } else {
            var startAngle = (360/n-90);
        }
        
        var outerPoints = [];
        var innerPoints = [];
        var vertices = [];
        var triangles = [];
        var x; var y; var angle;
    
        for (var i = 0; i < n; i++) {
            angle = (startAngle+(360/n)*i)*Math.PI/180;
            x = Math.cos(angle)*outerRadius;
            y = Math.sin(angle)*outerRadius;
            outerPoints.push([x, y]);
        }
    
        for (var i = 0; i < n; i++) {
            angle = (startAngle+(360/n)*i+(360/n)/2)*Math.PI/180;
            x = Math.cos(angle)*innerRadius;
            y = Math.sin(angle)*innerRadius;
            innerPoints.push([x, y]);
        }
        
        var points = outerPoints.concat(innerPoints);
        var indexes = [];
        var firstIndex = n-1;
        var secondIndex = points.length/2;
    
        for (var i = 0; i < points.length; i++) {
            if (i%2 === 0) {
                firstIndex = n-1;
                secondIndex--;
            } else {
                firstIndex = n-2;
            }
            secondIndex+1;
            if (i == points.length-1) {
                firstIndex = (points.length/2-1)*2;
            }
            indexes.push([firstIndex+1, secondIndex+1])
        }
    
        for (var i = 0; i < points.length; i++) {
            vertices.push(new Vertex(points[i][0], points[i][1], projectionPlaneZ, 1, [255,0,0,255]));
        }
        vertices.push(new Vertex(0,0,projectionPlaneZ+frontThickness,1,[255,0,0,255]), new Vertex(0,0,projectionPlaneZ-backThickness,1,[255,0,0,255]));
    
        var shifts = [0];
        var B = (vertices.length-2)/2;
        for (var i = 0; i < vertices.length-2; i++) {
            if (i%2 === 0) {
                B--;
            } else {
                B+1;
            }
            shifts.unshift(B);   
        }
        shifts.shift()
    
        for (var i = 0; i < indexes.length; i++) {
            triangles.push(
                new Triangle(
                    shifts[i],
                    shifts[i]+indexes[i][0],
                    shifts[i]+indexes[i][0]+indexes[i][1],
                    [255,0,0,255]
                )
            )
            triangles.push(
                new Triangle(
                    shifts[i],
                    shifts[i]+indexes[i][0],
                    shifts[i]+indexes[i][0]+indexes[i][1]+1,
                    [255,0,0,255]
                )
            )
        }
    
        return [vertices, triangles];
    }

    static Plane(width, height, widthSegments, heightSegments, color = [0,0,0,255]) {
        var xStart = -width/2;
        var yStart = height/2;
        var xStep = width/widthSegments;
        var yStep = height/heightSegments;
        var vertices = [];
        var triangles = [];

        for (var y = yStart; y >= -height/2; y -= yStep) {
            for (var x = xStart; x <= width/2; x += xStep) {
                vertices.push(new Vertex(x,y,0,color));
            }
        }
        
        for (var i = 0; i <= heightSegments-1; ++i) {
            var k1 = i*(widthSegments+1);
            var k2 = k1+widthSegments+1;
            for (var j = 0; j < widthSegments; ++j, ++k1, ++k2) {
                if (i%2 == 0 || i%2 == 1) {
                    triangles.push(new Triangle(k1,k2,k1+1));
                }
                triangles.push(new Triangle(k1+1,k2+1,k2));
            }
        }
        return [vertices, triangles];
    }
}


var canvas = new Canvas(document.getElementById('Canvas'));
var render = new Render(canvas);
var camera = new Camera(new Vertex(0,0,-15), Matrix.Rotation(0,0,0))

var cubeMesh = Mesh.Cube(4,4,2,[255,0,0,255]);
var cube = new Instance(new Model(cubeMesh[0],cubeMesh[1]), new Vertex(-4,3,0), Matrix.Rotation(0,20,90), Matrix.Scale(1,1,1));

var icosahedronMesh = Mesh.Icosahedron(1, [0,0,0,255]);
var icosahedron = new Instance(new Model(icosahedronMesh[0],icosahedronMesh[1]), new Vertex(4,3,0), Matrix.Rotation(0,0,50), Matrix.Scale(1,1,1));

var planeMesh = Mesh.Plane(4,4,10,10, [255,0,0,255]);
var planeXY = new Instance(new Model(planeMesh[0],planeMesh[1]), new Vertex(-4,-3,0), Matrix.Rotation(0,20,0), Matrix.Scale(1,1,1));

var starMesh = Mesh.Star(5,2,5,0.5,0.5);
var star = new Instance(new Model(starMesh[0],starMesh[1]), new Vertex(4,-3,0), Matrix.Rotation(0,20,180), Matrix.Scale(1,1,1));


var instances = [
    cube,
    icosahedron,
    planeXY,
    star
]

function RENDER() {
    canvas.Clear();
    render.Scene(camera, instances);
    canvas.Update();
}
RENDER();

var ax = 0;
var ay = 0;
var az = 0;
function animate() {
    ax += 1;
    ay += 1;
    az += 1;
    instances[0].transform = Matrix.Rotation(ax+90,ay,az);
    instances[1].transform = Matrix.Rotation(ax,ay+90,az);
    instances[2].transform = Matrix.Rotation(ax,ay,az+90);
    instances[3].transform = Matrix.Rotation(ax,ay,az+90);
    RENDER();
    requestAnimationFrame(animate);
}
//animate();



