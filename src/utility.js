var getUtil = (function () {

    this.cache = {};

    // Returns the combination of the two given two-dimensional arrays
    function concatColumns(array1, array2) {
        var copy1 = deepCopy(array1);
        var copy2 = deepCopy(array2);

        var concatInner = function (matrix1, matrix2) {
            while (matrix2.length) {
                matrix1.push(matrix2.shift());
            }
        };

        for (var i = 0; i < copy1.length; i++) {
            concatInner(copy1[i], copy2[i]);
        }
        return copy1;
    };

    // Perform a deep copy of all elements
    function deepCopy(matrix) {
        var copy = [];
        matrix.forEach(function (value, index, object) {
            if (value != null && typeof value == "object") {
                copy.push(deepCopy(value));
            } else {
                copy.push(value);
            }
        });
        return copy;
    };

    // Returns an array in which all the elements are positive and less then
    //  modulo
    function fix(array, modulo) {
        return map(array, function (value) {
            return ((Math.round(value) % modulo) + modulo) % modulo; // Ensures positive numbers
        });
    };

    // Performs Gauss reduction for a finite field, putting it into row-echelon
    //  form.
    function gaussReduction(matrix, modulo) {
        matrix = fix(matrix, modulo); // Note this makes a copy of the given matrix, and also ensures that all elements are < modulo
        if (!matrix[0]) {
            return matrix;
        }
        var numColumns = matrix[0].length;
        var numRows = matrix.length;
        var reducedRows = {};
        var scale = function (row, multiple) {
            for (var col = 0; col < row.length; col++) {
                row[col] *= multiple;
                row[col] = mod(row[col], modulo);
            }
        };
        var subtractRow = function (from, row, multiple) {
            //console.log("subtracting " + multiple + " multiples of " + row.toString() + " from " + from.toString());
            for (var col = 0; col < from.length; col++) {
                from[col] -= multiple * row[col];
                from[col] = mod(from[col], modulo);
            }
            //console.log("...got " + from.toString());
        };
        for (var col = 0; col < numColumns && col < numRows; col++) {
            //console.warn("working column " + col);
            // Get zeros in this column for all rows but one, which will become a
            //  reduced row
            var numZeros;
            var lastMinimumRow = -1;
            do {
                var minimumRow = -1;
                numZeros = 0;
                // Find the row with the minimum value for this column
                for (var row = 0; row < numRows; row++) {
                    // We don't want to scale a row that has already been reduced, nor do we want to consider it the minimum row
                    if (!reducedRows[row]) {
                        // Perhaps we can reduce this row to 1 in this column by using a multiplicative inverse
                        var multiplicativeInverseValue = multiplicativeInverse(matrix[row][col], modulo);
                        if (multiplicativeInverseValue) {
                            scale(matrix[row], multiplicativeInverseValue);
                            minimumRow = row;
                            //console.log("scaled row " + row + " " + multiplicativeInverseValue + " times to " + matrix[row].toString() + ", which made it the minimum row");
                            break; // Skips the rest of this for-loop
                        }
                        if ((minimumRow == -1 || matrix[row][col] < matrix[minimumRow][col]) && matrix[row][col] != 0) {
                            // Otherwise, see if this row should be considered the minimum row
                            minimumRow = row;
                        }
                    }
                }
                //console.log("minimum row is " + minimumRow);
                if (minimumRow == lastMinimumRow) {
                    // We're going around in circles. Time to stop
                    //console.log("we're being circular because we've already decided that " + minimumRow + " is the minimum row");
                    break;
                } else {
                    lastMinimumRow = minimumRow;
                }
                // Subtract multiples of that row from all the other rows
                for (var row = 0; row < numRows; row++) {
                    if (row != minimumRow && matrix[row][col]) {
                        subtractRow(matrix[row], matrix[minimumRow], Math.floor(matrix[row][col] / matrix[minimumRow][col]));
                    }
                    if (matrix[row][col] == 0) {
                        //console.log("row " + row + " has a zero");
                        numZeros++;
                    }
                }
                //console.log(matrix.toString());
            } while (numZeros < numRows - 1);
            for (var row = 0; row < numRows; row++) {
                if (matrix[row][col]) {
                    //console.log("row " + row + " is reduced");
                    reducedRows[row] = true;
                    break;
                }
            }
        }
        return matrix;
    };

    // Generates a random array of numbers between 0 and the given modulo
    function getRandomSeries(length, modulo) {
        var retVal = [];
        for (var i = 0; i < length; i++) {
            retVal.push(Math.floor(Math.random() * modulo));
        }
        return retVal;
    };

    // Generates an identity matrix
    function identity(dimension) {
        return make(dimension, dimension, function (row, column) {
            return row == column ? 1 : 0;
        });
    };

    // Finds the inverse of a 2x2 matrix in the given modulo
    function inverse2x2(matrix, modulo) {
        var a = matrix[0][0];
        var b = matrix[0][1];
        var c = matrix[1][0];
        var d = matrix[1][1];
        var determinantInv = multiplicativeInverse(a * d - b * c, modulo);
        return fix(
            [
                [determinantInv * d, determinantInv * -b],
                [determinantInv * -c, determinantInv * a]
            ],
            modulo
        );
    };

    // Returns a two-dimensional array with the given number of rows and
    //  columns. Each value is calculated based on the result of
    //  fn(row, column)
    function make(rows, columns, fn) {
        var array = [];
        for (var i = 0; i < rows; i++) {
            var row = [];
            for (var j = 0; j < columns; j++) {
                var element = fn(i, j);
                row.push(element);
            }
            array.push(row);
        }
        return array;
    };

    // Applies a function to each element of the given array (which can be
    //  multi-dimensional), returning this transformed version without
    //  affecting the given array
    function map(array, calculate) {
        // Recursive function for applying the calculation. Note that we
        //  recurse into nested arrays
        var mapInner = function (matrix, coordinates, calculate) {
            coordinates = coordinates || [];
            matrix.forEach(function (value, index, object) {
                var newValue;
                var newCoordinates = coordinates.slice();
                newCoordinates.push(index);
                if (typeof value == "object") {
                    newValue = mapInner(value, newCoordinates, calculate);
                } else {
                    newValue = calculate(value, newCoordinates);
                }
                if (typeof newValue != "undefined") {
                    matrix[index] = newValue;
                }
            });
        };

        var copy = deepCopy(array);
        mapInner(copy, [], calculate);
        return copy;
    };

    // Takes care of negative numbers as well
    function mod(x, modulo) {
        return ((x % modulo) + modulo) % modulo;
    };

    // Lazily-computes (and then caches) the multiplicative inverse of a number
    function multiplicativeInverse(x, modulo) {
        x = mod(x, modulo);
        var inverses;
        if (!(inverses = cache[modulo] || []).length) {
            for (var i = 1; i < modulo; i++) {
                if (inverses[i] != null) {
                    continue;
				} else {
                    for (var multiplier = 1; multiplier < modulo; multiplier++) {
                        if (mod(i * multiplier, modulo) == 1) {
                            // i and multiplier are multiplicative inverses of each other in mod modulo
                            inverses[i] = multiplier;
                            inverses[multiplier] = i;
                            break;
                        }
                    }
                }
            }
            cache[modulo] = inverses;
        }
        return inverses[x];
    };

    function clearCache() {
        this.cache = {};
    };

    // Multiplies two matrices in the given modulo
    function multiply(matrix1, matrix2, modulo) {
        var result = [[]];
        if (!matrix2.length) {
            return result;
        }
        for (var i = 0; i < matrix1.length; i++) {
            result[i] = [];
            for (var j = 0; j < matrix2[0].length; j++) {
                var sum = 0;
                for (var k = 0; k < matrix1[0].length; k++) {
                    sum += matrix1[i][k] * matrix2[k][j];
                    sum %= modulo;
                }
                result[i][j] = sum;
            }
        }
        return result;
    };

    // Returns a transposed copy of the given array
    function transpose(array) {
        var transposed = [];
        var rows = array.length;
        if (!rows) {
            return transposed;
        }
        var columns = array[0].length;
        for (var i = 0; i < columns; i++) {
            var row = [];
            for (var j = 0; j < rows; j++) {
                row.push(array[j][i]);
            }
            transposed.push(row);
        }
        return transposed;
    };

    return {
        cache: cache,
        clearCache: clearCache,
        concatColumns: concatColumns,
        deepCopy: deepCopy,
        fix: gaussReduction,
        gaussReduction: gaussReduction,
        getRandomSeries: getRandomSeries,
        identity: identity,
        inverse2x2: inverse2x2,
        make: make,
        map: map,
        mod: mod,
        multiplicativeInverse: multiplicativeInverse,
        multiply: multiply,
        transpose: transpose
    };
});

module.exports = getUtil;
