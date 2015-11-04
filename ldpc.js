function LDPC(options) {
	options	= options || {};
	if (typeof options.n == "undefined")
		options.n = 0;
	if (typeof options.k == "undefined")
		option.k = options.n;
	if (typeof options.modulo == "undefined")
		options.modulo = 2;
	if (typeof options.randomSeed == "undefined")
		options.randomSeed = Date.now();

	var random = new LDPC.Random(options.randomSeed);

	// See https://en.wikipedia.org/wiki/Low-density_parity-check_code
	//  for explanations of p, -p(transposed), the generator matrix, and
	//  the parity matrix
	var invertibleNumbers = [];
	var i = -1;
	for (var i = 1; i < options.modulo; i++) {
		if (LDPC.util.multiplicativeInverse(i, options.modulo))
			invertibleNumbers.push(i);
	}
	var p = LDPC.util.make(options.k, options.n - options.k, function(row, column) {
		if ((row % options.k) == (column % options.k) && column > 0)
			return 0;
		else
			return random.nextFromChoice(invertibleNumbers);
	});
	var negPT = LDPC.util.map(LDPC.util.transpose(p), function(value) {
		return LDPC.util.mod(-value, options.modulo);
	});

	var generator = LDPC.util.concatColumns(LDPC.util.identity(options.k), p);
	var parity = LDPC.util.concatColumns(negPT, LDPC.util.identity(options.n - options.k));

	// Wikipedia says that we can check that the row space of G is
	//  orthogonal to H by doing this:
	var test = LDPC.util.multiply(generator, LDPC.util.transpose(parity), options.modulo);
	// Every element of test should be zero
	LDPC.util.map(test, function(value) {
		if (value)
			throw "The generator and parity matrices don't have orthogonal row spaces";
	});

	/**
	 * Attempts to decode the given encoded message. Returns an array of
	 * options.k symbols. The returned array has a "decoded" property that tells
	 * if the encoded message was successfully decoded completely (the same as
	 * there being no erased symbols in the returned result)
	 */
	this.decode = function(encoded) {
		var encoded = LDPC.util.deepCopy(encoded);
		var matrix = [];
		for (var i = 0; i < parity.length; i++) {
			matrix.push([]);
		}
		var sums = [];
		var variables = [];

		for (var i = 0; i < encoded.length; i++) {
			var symbol = encoded[i];
			if (symbol == null || symbol < 0) {
				// This symbol is missing. Therefore we need to consider it
				//  as a variable to solve for
				variables.push(i);
				if (variables.length > parity.length)
					break; // No point in continuing because we won't be able to solve this system of equations completely
				for (var j = 0; j < parity.length; j++) {
					var parityRow = parity[j];
					matrix[j].push(parityRow[i]);
				}
			} else {
				// This symbol is present, so it should contribute toward
				//  the sum we're trying to solve for
				for (var j = 0; j < parity.length; j++) {
					var parityRow = parity[j];
					sums[j] = sums[j] || 0;
					sums[j] -= symbol * parityRow[i];
					//sums[j] = LDPC.util.mod(sums[j], options.modulo);
				}
			}
		}

		if (variables.length <= parity.length) { // Only continue if we've got a change of solving the system of equations
			for (var j = 0; j < parity.length; j++) {
				matrix[j].push(sums[j]); // The sums will fill the right-hand column
			}

			var reduced = LDPC.util.gaussReduction(matrix, options.modulo); // This tries to put the matrix into reduced row-echelon form

			for (var i = 0; i < reduced.length; i++) {
				var reducedRow = reduced[i];
				// If there's exactly one "1" and the rest "0" (except for the
				//  last element), then the last element is the answer for the
				//  variable indicated by the field with the "1" in it
				var variableIndex = -1;
				for (var j = 0; j < reducedRow.length - 1; j++) {
					var element = reducedRow[j];
					if (element == 1 && variableIndex == -1) {
						// We haven't yet found a "1" element
						variableIndex = j;
					} else if (element) {
						// This row doesn't tell us anything conclusively
						variableIndex = -1;
						break;
					}
				}
				if (variableIndex != -1) {
					// We found an answer: the variableIndex'th variable equals the last element in this row
					var answer = reducedRow[reducedRow.length - 1];
					//console.log("variable " + variableIndex + "=" + answer);
					encoded[variables[variableIndex]] = answer;
				}
			}
		}
		var retVal = encoded.slice(0, options.k);
		retVal.decoded = true;
		retVal.forEach(function(value) {
			if (value == null || value < 0) {
				retVal.decoded = false;
			}
		})
		retVal.all = encoded;

		return retVal;
	}

	/**
	 * Encodes the given message, returning an array options.n symbols
	 * long
	 */
	this.encode = function(message) {
		return LDPC.util.multiply([message], generator, options.modulo)[0];
	};

	/**
	 * Returns a copy of the generator matrix used by this instance of LDPC
	 */
	this.getGenerator = function() {
		return LDPC.util.deepCopy(generator);
	};

	/**
	 * Returns a copy of the parity matrix used by this instance of LDPC
	 */
	this.getParity = function() {
		return LDPC.util.deepCopy(parity);
	};
}

/**
 * A simple RNG per http://stackoverflow.com/questions/3062746/special-simple-random-number-generator
 * Note this will only generate numbers up to 2^32 - 1
 */
LDPC.Random = function(seed) {
	if (typeof seed == "undefined")
		seed = Date.now(); // Seed with current time if no seed given
	var a = 1103515245, c = 12345, m = Math.pow(2, 32);
	/**
	 * Generates a number between 0 and 1 (with 32-bit precision)
	 */
	this.next = function() {
		seed = (a * seed + c) % m;
		return seed / m;
	};

	/**
	 * Picks a random choice from the given choices
	 */
	this.nextFromChoice = function(choices) {
		return choices[Math.floor(this.next() * choices.length)];
	};

	/**
	 * Generates a normal number (mean of zero) with given standard
	 * deviation
	 */
	this.nextGaussian = function(standardDeviation) {
		var u1 = this.next();
		var u2 = this.next();
		var randStdNormal = Math.sqrt(-2 * Math.log(u1)) * Math.sin(2 * Math.PI * u2);
		return standardDeviation * randStdNormal;
	};

	/**
	 * http://wiki.q-researchsoftware.com/wiki/How_to_Generate_Random_Numbers:_Poisson_Distribution
	 */
	this.nextPoisson = function(lambda) {
		var L = Math.exp(-lambda);
		var p = 1.0;
		var k = 0;
		do {
			k++;
			p *= this.next();
		} while (p > L);
		return k - 1;
	};

	/**
	 * Scrambles the list of numbers from 0 to (length-1)
	 */
	this.nextSet = function(length) {
		var set = [];
		for (var i = 0; i < length; i++) {
			set.push(i);
		}
		for (var i = 0; i < set.length; i++) {
			var randomIndex = Math.floor(this.next() * set.length);
			var temp = set[i];
			set[i] = set[randomIndex];
			set[randomIndex] = temp;
		}
		return set;
	};
};

// Set up utility functions
LDPC.util = {};

// Returns the combination of the two given two-dimensional arrays
LDPC.util.concatColumns = function(array1, array2) {
	var copy1 = LDPC.util.deepCopy(array1);
	var copy2 = LDPC.util.deepCopy(array2);

	var concatInner = function(matrix1, matrix2) {
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
LDPC.util.deepCopy = function(matrix) {
	var copy = [];
	matrix.forEach(function(value, index, object) {
		if (value != null && typeof value == "object") {
			copy.push(LDPC.util.deepCopy(value));
		} else {
			copy.push(value);
		}
	});
	return copy;
};

// Returns an array in which all the elements are positive and less then
//  modulo
LDPC.util.fix = function(array, modulo) {
	return LDPC.util.map(array, function(value) {
		return ((Math.round(value) % modulo) + modulo) % modulo; // Ensures positive numbers
	});
};

// Performs Gauss reduction for a finite field, putting it into row-echelon
//  form.
LDPC.util.gaussReduction = function(matrix, modulo) {
	matrix = LDPC.util.fix(matrix, modulo); // Note this makes a copy of the given matrix, and also ensures that all elements are < modulo
	if (!matrix[0])
		return matrix;
	var numColumns = matrix[0].length;
	var numRows = matrix.length;
	var reducedRows = {};
	var scale = function(row, multiple) {
		for (var col = 0; col < row.length; col++) {
			row[col] *= multiple;
			row[col] = LDPC.util.mod(row[col], modulo);
		}
	};
	var subtractRow = function(from, row, multiple) {
		//console.log("subtracting " + multiple + " multiples of " + row.toString() + " from " + from.toString());
		for (var col = 0; col < from.length; col++) {
			from[col] -= multiple * row[col];
			from[col] = LDPC.util.mod(from[col], modulo);
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
				if (!reducedRows[row]) { // We don't want to scale a row that has already been reduced, nor do we want to consider it the minimum row
					// Perhaps we can reduce this row to 1 in this column by using a multiplicative inverse
					var multiplicativeInverse = LDPC.util.multiplicativeInverse(matrix[row][col], modulo);
					if (multiplicativeInverse) {
						scale(matrix[row], multiplicativeInverse);
						minimumRow = row;
						//console.log("scaled row " + row + " " + multiplicativeInverse + " times to " + matrix[row].toString() + ", which made it the minimum row");
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
LDPC.util.getRandomSeries = function(length, modulo) {
	var retVal = [];
	for (var i = 0; i < length; i++) {
		retVal.push(Math.floor(Math.random() * modulo));
	}
	return retVal;
};

// Generates an identity matrix
LDPC.util.identity = function(dimension) {
	return LDPC.util.make(dimension, dimension, function(row, column) {
		return row == column ? 1 : 0;
	});
};

// Finds the inverse of a 2x2 matrix in the given modulo
LDPC.util.inverse2x2 = function(matrix, modulo) {
	var a = matrix[0][0];
	var b = matrix[0][1];
	var c = matrix[1][0];
	var d = matrix[1][1];
	var determinantInv = LDPC.util.multiplicativeInverse(
		a * d - b * c,
		modulo);
	return LDPC.util.fix(
		[
			[determinantInv * d, determinantInv * -b],
			[determinantInv * -c, determinantInv * a]
		],
		modulo);
};

// Returns a two-dimensional array with the given number of rows and
//  columns. Each value is calculated based on the result of
//  fn(row, column)
LDPC.util.make = function(rows, columns, fn) {
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
LDPC.util.map = function(array, calculate) {
	// Recursive function for applying the calculation. Note that we
	//  recurse into nested arrays
	var mapInner = function(matrix, coordinates, calculate) {
		coordinates = coordinates || [];
		matrix.forEach(function(value, index, object) {
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

	var copy = LDPC.util.deepCopy(array);
	mapInner(copy, [], calculate);
	return copy;
};

// Takes care of negative numbers as well
LDPC.util.mod = function(x, modulo) {
	return ((x % modulo) + modulo) % modulo;
};

// Lazily-computes (and then caches) the multiplicative inverse of a number
LDPC.util.multiplicativeInverse = function() {
	var cache = {};
	return function(x, modulo) {
		x = LDPC.util.mod(x, modulo);
		var inverses;
		if (!(inverses = cache[modulo] || []).length) {
			for (var i = 1; i < modulo; i++) {
				if (inverses[i] != null){
					continue;
				} else {
					for (var multiplier = 1; multiplier < modulo; multiplier++) {
						if (LDPC.util.mod(i * multiplier, modulo) == 1) {
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
}();

// Multiplies two matrices in the given modulo
LDPC.util.multiply = function(matrix1, matrix2, modulo) {
    var result = [[]];
    if (!matrix2.length)
    	return result;
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
LDPC.util.transpose = function(array) {
	var transposed = [];
	var rows = array.length;
	if (!rows)
		return transposed;
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