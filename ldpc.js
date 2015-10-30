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
	var p = LDPC.util.make(options.k, options.n - options.k, function() {
		return Math.floor(random.next() * options.modulo);
	});
	var negPT = LDPC.util.fix(LDPC.util.map(LDPC.util.transpose(p), function(value) {
		return -1 * value;
	}), options.modulo);

	var generator = LDPC.util.concatColumns(LDPC.util.identity(options.k), p);
	var parity = LDPC.util.concatColumns(negPT, LDPC.util.identity(options.n - options.k));

	// Wikipedia says that we can check that the row space of G is
	//  orthogonal to H by doing this:
	var test = LDPC.util.fix(LDPC.util.multiply(generator, LDPC.util.transpose(parity)), options.modulo);
	// Every element of test should be zero
	LDPC.util.map(test, function(value) {
		if (value)
			throw "The generator and parity matrices don't have orthogonal row spaces";
	});

	/**
	 * Attempts to decode the given encoded message
	 */
	this.decode = function(encoded) {
		encoded = encoded.slice();
		var resolvedAtLeastOne = true;
		while (resolvedAtLeastOne) {
			resolvedAtLeastOne = false;

			// See if any of the parity constraints can be satisfied
			parity.forEach(function(constraints) {
				var missingSymbolIndices = [];
				var total = 0;
				for (var i = 0; i < constraints.length; i++) {
					var constraint = constraints[i];
					var symbol = encoded[i];
					if (constraint > 0) {
						if (symbol < 0 || symbol == null) {
							// This symbol has been erased / is missing
							missingSymbolIndices.push(i);
						}
						else {
							// This symbol hasn't been erased
							total += symbol * constraint;
						}
					}
				}
				if (missingSymbolIndices.length == 1) {
					var missingSymbolIndex = missingSymbolIndices.pop();
					var constraint = constraints[missingSymbolIndex];
					// We can figure out this one missing symbol. Its value
					//  will make the total add up to zero mod
					//  options.modulo
					var missingValue = (((-total * constraint) % options.modulo) + options.modulo) % options.modulo;
					encoded[missingSymbolIndex] = missingValue;
					resolvedAtLeastOne = true;
				} else if (missingSymbolIndices.length == 0) {
					// There are no missing symbols, so let's check that
					//  that this parity constraint holds true
					if (total % options.modulo)
						throw "Total isn't zero in mod options.modulo";
				}
			});
		}
		var retVal = encoded.slice(0, options.k); // Return the first k symbols
		// Determine if the result is completely decoded yet; .decoded will indicate that
		retVal.decoded = true;
		retVal.forEach(function(value) {
			retVal.decoded &= value > 0 || value === 0;
		});
		retVal.result = encoded; // .result will be the full array of symbols
		return retVal;
	};

	/**
	 * Encodes the given message, returning an array options.n symbols
	 * long
	 */
	this.encode = function(message) {
		return LDPC.util.fix(LDPC.util.multiply([message], generator), options.modulo)[0];
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
		if (typeof value == "object") {
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

// http://stackoverflow.com/questions/27205018/multiply-2-matrices-in-javascript
LDPC.util.multiply = function(matrix1, matrix2) {
    var result = [];
    for (var i = 0; i < matrix1.length; i++) {
        result[i] = [];
        for (var j = 0; j < matrix2[0].length; j++) {
            var sum = 0;
            for (var k = 0; k < matrix1[0].length; k++) {
                sum += matrix1[i][k] * matrix2[k][j];
            }
            result[i][j] = sum;
        }
    }
    return result;
};

// Returns a transposed copy of the given array
LDPC.util.transpose = function(array) {
	var rows = array.length;
	var columns = array[0].length;
	var transposed = [];
	for (var i = 0; i < columns; i++) {
		var row = [];
		for (var j = 0; j < rows; j++) {
			row.push(array[j][i]);
		}
		transposed.push(row);
	}
	return transposed;
};