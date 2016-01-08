var getLDPC = (function(options){

	var getUtility = require('./utility');
	var getRandomGenerator = require('./randomGenerator');

	options	= options || {};
	if (typeof options.n == "undefined") {
		options.n = 0;
	}
	if (typeof options.k == "undefined") {
		option.k = options.n;
	}
	if (typeof options.modulo == "undefined") {
		options.modulo = 2;
	}
	if (typeof options.randomSeed == "undefined") {
		options.randomSeed = Date.now();
	}

	var utility = getUtility();
	var random = getRandomGenerator(options.randomSeed);

	// See https://en.wikipedia.org/wiki/Low-density_parity-check_code
	// for explanations of p, -p(transposed), the generator matrix, and
	// the parity matrix
	var invertibleNumbers = [];
	var i = -1;
	for (var i = 1; i < options.modulo; i++) {
		if (utility.multiplicativeInverse(i, options.modulo)) {
			invertibleNumbers.push(i);
		}
	}
	var p = utility.make(options.k, options.n - options.k, function(row, column) {
		if ((row % options.k) == (column % options.k) && column > 0)
			return 0;
		else
			return random.nextFromChoice(invertibleNumbers);
	});
	var negPT = utility.map(utility.transpose(p), function(value) {
		return utility.mod(-value, options.modulo);
	});

	var generator = utility.concatColumns(utility.identity(options.k), p);
	var parity = utility.concatColumns(negPT, utility.identity(options.n - options.k));

	// Wikipedia says that we can check that the row space of G is
	// orthogonal to H by doing this:
	var test = utility.multiply(generator, utility.transpose(parity), options.modulo);
	// Every element of test should be zero
	utility.map(test, function(value) {
		if (value) {
			throw "The generator and parity matrices don't have orthogonal row spaces";
		}
	});

	/**
	 * Attempts to decode the given encoded message. Returns an array of
	 * options.k symbols. The returned array has a "decoded" property that tells
	 * if the encoded message was successfully decoded completely (the same as
	 * there being no erased symbols in the returned result)
	 */
	function decode(encoded) {
		var encoded = utility.deepCopy(encoded);
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
					//sums[j] = utility.mod(sums[j], options.modulo);
				}
			}
		}

		if (variables.length <= parity.length) { // Only continue if we've got a change of solving the system of equations
			for (var j = 0; j < parity.length; j++) {
				matrix[j].push(sums[j]); // The sums will fill the right-hand column
			}

			var reduced = utility.gaussReduction(matrix, options.modulo); // This tries to put the matrix into reduced row-echelon form

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
	function encode(message) {
		return utility.multiply([message], generator, options.modulo)[0];
	};

	/**
	 * Returns a copy of the generator matrix used by this instance of LDPC
	 */
	function getGenerator() {
		return utility.deepCopy(generator);
	};

	/**
	 * Returns a copy of the parity matrix used by this instance of LDPC
	 */
	function getParity() {
		return utility.deepCopy(parity);
	};

	return {
		options: options,
		decode: decode,
		encode: encode,
		getGenerator: getGenerator,
		getParity: getParity
	};

});

module.exports = getLDPC;
