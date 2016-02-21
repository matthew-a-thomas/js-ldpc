var LDPC = (function (options) {

    var getUtility = require('./utility');
    var getRandomGenerator = require('./randomGenerator');

    var Options = options || {};
    if (Options.n === undefined) {
        Options.n = 0;
    }
    if (Options.k === undefined) {
        Options.k = Options.n;
    }
    if (Options.modulo === undefined) {
        Options.modulo = 2;
    }
    if (Options.randomSeed === undefined) {
        Options.randomSeed = Date.now();
    }

    var Utility = getUtility();
    var RandomGenerator = getRandomGenerator(Options.randomSeed);

    /*
      See https://en.wikipedia.org/wiki/Low-density_parity-check_code
      for explanations of p, -p(transposed), the generator matrix, and
      the parity matrix
    */
    var invertibleNumbers = [];
    var i = -1;
    for (var i = 1; i < Options.modulo; i++) {
        if (Utility.multiplicativeInverse(i, Options.modulo)) {
            invertibleNumbers.push(i);
        }
    }

    var p = Utility.make(Options.k, Options.n - Options.k, function (row, column) {
        if ((row % Options.k) == (column % Options.k) && column > 0) {
            return 0;
        } else {
            return RandomGenerator.nextFromChoice(invertibleNumbers);
        }
    });

    var negPT = Utility.map(Utility.transpose(p), function (value) {
        return Utility.mod(-value, Options.modulo);
    });

    var generatorMatrix = Utility.concatColumns(Utility.identity(Options.k), p);
    var parityMatrix = Utility.concatColumns(negPT, Utility.identity(Options.n - Options.k));

    // Wikipedia says that we can check that the row space of G is
    // orthogonal to H by doing this:
    var test = Utility.multiply(generatorMatrix, Utility.transpose(parityMatrix), Options.modulo);

    // Every element of test should be zero
    Utility.map(test, function (value) {
        if (value) {
            throw "The generator and parity matrices don't have orthogonal row spaces";
        }
    });

    /**
    * Attempts to decode the given encoded message. Returns an array of
    * Options.k symbols. The returned array has a "decoded" property that tells
    * if the encoded message was successfully decoded completely (the same as
    * there being no erased symbols in the returned result)
    */
    function decode(encoded) {
        var encoded = Utility.deepCopy(encoded);
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
                if (variables.length > parity.length) {
                    break; // No point in continuing because we won't be able to solve this system of equations completely
                }
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
                    //sums[j] = mod(sums[j], Options.modulo);
                }
            }
        }

        // Only continue if we've got a change of solving the system of equations
        if (variables.length <= parity.length) {
            for (var j = 0; j < parity.length; j++) {
                // The sums will fill the right-hand column
                matrix[j].push(sums[j]);
            }
            // This tries to put the matrix into reduced row-echelon form
            var reduced = Utility.gaussReduction(matrix, Options.modulo);

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
        var retVal = encoded.slice(0, Options.k);
        retVal.decoded = true;
        retVal.forEach(function (value) {
            if (value === null || value < 0) {
                retVal.decoded = false;
            }
        });
        retVal.all = encoded;

        return retVal;
    }

    /**
    * Encodes the given message, returning an array Options.n symbols
    * long
    */
    function encode(message) {
        return Utility.multiply([message], RandomGenerator, Options.modulo)[0];
    };

    /**
    * Returns a copy of the generator matrix used by this instance of LDPC
    */
    function getGenerator() {
        return Utility.deepCopy(generatorMatrix);
    };

    /**
    * Returns a copy of the parity matrix used by this instance of LDPC
    */
    function getParity() {
        return Utility.deepCopy(parityMatrix);
    };

    return {
        options: Options,
        decode: decode,
        encode: encode,
        getGenerator: getGenerator,
        getParity: getParity
    };
});

module.exports = LDPC;
