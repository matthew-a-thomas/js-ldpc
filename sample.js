var getLDPC = require('./src/ldpc');
var getUtility = require('./src/utility');
var getRandomGenerator = require('./src/randomGenerator');

// Configure this test
var testParameters = {
	iterations: 1,
	k: {
		min: 1,
		max: 16
	},
	modulo: {
		min: 2,
		max: 16
	},
	n: {
		min: 1,
		max: 16
	}
};

function RunningAverage() {
	var sum = 0;
	var num = 0;
	this.getAverage = function() {
		return sum / num;
	};
	this.update = function(value) {
		sum += value;
		num++;
		return this.getAverage();
	};
}

// Override arrays' toString function for this test
Array.prototype.toString = function() {
	return JSON.stringify(this);
};

// An array of functions to execute in the course of testing
var testSteps = [
	function() {
		for (var modulo = testParameters.modulo.min; modulo <= testParameters.modulo.max; modulo++) {
			for (var n = testParameters.n.min; n <= testParameters.n.max; n++) {
				for (var k = testParameters.k.min; k < n && k <= testParameters.k.max; k++) {
					for (var iteration = 0; iteration < testParameters.iterations; iteration++) {
						// Set up options for this iteration
						var options = {
							n: n,
							k: k,
							modulo: modulo,
							randomSeed: Math.floor(Math.random() * Math.pow(2, 31))
						}

						// Queue up a function to execute this iteration
						testSteps.push(function(options) {
							return function() {
								// Run a test iteration with these options
								var report = testIteration(options);
								// Process the test report
								//processTestReport(report);
                                console.log(report);
							};
						}(options))
					}
				}
			}
		}
	}
];

// Executes tests
var go = function() {
	if (testSteps.length > 0) {
		var fn = testSteps.shift();
		fn();
		setTimeout(go, 0); // Go again (asynchronously, so that any other JS events can happen)
	}
};

// Respond to test reports
var passed = 0;
var failed = 0;
var startTime = 0;
var numTests = 0;


// Perform a single test
function testIteration(options) {
	// Set up a report object, which is what this function will return
	var report = {
		error: null,
		ldpc: {
			generator: [[]],
			parity: [[]]
		},
		message: {
			raw: [],
			encoded: [],
			erased: [],
			decoded: []
		},
		options: options
	};

	try {
		// Make our encoder/decoder
		var ldpc = getLDPC(options);
        var ldpcUtility = getUtility();
        var random = getRandomGenerator();

		report.ldpc.generator = ldpc.getGenerator();
		report.ldpc.parity = ldpc.getParity();

		// Make a random message
		var message = ldpcUtility.getRandomSeries(options.k, options.modulo);
		report.message.raw = message;

		// Encode it
		var encoded = ldpc.encode(message);
		report.message.encoded = encoded;

		// Simulate erasure channel by nulling some symbols
		encoded = ldpcUtility.deepCopy(encoded);
		var lostSymbols = random.nextSet(encoded.length);
		for (var i = 0; i < options.n - options.k; i++) { // We'll null n-k symbols
			encoded[lostSymbols.pop()] = null;
		}
		report.message.erased = encoded;

		// Try decoding from the encoded message
		var decoded = ldpc.decode(encoded);
		report.message.decoded = decoded;
	} catch (e) {
		report.error = e;
	} finally {
		return report;
	}
}

go();