/**
 * Creates a new Low-Density Parity-Check encoder/decoder
 * @param numSymbols The number of symbols which will be encoded/decoded
 * @param overhead The number of additional symbols used as overhead
 * @param modulo The number space we're working in
 * @param mixing The percent of constraints that factor in each node
 *  (expressed as a number between 0 and 1 / a percentage)
 */
function LDPC(options) {
	options = options || {};
	options.numSymbols = options.numSymbols || 0; // This isn't a useful value, so you'll want to provide your own
	options.overhead = options.overhead || 0;
	options.modulo = options.modulo || 2; // Default to base 2
	options.mixing = options.mixing || 0.1; // Default to 10% mixing
	options.randomSeed = options.randomSeed || undefined; // Default to current time (see Random class's constructor)

	/**
	 * A constraint defines a series of nodes which must add to zero (within the
	 * given modulo space)
	 */
	function Constraint() {
		this.nodes = [];
		this.tryToSatisfy = function() {
			var total = 0;
			var singleErasedNode = null;
			for (var i = this.nodes.length - 1; i >= 0; i--) {
				var node = this.nodes[i];
				if (!node.getErased()) {
					total += node.getValue();
					total %= options.modulo;
				} else {
					if (!singleErasedNode) {
						singleErasedNode = node;
					} else {
						// There's already an erased node, so having another one
						//  means that we won't be able to deduce the right answer
						return false;
					}
				}
			};
			if (singleErasedNode) {
				// Only one of our nodes is erased, so we can deduce its value from
				//  the rest of the nodes (its value must add to the total to make
				//  it zero)
				singleErasedNode.setValue((options.modulo - total) % options.modulo);
				return true; // We've resolved a previously-erased node
			} else {
				return total == 0; // This constraint is satisfied
			}
		};
	}

	/**
	 * nextID returns monotomically-increasing numbers
	 */
	var nextID = function() {
		var lastID = -1;
		return function() {
			return ++lastID;
		};
	}();

	/**
	 * A node is just a value and an indication of it being erased or not
	 */
	function Node() {
		var value = 0;
		this.constraints = [];
		/**
		 * Returns true if the value is < 0
		 */
		this.getErased = function() {
			return value < 0;
		};
		this.getValue = function() {
			return value;
		};
		/**
		 * Use newValue < 0 to indicate that it's missing/erased. Returns this
		 * node
		 */
		this.setValue = function(newValue) {
			value = newValue;
			return this;
		};
	}

	/**
	 * A simple RNG per http://stackoverflow.com/questions/3062746/special-simple-random-number-generator
	 * Note this will only generate numbers up to 2^32 - 1
	 */
	function Random(seed) {
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
	}

	/**
	 * Returns a shuffled set of indices from the given array, sorted
	 * randomly but in preference of the lesser-used elements of the array.
	 * arrayUsage is updated with new values. rand is the random number
	 * generator, and perturbation is the percent to mix up the results
	 * (ranging from 0 to 1). Perturbation of 1 means there's a 100% chance
	 * that the result will be out of order
	 */
	function sortStocastic(arrayUsage, rand, perturbation) {
		// Make an array of objects that link index to usage
		var detailed = [];
		for (var i = 0; i < arrayUsage.length; i++) {
			detailed.push({
				index: i,
				usage: arrayUsage[i]
			});
		}
		detailed = detailed.sort(function(a, b) {
			return a.usage - b.usage;
		});
		// detailed is now sorted from least-used to most-used. Time to mix
		//  it up
		for (var i = 0; i < detailed.length - 1; i++) {
			if (rand.next() < perturbation) {
				// Swap this element with the next one
				perturbedIndex = i + 1;
				var temp = detailed[perturbedIndex];
				detailed[perturbedIndex] = detailed[i];
				detailed[i] = temp;
			}
		}
		// Now pull out the list of indices from detailed and return that
		var retVal = [];
		for (var i = 0; i < detailed.length; i++) {
			var detail = detailed[i];
			retVal.push(detail.index);
		}
		return retVal;
	}

	/**
	 * Hooks a node and constraint together
	 */
	function link(node, constraint) {
		node.constraints.push(constraint);
		constraint.nodes.push(node);
	}

	/**
	 * Returns the given array, shuffled using the given random number generator
	 */
	function shuffle(array, rng) {
		array = array.slice();
		for (var i = 0; i < array.length; i++) {
			var swapIndex = Math.floor(rng.next() * array.length);
			var temp = array[swapIndex];
			array[swapIndex] = array[i];
			array[i] = temp;
		}
		return array;
	}

	// Set up nodes and constraints. This is basically us creating the sparse
	//  parity check matrix
	var nodes = [];
	var constraints = [];
	var constraintUsage = [];
	for (var i = 0; i < options.overhead; i++) {
		constraints.push(new Constraint());
		constraintUsage.push(0);
	}
	for (var i = 0; i < options.numSymbols + options.overhead; i++) {
		var node = new Node();
		nodes.push(node);
	}
	// Seed the RNG in preparation of linking nodes and constraints
	var rand = new Random(options.randomSeed);
	// Note that we don't want all of the constraints to have an even number of
	//  corresponding nodes.
	// Note that as we go along, we want to make sure that each constraint has
	//  about the same number of nodes, so we'll give statistical preference to
	//  constraints which don't have as many nodes
	for (var i = 0; i < nodes.length; i++) {
		var node = nodes[i];
		var constraintIndices = sortStocastic(constraintUsage, rand, options.mixing);
		for (var j = constraintIndices.length * options.mixing + rand.nextGaussian(options.mixing); j > 0 && constraintIndices.length > 0; j--) {
			var constraintIndex = constraintIndices.shift();
			var constraint = constraints[constraintIndex];
			link(node, constraint);
			constraintUsage[constraintIndex]++;
		}
	}

	this.decode = function(symbolArray) {
		if (symbolArray.length != nodes.length) {
			throw "Wrong number of symbols. This LDPC needs " + nodes.length;
		}

		// Set the nodes' values
		for (var i = 0; i < symbolArray.length; i++) {
			var symbol = symbolArray[i];
			var node = nodes[i];
			node.setValue(symbol);
		}

		// Continue trying to resolve constraints until we've resolved as many as we can
		var lastNumConstraintsSatisfied = 0;
		while (true) {
			var numConstraintsSatisfied = 0;
			for (var i = 0; i < constraints.length; i++) {
				var constraint = constraints[i];
				numConstraintsSatisfied += constraint.tryToSatisfy() ? 1 : 0;
			}
			if (numConstraintsSatisfied <= lastNumConstraintsSatisfied)
				break;
			lastNumConstraintsSatisfied = numConstraintsSatisfied;
		}

		// Gather results together to return
		var result = [];
		var decoded = true;
		for (var i = 0; i < nodes.length; i++) {
			var node = nodes[i];
			decoded &= !node.getErased();
			result.push(node.getValue());
		}
		return {
			decoded: decoded ? true : false,
			result: result
		};
	};
}