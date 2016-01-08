var getRandomGenerator = (function (seed) {
    /**
    * A simple RNG per http://stackoverflow.com/questions/3062746/special-simple-random-number-generator
    * Note this will only generate numbers up to 2^32 - 1
    */
    if (typeof seed === undefined) {
        // Seed with current time if no seed given
        seed = Date.now();
    }

    var a = 1103515245;
    var c = 12345;
    var m = Math.pow(2, 32);
    /**
    * Generates a number between 0 and 1 (with 32-bit precision)
    */
    function next() {
        seed = (a * seed + c) % m;
        return seed / m;
    };

    /**
    * Picks a random choice from the given choices
    */
    function nextFromChoice(choices) {
        return choices[Math.floor(this.next() * choices.length)];
    };

    /**
    * Generates a normal number (mean of zero) with given standard
    * deviation
    */
    function nextGaussian(standardDeviation) {
        var u1 = this.next();
        var u2 = this.next();
        var randStdNormal = Math.sqrt(-2 * Math.log(u1)) * Math.sin(2 * Math.PI * u2);
        return standardDeviation * randStdNormal;
    };

    /**
    * http://wiki.q-researchsoftware.com/wiki/How_to_Generate_Random_Numbers:_Poisson_Distribution
    */
    function nextPoisson(lambda) {
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
    function nextSet(length) {
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

    return {
        next: next,
        nextFromChoice: nextFromChoice,
        nextGaussian: nextGaussian,
        nextPoisson: nextPoisson,
        nextSet: nextSet
    }
});

module.exports = getRandomGenerator;
