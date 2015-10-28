# js-ldpc
A JavaScript implementation of a Low-Density Parity-Check code.

The goal of this project is to create an implementation of an LDPC (https://en.wikipedia.org/wiki/Low-density_parity-check_code) because an existing one doesn't seem readily available.

Robert Gallager's original paper on the topic can be found here: http://www.rle.mit.edu/rgallager/documents/ldpc.pdf

##API
These methods are planned:

###Constructor
```JavaScript
var ldpc = new LDPC(options);
```
`options` is an object that can have these properties:
* `numSymbols` the number of symbols in an unencoded message
* `overhead` the number of symbols to add for redundancy
* `modulo` currently only base-2 is planned. Defaults to 2
* `mixing` a parameter that tries to capture the essence of tuning a LDPC code. This affects how dense the parity-check matrix is among other things. Defaults to 10% (which is probably too low)
* `randomSeed` value to initialize the random number generator to. Defaults to `Date.now()`

###Decoding
```JavaScript
var result = ldpc.decode(array);
if (result.decoded) {
  // result.result should have all-positive numbers
} else {
  // result.result might have more positive numbers than you gave the decode function
}
```
`array` is an array `options.numSymbols + options.overhead` in length of numbers that are less than `options.modulo`. Use negative numbers to indicate a missing symbol (we're targeting an erasure channel)
