# js-ldpc
A JavaScript implementation of a Low-Density Parity-Check code.

The goal of this project is to create an implementation of an LDPC (https://en.wikipedia.org/wiki/Low-density_parity-check_code) because an existing one doesn't seem readily available.

Robert Gallager's original paper on the topic can be found here: http://www.rle.mit.edu/rgallager/documents/ldpc.pdf

##Browser support
Currently this has only been tested in Chrome >= 46

##Examples
* `example.html` - shows decoding in the face of symbols being lost
* `example2.html` - shows how many symbols are necessary from an encoded message to successfully decode

##Documentation

###Constructor
```JavaScript
var ldpc = new LDPC(options);
```
`options` is an object that can have these properties:
* `n` the number of symbols in an encoded message
* `k` the number of non-redundant symbols (number of symbols in unencoded message)
* `modulo` e.g. `2` for binary digits, `256` for bytes, etc
* `randomSeed` value to which to initialize an internal random number generator, which is used to build the parity-check and generator matrices. Defaults to `Date.now()`

###Encoding
```JavaScript
var encoded = ldpc.encode(message);
console.log(encoded.join()); // Prints the encoded message
```
`message` is an array `options.k` in length of numbers that are less than `options.modulo`

###Decoding
```JavaScript
var result = ldpc.decode(encoded);
if (result.decoded) {
  // result should have all-positive and non-null numbers
} else {
  // result might have more positive numbers than you gave the decode function
  // Also see result.result, which has all the encoded symbols including any additional which were deduced. Note that the first k symbols are the same as result
}
```
`encoded` is an array `options.n` in length of numbers that are less than `options.modulo`. Use negative numbers or `null` to indicate a missing symbol (we're targeting an erasure channel)
