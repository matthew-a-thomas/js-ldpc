# js-ldpc
A JavaScript implementation of a Low-Density Parity-Check code.

The goal of this project is to create an implementation of an LDPC (https://en.wikipedia.org/wiki/Low-density_parity-check_code) because an existing one doesn't seem readily available.

##API
These methods are planned:

###Create buffer
```JavaScript
var buffer = LDPC.createBuffer(message);
```
The type of `message` is not yet determined.

`returns` most likely an instance of (https://nodejs.org/api/buffer.html)
