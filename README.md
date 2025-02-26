# Batch-PoSO Pseudocode: Prover and Verifier

Implementation in C (and PARI/GP for the 3square decomposition) of Batch-PoSO from https://eprint.iacr.org/2022/1153.pdf (CCS 2022). This implementation uses OpenSSL for portability and generic purpose (more diversity on the usable curves for instance). The paper timing were done by adapting an old version of bitcoin-core [secp256k1](https://github.com/bitcoin-core/secp256k1) and is not part of the scope of this repository.

The code was tested on a 2,3 GHz Intel Core i7 MacBook Pro.

## Requirements

Requires libpari (https://pari.math.u-bordeaux.fr/pub/pari/manuals/2.13.3/INSTALL.pdf), OpenSSL's libcrypto and libssl (https://github.com/openssl/openssl/blob/master/INSTALL.md).
The version used when testing the codes where the following:

```bash
> otool -L main-dyn 
main-dyn:
	/usr/local/lib/libpari-gmp.dylib (compatibility version 2.13.0, current version 2.13.2)
	/usr/local/opt/openssl@3/lib/libssl.3.dylib (compatibility version 3.0.0, current version 3.0.0)
	/usr/local/opt/openssl@3/lib/libcrypto.3.dylib (compatibility version 3.0.0, current version 3.0.0)
	/usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 1351.0.0)
```

## Test

The code can be simply tested using the following commands:

```bash
> make clean; make
> ./main-dyn
```

The outputs provides the number of protocols aborted due to rejection sampling and the average time for the prover part and the verifier part.

```bash
> ./main-dyn
Total number of rejected protocol: 215 / 1000
Total time taken by CPU for PROVE: 0.005061s
Total time taken by CPU for VERIF: 0.003228s
```

## Notes

_OpenSSL_: We used OpenSSL to perform the arithmetic operations for portability and simplicity in order to produce a reference implementation loyal to the specifications. This allowed also to do tests on order curves by simply changing the NID in the main. In order to achieve the best performance, one can make used of libraries such as bitcoin-core', which implements seckp256k1 operations with higly optimized code (https://github.com/bitcoin-core/secp256k1).

_Clean Up_: Due to the nature of this implementation being a reference implementation, numerous variables are instantiated and then cleaned/freed, which makes the clean up a bit tedious/ugly. The number of variables used, and their cleaning could be easily reduced and made compact for further implementations.

