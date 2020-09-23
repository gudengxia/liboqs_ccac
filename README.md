[AppVeyor](https://ci.appveyor.com/project/dstebila/liboqs): ![Build status image](https://ci.appveyor.com/api/projects/status/9d2ts78x88r8wnii/branch/master?svg=true), [CircleCI](https://circleci.com/gh/open-quantum-safe/liboqs/tree/master): ![Build status image](https://circleci.com/gh/open-quantum-safe/liboqs/tree/master.svg?style=svg)

liboqs-ccac
======================

liboqs is an open source C library for quantum-safe cryptographic algorithms. liboqs-ccac comes from liboqs, and add mulan and aigis signature algorithm to liboqs.


## Status

### Supported Algorithms

Details on each supported algorithm can be found in the [docs/algorithms folder](https://github.com/open-quantum-safe/liboqs/tree/master/docs/algorithms).

#### Key encapsulation mechanisms

- **BIKE**: BIKE1-L1-CPA, BIKE1-L3-CPA, BIKE1-L1-FO, BIKE1-L3-FO
- **Classic McEliece**: Classic-McEliece-348864†, Classic-McEliece-348864f†, Classic-McEliece-460896†, Classic-McEliece-460896f†, Classic-McEliece-6688128†, Classic-McEliece-6688128f†, Classic-McEliece-6960119†, Classic-McEliece-6960119f†, Classic-McEliece-8192128†, Classic-McEliece-8192128f†
- **FrodoKEM**: FrodoKEM-640-AES, FrodoKEM-640-SHAKE, FrodoKEM-976-AES, FrodoKEM-976-SHAKE, FrodoKEM-1344-AES, FrodoKEM-1344-SHAKE
- **HQC**: HQC-128-1-CCA2, HQC-192-1-CCA2, HQC-192-2-CCA2, HQC-256-1-CCA2†, HQC-256-2-CCA2†, HQC-256-3-CCA2†
- **Kyber**: Kyber512, Kyber768, Kyber1024, Kyber512-90s, Kyber768-90s, Kyber1024-90s
- **NTRU**: NTRU-HPS-2048-509, NTRU-HPS-2048-677, NTRU-HPS-4096-821, NTRU-HRSS-701
- **SABER**: LightSaber-KEM, Saber-KEM, FireSaber-KEM
- **SIKE**: SIDH-p434, SIDH-p503, SIDH-p610, SIDH-p751, SIKE-p434, SIKE-p503, SIKE-p610, SIKE-p751, SIDH-p434-compressed, SIDH-p503-compressed, SIDH-p610-compressed, SIDH-p751-compressed, SIKE-p434-compressed, SIKE-p503-compressed, SIKE-p610-compressed, SIKE-p751-compressed

#### Signature schemes

- **Dilithium**: Dilithium2, Dilithium3, Dilithium4
- **Falcon**: Falcon-512, Falcon-1024
- **Picnic**: Picnic-L1-FS, Picnic-L1-UR, Picnic-L1-full, Picnic-L3-FS, Picnic-L3-UR, Picnic-L3-full, Picnic-L5-FS, Picnic-L5-UR, Picnic-L5-full, Picnic3-L1, Picnic3-L3, Picnic3-L5
- **Rainbow**: Rainbow-Ia-Classic, Rainbow-Ia-Cyclic, Rainbow-Ia-Cyclic-Compressed, Rainbow-IIIc-Classic†, Rainbow-IIIc-Cyclic†, Rainbow-IIIc-Cyclic-Compressed†, Rainbow-Vc-Classic†, Rainbow-Vc-Cyclic†, Rainbow-Vc-Cyclic-Compressed†
- **SPHINCS+-Haraka**: SPHINCS+-Haraka-128f-robust, SPHINCS+-Haraka-128f-simple, SPHINCS+-Haraka-128s-robust, SPHINCS+-Haraka-128s-simple, SPHINCS+-Haraka-192f-robust, SPHINCS+-Haraka-192f-simple, SPHINCS+-Haraka-192s-robust, SPHINCS+-Haraka-192s-simple, SPHINCS+-Haraka-256f-robust, SPHINCS+-Haraka-256f-simple, SPHINCS+-Haraka-256s-robust, SPHINCS+-Haraka-256s-simple
- **SPHINCS+-SHA256**: SPHINCS+-SHA256-128f-robust, SPHINCS+-SHA256-128f-simple, SPHINCS+-SHA256-128s-robust, SPHINCS+-SHA256-128s-simple, SPHINCS+-SHA256-192f-robust, SPHINCS+-SHA256-192f-simple, SPHINCS+-SHA256-192s-robust, SPHINCS+-SHA256-192s-simple, SPHINCS+-SHA256-256f-robust, SPHINCS+-SHA256-256f-simple, SPHINCS+-SHA256-256s-robust, SPHINCS+-SHA256-256s-simple
- **SPHINCS+-SHAKE256**: SPHINCS+-SHAKE256-128f-robust, SPHINCS+-SHAKE256-128f-simple, SPHINCS+-SHAKE256-128s-robust, SPHINCS+-SHAKE256-128s-simple, SPHINCS+-SHAKE256-192f-robust, SPHINCS+-SHAKE256-192f-simple, SPHINCS+-SHAKE256-192s-robust, SPHINCS+-SHAKE256-192s-simple, SPHINCS+-SHAKE256-256f-robust, SPHINCS+-SHAKE256-256f-simple, SPHINCS+-SHAKE256-256s-robust, SPHINCS+-SHAKE256-256s-simple
- **Mulan**: mulan
- **Aigis**: aigis
Note that algorithms marked with a dagger (†) have large stack usage and may cause failures when run on threads or in constrained environments.

### Limitations and Security

As research advances, the supported algorithms may see rapid changes in their security, and may even prove insecure against both classical and quantum computers.

liboqs does not intend to "pick winners": algorithm support is informed by the NIST [Post-Quantum Cryptography Standardization](https://csrc.nist.gov/Projects/Post-Quantum-Cryptography/Post-Quantum-Cryptography-Standardization) project. We strongly recommend that applications and protocols rely on the outcomes of ths effort when deploying post-quantum cryptography.

We realize some parties may want to deploy quantum-safe cryptography prior to the conclusion of the NIST standardization project.  We strongly recommend such attempts make use of so-called **hybrid cryptography**, in which quantum-safe public-key algorithms are used alongside traditional public key algorithms (like RSA or elliptic curves) so that the solution is at least no less secure than existing traditional cryptography.

## Quickstart

### Linux/macOS

1. Install dependencies:

	On Ubuntu:

		 sudo apt install cmake gcc ninja-build libssl-dev python3-pytest python3-pytest-xdist unzip xsltproc doxygen graphviz

	On macOS, using a package manager of your choice (we've picked Homebrew):

		brew install cmake ninja openssl@1.1 wget doxygen graphviz
		pip3 install pytest pytest-xdist
	
	Note that, if you want liboqs to use OpenSSL for various symmetric crypto algorithms (AES, SHA-2, etc.) then you must have OpenSSL version 1.1.1 or higher.

2. Get the source:

		git clone -b master https://github.com/gudengxia/liboqs_ccac.git
		cd liboqs

	and build:

		mkdir build && cd build
		cmake .. -GNinja -DOQS_USE_OPENSSL=OFF -DBUILD_SHARED_LIBS=OFF
		ninja

Various options can be passed to `cmake` to customize the build. Some of them include:

- `-DOQS_USE_OPENSSL=<val>`: `<val>` can be `ON` or `OFF`; when `ON`, liboqs uses OpenSSL's AES, SHA-2, and SHA-3 implementations.
- `-DBUILD_SHARED_LIBS=<val>`: `<val>` can be `ON` or `OFF`; when `ON`, CMake generates instructions for building a shared library, otherwise it generates instructions for building a static library.
- `-DOPENSSL_ROOT_DIR=<dir>`: `<dir>` specifies the directory in which CMake will look for OpenSSL.

All supported options are listed in the `.CMake/alg-support.cmake` file, and can be viewed by running `cmake -LAH ..` in the `build` directory. They are also listed and explained in [the wiki](https://github.com/open-quantum-safe/liboqs/wiki/Customizing-liboqs).

The following instructions assume we are in `build`.

3. The main build result is `lib/liboqs.a`, a static library. The public headers are located in the `include` directory. There are also a variety of programs built under the `tests` directory:

	- `test_kem`: Simple test harness for key encapsulation mechanisms
	- `test_sig`: Simple test harness for key signature schemes
	- `kat_kem`: Program that generates known answer test (KAT) values for key encapsulation mechanisms using the same procedure as the NIST submission requirements, for checking against submitted KAT values using `tests/test_kat.py`
	- `kat_sig`: Program that generates known answer test (KAT) values for signature schemes using the same procedure as the NIST submission requirements, for checking against submitted KAT values using `tests/test_kat.py`
	- `speed_kem`: Benchmarking program for key encapsulation mechanisms; see `./speed_kem --help` for usage instructions
	- `speed_sig`: Benchmarking program for signature mechanisms; see `./speed_sig --help` for usage instructions
	- `example_kem`: Minimal runnable example showing the usage of the KEM API
	- `example_sig`: Minimal runnable example showing the usage of the signature API
	- `test_aes`, `test_sha3`: Simple test harnesses for crypto sub-components

	The test suite can be run using

		ninja run_tests

4. To generate HTML documentation of the API, run:

		ninja gen_docs

	Then open `docs/doxygen/html/index.html` in your web browser.

4. Finally, `ninja install` can be run to install the built library and `include` files to a location of choice, which can be specified by passing the `-DCMAKE_INSTALL_PREFIX=<dir>` option to `cmake` at configure time.

