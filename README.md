`dint` - Dictionary of INTeger sequences
-------

More informative README is coming soon.

Benchmark
---------

A comparison between the space of `single_rect`, `single_packed` and `multi_packed` on the provided `test_collection` is shown below (`bpi` stands for "bits per integer"). Results have been collected on a machine with an Intel i7-7700 processor clocked at 3.6 GHz and running Linux 4.4.0, 64 bits. The code was compiled using the highest optimization setting (see CMakeLists.txt).

|     **Index**     |**docs [bpi]**  |**freqs [bpi]**  |
|-------------------|---------------:|----------------:|
|`single_rect`      | 5.939          | 3.047           |
|`single_packed`    | 5.939          | 3.047           |
|`multi_packed`     | 4.766          | 2.455           |
|`PEF eps-opt`      | 6.369          | 3.479           |


Authors
-------
* Giulio Ermanno Pibiri, <giulio.pibiri@di.unipi.it>
* Matthias Petri, <matthias.petri@gmail.com>
* Alistair Moffat, <ammoffat@unimelb.edu.au>

Bibliography
------------
* [1] Giulio Ermanno Pibiri, Matthias Petri and Alistair Moffat, *Fast Dictionary-based Compression for Inverted Indexes*. In the Proceedings of the 12-th ACM Conference on Web Search and Data Mining (WSDM 2019).