DINT
====

Gov2 Space (bpi)
----------------

| Encoder           | docs              | freqs            |
|:------------------|------------------:|-----------------:|
| Interpolative     |  2.50672          |  2.09706         |
| OptPFOR           |  2.99748          |  2.29313         |
| VarintG8IU        |  9.21201          |  9.00871         |
| QMX               |  3.76381          |  3.016           |
| VByte             |  8.26401          |  8.02            |
| Uncompressed      | 32.0004           | 32.0003          |
| Simple16          |  3.36379          |  2.62227         |
| StreamVByte       | 10.1672           | 10.0069          |
| MaskedVByte       |  8.26401          |  8.02            |
| VarintGB          | 10.1672           | 10.0069          |
| DINT*  (actual)   |  3.2332           |  2.2739          |
| DINT*  (estimated)|  3.7467 (+15.88%) |  2.6341 (+15.84%)|
| DINT** (actual)   |  3.5112           |  2.5285          |
| DINT** (estimated)|  4.4372 (+26.37%) |  3.1382 (+24.11%)|

CC_News Space (bpi)
-------------------

| Encoder           | docs              | freqs            |
|:------------------|------------------:|-----------------:|
| Interpolative     |  4.3954           |  3.2474          |
| OptPFOR           |  4.1970           |  2.7653          |
| VarintG8IU        |  9.2909           |  9.0007          |
| QMX               |  4.3301           |  3.0798          |
| VByte             |  8.3719           |  8.0026          |
| Uncompressed      | 32.0002           | 32.0002          |
| Simple16          |  4.1714           |  2.7591          |
| StreamVByte       | 10.2307           | 10.0006          |
| MaskedVByte       |  8.3719           |  8.0026          |
| VarintGB          | 10.2307           | 10.0006          |
| DINT*  (actual)   |  4.0660           |  2.0655          |
| DINT*  (estimated)|  4.4729 (+10.00%) |  2.2752 (+10.15%)|
| DINT** (actual)   |  4.3817           |  2.3267          |
| DINT** (estimated)|  5.0632 (+15.55%) |  2.7294 (+14.75%)|

Gov2 Decoding Time (ns x int)
-----------------------------

| Encoder           | docs        | freqs       |
|:------------------|------------:|------------:|
| Interpolative     | 8.32194     | 8.08594     |
| OptPFOR           | 0.88902     | 0.74642     |
| VarintG8IU        | 0.43187     | 0.43020     |
| QMX               | 0.75280     | 0.74732     |
| VByte             | 0.84330     | 0.71624     |
| Uncompressed      | 1.21538     | 1.01599     |
| Simple16          | 1.08865     | 0.88144     |
| StreamVByte       | 0.41451     | 0.41120     |
| MaskedVByte       | 0.44983     | 0.40227     |
| VarintGB          | 0.55312     | 0.51268     |
| DINT*             | 0.70945     | 0.42871     |
| DINT**            | 0.58482     | 0.39663     |
| DINT***           | 1.33705     | 0.94880     |

CC_News Decoding Time (ns x int)
--------------------------------

| Encoder           | docs        | freqs       |
|:------------------|------------:|------------:|
| Interpolative     | 9.28114     | 8.52188     |
| OptPFOR           | 0.98962     | 0.88705     |
| VarintG8IU        | 0.56313     | 0.44866     |
| QMX               | 1.38498     | 1.28182     |
| VByte             | 0.95535     | 0.81070     |
| Uncompressed      | 0.90068     | 0.89023     |
| Simple16          | 1.33971     | 1.02731     |
| StreamVByte       | 0.53137     | 0.49396     |
| MaskedVByte       | 0.53541     | 0.44222     |
| VarintGB          | 0.61775     | 0.58112     |
| DINT*             | 0.87670     | 0.64105     |
| DINT**            | 0.75852     | 0.47938     |

- DINT*   uses rectangular dictionaries of 2^16 x 16 x 4 bytes = 4 MiB
- DINT**  uses rectangular dictionaries of 2^16 x  8 x 4 bytes = 2 MiB
- DINT*** uses packed dictionaries: 1.371 MiB for the docs; 0.591 MiB for the freqs.

Gov2.docs DINT Statistics
-------------------------

-  4847 16-int entries out of 156569293 (0.00309575%)
-  6456  8-int entries out of 200793248 (0.00321525%)
- 12654  4-int entries out of 121968050 (0.0103748%)
-  7087  2-int entries out of  21792274 (0.0325207%)
- 34492  1-int entries out of    695841 (4.95688%)

| total ints      | ints covered by runs | ints covered by table | ints covered by exceptions |
|----------------:|---------------------:|----------------------:|---------------------------:|
| 5406586692      | 2229482656 (41.23%)  | 3173457461 (58.69%)   | 3646575  (0.067%)          |

| total codewords | codewords for runs   | codewords for table   | codewords for exceptions   |
|----------------:|---------------------:|----------------------:|---------------------------:|
| 1092408356      | 14132727 (1.293%)    | 1067335904 (97.70%)   | 10939725  (1.001%)         |

| total codewords         | 1092408356           |
|-------------------------|---------------------:|
| codewords for 16+ ints  |  39752184 ( 3.638%)  |
| codewords for  8  ints  |  69770445 ( 6.386%)  |
| codewords for  4  ints  | 265554905 (24.309%)  |
| codewords for  2  ints  | 436771872 (39.982%)  |
| codewords for  1  ints  | 280558950 (25.682%)  |


CC_News.docs DINT Statistics
----------------------------

-  4847 16-int entries out of 680415119 (0.000712359%)
- 10423  8-int entries out of 870845771 (0.00119688%)
- 14955  4-int entries out of 490138856 (0.00305118%)
-  9135  2-int entries out of  67162184 (0.0136014%)
- 26176  1-int entries out of    833660 (3.13989%)

| total ints      | ints covered by runs | ints covered by table | ints covered by exceptions |
|----------------:|---------------------:|----------------------:|---------------------------:|
| 19691599096     | 5385274048 (27.34%)  | 14288693355 (72.56%)  | 17631693  (0.089%)         |

| total codewords | codewords for runs   | codewords for table   | codewords for exceptions   |
|----------------:|---------------------:|----------------------:|---------------------------:|
| 5003968348      | 46411036 (0.927%)    | 4904662233 (98.01%)   | 52895079  (1.057%)         |

| total codewords         | 5003968348           |
|-------------------------|---------------------:|
| codewords for 16+ ints  |  158635935 ( 3.17%)  |
| codewords for  8  ints  |  315066971 ( 6.29%)  |
| codewords for  4  ints  | 1141219579 (22.80%)  |
| codewords for  2  ints  | 2071530103 (41.39%)  |
| codewords for  1  ints  | 1317515760 (26.32%)  |

- A codeword is a short, namely 16 bits.
- Exceptions are encoded using 3 shorts.

