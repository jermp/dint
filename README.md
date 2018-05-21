DINT
====

Gov2 Space (bpi)
----------------

| Encoder           | docs              | freqs            |
|:------------------|------------------:|-----------------:|
| Interpolative     |  3.9932           |  3.4735          |
| OptPFOR           |  3.5859           |  2.9141          |
| VarintG8IU        |  9.2124           |  9.0087          |
| QMX               |  3.5031           |  3.0636          |
| VByte             |  8.2653           |  8.0202          |
| Uncompressed      | 32.0004           | 32.0003          |
| Simple16          |  3.5154           |  2.9305          |
| StreamVByte       | 10.1676           | 10.0070          |
| MaskedVByte       |  8.2653           |  8.0202          |
| VarintGB          | 10.1676           | 10.0070          |
| DINT (actual)     |  3.2332           |  2.2739          |
| DINT (estimated)  |  3.7467 (+15.88%) |  2.6341 (+15.84%)|

CC_News Space (bpi)
-------------------

| Encoder           | docs              | freqs            |
|:------------------|------------------:|-----------------:|
| Interpolative     |  4.3954           |                  |
| OptPFOR           |  4.1970           |                  |
| VarintG8IU        |  9.2909           |                  |
| QMX               |  4.3301           |                  |
| VByte             |  8.3719           |                  |
| Uncompressed      | 32.0002           |                  |
| Simple16          |  4.1714           |                  |
| StreamVByte       | 10.2307           |                  |
| MaskedVByte       |  8.3719           |                  |
| VarintGB          | 10.2307           |                  |
| DINT (actual)     |  4.0660           |  2.0655          |
| DINT (estimated)  |  4.4729 (+10.00%) |  2.2752 (+10.15%)|

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

CC_News Decoding Time (ns x int)
--------------------------------

| Encoder           | docs        | freqs       |
|:------------------|------------:|------------:|
| Interpolative     | 9.28114     |             |
| OptPFOR           | 0.98962     |             |
| VarintG8IU        | 0.56313     |             |
| QMX               | 1.38498     |             |
| VByte             | 0.95535     |             |
| Uncompressed      | 0.90068     |             |
| Simple16          | 1.33971     |             |
| StreamVByte       | 0.53137     |             |
| MaskedVByte       | 0.53541     |             |
| VarintGB          | 0.61775     |             |
| DINT*             | 0.87670     |             |


* DINT uses rectangular dictionaries of 2^16 x 16 x 4 bytes

Gov2 DINT Statistics
--------------------

| total ints      | ints covered by runs | ints covered by table | ints covered by exceptions |
|----------------:|---------------------:|----------------------:|---------------------------:|
| 5406586692      | 2229482656 (41.23%)  | 3173457461 (58.69%)   | 3646575  (0.067%)          |

| total codewords | codewords for runs   | codewords for table   | codewords for exceptions   |
|----------------:|---------------------:|----------------------:|---------------------------:|
| 1092408356      | 14132727 (1.293%)    | 1067335904 (97.70%)   | 10939725  (1.001%)         |

| total codewords         | 1092408356           |
|-------------------------|---------------------:|
| codewords for 16+ ints  | 39752184 (3.638%)    |
| codewords for 8 ints    | 69770445 (6.386%)    |
| codewords for 4 ints    | 265554905 (24.309%)  |
| codewords for 2 ints    | 436771872 (39.982%)  |
| codewords for 1 ints    | 280558950 (25.682%)  |


CC_News DINT Statistics
-----------------------

| total ints      | ints covered by runs | ints covered by table | ints covered by exceptions |
|----------------:|---------------------:|----------------------:|---------------------------:|
| 19691599096     | 5385274048 (27.34%)  | 14288693355 (72.56%)  | 17631693  (0.089%)         |

| total codewords | codewords for runs   | codewords for table   | codewords for exceptions   |
|----------------:|---------------------:|----------------------:|---------------------------:|
| 5003968348      | 46411036 (0.927%)    | 4904662233 (98.01%)   | 52895079  (1.057%)         |

| total codewords         | 5003968348           |
|-------------------------|---------------------:|
| codewords for 16+ ints  | 158635935 (3.170%)   |
| codewords for 8 ints    | 315066971 (6.296%)   |
| codewords for 4 ints    | 1141219579 (22.80%)  |
| codewords for 2 ints    | 2071530103 (41.39%)  |
| codewords for 1 ints    | 1317515760 (26.32%)  |

- A codeword is a short, namely 16 bits.
- Exceptions are encoded using 3 shorts.

