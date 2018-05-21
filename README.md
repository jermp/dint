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


* DINT uses rectangular dictionaries of 2^16 x 16 x 4 bytes
