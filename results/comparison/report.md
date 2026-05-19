# Agent vs Manual Spidroin Typing Comparison

Reciprocal overlap threshold: **50%**

## Locus-level metrics (overall)

- True positives (matched): **1309**
- False negatives (manual_only): **435**
- False positives (agent_only): **861**
- Precision: **0.603**
- Recall: **0.751**
- F1: **0.669**

## Classification agreement on matched pairs

- Spidroin_type identical: 1103 / 1309 (84.3%)
- Full_length identical: 1262 / 1309 (96.4%)
- Hint_type identical: 1291 / 1309 (98.6%)

## Species coverage

| side        |   len |
|:------------|------:|
| agent_only  |    18 |
| both        |   109 |
| manual_only |     4 |

## Recall stratified by manual Scoring

|   manual_scoring |   total |   matched |   recall |
|-----------------:|--------:|----------:|---------:|
|          nan     |   2.000 |     1.000 |    0.500 |
|            1.000 |  69.000 |    33.000 |    0.478 |
|            2.000 | 277.000 |   146.000 |    0.527 |
|            3.000 | 273.000 |   120.000 |    0.440 |
|            4.000 | 697.000 |   613.000 |    0.879 |
|            5.000 | 426.000 |   396.000 |    0.930 |

## Boundary offsets (matched pairs, agent - manual)

| statistic   |   start_diff |    end_diff |
|:------------|-------------:|------------:|
| count       |    1309      |   1309      |
| null_count  |       0      |      0      |
| mean        |      33.5523 |    -47.0519 |
| std         |    1454.87   |   1311.41   |
| min         |   -9527      | -46111      |
| 25%         |       0      |      0      |
| 50%         |       0      |      0      |
| 75%         |       0      |      0      |
| max         |   50419      |   2481      |

## Agent confidence vs manual Scoring crosstab

| agent_confidence   |   manual_scoring |   count |
|:-------------------|-----------------:|--------:|
| high               |              nan |       1 |
| high               |                1 |       1 |
| high               |                2 |      51 |
| high               |                3 |      59 |
| high               |                4 |     578 |
| high               |                5 |     384 |
| medium             |                1 |      32 |
| medium             |                2 |      95 |
| medium             |                3 |      61 |
| medium             |                4 |      35 |
| medium             |                5 |      12 |

## Bottom 10 species by Recall

| species                         |   tp |   fn |   fp |   precision |   recall |   f1 |
|:--------------------------------|-----:|-----:|-----:|------------:|---------:|-----:|
| 062.Anelosimus studiosus        |    0 |    0 |   23 |       0.000 |      nan |  nan |
| 061.Amaurobius ferox            |    0 |    0 |   20 |       0.000 |      nan |  nan |
| 117.Tetragnatha versicolor      |    0 |    0 |   19 |       0.000 |      nan |  nan |
| 125.Aptostichus stephencolberti |    0 |    0 |    7 |       0.000 |      nan |  nan |
| 124.Uloborus plumipes           |    0 |    0 |   11 |       0.000 |      nan |  nan |
| 123.Uloborus diversus           |    0 |    0 |   14 |       0.000 |      nan |  nan |
| 046.Scytodidae sp               |    0 |    0 |    1 |       0.000 |      nan |  nan |
| 066.Argiope aurantia            |    0 |    0 |   56 |       0.000 |      nan |  nan |
| 120.Trichonephila clavipes      |    0 |    0 |   34 |       0.000 |      nan |  nan |
| 122.Troglohyphantes excavatus   |    0 |    0 |   10 |       0.000 |      nan |  nan |

## Top Spidroin_type confusions

| manual   | agent     |   count |
|:---------|:----------|--------:|
| acsp     | acsp      |     182 |
| acsp     | masp      |       3 |
| acsp     | cysp      |       2 |
| acsp     | masp/acsp |       1 |
| acsp     | pysp      |       1 |
| agsp1    | agsp1     |      24 |
| agsp1    | unknown   |       1 |
| agsp1    | acsp      |       1 |
| agsp2    | agsp2     |      29 |
| crsp     | crsp      |       7 |
| cysp     | cysp      |      72 |
| cysp     | crsp      |       2 |
| flag     | flag      |      65 |
| masp     | masp      |     120 |
| masp     | masp1     |      66 |
| masp     | misp      |      36 |
| masp     | masp2     |       8 |
| masp     | masp/acsp |       3 |
| masp     | acsp      |       2 |
| masp     | masp/misp |       1 |

## Full_length confusion matrix

| manual   | agent   |   count |
|:---------|:--------|--------:|
| False    | False   |     132 |
| False    | True    |      41 |
| True     | True    |    1130 |
| True     | False   |       6 |

## Hint_type confusion matrix

| manual      | agent       |   count |
|:------------|:------------|--------:|
| C-terminal  | C-terminal  |      56 |
| C-terminal  | Full_length |       1 |
| Full_length | Full_length |    1166 |
| Full_length | C-terminal  |       8 |
| Full_length | N-terminal  |       4 |
| N-terminal  | N-terminal  |      69 |
| N-terminal  | Full_length |       4 |
| N-terminal  | C-terminal  |       1 |
