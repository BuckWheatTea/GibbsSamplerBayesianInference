# GibbsSamplerBayesianInference
## main function
1. GibbsSampler
2. GibbsSamplerOrdinary

### GibbsSampler
GibbsSampler using Bayesian Inference
#### parameters:
Sequences
MotifLength
iterT: iteriate times, the higher, the more accurate, the slower
CtrlN: Whether or not use strategy 1,
       set 0 means DO NOT apply strategy 1,
       set -1 means apply strategy 1, but the real value of CtrlN is given by program automaticall
       set a positive number as real value of CtrlN, notice that it must be less than iterT
improve: Whether or not use strategy 2
       set false means DO NOT apply strategy 2
       set true means apply strategy 2
alpha: DO NOT recommend that set it by yourself
beta:  DO NOT recommend that set it by yourself

#### Example 1:
Motifs = GibbsSampler(Sequences, 10, 3000, 0, true)

MotifLength = 10,
iterT=3000,
CtrlN=0 means DO NOT apply strategy 1,
improve = true means apply strategy 2

## Example 2:
Motifs = GibbsSampler(Sequences, 12, 5000, -1, false)

MotifLength = 12,
iterT=5000,
CtrlN=-1 means apply strategy 1 but the real value is given by program,
improve = flase means DO NOT apply strategy 2 

### GibbsSamplerOrdinary
Orindar GibbsSampler method
#### parameters:
Sequences
MotifLength
iterT: iteriate times, the higher, the more accurate, the slower


#### Example:
Motifs, score, _, _ = GibbsSamplerOrdinary(Sequences, 10, 10000)

MotifLength = 10,
iterT=10000,
