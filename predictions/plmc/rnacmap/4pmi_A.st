#Name: 4pmi_A
#Length:  40 
#PageNumber: 2
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GGGAGUAUAUGGGCGCACUUCGGUGACGGUACAGGCUCCU
..[........(....]....)........(......)..
EEEEEEEEEEESHHHHHHHHHSXXXXXXXXSHHHHHHSEE
NNKNNNNNNNNNNNNNKNNNNNNNNNNNNNNNNNNNNNNN
S1 12..12 "G" 22..22 "G"
S2 31..31 "A" 38..38 "C"
H1 13..21 "GCGCACUUC" (12,22) G:G PK{1}
H2 32..37 "CAGGCU" (31,38) A:C 
X1 23..30 "GUGACGGU" (22,12) G:G (31,38) A:C 
E1 1..11 "GGGAGUAUAUG" PK{1}
E2 39..40 "CU" 
PK1 1bp 3..3 17..17 E1 1..11 H1 13..21
PK1.1 3 G 17 A
NCBP1 31 A 38 C S2
NCBP2 12 G 22 G S1
NCBP3 3 G 17 A PK1.1
segment1 1bp 12..12 G 22..22 G
segment2 1bp 31..31 A 38..38 C
