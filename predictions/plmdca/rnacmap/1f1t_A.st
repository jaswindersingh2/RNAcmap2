#Name: 1f1t_A
#Length: ,38,
#PageNumber: 2
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GGAUCCCGACUGGCGAGAGCCAGGUAACGAAUGGAUCC
(..(.[(...].................).....)..)
SIISIISHHHHHHHHHHHHHHHHHHHHHSIIIIISIIS
NNNNNKNNNNKNNNNNNNNNNNNNNNNNNNNNNNNNNN
S1 1..1 "G" 38..38 "C"
S2 4..4 "U" 35..35 "A"
S3 7..7 "C" 29..29 "G"
H1 8..28 "GACUGGCGAGAGCCAGGUAAC" (7,29) C:G PK{1}
I1.1 2..3 "GA" (1,38) G:C 
I1.2 36..37 "UC" (35,4) A:U 
I2.1 5..6 "CC" (4,35) U:A PK{1}
I2.2 30..34 "AAUGG" (29,7) G:C 
PK1 1bp 6..6 11..11 I2.1 5..6 H1 8..28
PK1.1 6 C 11 U
NCBP1 6 C 11 U PK1.1
segment1 3bp 1..7 GGAUCCC 29..38 GAAUGGAUCC
