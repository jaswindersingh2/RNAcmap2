#Name: PKB246
#Length:  34 
#PageNumber: 2
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GACCGTGCACCATGAGCACGAATCCTAAACCTCA
(...)...(.(...)....)........[(..])
SHHHSMMMSISHHHSIIIISMMMMMMMMMSHHHS
NNNNNNNNNNNNNNNNNNNNNNNNNNNNKNNNKN
S1 1..1 "G" 5..5 "G"
S2 9..9 "A" 20..20 "G"
S3 11..11 "C" 15..15 "A"
S4 30..30 "C" 34..34 "A"
H1 2..4 "ACC" (1,5) G:G 
H2 12..14 "ATG" (11,15) C:A 
H3 31..33 "CTC" (30,34) C:A PK{1}
I1.1 10..10 "C" (9,20) A:G 
I1.2 16..19 "GCAC" (15,11) A:C 
M1.1 6..8 "TGC" (5,1) G:G (9,20) A:G 
M1.2 21..29 "AATCCTAAA" (20,9) G:A (30,34) C:A PK{1}
PK1 1bp 29..29 33..33 M1.2 21..29 H3 31..33
PK1.1 29 A 33 C
NCBP1 11 C 15 A S3
NCBP2 30 C 34 A S4
NCBP3 1 G 5 G S1
NCBP4 9 A 20 G S2
NCBP5 29 A 33 C PK1.1
segment1 1bp 1..1 G 5..5 G
segment2 2bp 9..11 ACC 15..20 AGCACG
segment3 1bp 30..30 C 34..34 A
