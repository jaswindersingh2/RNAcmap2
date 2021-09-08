#Name: PKB246
#Length:  34 
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GACCGTGCACCATGAGCACGAATCCTAAACCTCA
...[({.])}.............(...)(....)
EEEESHHHSMMMMMMMMMMMMMMSHHHSSHHHHS
NNNKNKNKNKNNNNNNNNNNNNNNNNNNNNNNNN
S1 5..5 "G" 9..9 "A"
S2 24..24 "C" 28..28 "A"
S3 29..29 "A" 34..34 "A"
H1 6..8 "TGC" (5,9) G:A PK{1,2}
H2 25..27 "CTA" (24,28) C:A 
H3 30..33 "CCTC" (29,34) A:A 
M1.1 10..23 "CCATGAGCACGAAT" (9,5) A:G (24,28) C:A PK{2}
M1.2 29..28 "" (28,24) A:C (29,34) A:A 
E1 1..4 "GACC" PK{1}
PK1 1bp 4..4 8..8 E1 1..4 H1 6..8
PK2 1bp 6..6 10..10 H1 6..8 M1.1 10..23
PK1.1 4 C 8 C
PK2.1 6 T 10 C
NCBP1 24 C 28 A S2
NCBP2 29 A 34 A S3
NCBP3 5 G 9 A S1
NCBP4 4 C 8 C PK1.1
NCBP5 6 T 10 C PK2.1
segment1 1bp 5..5 G 9..9 A
segment2 1bp 24..24 C 28..28 A
segment3 1bp 29..29 A 34..34 A
