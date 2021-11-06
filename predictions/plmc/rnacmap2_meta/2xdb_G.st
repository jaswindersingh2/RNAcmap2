#Name: 2xdb_G
#Length: ,40,
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
AUUCAGGUGAUUUGCUACCUUUAAGUGCAGCUAGAAAUUC
..([{(...([(..].....]..)[)}...)..]).....
EESIISIIISISHHHHHHHHHHHSISIIIISIIISEEEEE
NNNKKNNNNNKNNNKNNNNNKNNNKNKNNNNNNKNNNNNN
S1 3..3 "U" 35..35 "A"
S2 6..6 "G" 31..31 "C"
S3 10..10 "A" 26..26 "U"
S4 12..12 "U" 24..24 "A"
H1 13..23 "UGCUACCUUUA" (12,24) U:A PK{1,3}
I1.1 4..5 "CA" (3,35) U:A PK{1,2}
I1.2 32..34 "UAG" (31,6) C:G PK{4}
I2.1 7..9 "GUG" (6,31) G:C 
I2.2 27..30 "GCAG" (26,10) U:A PK{2}
I3.1 11..11 "U" (10,26) A:U PK{3}
I3.2 25..25 "G" (24,12) A:U PK{4}
E1 1..2 "AU" 
E2 36..40 "AAUUC" 
PK1 1bp 4..4 21..21 I1.1 4..5 H1 13..23
PK2 1bp 5..5 27..27 I1.1 4..5 I2.2 27..30
PK3 1bp 11..11 15..15 I3.1 11..11 H1 13..23
PK4 1bp 25..25 34..34 I3.2 25..25 I1.2 32..34
PK1.1 4 C 21 U
PK2.1 5 A 27 G
PK3.1 11 U 15 C
PK4.1 25 G 34 G
NCBP1 4 C 21 U PK1.1
NCBP2 5 A 27 G PK2.1
NCBP3 11 U 15 C PK3.1
NCBP4 25 G 34 G PK4.1
segment1 4bp 3..12 UCAGGUGAUU 24..35 AGUGCAGCUAGA
