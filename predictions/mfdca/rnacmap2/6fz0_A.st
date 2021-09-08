#Name: 6fz0_A
#Length:  53 
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
AGGCGCAUUUGAACUGUAUUGUACGCCUUGCAUAAAGCAAAAGUACUAAAAAA
.[.(..................{.)....(..]....).[...}(......])
EEESHHHHHHHHHHHHHHHHHHHHSMMMMSHHHHHHHSMMMMMMSHHHHHHHS
NKNNNNNNNNNNNNNNNNNNNNKNNNNNNNNNKNNNNNNKNNNKNNNNNNNKN
S1 4..4 "C" 25..25 "G"
S2 30..30 "G" 38..38 "C"
S3 45..45 "A" 53..53 "A"
H1 5..24 "GCAUUUGAACUGUAUUGUAC" (4,25) C:G PK{2}
H2 31..37 "CAUAAAG" (30,38) G:C PK{1}
H3 46..52 "CUAAAAA" (45,53) A:A PK{3}
M1.1 26..29 "CCUU" (25,4) G:C (30,38) G:C 
M1.2 39..44 "AAAAGU" (38,30) C:G (45,53) A:A PK{2,3}
E1 1..3 "AGG" PK{1}
PK1 1bp 2..2 33..33 E1 1..3 H2 31..37
PK2 1bp 23..23 44..44 H1 5..24 M1.2 39..44
PK3 1bp 40..40 52..52 M1.2 39..44 H3 46..52
PK1.1 2 G 33 U
PK2.1 23 A 44 U
PK3.1 40 A 52 A
NCBP1 45 A 53 A S3
NCBP2 40 A 52 A PK3.1
segment1 1bp 4..4 C 25..25 G
segment2 1bp 30..30 G 38..38 C
segment3 1bp 45..45 A 53..53 A
