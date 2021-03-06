#Name: 4pkd_V
#Length: ,55,
#PageNumber: 2
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GGAUCCAUUGCACUCCGGAUCCAGGAGAUACCAUGAUCACGAAGGUGGUUUUCCU
.......[........[.....(](.(.[(..(......].....)].).).).)
EEEEEEEEEEEEEEEEEEEEEESISISIISIISHHHHHHHHHHHHSIISISISIS
NNNNNNNKNNNNNNNNKNNNNNNKNNNNKNNNNNNNNNNKNNNNNNKNNNNNNNN
S1 23..23 "A" 55..55 "U"
S2 25..25 "G" 53..53 "C"
S3 27..27 "G" 51..51 "U"
S4 30..30 "A" 49..49 "U"
S5 33..33 "A" 46..46 "U"
H1 34..45 "UGAUCACGAAGG" (33,46) A:U PK{3}
I1.1 24..24 "G" (23,55) A:U PK{2}
I1.2 54..54 "C" (53,25) C:G 
I2.1 26..26 "A" (25,53) G:C 
I2.2 52..52 "U" (51,27) U:G 
I3.1 28..29 "AU" (27,51) G:U PK{3}
I3.2 50..50 "U" (49,30) U:A 
I4.1 31..32 "CC" (30,49) A:U 
I4.2 47..48 "GG" (46,33) U:A PK{1}
E1 1..22 "GGAUCCAUUGCACUCCGGAUCC" PK{1,2}
PK1 1bp 8..8 47..47 E1 1..22 I4.2 47..48
PK2 1bp 17..17 24..24 E1 1..22 I1.1 24..24
PK3 1bp 29..29 40..40 I3.1 28..29 H1 34..45
PK1.1 8 U 47 G
PK2.1 17 G 24 G
PK3.1 29 U 40 C
NCBP1 17 G 24 G PK2.1
NCBP2 29 U 40 C PK3.1
segment1 5bp 23..33 AGGAGAUACCA 46..55 UGGUUUUCCU
