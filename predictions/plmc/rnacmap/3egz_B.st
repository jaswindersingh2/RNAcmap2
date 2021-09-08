#Name: 3egz_B
#Length:  65 
#PageNumber: 4
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GAGGGAGAGGUGAAGAAUACGACCACCUAGGUACCAUUGCACUCCGGUACCUAAAACAUACCCUC
[..[..(.....(.........)[..(.{....(.<.....]]...).(..}....)..)].)>.
EEEEEESMMMMMSHHHHHHHHHSMMMSMMMMMMSHHHHHHHHHHHHSMSHHHHHHHSMMSMMSEE
KNNKNNNNNNNNNNNNNNNNNNNKNNNNKNNNNNNKNNNNNKKNNNNNNNNKNNNNNNNNKNNKN
S1 7..7 "G" 63..63 "C"
S2 13..13 "A" 23..23 "C"
S3 27..27 "C" 60..60 "A"
S4 34..34 "C" 47..47 "G"
S5 49..49 "A" 57..57 "C"
H1 14..22 "AGAAUACGA" (13,23) A:C 
H2 35..46 "CAUUGCACUCCG" (34,47) C:G PK{2,3,5}
H3 50..56 "CCUAAAA" (49,57) A:C PK{4}
M1.1 8..12 "AGGUG" (7,63) G:C (13,23) A:C 
M1.2 24..26 "CAC" (23,13) C:A (27,60) C:A PK{3}
M1.3 61..62 "CC" (60,27) A:C (63,7) C:G PK{1}
M2.1 28..33 "UAGGUA" (27,60) C:A (34,47) C:G PK{4}
M2.2 48..48 "U" (47,34) G:C (49,57) A:C 
M2.3 58..59 "AU" (57,49) C:A (60,27) A:C 
E1 1..6 "GAGGGA" PK{1,2}
E2 64..65 "UC" PK{5}
PK1 1bp 1..1 61..61 E1 1..6 M1.3 61..62
PK2 1bp 4..4 43..43 E1 1..6 H2 35..46
PK3 1bp 24..24 42..42 M1.2 24..26 H2 35..46
PK4 1bp 29..29 52..52 M2.1 28..33 H3 50..56
PK5 1bp 36..36 64..64 H2 35..46 E2 64..65
PK1.1 1 G 61 C
PK2.1 4 G 43 U
PK3.1 24 C 42 C
PK4.1 29 A 52 U
PK5.1 36 A 64 U
NCBP1 27 C 60 A S3
NCBP2 13 A 23 C S2
NCBP3 49 A 57 C S5
NCBP4 24 C 42 C PK3.1
segment1 1bp 7..7 G 63..63 C
segment2 1bp 13..13 A 23..23 C
segment3 1bp 27..27 C 60..60 A
segment4 1bp 34..34 C 47..47 G
segment5 1bp 49..49 A 57..57 C
