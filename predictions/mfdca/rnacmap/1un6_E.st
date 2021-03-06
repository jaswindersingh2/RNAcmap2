#Name: 1un6_E
#Length:  61 
#PageNumber: 2
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GCCGGCCACACCUACGGGGCCUGGUUAGUACCUGGGAAACCUGGGAAUACCAGGUGCCGGC
([(.....[..(][..)....(............((...)).]........).....])).
SBSMMMMMMMMSHHHHSMMMMSIIIIIIIIIIIISSHHHSSIIIIIIIIIISMMMMMMSSE
NKNNNNNNKNNNKKNNNNNNNNNNNNNNNNNNNNNNNNNNNNKNNNNNNNNNNNNNNKNNN
S1 1..1 "G" 60..60 "G"
S2 3..3 "C" 59..59 "G"
S3 12..12 "C" 17..17 "G"
S4 22..22 "U" 52..52 "A"
S5 35..36 "GG" 40..41 "CC"
H1 13..16 "UACG" (12,17) C:G PK{2,3}
H2 37..39 "AAA" (36,40) G:C 
B1 2..2 "C" (1,60) G:G (3,59) C:G PK{1}
I1.1 23..34 "GGUUAGUACCUG" (22,52) U:A 
I1.2 42..51 "UGGGAAUACC" (41,35) C:G PK{3}
M1.1 4..11 "GGCCACAC" (3,59) C:G (12,17) C:G PK{2}
M1.2 18..21 "GGCC" (17,12) G:C (22,52) U:A 
M1.3 53..58 "GGUGCC" (52,22) A:U (59,3) G:C PK{1}
E1 61..61 "C" 
PK1 1bp 2..2 58..58 B1 2..2 M1.3 53..58
PK2 1bp 9..9 13..13 M1.1 4..11 H1 13..16
PK3 1bp 14..14 43..43 H1 13..16 I1.2 42..51
PK1.1 2 C 58 C
PK2.1 9 C 13 U
PK3.1 14 A 43 G
NCBP1 1 G 60 G S1
NCBP2 2 C 58 C PK1.1
NCBP3 9 C 13 U PK2.1
NCBP4 14 A 43 G PK3.1
segment1 2bp 1..3 GCC 59..60 GG
segment2 1bp 12..12 C 17..17 G
segment3 3bp 22..36 UGGUUAGUACCUGGG 40..52 CCUGGGAAUACCA
