#Name: PKB48
#Length:  61 
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
UAGGGCCGUCAGGGUUCUCCCGCCCUCCGUGAACCCAGGCUAAAAGUUAAGGUAGGGGGGC
[.(..{.........]).[}.(((.(.(........(......]...[)...)]..).)))
EESHHHHHHHHHHHHHSXXXXSSSISISIIIIIIIISHHHHHHHHHHHSIIISIIISISSS
KNNNNKNNNNNNNNNKNNKKNNNNNNNNNNNNNNNNNNNNNNNKNNNKNNNNNKNNNNNNN
S1 3..3 "G" 17..17 "C"
S2 22..24 "GCC" 59..61 "GGC"
S3 26..26 "U" 57..57 "G"
S4 28..28 "C" 53..53 "U"
S5 37..37 "A" 49..49 "A"
H1 4..16 "GGCCGUCAGGGUU" (3,17) G:C PK{1,2}
H2 38..48 "GGCUAAAAGUU" (37,49) A:A PK{3,4}
I1.1 25..25 "C" (24,59) C:G 
I1.2 58..58 "G" (57,26) G:U 
I2.1 27..27 "C" (26,57) U:G 
I2.2 54..56 "AGG" (53,28) U:C PK{4}
I3.1 29..36 "GUGAACCC" (28,53) C:U 
I3.2 50..52 "AGG" (49,37) A:A 
X1 18..21 "UCCC" (17,3) C:G (22,61) G:C PK{2,3}
E1 1..2 "UA" PK{1}
PK1 1bp 1..1 16..16 E1 1..2 H1 4..16
PK2 1bp 6..6 20..20 H1 4..16 X1 18..21
PK3 1bp 19..19 44..44 X1 18..21 H2 38..48
PK4 1bp 48..48 54..54 H2 38..48 I2.2 54..56
PK1.1 1 U 16 U
PK2.1 6 C 20 C
PK3.1 19 C 44 A
PK4.1 48 U 54 A
NCBP1 37 A 49 A S5
NCBP2 28 C 53 U S4
NCBP3 1 U 16 U PK1.1
NCBP4 6 C 20 C PK2.1
NCBP5 19 C 44 A PK3.1
segment1 1bp 3..3 G 17..17 C
segment2 6bp 22..37 GCCCUCCGUGAACCCA 49..61 AAGGUAGGGGGGC
