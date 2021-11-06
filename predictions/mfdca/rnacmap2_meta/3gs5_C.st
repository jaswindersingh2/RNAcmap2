#Name: 3gs5_C
#Length:  36 
#PageNumber: 2
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GGCAGAGAAACACACGAUCGUGGUACAUUACCUGCC
((.(...[..(.(((....))))(]..)....).))
SSISMMMMMMSBSSSHHHHSSSSSHHHSMMMMSISS
NNNNNNNKNNNNNNNNNNNNNNNNKNNNNNNNNNNN
S1 1..2 "GG" 35..36 "CC"
S2 4..4 "A" 33..33 "U"
S3 11..11 "C" 23..23 "G"
S4 13..15 "CAC" 20..22 "GUG"
S5 24..24 "U" 28..28 "U"
H1 16..19 "GAUC" (15,20) C:G 
H2 25..27 "ACA" (24,28) U:U PK{1}
B1 12..12 "A" (11,23) C:G (13,22) C:G 
I1.1 3..3 "C" (2,35) G:C 
I1.2 34..34 "G" (33,4) U:A 
M1.1 5..10 "GAGAAA" (4,33) A:U (11,23) C:G PK{1}
M1.2 24..23 "" (23,11) G:C (24,28) U:U 
M1.3 29..32 "UACC" (28,24) U:U (33,4) U:A 
PK1 1bp 8..8 25..25 M1.1 5..10 H2 25..27
PK1.1 8 A 25 A
NCBP1 24 U 28 U S5
NCBP2 8 A 25 A PK1.1
segment1 3bp 1..4 GGCA 33..36 UGCC
segment2 4bp 11..15 CACAC 20..23 GUGG
segment3 1bp 24..24 U 28..28 U
