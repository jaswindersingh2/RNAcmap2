#Name: 5dcv_B
#Length:  51 
#PageNumber: 2
GGUAUGGUGGAUGAAAGCGGUGAAGGGAAACCGAGUUAACCCGCCUAAGCC
(((..(.(..[....)..(](..[..)....).....]......)...)))
SSSIISMSHHHHHHHSMMSISHHHHHSIIIISMMMMMMMMMMMMSIIISSS
NNNNNNNNNNKNNNNNNNNKNNNKNNNNNNNNNNNNNKNNNNNNNNNNNNN
S1 1..3 "GGU" 49..51 "GCC"
S2 6..6 "G" 45..45 "C"
S3 8..8 "U" 16..16 "A"
S4 19..19 "G" 32..32 "C"
S5 21..21 "U" 27..27 "G"
H1 9..15 "GGAUGAA" (8,16) U:A PK{1}
H2 22..26 "GAAGG" (21,27) U:G PK{2}
I1.1 4..5 "AU" (3,49) U:G 
I1.2 46..48 "UAA" (45,6) C:G 
I2.1 20..20 "G" (19,32) G:C PK{1}
I2.2 28..31 "AAAC" (27,21) G:U 
M1.1 7..7 "G" (6,45) G:C (8,16) U:A 
M1.2 17..18 "GC" (16,8) A:U (19,32) G:C 
M1.3 33..44 "GAGUUAACCCGC" (32,19) C:G (45,6) C:G PK{2}
PK1 1bp 11..11 20..20 H1 9..15 I2.1 20..20
PK2 1bp 24..24 38..38 H2 22..26 M1.3 33..44
PK1.1 11 A 20 G
PK2.1 24 A 38 A
NCBP1 11 A 20 G PK1.1
NCBP2 24 A 38 A PK2.1
segment1 4bp 1..6 GGUAUG 45..51 CUAAGCC
segment2 1bp 8..8 U 16..16 A
segment3 2bp 19..21 GGU 27..32 GAAACC
