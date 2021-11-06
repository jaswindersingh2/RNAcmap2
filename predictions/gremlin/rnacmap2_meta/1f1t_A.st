#Name: 1f1t_A
#Length:  38 
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GGAUCCCGACUGGCGAGAGCCAGGUAACGAAUGGAUCC
.(.((([...(....{.....)....(.]..))))}).
ESISSSMMMMSHHHHHHHHHHSMMMMSHHHHSSSSISE
NNNNNNKNNNNNNNNKNNNNNNNNNNNNKNNNNNNKNN
S1 2..2 "G" 37..37 "C"
S2 4..6 "UCC" 33..35 "GGA"
S3 11..11 "U" 22..22 "A"
S4 27..27 "A" 32..32 "U"
H1 12..21 "GGCGAGAGCC" (11,22) U:A PK{2}
H2 28..31 "CGAA" (27,32) A:U PK{1}
I1.1 3..3 "A" (2,37) G:C 
I1.2 36..36 "U" (35,4) A:U PK{2}
M1.1 7..10 "CGAC" (6,33) C:G (11,22) U:A PK{1}
M1.2 23..26 "GGUA" (22,11) A:U (27,32) A:U 
M1.3 33..32 "" (32,27) U:A (33,6) G:C 
E1 1..1 "G" 
E2 38..38 "C" 
PK1 1bp 7..7 29..29 M1.1 7..10 H2 28..31
PK2 1bp 16..16 36..36 H1 12..21 I1.2 36..36
PK1.1 7 C 29 G
PK2.1 16 A 36 U
segment1 4bp 2..6 GAUCC 33..37 GGAUC
segment2 1bp 11..11 U 22..22 A
segment3 1bp 27..27 A 32..32 U
