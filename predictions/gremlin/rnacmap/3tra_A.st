#Name: 3tra_A
#Length:  75 
#PageNumber: 2
#Warning: Structure contains linked PK-segments of same sizes 2 and 2. Using PK brackets for the more 5' segment
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
UCCGUGAUAGUUUAAUGGUCAGAAUGGGCGCUUGUCGCGUGCCAGAUCGGGGUUCAAUUCCCCGUCGCGGAGCCA
(((((((...[[......((..]](.([(((.....)).]).)))).((((.........)))))))))))....
SSSSSSSMMMMMMMMMMMSSBBBBSBSISSSHHHHHSSBBSISSSSMSSSSHHHHHHHHHSSSSSSSSSSSEEEE
NNNNNNNNNNKKNNNNNNNNNNKKNNNKNNNNNNNNNNNKNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
S1 1..7 "UCCGUGA" 65..71 "UCGCGGA"
S2 19..20 "UC" 45..46 "GA"
S3 25..25 "U" 44..44 "A"
S4 27..27 "G" 43..43 "C"
S5 29..29 "C" 41..41 "G"
S6 30..31 "GC" 37..38 "GC"
S7 48..51 "CGGG" 61..64 "CCCG"
H1 32..36 "UUGUC" (31,37) C:G 
H2 52..60 "GUUCAAUUC" (51,61) G:C 
B1 21..24 "AGAA" (20,45) C:G (25,44) U:A PK{1}
B2 26..26 "G" (25,44) U:A (27,43) G:C 
B3 39..40 "GU" (38,30) C:G (41,29) G:C PK{2}
I1.1 28..28 "G" (27,43) G:C PK{2}
I1.2 42..42 "C" (41,29) G:C 
M1.1 8..18 "UAGUUUAAUGG" (7,65) A:U (19,46) U:A PK{1}
M1.2 47..47 "U" (46,19) A:U (48,64) C:G 
M1.3 65..64 "" (64,48) G:C (65,7) U:A 
E1 72..75 "GCCA" 
PK1 2bp 11..12 23..24 M1.1 8..18 B1 21..24
PK2 1bp 28..28 40..40 I1.1 28..28 B3 39..40
PK1.1 11 U 24 A
PK1.2 12 U 23 A
PK2.1 28 G 40 U
segment1 7bp 1..7 UCCGUGA 65..71 UCGCGGA
segment2 7bp 19..31 UCAGAAUGGGCGC 37..46 GCGUGCCAGA
segment3 4bp 48..51 CGGG 61..64 CCCG
