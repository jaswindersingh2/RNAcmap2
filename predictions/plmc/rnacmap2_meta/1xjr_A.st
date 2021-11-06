#Name: 1xjr_A
#Length:  47 
#PageNumber: 2
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GGAGUUCACCGAGGCCACGCGGAGUACGAUCGAGGGUACAGUGAAUU
...(((((..(....)...(((...))).....[.(.].).))))).
EEESSSSSMMSHHHHSMMMSSSHHHSSSMMMMMMMSHHHSMSSSSSE
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNKNNNKNNNNNNNNN
S1 4..8 "GUUCA" 42..46 "UGAAU"
S2 11..11 "G" 16..16 "C"
S3 20..22 "CGG" 26..28 "ACG"
S4 36..36 "G" 40..40 "A"
H1 12..15 "AGGC" (11,16) G:C 
H2 23..25 "AGU" (22,26) G:A 
H3 37..39 "UAC" (36,40) G:A PK{1}
M1.1 9..10 "CC" (8,42) A:U (11,16) G:C 
M1.2 17..19 "ACG" (16,11) C:G (20,28) C:G 
M1.3 29..35 "AUCGAGG" (28,20) G:C (36,40) G:A PK{1}
M1.4 41..41 "G" (40,36) A:G (42,8) U:A 
E1 1..3 "GGA" 
E2 47..47 "U" 
PK1 1bp 34..34 38..38 M1.3 29..35 H3 37..39
PK1.1 34 G 38 A
NCBP1 36 G 40 A S4
NCBP2 22 G 26 A S3
NCBP3 34 G 38 A PK1.1
segment1 5bp 4..8 GUUCA 42..46 UGAAU
segment2 1bp 11..11 G 16..16 C
segment3 3bp 20..22 CGG 26..28 ACG
segment4 1bp 36..36 G 40..40 A
