#Name: 7jju_A
#Length: ,51,
#PageNumber: 4
GGAACUAGCUCAGCAUACACGAGUUGCAAGCAUGGGAAGUUCAAGCCUCGU
.[{((]....((((.[....<<.......)).))....).}.]....>.>)
EEESSIIIIISSSSHHHHHHHHHHHHHHHSSBSSIIIISBBBBBBBBBBBS
NKKNNKNNNNNNNNNKNNNNKKNNNNNNNNNNNNNNNNNNKNKNNNNKNKN
S1 4..4 "A" 51..51 "U"
S2 5..5 "C" 39..39 "G"
S3 11..12 "CA" 33..34 "UG"
S4 13..14 "GC" 30..31 "GC"
H1 15..29 "AUACACGAGUUGCAA" (14,30) C:G PK{3,4}
B1 32..32 "A" (31,13) C:G (33,12) U:A 
B2 40..50 "UUCAAGCCUCG" (39,5) G:C (51,4) U:A PK{2,3,4}
I1.1 6..10 "UAGCU" (5,39) C:G PK{1}
I1.2 35..38 "GGAA" (34,11) G:C 
E1 1..3 "GGA" PK{1,2}
PK1 1bp 2..2 6..6 E1 1..3 I1.1 6..10
PK2 1bp 3..3 41..41 E1 1..3 B2 40..50
PK3 1bp 16..16 43..43 H1 15..29 B2 40..50
PK4 2bp 21..22 48..50 H1 15..29 B2 40..50
PK1.1 2 G 6 U
PK2.1 3 A 41 U
PK3.1 16 U 43 A
PK4.1 21 G 50 G
PK4.2 22 A 48 U
NCBP1 21 G 50 G PK4.1
segment1 6bp 4..14 ACUAGCUCAGC 30..51 GCAUGGGAAGUUCAAGCCUCGU
