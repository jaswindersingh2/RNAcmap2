#Name: 6ck5_A
#Length:  117 
#PageNumber: 3
GUGAAAGUGUACCUAGGGUUCCAGCCUAUUUGUAGGUGUUCGGACCGAGCGGUACAGGUAUAUUUUUAUAUACCACACCUUAGGGACAAAAGCCCGGGAGGAUAGGUUUCACUCGUA
.[(((((((((((.........(((((({...)))))...)........{))))))(((((((]...)))))))...[}...(.(.......).)..(..].).}))))).......
EESSSSSSSSSSSIIIIIIIIISSSSSSHHHHSSSSSBBBSIIIIIIIIISSSSSSSSSSSSSHHHHSSSSSSSMMMMMMMMSISHHHHHHHSISMMSHHHHSMMSSSSSEEEEEEE
NKNNNNNNNNNNNNNNNNNNNNNNNNNNKNNNNNNNNNNNNNNNNNNNNKNNNNNNNNNNNNNKNNNNNNNNNNNNNKKNNNNNNNNNNNNNNNNNNNNNKNNNKNNNNNNNNNNNN
S1 3..7 "GAAAG" 106..110 "GUUUC"
S2 8..13 "UGUACC" 51..56 "GGUACA"
S3 23..23 "A" 41..41 "C"
S4 24..28 "GCCUA" 33..37 "UAGGU"
S5 57..63 "GGUAUAU" 68..74 "AUAUACC"
S6 83..83 "G" 95..95 "C"
S7 85..85 "G" 93..93 "C"
S8 98..98 "G" 103..103 "U"
H1 29..32 "UUUG" (28,33) A:U PK{2}
H2 64..67 "UUUU" (63,68) U:A PK{1}
H3 86..92 "ACAAAAG" (85,93) G:C 
H4 99..102 "AGGA" (98,103) G:U PK{4}
B1 38..40 "GUU" (37,24) U:G (41,23) C:A 
I1.1 14..22 "UAGGGUUCC" (13,51) C:G 
I1.2 42..50 "GGACCGAGC" (41,23) C:A PK{3}
I2.1 84..84 "G" (83,95) G:C 
I2.2 94..94 "C" (93,85) C:G 
M1.1 8..7 "" (7,106) G:G (8,56) U:A 
M1.2 57..56 "" (56,8) A:U (57,74) G:C 
M1.3 75..82 "ACACCUUA" (74,57) C:G (83,95) G:C PK{3,4}
M1.4 96..97 "GG" (95,83) C:G (98,103) G:U 
M1.5 104..105 "AG" (103,98) U:G (106,7) G:G PK{2}
E1 1..2 "GU" PK{1}
E2 111..117 "ACUCGUA" 
PK1 1bp 2..2 64..64 E1 1..2 H2 64..67
PK2 1bp 29..29 105..105 H1 29..32 M1.5 104..105
PK3 1bp 50..50 79..79 I1.2 42..50 M1.3 75..82
PK4 1bp 78..78 101..101 M1.3 75..82 H4 99..102
PK1.1 2 U 64 U
PK2.1 29 U 105 G
PK3.1 50 C 79 C
PK4.1 78 C 101 G
NCBP1 7 G 106 G S1
NCBP2 23 A 41 C S3
NCBP3 2 U 64 U PK1.1
NCBP4 50 C 79 C PK3.1
segment1 5bp 3..7 GAAAG 106..110 GUUUC
segment2 12bp 8..28 UGUACCUAGGGUUCCAGCCUA 33..56 UAGGUGUUCGGACCGAGCGGUACA
segment3 7bp 57..63 GGUAUAU 68..74 AUAUACC
segment4 2bp 83..85 GGG 93..95 CCC
segment5 1bp 98..98 G 103..103 U
