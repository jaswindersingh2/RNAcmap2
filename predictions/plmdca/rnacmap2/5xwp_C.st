#Name: 5xwp_C
#Length:  59 
#PageNumber: 4
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GGACCACCCCAAAAAUGAAGGGGACUAAAACACAAAUCUAUCUGAAUAAACUCUUCUUC
....[..[......(.]..[{...(]..)...)(...<....).(......].}...>)
EEEEEEEEEEEEEESIIIIIIIIISHHHSIIISSHHHHHHHHSMSHHHHHHHHHHHHHS
NNNNKNNKNNNNNNNNKNNKKNNNNKNNNNNNNNNNNKNNNNNNNNNNNNNKNKNNNKN
S1 15..15 "A" 33..33 "C"
S2 25..25 "C" 29..29 "A"
S3 34..34 "A" 43..43 "U"
S4 45..45 "A" 59..59 "C"
H1 26..28 "UAA" (25,29) C:A PK{3}
H2 35..42 "AAUCUAUC" (34,43) A:U PK{5}
H3 46..58 "AUAAACUCUUCUU" (45,59) A:C PK{1,4,5}
I1.1 16..24 "UGAAGGGGA" (15,33) A:C PK{2,3,4}
I1.2 30..32 "ACA" (29,25) A:C 
M1.1 34..33 "" (33,15) C:A (34,43) A:U 
M1.2 44..44 "G" (43,34) U:A (45,59) A:C 
E1 1..14 "GGACCACCCCAAAA" PK{1,2}
PK1 1bp 5..5 52..52 E1 1..14 H3 46..58
PK2 1bp 8..8 17..17 E1 1..14 I1.1 16..24
PK3 1bp 20..20 26..26 I1.1 16..24 H1 26..28
PK4 1bp 21..21 54..54 I1.1 16..24 H3 46..58
PK5 1bp 38..38 58..58 H2 35..42 H3 46..58
PK1.1 5 C 52 U
PK2.1 8 C 17 G
PK3.1 20 G 26 U
PK4.1 21 G 54 U
PK5.1 38 C 58 U
NCBP1 45 A 59 C S4
NCBP2 25 C 29 A S2
NCBP3 15 A 33 C S1
NCBP4 5 C 52 U PK1.1
NCBP5 38 C 58 U PK5.1
segment1 2bp 15..25 AUGAAGGGGAC 29..33 AACAC
segment2 1bp 34..34 A 43..43 U
segment3 1bp 45..45 A 59..59 C
