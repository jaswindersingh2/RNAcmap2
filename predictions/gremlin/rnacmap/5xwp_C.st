#Name: 5xwp_C
#Length: ,59,
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GGACCACCCCAAAAAUGAAGGGGACUAAAACACAAAUCUAUCUGAAUAAACUCUUCUUC
...................(.(...........[(..(..{)])..[.....}...))]
EEEEEEEEEEEEEEEEEEESBSIIIIIIIIIIIISIISHHHSISIIIIIIIIIIIISSE
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNKNNNNNNKNKNNNKNNNNNKNNNNNK
S1 20..20 "G" 58..58 "U"
S2 22..22 "G" 57..57 "U"
S3 35..35 "A" 44..44 "G"
S4 38..38 "C" 42..42 "C"
H1 39..41 "UAU" (38,42) C:C PK{2}
B1 21..21 "G" (20,58) G:U (22,57) G:U 
I1.1 23..34 "GACUAAAACACA" (22,57) G:U PK{1}
I1.2 45..56 "AAUAAACUCUUC" (44,35) G:A PK{2,3}
I2.1 36..37 "AU" (35,44) A:G 
I2.2 43..43 "U" (42,38) C:C PK{1}
E1 1..19 "GGACCACCCCAAAAAUGAA" 
E2 59..59 "C" PK{3}
PK1 1bp 34..34 43..43 I1.1 23..34 I2.2 43..43
PK2 1bp 41..41 53..53 H1 39..41 I1.2 45..56
PK3 1bp 47..47 59..59 I1.2 45..56 E2 59..59
PK1.1 34 A 43 U
PK2.1 41 U 53 C
PK3.1 47 U 59 C
NCBP1 38 C 42 C S4
NCBP2 35 A 44 G S3
NCBP3 41 U 53 C PK2.1
NCBP4 47 U 59 C PK3.1
segment1 4bp 20..38 GGGGACUAAAACACAAAUC 42..58 CUGAAUAAACUCUUCUU
