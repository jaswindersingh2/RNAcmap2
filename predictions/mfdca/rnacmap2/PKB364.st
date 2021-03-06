#Name: PKB364
#Length:  80 
#PageNumber: 4
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
UUUUCCGGUGGCCACCCGCCCGGGCCACCGUUACUCCACUCCACUCCUUCGGGACUGGUUUGGAGGAACAUAACAGGGCU
...........[.(.......)(].....)....((.(..[.{..(......).(.<]).).))....}..........>
EEEEEEEEEEEEESHHHHHHHSSHHHHHHSMMMMSSISMMMMMMMSHHHHHHSMSHHHSMSISSEEEEEEEEEEEEEEEE
NNNNNNNNNNNKNNNNNNNNNNNKNNNNNNNNNNNNNNNNKNKNNNNNNNNNNNNNKKNNNNNNNNNNKNNNNNNNNNNK
S1 14..14 "A" 22..22 "G"
S2 23..23 "G" 30..30 "G"
S3 35..36 "UC" 63..64 "GA"
S4 38..38 "A" 61..61 "U"
S5 46..46 "C" 53..53 "G"
S6 55..55 "C" 59..59 "U"
H1 15..21 "CCCGCCC" (14,22) A:G 
H2 24..29 "GCCACC" (23,30) G:G PK{1}
H3 47..52 "CUUCGG" (46,53) C:G 
H4 56..58 "UGG" (55,59) C:U PK{2,4}
I1.1 37..37 "C" (36,63) C:G 
I1.2 62..62 "G" (61,38) U:A 
M1.1 23..22 "" (22,14) G:A (23,30) G:G 
M1.2 31..34 "UUAC" (30,23) G:G (35,64) U:A 
M2.1 39..45 "CUCCACU" (38,61) A:U (46,53) C:G PK{2,3}
M2.2 54..54 "A" (53,46) G:C (55,59) C:U 
M2.3 60..60 "U" (59,55) U:C (61,38) U:A 
E1 1..13 "UUUUCCGGUGGCC" PK{1}
E2 65..80 "GGAACAUAACAGGGCU" PK{3,4}
PK1 1bp 12..12 24..24 E1 1..13 H2 24..29
PK2 1bp 41..41 58..58 M2.1 39..45 H4 56..58
PK3 1bp 43..43 69..69 M2.1 39..45 E2 65..80
PK4 1bp 57..57 80..80 H4 56..58 E2 65..80
PK1.1 12 C 24 G
PK2.1 41 C 58 G
PK3.1 43 A 69 C
PK4.1 57 G 80 U
NCBP1 14 A 22 G S1
NCBP2 55 C 59 U S6
NCBP3 23 G 30 G S2
NCBP4 43 A 69 C PK3.1
segment1 1bp 14..14 A 22..22 G
segment2 1bp 23..23 G 30..30 G
segment3 3bp 35..38 UCCA 61..64 UGGA
segment4 1bp 46..46 C 53..53 G
segment5 1bp 55..55 C 59..59 U
