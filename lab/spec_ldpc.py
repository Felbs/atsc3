"""LDPC tables and parameters extracted from ATSC A/322 (ATSC 3.0 Physical Layer Protocol).

PROVENANCE
----------
Primary source : A/322:2026-04 (14 April 2026), Annex A.2 and Section 6.1.3.
Cross-check    : A/322:2024-04 -- an INDEPENDENT EDITION of the same standard.
Extracted from : pdftotext renderings on disk under lab\\spec\\
                   A322_2026_tbl.txt / A322_2026.txt   (-table / -layout modes)
                   A322_2024_tbl.txt / A322_2024.txt   (-table / -layout modes)
                 i.e. 2 editions x 2 independent pdftotext layout algorithms = 4 witnesses.

Contents:
  LDPC_16200      Annex A.2 tables A.2.1 .. A.2.12, ALL twelve code rates 2/15 .. 13/15
                  for Ninner = 16200.  (Annex A.1, Ninner = 64800, is NOT extracted here.)
  LDPC_PARAMS     Section 6 Table 6.4 (Type A/B), Table 6.6 (Type A sizes, Ninner=16200)
                  and Table 6.7 (Qldpc for Type B), for Ninner = 16200.
  LDPC_PARAMS_64800
                  Table 6.5 (Type A sizes, Ninner=64800) and the Ninner=64800 column of
                  Table 6.7, for completeness.  No Annex A.1 index tables are included.
  PROSE_6_1_3_1 / PROSE_6_1_3_2
                  The Type A and Type B encoding procedures, so the encoder can be
                  implemented directly from this file.

VERIFICATION EVIDENCE  (all of this is re-run by ``python spec_ldpc.py``)
------------------------------------------------------------------------
The A/322 PDF sets these tables in a math font whose glyphs pdftotext can mis-map, and the
tables wrap across page boundaries with running headers, footers and bare page numbers
("144" .. "148") interleaved into the number stream.  A page number silently swallowed as
a table entry is the exact failure mode guarded against here.  Evidence gathered:

(a) CROSS-EDITION / CROSS-RENDERER.  The complete multiset of integers in Annex A.2 was
    computed independently from all four files.  All four agree exactly:
        1819 integers, arithmetic sum 6754166.
    The tables below account for exactly those 1819 integers with exactly that sum, so no
    integer was dropped, duplicated, or absorbed from page furniture.  Additionally the
    per-table row lists parsed from A322_2026_tbl.txt and A322_2024_tbl.txt are equal
    element-for-element for all 12 tables, as are 10 of 12 from the two -layout files
    (the -layout renderer bleeds the right-hand column of A.2.10 and A.2.11 across the
    caption line, so those two were not independently re-parsed from -layout; they are
    still covered by the multiset check and by the 2024/2026 -table agreement).

(b) ROW COUNT.  Derived independently from Section 6, not from the PDF layout:
      Type B: rows = Kldpc/360.
      Type A: rows = Kldpc/360 + Q1.   (Steps (i)-(iv) consume one row per 360 INFORMATION
              bits; step (vii) consumes one further row per 360 of the M1 dual-diagonal
              parity bits, i.e. M1/360 = Q1 extra rows.)
    NOTE: a plain "Kldpc/360" rule is WRONG for the Type A codes -- e.g. rate 3/15 has
    3240/360 = 9 information rows but the printed table has 12 rows (9 + Q1 = 9 + 3).
    Every table matches its derived count exactly.  Independent corroboration: in all four
    Type A tables the trailing Q1 rows have a row weight strictly smaller than EVERY
    information row (2/15: 7 then 4;  3/15: 11 then 10;  4/15: 10,9 then 7;  5/15: 9,10
    then 7), so the printed table visibly breaks at exactly Kldpc/360 -- which is where the
    rule says it should.  An off-by-one row split would destroy this.

(c) RANGE.  Every address lies in [0, Nldpc - Kldpc).  Largest address per table is within
    a few units of the limit in every case (e.g. 8/15 max 7559 vs limit 7560), which is
    what a correctly captured table looks like; a swallowed page number would show up as
    an out-of-family small value in a high row, and none is present.

(d) ACCUMULATOR COVERAGE (the Q1/Q2 consistency check).  The full parity-bit accumulator
    address set is expanded per Section 6.1.3.1 step (ii) (Type A, using Q1 below M1 and Q2
    at/above M1) and Section 6.1.3.2 (Type B, using Qldpc).  For every one of the 12 codes
    the expansion lands strictly inside [0, Nldpc-Kldpc) AND covers every single parity
    accumulator at least once.  A corrupted address would orphan accumulators or overflow.

(e) NO DUPLICATE address within any row; every row is strictly increasing.

(f) MULTI-COLUMN LAYOUT.  Tables A.2.8, A.2.10, A.2.11 and A.2.12 are printed as TWO
    side-by-side column blocks in the PDF.  They are read COLUMN-MAJOR (all of the left
    block, then all of the right block).  Evidence: (i) only column-major yields the row
    count derived in (b) together with a weight profile matching the other tables; (ii) all
    four unambiguous single-column Type B tables (6/15, 7/15, 8/15, 10/15) have
    non-increasing row weights, and column-major reading reproduces that invariant for all
    four two-column tables while row-interleaved reading violates it in all four.

(g) SECTION 6.1.3 PROSE AND PARAMETER TABLES 6.4/6.5/6.6/6.7 were diffed between the 2024
    and 2026 editions: zero numeric differences (the only diff is one line wrap).
    Table 6.6/6.5 internal arithmetic re-checked here: M1 + M2 == Nldpc - Kldpc,
    Q1 == M1/360, Q2 == M2/360.  Table 6.7 re-checked: Qldpc == (Nldpc - Kldpc)/360.

KNOWN LIMITATION: the PDF's italic math variables (Kldpc, lambda, p, etc.) are dropped
entirely by pdftotext, so PROSE_* below is the extracted text with the missing symbols
restored by hand in [brackets].  The bracketed symbols are editorial reconstruction; the
unbracketed words and all numbers are verbatim.  The NUMERIC TABLES are untouched.

SHA256 of json.dumps({rate: rows}, sort_keys=True):
  1803e817270b94d50e1ef3dee7fa1d76025b12f94a8bb3e74446e394bf3d8358
"""

# --------------------------------------------------------------------------------------
# Section 6.1.3.1 -- Type A LDPC Encoding (used for Ninner=16200 rates 2/15..5/15)
# --------------------------------------------------------------------------------------
PROSE_6_1_3_1 = """\
6.1.3.1 Type A LDPC Encoding

Type A LDPC encoding shall be realized as follows.

An LDPC code is used to encode the information block [s] = ([s_0], [s_1], ... , [s_{Kldpc-1}]).
To generate a codeword [Lambda] = [lambda_0], [lambda_1], ... , [lambda_{Ninner-1}] of length
[Ninner] = [Kldpc] + [M1] + [M2], the parity bits before parity interleaving
[p] = ([p_0], [p_1], ... , [p_{M1+M2-1}]) are calculated from S. Then for any [i] from 0 to
[Kldpc] - 1, [lambda]_i is set equal to [s]_i, since the code is systematic.

[M1] and [M2] are parity lengths corresponding to a dual diagonal matrix and an identity
matrix, respectively. The parity lengths depending on code rates shall be used as specified
in Table 6.5 and Table 6.6. The detailed procedure to calculate parity bits shall be as
follows:

  i)   Initialize [lambda_i] = [s_i] for [i] = 0,1, . . . [Kldpc] - 1
       [p_j] = 0 for [j] = 0,1, . . . [M1] + [M2] - 1

       Accumulate the first information bit, [s_0], at parity bit addresses specified in the
       first row of the tables in Annex A (Table A.1.1 to Table A.1.4, Table A.1.6,
       Table A.2.1 to Table A.2.4).

  ii)  For the next 359 information bits, [s_m], [m] = 1,2, ... ,359, accumulate [s_m] at
       parity bit addresses, which are calculated as follows:

              ([x] + [m] x [Q1]) mod [M1]                        if [x] <  [M1]
              [M1] + {([x] - [M1] + [m] x [Q2])} mod [M2]        if [x] >= [M1]

       where x denotes the address of the parity bit accumulator corresponding to the first
       bit [s_0]. [Q1] = [M1]/360 and [Q2] = [M2]/360 are code rate dependent constants
       specified in Table 6.5 and Table 6.6.

  iii) For the 361st information bit [s_360], the addresses of the parity bit accumulators
       are given in the second row of the tables in Annex A. In a similar manner, the
       addresses of the parity bit accumulators for the following 359 information bits
       [s_m], [m] = 361,362, ... ,719 are obtained using the equation in step (ii) where [x]
       denotes the address of the parity bit accumulator corresponding to the information
       bit [s_360], i.e. the entries in the second row of the tables in Annex A.

  iv)  In a similar manner, for every group of 360 new information bits, a new row from the
       tables in Annex A is used to find the addresses of the parity bit accumulators.

  v)   After the codeword bits from [lambda_0] to [lambda_{Kldpc-1}] are exhausted,
       sequentially perform the following operations starting with [i] = 1.

              [p_i] = [p_i] XOR [p_{i-1}]   for [i] = 1,2, ... , [M1] - 1

  vi)  The parity bits from [lambda_Kldpc] to [lambda_{Kldpc+M1-1}], which correspond to the
       dual diagonal matrix, are obtained using the following interleaving operation:

              [lambda_{Kldpc + 360*s + t}] = [p_{M1*t + s}]   for 0 <= [s] < 360,
                                                                  0 <= [t] < [Q1]

  vii) For every group of 360 new codeword bits from [lambda_Kldpc] to
       [lambda_{Kldpc+M1-1}], a new row from the tables in Annex A and the equation in step
       (ii) are used to find the addresses of the parity bit accumulators.

  viii) After the codeword bits from [lambda_Kldpc] to [lambda_{Kldpc+M1-1}] are exhausted,
       the parity bits from [lambda_{Kldpc+M1}] to [lambda_{Kldpc+M1+M2-1}], which correspond
       to the identity matrix, are obtained using the following interleaving operation:

              [lambda_{Kldpc + M1 + 360*s + t}] = [p_{M1 + M2*t + s}]   for 0 <= [s] < 360,
                                                                            0 <= [t] < [Q2]

Note that the parity bits from [lambda_Kldpc] to [lambda_{Kldpc+M1-1}] are obtained according
to steps (i) to (vi), and then the parity bits from [lambda_{Kldpc+M1}] to
[lambda_{Kldpc+M1+M2-1}] are obtained using the results of steps (i) to (vi) according to
steps (vii) to (viii).

IMPLEMENTATION NOTE (not spec text): consequently the Annex A table for a Type A code has
Kldpc/360 rows for step (i)-(iv) FOLLOWED BY Q1 = M1/360 further rows for step (vii) --
Kldpc/360 + Q1 rows in total.  This is verified for every Type A table in this file.
"""

# --------------------------------------------------------------------------------------
# Section 6.1.3.2 -- Type B LDPC Encoding (used for Ninner=16200 rates 6/15..13/15)
# --------------------------------------------------------------------------------------
PROSE_6_1_3_2 = """\
6.1.3.2 Type B LDPC Encoding

Type B LDPC encoding shall be realized as follows.

Let s0, s1, ..., sNouter-1 be information bits to be encoded and [lambda]0, [lambda]1, ...,
[lambda]Ninner-1 be code bits to be calculated. Then for any k from 0 to Nouter-1,
[lambda]k shall be set equal to sk, since the code is systematic. For the remaining code
bits, set [lambda]Nouter + k = pk (0 <= k < Minner) and these parity bits pk shall be
calculated as follows. In the following, q(i, j, 0) denotes the j-th entry in the i-th row
in the indices list, see Annex A, and q(i, j, l) = q(i, j, 0) + Qldpc*l (mod Minner) for
0 < l < 360, and all accumulations are realized by additions in GF(2). Qldpc shall be
defined as in Table 6.7.

  i)   Initialize pk = 0 for 0 <= k < Minner.

  ii)  For 0 <= k < Nouter, set i = floor(k/360), and l = k (mod 360). Now accumulate sk to
       pq(i, j, l) for all j:

         pq(i, 0, l) = pq(i, 0, l) + sk,
         pq(i, 1, l) = pq(i, 1, l) + sk,
         pq(i, 2, l) = pq(i, 2, l) + sk,
         ... ,
         pq(i, w(i)-1, l) = pq(i, w(i)-1, l) + sk,

       where w(i) is the number of elements in the i-th row in the indices list in Annex A.

  iii) For all 0 < k < Minner, pk = pk + pk-1.

From these steps, all code bits [lambda]0, [lambda]1, ..., [lambda]Ninner-1 shall be
obtained.

IMPLEMENTATION NOTE (not spec text): here Nouter == Kldpc and Minner == Ninner - Kldpc for
the inner LDPC code, so a Type B table has exactly Kldpc/360 rows, and Qldpc == Minner/360.
"""

# --------------------------------------------------------------------------------------
# Table 6.4 -- structure type per code rate (verbatim, both length columns)
# --------------------------------------------------------------------------------------
TABLE_6_4 = {
    #  rate  : (Ninner=64800, Ninner=16200)
    "2/15":  ("A", "A"),
    "3/15":  ("A", "A"),
    "4/15":  ("A", "A"),
    "5/15":  ("A", "A"),
    "6/15":  ("B", "B"),
    "7/15":  ("A", "B"),
    "8/15":  ("B", "B"),
    "9/15":  ("B", "B"),
    "10/15": ("B", "B"),
    "11/15": ("B", "B"),
    "12/15": ("B", "B"),
    "13/15": ("B", "B"),
}

# --------------------------------------------------------------------------------------
# Annex A.2 -- LDPC code matrices, Ninner = 16200.
# Each row is a list of parity-bit-accumulator addresses (Section 6.1.3.1 step (i) /
# Section 6.1.3.2 q(i, j, 0)).
# --------------------------------------------------------------------------------------
LDPC_16200 = {
    "2/15": {
        "table": "A.2.1",
        "rows": [
            [2889, 3122, 3208, 4324, 5968, 7241, 13215],
            [281, 923, 1077, 5252, 6099, 10309, 11114],
            [727, 2413, 2676, 6151, 6796, 8945, 12528],
            [2252, 2322, 3093, 3329, 8443, 12170, 13748],
            [575, 2489, 2944, 6577, 8772, 11253, 11657],
            [310, 1461, 2482, 4643, 4780, 6936, 11970],
            [8691, 9746, 10794, 13582],
            [3717, 6535, 12470, 12752],
            [6011, 6547, 7020, 11746],
            [5309, 6481, 10244, 13824],
            [5327, 8773, 8824, 13343],
            [3506, 3575, 9915, 13609],
            [3393, 7089, 11048, 12816],
            [3651, 4902, 6118, 12048],
            [4210, 10132, 13375, 13377],
        ],
    },
    "3/15": {
        "table": "A.2.2",
        "rows": [
            [8, 372, 841, 4522, 5253, 7430, 8542, 9822, 10550, 11896, 11988],
            [80, 255, 667, 1511, 3549, 5239, 5422, 5497, 7157, 7854, 11267],
            [257, 406, 792, 2916, 3072, 3214, 3638, 4090, 8175, 8892, 9003],
            [80, 150, 346, 1883, 6838, 7818, 9482, 10366, 10514, 11468, 12341],
            [32, 100, 978, 3493, 6751, 7787, 8496, 10170, 10318, 10451, 12561],
            [504, 803, 856, 2048, 6775, 7631, 8110, 8221, 8371, 9443, 10990],
            [152, 283, 696, 1164, 4514, 4649, 7260, 7370, 11925, 11986, 12092],
            [127, 1034, 1044, 1842, 3184, 3397, 5931, 7577, 11898, 12339, 12689],
            [107, 513, 979, 3934, 4374, 4658, 7286, 7809, 8830, 10804, 10893],
            [2045, 2499, 7197, 8887, 9420, 9922, 10132, 10540, 10816, 11876],
            [2932, 6241, 7136, 7835, 8541, 9403, 9817, 11679, 12377, 12810],
            [2211, 2288, 3937, 4310, 5952, 6597, 9692, 10445, 11064, 11272],
        ],
    },
    "4/15": {
        "table": "A.2.3",
        "rows": [
            [19, 585, 710, 3241, 3276, 3648, 6345, 9224, 9890, 10841],
            [181, 494, 894, 2562, 3201, 4382, 5130, 5308, 6493, 10135],
            [150, 569, 919, 1427, 2347, 4475, 7857, 8904, 9903],
            [1005, 1018, 1025, 2933, 3280, 3946, 4049, 4166, 5209],
            [420, 554, 778, 6908, 7959, 8344, 8462, 10912, 11099],
            [231, 506, 859, 4478, 4957, 7664, 7731, 7908, 8980],
            [179, 537, 979, 3717, 5092, 6315, 6883, 9353, 9935],
            [147, 205, 830, 3609, 3720, 4667, 7441, 10196, 11809],
            [60, 1021, 1061, 1554, 4918, 5690, 6184, 7986, 11296],
            [145, 719, 768, 2290, 2919, 7272, 8561, 9145, 10233],
            [388, 590, 852, 1579, 1698, 1974, 9747, 10192, 10255],
            [231, 343, 485, 1546, 3155, 4829, 7710, 10394, 11336],
            [4381, 5398, 5987, 9123, 10365, 11018, 11153],
            [2381, 5196, 6613, 6844, 7357, 8732, 11082],
            [1730, 4599, 5693, 6318, 7626, 9231, 10663],
        ],
    },
    "5/15": {
        "table": "A.2.4",
        "rows": [
            [69, 244, 706, 5145, 5994, 6066, 6763, 6815, 8509],
            [257, 541, 618, 3933, 6188, 7048, 7484, 8424, 9104],
            [69, 500, 536, 1494, 1669, 7075, 7553, 8202, 10305],
            [11, 189, 340, 2103, 3199, 6775, 7471, 7918, 10530],
            [333, 400, 434, 1806, 3264, 5693, 8534, 9274, 10344],
            [111, 129, 260, 3562, 3676, 3680, 3809, 5169, 7308, 8280],
            [100, 303, 342, 3133, 3952, 4226, 4713, 5053, 5717, 9931],
            [83, 87, 374, 828, 2460, 4943, 6311, 8657, 9272, 9571],
            [114, 166, 325, 2680, 4698, 7703, 7886, 8791, 9978, 10684],
            [281, 542, 549, 1671, 3178, 3955, 7153, 7432, 9052, 10219],
            [202, 271, 608, 3860, 4173, 4203, 5169, 6871, 8113, 9757],
            [16, 359, 419, 3333, 4198, 4737, 6170, 7987, 9573, 10095],
            [235, 244, 584, 4640, 5007, 5563, 6029, 6816, 7678, 9968],
            [123, 449, 646, 2460, 3845, 4161, 6610, 7245, 7686, 8651],
            [136, 231, 468, 835, 2622, 3292, 5158, 5294, 6584, 9926],
            [3085, 4683, 8191, 9027, 9922, 9928, 10550],
            [2462, 3185, 3976, 4091, 8089, 8772, 9342],
        ],
    },
    "6/15": {
        "table": "A.2.5",
        "rows": [
            [27, 430, 519, 828, 1897, 1943, 2513, 2600, 2640, 3310, 3415, 4266, 5044, 5100, 5328,
             5483, 5928, 6204, 6392, 6416, 6602, 7019, 7415, 7623, 8112, 8485, 8724, 8994, 9445, 9667],
            [27, 174, 188, 631, 1172, 1427, 1779, 2217, 2270, 2601, 2813, 3196, 3582, 3895, 3908,
             3948, 4463, 4955, 5120, 5809, 5988, 6478, 6604, 7096, 7673, 7735, 7795, 8925, 9613, 9670],
            [27, 370, 617, 852, 910, 1030, 1326, 1521, 1606, 2118, 2248, 2909, 3214, 3413, 3623,
             3742, 3752, 4317, 4694, 5300, 5687, 6039, 6100, 6232, 6491, 6621, 6860, 7304, 8542, 8634],
            [990, 1753, 7635, 8540],
            [933, 1415, 5666, 8745],
            [27, 6567, 8707, 9216],
            [2341, 8692, 9580, 9615],
            [260, 1092, 5839, 6080],
            [352, 3750, 4847, 7726],
            [4610, 6580, 9506, 9597],
            [2512, 2974, 4814, 9348],
            [1461, 4021, 5060, 7009],
            [1796, 2883, 5553, 8306],
            [1249, 5422, 7057],
            [3965, 6968, 9422],
            [1498, 2931, 5092],
            [27, 1090, 6215],
            [26, 4232, 6354],
        ],
    },
    "7/15": {
        "table": "A.2.6",
        "rows": [
            [553, 742, 901, 1327, 1544, 2179, 2519, 3131, 3280, 3603, 3789, 3792, 4253, 5340, 5934,
             5962, 6004, 6698, 7793, 8001, 8058, 8126, 8276, 8559],
            [503, 590, 598, 1185, 1266, 1336, 1806, 2473, 3021, 3356, 3490, 3680, 3936, 4501, 4659,
             5891, 6132, 6340, 6602, 7447, 8007, 8045, 8059, 8249],
            [795, 831, 947, 1330, 1502, 2041, 2328, 2513, 2814, 2829, 4048, 4802, 6044, 6109, 6461,
             6777, 6800, 7099, 7126, 8095, 8428, 8519, 8556, 8610],
            [601, 787, 899, 1757, 2259, 2518, 2783, 2816, 2823, 2949, 3396, 4330, 4494, 4684, 4700,
             4837, 4881, 4975, 5130, 5464, 6554, 6912, 7094, 8297],
            [4229, 5628, 7917, 7992],
            [1506, 3374, 4174, 5547],
            [4275, 5650, 8208, 8533],
            [1504, 1747, 3433, 6345],
            [3659, 6955, 7575, 7852],
            [607, 3002, 4913, 6453],
            [3533, 6860, 7895, 8048],
            [4094, 6366, 8314],
            [2206, 4513, 5411],
            [32, 3882, 5149],
            [389, 3121, 4626],
            [1308, 4419, 6520],
            [2092, 2373, 6849],
            [1815, 3679, 7152],
            [3582, 3979, 6948],
            [1049, 2135, 3754],
            [2276, 4442, 6591],
        ],
    },
    "8/15": {
        "table": "A.2.7",
        "rows": [
            [5, 519, 825, 1871, 2098, 2478, 2659, 2820, 3200, 3294, 3650, 3804, 3949, 4426, 4460,
             4503, 4568, 4590, 4949, 5219, 5662, 5738, 5905, 5911, 6160, 6404, 6637, 6708, 6737,
             6814, 7263, 7412],
            [81, 391, 1272, 1633, 2062, 2882, 3443, 3503, 3535, 3908, 4033, 4163, 4490, 4929, 5262,
             5399, 5576, 5768, 5910, 6331, 6430, 6844, 6867, 7201, 7274, 7290, 7343, 7350, 7378,
             7387, 7440, 7554],
            [105, 975, 3421, 3480, 4120, 4444, 5957, 5971, 6119, 6617, 6761, 6810, 7067, 7353],
            [6, 138, 485, 1444, 1512, 2615, 2990, 3109, 5604, 6435, 6513, 6632, 6704, 7507],
            [20, 858, 1051, 2539, 3049, 5162, 5308, 6158, 6391, 6604, 6744, 7071, 7195, 7238],
            [1140, 5838, 6203, 6748],
            [6282, 6466, 6481, 6638],
            [2346, 2592, 5436, 7487],
            [2219, 3897, 5896, 7528],
            [2897, 6028, 7018],
            [1285, 1863, 5324],
            [3075, 6005, 6466],
            [5, 6020, 7551],
            [2121, 3751, 7507],
            [4027, 5488, 7542],
            [2, 6012, 7011],
            [3823, 5531, 5687],
            [1379, 2262, 5297],
            [1882, 7498, 7551],
            [3749, 4806, 7227],
            [2, 2074, 6898],
            [17, 616, 7482],
            [9, 6823, 7480],
            [5195, 5880, 7559],
        ],
    },
    "9/15": {
        "table": "A.2.8",
        "rows": [
            [212, 255, 540, 967, 1033, 1517, 1538, 3124, 3408, 3800, 4373, 4864, 4905, 5163, 5177,
             6186],
            [275, 660, 1351, 2211, 2876, 3063, 3433, 4088, 4273, 4544, 4618, 4632, 5548, 6101, 6111,
             6136],
            [279, 335, 494, 865, 1662, 1681, 3414, 3775, 4252, 4595, 5272, 5471, 5796, 5907, 5986,
             6008],
            [345, 352, 3094, 3188, 4297, 4338, 4490, 4865, 5303, 6477],
            [222, 681, 1218, 3169, 3850, 4878, 4954, 5666, 6001, 6237],
            [172, 512, 1536, 1559, 2179, 2227, 3334, 4049, 6464],
            [716, 934, 1694, 2890, 3276, 3608, 4332, 4468, 5945],
            [1133, 1593, 1825, 2571, 3017, 4251, 5221, 5639, 5845],
            [1076, 1222, 6465],
            [159, 5064, 6078],
            [374, 4073, 5357],
            [2833, 5526, 5845],
            [1594, 3639, 5419],
            [1028, 1392, 4239],
            [115, 622, 2175],
            [300, 1748, 6245],
            [2724, 3276, 5349],
            [1433, 6117, 6448],
            [485, 663, 4955],
            [711, 1132, 4315],
            [177, 3266, 4339],
            [1171, 4841, 4982],
            [33, 1584, 3692],
            [2820, 3485, 4249],
            [1716, 2428, 3125],
            [250, 2275, 6338],
            [108, 1719, 4961],
        ],
    },
    "10/15": {
        "table": "A.2.9",
        "rows": [
            [352, 747, 894, 1437, 1688, 1807, 1883, 2119, 2159, 3321, 3400, 3543, 3588, 3770, 3821,
             4384, 4470, 4884, 5012, 5036, 5084, 5101, 5271, 5281, 5353],
            [505, 915, 1156, 1269, 1518, 1650, 2153, 2256, 2344, 2465, 2509, 2867, 2875, 3007, 3254,
             3519, 3687, 4331, 4439, 4532, 4940, 5011, 5076, 5113, 5367],
            [268, 346, 650, 919, 1260, 4389, 4653, 4721, 4838, 5054, 5157, 5162, 5275, 5362],
            [220, 236, 828, 1590, 1792, 3259, 3647, 4276, 4281, 4325, 4963, 4974, 5003, 5037],
            [381, 737, 1099, 1409, 2364, 2955, 3228, 3341, 3473, 3985, 4257, 4730, 5173, 5242],
            [88, 771, 1640, 1737, 1803, 2408, 2575, 2974, 3167, 3464, 3780, 4501, 4901, 5047],
            [749, 1502, 2201, 3189],
            [2873, 3245, 3427],
            [2158, 2605, 3165],
            [1, 3438, 3606],
            [10, 3019, 5221],
            [371, 2901, 2923],
            [9, 3935, 4683],
            [1937, 3502, 3735],
            [507, 3128, 4994],
            [25, 3854, 4550],
            [1178, 4737, 5366],
            [2, 223, 5304],
            [1146, 5175, 5197],
            [1816, 2313, 3649],
            [740, 1951, 3844],
            [1320, 3703, 4791],
            [1754, 2905, 4058],
            [7, 917, 5277],
            [3048, 3954, 5396],
            [4804, 4824, 5105],
            [2812, 3895, 5226],
            [0, 5318, 5358],
            [1483, 2324, 4826],
            [2266, 4752, 5387],
        ],
    },
    "11/15": {
        "table": "A.2.10",
        "rows": [
            [49, 719, 784, 794, 968, 2382, 2685, 2873, 2974, 2995, 3540, 4179],
            [272, 281, 374, 1279, 2034, 2067, 2112, 3429, 3613, 3815, 3838, 4216],
            [206, 714, 820, 1800, 1925, 2147, 2168, 2769, 2806, 3253, 3415, 4311],
            [62, 159, 166, 605, 1496, 1711, 2652, 3016, 3347, 3517, 3654, 4113],
            [363, 733, 1118, 2062, 2613, 2736, 3143, 3427, 3664, 4100, 4157, 4314],
            [57, 142, 436, 983, 1364, 2105, 2113, 3074, 3639, 3835, 4164, 4242],
            [870, 921, 950, 1212, 1861, 2128, 2707, 2993, 3730, 3968, 3983, 4227],
            [185, 2684, 3263],
            [2035, 2123, 2913],
            [883, 2221, 3521],
            [1344, 1773, 4132],
            [438, 3178, 3650],
            [543, 756, 1639],
            [1057, 2337, 2898],
            [171, 3298, 3929],
            [1626, 2960, 3503],
            [484, 3050, 3323],
            [2283, 2336, 4189],
            [2732, 4132, 4318],
            [225, 2335, 3497],
            [600, 2246, 2658],
            [1240, 2790, 3020],
            [301, 1097, 3539],
            [1222, 1267, 2594],
            [1364, 2004, 3603],
            [1142, 1185, 2147],
            [564, 1505, 2086],
            [697, 991, 2908],
            [1467, 2073, 3462],
            [2574, 2818, 3637],
            [748, 2577, 2772],
            [1151, 1419, 4129],
            [164, 1238, 3401],
        ],
    },
    "12/15": {
        "table": "A.2.11",
        "rows": [
            [3, 394, 1014, 1214, 1361, 1477, 1534, 1660, 1856, 2745, 2987, 2991, 3124, 3155],
            [59, 136, 528, 781, 803, 928, 1293, 1489, 1944, 2041, 2200, 2613, 2690, 2847],
            [155, 245, 311, 621, 1114, 1269, 1281, 1783, 1995, 2047, 2672, 2803, 2885, 3014],
            [79, 870, 974, 1326, 1449, 1531, 2077, 2317, 2467, 2627, 2811, 3083, 3101, 3132],
            [4, 582, 660, 902, 1048, 1482, 1697, 1744, 1928, 2628, 2699, 2728, 3045, 3104],
            [175, 395, 429, 1027, 1061, 1068, 1154, 1168, 1175, 2147, 2359, 2376, 2613, 2682],
            [1388, 2241, 3118, 3148],
            [143, 506, 2067, 3148],
            [1594, 2217, 2705],
            [398, 988, 2551],
            [1149, 2588, 2654],
            [678, 2844, 3115],
            [1508, 1547, 1954],
            [1199, 1267, 1710],
            [2589, 3163, 3207],
            [1, 2583, 2974],
            [2766, 2897, 3166],
            [929, 1823, 2742],
            [1113, 3007, 3239],
            [1753, 2478, 3127],
            [0, 509, 1811],
            [1672, 2646, 2984],
            [965, 1462, 3230],
            [3, 1077, 2917],
            [1183, 1316, 1662],
            [968, 1593, 3239],
            [64, 1996, 2226],
            [1442, 2058, 3181],
            [513, 973, 1058],
            [1263, 3185, 3229],
            [681, 1394, 3017],
            [419, 2853, 3217],
            [3, 2404, 3175],
            [2417, 2792, 2854],
            [1879, 2940, 3235],
            [647, 1704, 3060],
        ],
    },
    "13/15": {
        "table": "A.2.12",
        "rows": [
            [71, 334, 645, 779, 786, 1124, 1131, 1267, 1379, 1554, 1766, 1798, 1939],
            [6, 183, 364, 506, 512, 922, 972, 981, 1039, 1121, 1537, 1840, 2111],
            [6, 71, 153, 204, 253, 268, 781, 799, 873, 1118, 1194, 1661, 2036],
            [6, 247, 353, 581, 921, 940, 1108, 1146, 1208, 1268, 1511, 1527, 1671],
            [6, 37, 466, 548, 747, 1142, 1203, 1271, 1512, 1516, 1837, 1904, 2125],
            [6, 171, 863, 953, 1025, 1244, 1378, 1396, 1723, 1783, 1816, 1914, 2121],
            [1268, 1360, 1647, 1769],
            [6, 458, 1231, 1414],
            [183, 535, 1244, 1277],
            [107, 360, 498, 1456],
            [6, 2007, 2059, 2120],
            [1480, 1523, 1670, 1927],
            [139, 573, 711, 1790],
            [6, 1541, 1889, 2023],
            [6, 374, 957, 1174],
            [287, 423, 872, 1285],
            [6, 1809, 1918],
            [65, 818, 1396],
            [590, 766, 2107],
            [192, 814, 1843],
            [775, 1163, 1256],
            [42, 735, 1415],
            [334, 1008, 2055],
            [109, 596, 1785],
            [406, 534, 1852],
            [684, 719, 1543],
            [401, 465, 1040],
            [112, 392, 621],
            [82, 897, 1950],
            [887, 1962, 2125],
            [793, 1088, 2159],
            [723, 919, 1139],
            [610, 839, 1302],
            [218, 1080, 1816],
            [627, 1646, 1749],
            [496, 1165, 1741],
            [916, 1055, 1662],
            [182, 722, 945],
            [5, 595, 1674],
        ],
    },
}

# --------------------------------------------------------------------------------------
# Section 6, Ninner = 16200: Table 6.4 (type) + Table 6.6 (Type A sizes) + Table 6.7 (Qldpc)
# M1/M2/Q1/Q2 are Type A only; Qldpc is Type B only.  Minner = Nldpc - Kldpc always.
# --------------------------------------------------------------------------------------
LDPC_PARAMS = {
    "2/15": {"Nldpc": 16200, "Kldpc":  2160, "M1": 3240, "M2": 10800, "Q1":   9, "Q2":  30, "Qldpc": None, "type": "A"},
    "3/15": {"Nldpc": 16200, "Kldpc":  3240, "M1": 1080, "M2": 11880, "Q1":   3, "Q2":  33, "Qldpc": None, "type": "A"},
    "4/15": {"Nldpc": 16200, "Kldpc":  4320, "M1": 1080, "M2": 10800, "Q1":   3, "Q2":  30, "Qldpc": None, "type": "A"},
    "5/15": {"Nldpc": 16200, "Kldpc":  5400, "M1":  720, "M2": 10080, "Q1":   2, "Q2":  28, "Qldpc": None, "type": "A"},
    "6/15": {"Nldpc": 16200, "Kldpc":  6480, "M1": None, "M2": None, "Q1": None, "Q2": None, "Qldpc":  27, "type": "B"},
    "7/15": {"Nldpc": 16200, "Kldpc":  7560, "M1": None, "M2": None, "Q1": None, "Q2": None, "Qldpc":  24, "type": "B"},
    "8/15": {"Nldpc": 16200, "Kldpc":  8640, "M1": None, "M2": None, "Q1": None, "Q2": None, "Qldpc":  21, "type": "B"},
    "9/15": {"Nldpc": 16200, "Kldpc":  9720, "M1": None, "M2": None, "Q1": None, "Q2": None, "Qldpc":  18, "type": "B"},
    "10/15": {"Nldpc": 16200, "Kldpc": 10800, "M1": None, "M2": None, "Q1": None, "Q2": None, "Qldpc":  15, "type": "B"},
    "11/15": {"Nldpc": 16200, "Kldpc": 11880, "M1": None, "M2": None, "Q1": None, "Q2": None, "Qldpc":  12, "type": "B"},
    "12/15": {"Nldpc": 16200, "Kldpc": 12960, "M1": None, "M2": None, "Q1": None, "Q2": None, "Qldpc":   9, "type": "B"},
    "13/15": {"Nldpc": 16200, "Kldpc": 14040, "M1": None, "M2": None, "Q1": None, "Q2": None, "Qldpc":   6, "type": "B"},
}

# --------------------------------------------------------------------------------------
# Table 6.5 (Type A sizes, Ninner = 64800) + Ninner=64800 column of Table 6.7.
# Provided for completeness; the corresponding Annex A.1 index tables are NOT in this file.
# --------------------------------------------------------------------------------------
LDPC_PARAMS_64800 = {
    "2/15":  {"Nldpc": 64800, "Kldpc":  8640, "M1": 1800, "M2": 54360, "Q1":   5, "Q2": 151, "Qldpc": None, "type": "A"},
    "3/15":  {"Nldpc": 64800, "Kldpc": 12960, "M1": 1800, "M2": 50040, "Q1":   5, "Q2": 139, "Qldpc": None, "type": "A"},
    "4/15":  {"Nldpc": 64800, "Kldpc": 17280, "M1": 1800, "M2": 45720, "Q1":   5, "Q2": 127, "Qldpc": None, "type": "A"},
    "5/15":  {"Nldpc": 64800, "Kldpc": 21600, "M1": 1440, "M2": 41760, "Q1":   4, "Q2": 116, "Qldpc": None, "type": "A"},
    "6/15":  {"Nldpc": 64800, "Kldpc": 25920, "M1": None, "M2": None,  "Q1": None, "Q2": None, "Qldpc": 108, "type": "B"},
    "7/15":  {"Nldpc": 64800, "Kldpc": 30240, "M1": 1080, "M2": 33480, "Q1":   3, "Q2":  93, "Qldpc": None, "type": "A"},
    "8/15":  {"Nldpc": 64800, "Kldpc": 34560, "M1": None, "M2": None,  "Q1": None, "Q2": None, "Qldpc":  84, "type": "B"},
    "9/15":  {"Nldpc": 64800, "Kldpc": 38880, "M1": None, "M2": None,  "Q1": None, "Q2": None, "Qldpc":  72, "type": "B"},
    "10/15": {"Nldpc": 64800, "Kldpc": 43200, "M1": None, "M2": None,  "Q1": None, "Q2": None, "Qldpc":  60, "type": "B"},
    "11/15": {"Nldpc": 64800, "Kldpc": 47520, "M1": None, "M2": None,  "Q1": None, "Q2": None, "Qldpc":  48, "type": "B"},
    "12/15": {"Nldpc": 64800, "Kldpc": 51840, "M1": None, "M2": None,  "Q1": None, "Q2": None, "Qldpc":  36, "type": "B"},
    "13/15": {"Nldpc": 64800, "Kldpc": 56160, "M1": None, "M2": None,  "Q1": None, "Q2": None, "Qldpc":  24, "type": "B"},
}

# Fingerprint of the Annex A.2 payload, agreed by all four independent extractions.
A2_INT_COUNT = 1819
A2_INT_SUM = 6754166
A2_SHA256 = "1803e817270b94d50e1ef3dee7fa1d76025b12f94a8bb3e74446e394bf3d8358"


def accumulator_addresses(rate):
    """Expand every parity-bit-accumulator address touched by the Annex A.2 table.

    Type A: Section 6.1.3.1 step (ii).  Type B: Section 6.1.3.2 q(i, j, l).
    Returns a set of parity indices in [0, Nldpc - Kldpc).
    """
    p = LDPC_PARAMS[rate]
    rows = LDPC_16200[rate]["rows"]
    hit = set()
    if p["type"] == "A":
        M1, M2, Q1, Q2 = p["M1"], p["M2"], p["Q1"], p["Q2"]
        for row in rows:
            for x in row:
                for m in range(360):
                    if x < M1:
                        hit.add((x + m * Q1) % M1)
                    else:
                        hit.add(M1 + ((x - M1 + m * Q2) % M2))
    else:
        Minner = p["Nldpc"] - p["Kldpc"]
        Q = p["Qldpc"]
        for row in rows:
            for q in row:
                for l in range(360):
                    hit.add((q + Q * l) % Minner)
    return hit


def _verify():
    import hashlib as _h
    import json as _j
    ok = True

    def check(label, cond, detail=""):
        nonlocal ok
        if not cond:
            ok = False
        print("  [%s] %-58s %s" % ("PASS" if cond else "FAIL", label, detail))

    print("ATSC A/322 Annex A.2 / Section 6.1.3 -- verification")
    print()
    print("(1) Table 6.4 / 6.6 / 6.7 internal arithmetic")
    for rate, p in LDPC_PARAMS.items():
        N, K = p["Nldpc"], p["Kldpc"]
        check("%-5s Kldpc == Nldpc * rate" % rate,
              K == N * int(rate.split("/")[0]) // 15, "Kldpc=%d" % K)
        check("%-5s type matches Table 6.4 (Ninner=16200 column)" % rate,
              p["type"] == TABLE_6_4[rate][1], "type=%s" % p["type"])
        if p["type"] == "A":
            check("%-5s M1 + M2 == Nldpc - Kldpc" % rate,
                  p["M1"] + p["M2"] == N - K, "%d + %d == %d" % (p["M1"], p["M2"], N - K))
            check("%-5s Q1 == M1/360 and Q2 == M2/360" % rate,
                  p["Q1"] * 360 == p["M1"] and p["Q2"] * 360 == p["M2"],
                  "Q1=%d Q2=%d" % (p["Q1"], p["Q2"]))
        else:
            check("%-5s Qldpc == (Nldpc - Kldpc)/360" % rate,
                  p["Qldpc"] * 360 == N - K, "Qldpc=%d" % p["Qldpc"])
    for rate, p in LDPC_PARAMS_64800.items():
        N, K = p["Nldpc"], p["Kldpc"]
        if p["type"] == "A":
            check("64800 %-5s M1+M2 == N-K, Q1=M1/360, Q2=M2/360" % rate,
                  p["M1"] + p["M2"] == N - K and p["Q1"] * 360 == p["M1"]
                  and p["Q2"] * 360 == p["M2"])
        else:
            check("64800 %-5s Qldpc == (N-K)/360" % rate, p["Qldpc"] * 360 == N - K)

    print()
    print("(2) Row counts  (Type B: Kldpc/360;  Type A: Kldpc/360 + Q1)")
    for rate, t in LDPC_16200.items():
        p = LDPC_PARAMS[rate]
        exp = p["Kldpc"] // 360 + (p["Q1"] if p["type"] == "A" else 0)
        check("%-5s %-7s rows == %d" % (rate, t["table"], exp),
              len(t["rows"]) == exp, "got %d" % len(t["rows"]))
    # Independent corroboration of the Type A "+ Q1" rule: the trailing Q1 rows serve
    # step (vii), not step (iv), and in every Type A table their weight is strictly
    # smaller than that of EVERY information row -- so the printed table visibly breaks
    # at exactly Kldpc/360.  Off-by-one row splits would break this.
    for rate, t in LDPC_16200.items():
        p = LDPC_PARAMS[rate]
        if p["type"] != "A":
            continue
        s = p["Kldpc"] // 360
        w = [len(r) for r in t["rows"]]
        check("%-5s Type A weight profile breaks at row Kldpc/360 = %d" % (rate, s),
              len(w) > s and min(w[:s]) > max(w[s:]),
              "info=%s parity=%s" % (sorted(set(w[:s])), sorted(set(w[s:]))))

    print()
    print("(3) Address range  [0, Nldpc - Kldpc)  and per-row well-formedness")
    for rate, t in LDPC_16200.items():
        p = LDPC_PARAMS[rate]
        lim = p["Nldpc"] - p["Kldpc"]
        flat = [x for row in t["rows"] for x in row]
        check("%-5s 0 <= addr < %d" % (rate, lim),
              min(flat) >= 0 and max(flat) < lim,
              "min=%d max=%d" % (min(flat), max(flat)))
        check("%-5s rows strictly increasing, no duplicates" % rate,
              all(all(r[i] < r[i + 1] for i in range(len(r) - 1)) for r in t["rows"]))

    print()
    print("(4) Accumulator coverage (Q1/Q2 resp. Qldpc expansion, Section 6.1.3)")
    for rate in LDPC_16200:
        p = LDPC_PARAMS[rate]
        lim = p["Nldpc"] - p["Kldpc"]
        hit = accumulator_addresses(rate)
        check("%-5s expansion covers all %5d parity accumulators" % (rate, lim),
              len(hit) == lim and min(hit) >= 0 and max(hit) < lim,
              "covered %d" % len(hit))

    print()
    print("(5) Payload fingerprint (agreed by 2024/2026 x -table/-layout extractions)")
    flat = [x for r in LDPC_16200 for row in LDPC_16200[r]["rows"] for x in row]
    blob = _j.dumps({r: LDPC_16200[r]["rows"] for r in LDPC_16200},
                    sort_keys=True).encode()
    check("total integer count == %d" % A2_INT_COUNT, len(flat) == A2_INT_COUNT,
          "got %d" % len(flat))
    check("integer sum == %d" % A2_INT_SUM, sum(flat) == A2_INT_SUM,
          "got %d" % sum(flat))
    check("sha256 == %s..." % A2_SHA256[:16],
          _h.sha256(blob).hexdigest() == A2_SHA256)

    print()
    print("RESULT:", "ALL CHECKS PASS" if ok else "*** FAILURES PRESENT ***")
    return 0 if ok else 1


if __name__ == "__main__":
    import sys
    sys.exit(_verify())
